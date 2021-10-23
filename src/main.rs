//! This crate serves as a command line utility to handle PDB files for use with the Quantum
//! Chemistry package ORCA. Its capabilities include querying the file for residues and atoms,
//! analyzing the defined QM and Active regions (according to ORCA nomenclature using the
//! occupancy and B values) and editing these same regions so the PDB file can be used as input
//! in a QM/MM calculation.
//!
//! # Example usage:
//! Query for atoms with the name 'Cu':
//!
//! `pdbman myfile.pdb Query -tl Cu`
//!
//! Analyze QM residues:
//!
//! `pdbman myfile.pdb Analyze -rq`
//!
//! Remove all QM and Active atoms and overwrite input file:
//!
//! `pdbman myfile.pdb Remove -w`
//!
//! Add atoms in a sphere around a given atom to Active region:
//!
//! `pdbman myfile.pdb Add -raws 2589 10`
//!
//! The `--sphere` or `-s` flag takes an atom ID and a radius in Angstrom as arguments.

#![allow(clippy::float_cmp)]

#[macro_use]
extern crate clap;

#[macro_use]
extern crate prettytable;

#[macro_use]
extern crate lazy_static;

mod dispatch;
mod functions;
mod options;
mod residue_ascii;
mod shell;

use std::error::Error;
use std::process;

use pdbtbx::StrictnessLevel;
use rustyline::highlight::MatchingBracketHighlighter;
use rustyline::validate::MatchingBracketValidator;
use rustyline::completion::FilenameCompleter;
use rustyline::hint::HistoryHinter;
use rustyline::config::OutputStreamType;
use rustyline::{Cmd, CompletionType, Config, EditMode, Editor, KeyEvent};

use options::*;
use dispatch::dispatch;
use shell::*;

fn run() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().collect();

    if args.len() == 1 {
        return Err("Please give a PDB file as argument!".into());
    } else if args.len() > 2 {
        return Err("Please give only one PDB file as argument!".into());
    }

    let filename = args[1].as_str();

    let mut pdb;

    match pdbtbx::open_pdb(filename, StrictnessLevel::Strict) {
        Ok((pdb_read, errors)) => {
            pdb = pdb_read;
            errors.iter().for_each(|x| println!("{}", x))
        }
        Err(errors) => {
            errors.iter().for_each(|x| println!("{}", x));
            return Err("Exiting".into());
        }
    }

    // if check_residue_overflow(&pdb) {
    //     let filename_insert = filename.to_string() + "_insert";
    //     if !std::path::Path::new(&filename_insert).exists() {
    //         eprintln!(
    //             "WARNING: More than 9999 residues present and not all of them have insertion codes.\n\
    //         PDB file with appropriate insertion code '{}_insert' will now be created.\n\
    //         Please use this file, otherwise the correctness of the results cannot be guaranteed.\n",
    //             filename
    //         );
    //         add_insertion_codes(&mut pdb)?;
    //         if let Err(e) = pdbtbx::save_pdb(pdb.clone(), &filename_insert, StrictnessLevel::Loose)
    //         {
    //             e.iter().for_each(|x| println!("{}", x))
    //         };
    //     } else {
    //         eprintln!(
    //             "WARNING: More than 9999 residues present and not all of them have insertion codes.\n\
    //         All output generated with this file is untrustworthy!\n\
    //         Please use existing PDB file '{}_insert' instead!\n", filename
    //         );
    //     }
    // }

    let mut parse = false;
    env_logger::init();
    let config = Config::builder()
        .history_ignore_space(true)
        .history_ignore_dups(true)
        .completion_type(CompletionType::List)
        .edit_mode(EditMode::Emacs)
        .output_stream(OutputStreamType::Stdout)
        .build();
    
    let helper = ShellHelper {
        completer: FilenameCompleter::new(),
        highlighter: MatchingBracketHighlighter::new(),
        hinter: HistoryHinter {},
        colored_prompt: "".to_owned(),
        validator: MatchingBracketValidator::new(),
    };

    let mut rl = Editor::with_config(config);
    rl.set_helper(Some(helper));
    rl.bind_sequence(KeyEvent::alt('n'), Cmd::HistorySearchForward);
    rl.bind_sequence(KeyEvent::alt('p'), Cmd::HistorySearchBackward);

    // Be careful not to return any error unnecessarily because they would break the loop
    loop {
        let p = "\npdbman> ";
        rl.helper_mut().expect("No helper").colored_prompt = format!("\x1b[1;32m{}\x1b[0m", p);
        let command = rl.readline(&p)?;
        rl.add_history_entry(command.as_str());

        if command == "exit" || command == "quit" || command == "q" {
            break;
        }

        let args = parse_args();

        // Don't return when an error occurs because it would break the loop and disrupt the workflow.
        // Errors returned from here mostly come from parsing the in-shell command line options.
        let matches = match args.try_get_matches_from(command.split_ascii_whitespace()) {
            Ok(m) => m,
            Err(e) => {
                println!("{}", e);
                continue;
            }
        };

        // If something goes wrong here it should actually return
        let mode = Mode::new(&matches)?;

        // Print errors instead of returning them
        if let Err(e) = dispatch(matches, mode, &mut pdb) {
            println! {"{}", e};
        }
        if mode.to_string() == "Add" || mode.to_string() == "Remove" {
            parse = true;
        }
    }

    // Only print to file if any changes were made and only once the loop was broken
    if parse {
        let filename_new = filename.to_string() + "_new";

        println!("Saving changes to {}", filename_new);
        if let Err(e) = pdbtbx::save_pdb(pdb, &filename_new, StrictnessLevel::Loose) {
            e.iter().for_each(|x| println!("{}", x))
        };

        // fs::remove_file(filename)?;
        // fs::rename(filename_new, filename)?;
    }

    // Figure out a way to print to another file on demand and not
    // necessarily overwrite the existing one
    // if mode.to_string() == "Add" || mode.to_string() == "Remove" {
    //     if matches
    //         .subcommand_matches(mode.to_string())
    //         .ok_or("Something wrong with mode 'Add' or 'Remove'")?
    //         .is_present("Overwrite")
    //     {
    //         let filename_new = filename.to_string() + "_new";

    //         if let Err(e) = pdbtbx::save_pdb(pdb.clone(), &filename_new, StrictnessLevel::Loose) {
    //             e.iter().for_each(|x| println!("{}", x))
    //         };

    //         fs::remove_file(filename)?;
    //         fs::rename(filename_new, filename)?;
    //     } else if matches
    //         .subcommand_matches(mode.to_string())
    //         .ok_or("Something wrong with mode 'Add' or 'Remove'")?
    //         .is_present("Outfile")
    //     {
    //         if let Err(e) = pdbtbx::save_pdb(
    //             pdb.clone(),
    //             matches
    //                 .subcommand_matches(mode.to_string())
    //                 .ok_or("Something wrong with mode 'Add' or 'Remove'")?
    //                 .value_of("Outfile")
    //                 .ok_or("Value for Outfile could not be parsed")?,
    //             StrictnessLevel::Loose,
    //         ) {
    //             e.iter().for_each(|x| println!("{}", x))
    //         }
    //     } else {
    //         print_to_stdout(&pdb)?;
    //     }
    // }
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    if let Err(e) = run() {
        eprintln!("{}", e);
        process::exit(1);
    }
    Ok(())
}
