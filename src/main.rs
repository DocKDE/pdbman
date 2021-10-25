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

#[macro_use]
extern crate anyhow;

mod dispatch;
mod functions;
mod options;
mod residue_ascii;
mod shell;

use std::error::Error;
use std::process;

use clap::{App, Arg};
use pdbtbx::StrictnessLevel;
use rustyline::completion::FilenameCompleter;
use rustyline::config::OutputStreamType;
use rustyline::error::ReadlineError;
use rustyline::highlight::MatchingBracketHighlighter;
use rustyline::hint::HistoryHinter;
use rustyline::validate::MatchingBracketValidator;
use rustyline::{Cmd, ColorMode, CompletionType, Config, EditMode, Editor, KeyEvent};

use dispatch::dispatch;
use options::*;
use shell::*;

fn run() -> Result<(), Box<dyn Error>> {
    let pdbman_match = App::new(crate_name!())
        .about(crate_description!())
        .version(crate_version!())
        .author(crate_authors!())
        .arg(Arg::new("PDBFILE").about("Path to PDB file").required(true))
        .get_matches();

    // Must be present, otherwise clap would complain, thus `unwrap()` is ok.
    let filename = pdbman_match.value_of("PDBFILE").unwrap();

    let mut pdb = match pdbtbx::open_pdb(filename, StrictnessLevel::Strict) {
        Ok((pdb_read, errors)) => {
            errors.iter().for_each(|x| println!("{}", x));
            pdb_read
        }
        Err(errors) => {
            errors.iter().for_each(|x| println!("{}", x));
            return Err("Exiting".into());
        }
    };

    let mut parse = false;
    env_logger::init();
    let config = Config::builder()
        .history_ignore_space(true)
        .history_ignore_dups(true)
        .completion_type(CompletionType::List)
        .edit_mode(EditMode::Emacs)
        .output_stream(OutputStreamType::Stdout)
        .color_mode(ColorMode::Forced)
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

        let command = match rl.readline(p) {
            Ok(c) => c,
            Err(ReadlineError::Interrupted) => {
                println!("CTRL-C");
                break;
            }
            Err(ReadlineError::Eof) => {
                println!("CTRL-D");
                break;
            }
            Err(e) => return Err(Box::new(e)),
        };

        if command == "exit" || command == "quit" || command == "e" {
            break;
        }

        let args = parse_args();
        // Don't return when an error occurs because it would break the loop and disrupt the workflow.
        // Errors returned from here mostly come from parsing the in-shell command line options.
        let matches = match args.try_get_matches_from(command.split_ascii_whitespace()) {
            Ok(m) => {
                rl.add_history_entry(command);
                m
            }
            Err(e) => {
                println!("{}", e);
                continue;
            }
        };

        // If something goes wrong here, it should actually return
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
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    if let Err(e) = run() {
        eprintln!("{}", e);
        process::exit(1);
    }
    Ok(())
}
