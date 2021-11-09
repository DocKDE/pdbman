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

use std::env;
use std::fs;
use std::process;

use anyhow::Context;
use anyhow::Result;
use clap::{AppSettings, Arg};
use colored::*;
use itertools::Itertools;
use pdbtbx::StrictnessLevel;
use rustyline::completion::FilenameCompleter;
use rustyline::config::OutputStreamType;
use rustyline::error::ReadlineError;
use rustyline::highlight::MatchingBracketHighlighter;
use rustyline::hint::HistoryHinter;
use rustyline::validate::MatchingBracketValidator;
use rustyline::{Cmd, ColorMode, CompletionType, Config, EditMode, Editor, KeyEvent};

use dispatch::dispatch;
use options::{parse_args, Mode};
use shell::ShellHelper;

const HELP_LONG: &str = "
pdbman 0.8.0
    
Benedikt M. Flöser <benedikt.floeser@cec.mpg.de>
    
Analyzes and edits PDB files for usage in QM/MM calculations with the ORCA Quantum Chemistry package
    
USAGE:
    pdbman <PDBFILE> <[OPTIONS]|[SUBCOMMAND]>
    
ARGS:
    <PDBFILE>    Path to PDB file
    \
OPTIONS:
    -f, --file <File>    Read commands from file
    -h, --help           Display help message
    -i, --interactive    Enter interactive mode

SUBCOMMANDS:
    Add                  Add atoms or residues to QM1/QM2/Active region
    Remove               Remove atoms or residues from QM1/QM2/Active region
    Query                Query atoms or residues
    Analyse              Analyze PDB file and QM1/QM2/Active region
    Write                Write changes to file or stdout

Calling a subcommand with the '--help/-h' flag will display a help message for it
";

const HELP_SHORT: &str = "
USAGE:
    pdbman <PDBFILE> <[OPTIONS]|[SUBCOMMAND]>
    
ARGS:
    <PDBFILE>    Path to PDB file
    \
OPTIONS:
    -f, --file <File>    Read commands from file
    -h, --help           Display help message
    -i, --interactive    Enter interactive mode

SUBCOMMANDS:
    Add                  Add atoms or residues to QM1/QM2/Active region
    Remove               Remove atoms or residues from QM1/QM2/Active region
    Query                Query atoms or residues
    Analyse              Analyze PDB file and QM1/QM2/Active region
    Write                Write changes to file or stdout

Calling a subcommand with the '--help/-h' flag will display a help message for it
";

fn run() -> Result<(), anyhow::Error> {
    let pdbman_match = app_from_crate!()
        .setting(AppSettings::DisableVersionFlag)
        .setting(AppSettings::IgnoreErrors)
        // .setting(AppSettings::AllowExternalSubcommands)
        .setting(AppSettings::NoAutoHelp)
        .arg(Arg::new("PDBFILE").about("Path to PDB file").required(true))
        .arg(
            Arg::new("Interactive")
                .about("Interactive Mode")
                .long("interactive")
                .short('i')
                .conflicts_with("File"),
        )
        .arg(
            Arg::new("File")
                .about("Read Commands from file")
                .long("file")
                .short('f')
                .takes_value(true)
                .conflicts_with("Interactive"),
        )
        .arg(
            Arg::new("Help")
                .about("Display help message")
                .long("help")
                .short('h'),
        )
        .get_matches();

    if pdbman_match.is_present("Help") {
        println!("{}", HELP_LONG);
        return Ok(());
    }

    // Must be present, otherwise clap would complain, thus `unwrap()` is ok.
    // let filename = pdbman_match.value_of("PDBFILE").unwrap();
    let filename = match pdbman_match.value_of("PDBFILE") {
        Some(s) => s,
        None => bail!("No PDB file path was given!".red()),
    };

    let read_pdb = || -> Result<pdbtbx::PDB, anyhow::Error> {
        match pdbtbx::open_pdb(filename, StrictnessLevel::Strict) {
            Ok((pdb_read, errors)) => {
                errors.iter().for_each(|x| println!("{}", x));
                Ok(pdb_read)
            }
            Err(errors) => {
                errors.iter().for_each(|x| println!("{}", x));
                bail!("Exiting...".red())
            }
        }
    };

    if pdbman_match.is_present("Interactive") {
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

        let mut pdb = read_pdb()?;

        // Be careful not to return any error unnecessarily because they would break the loop
        loop {
            let p = "\npdbman> ";
            rl.helper_mut().expect("No helper").colored_prompt = format!("\x1b[1;32m{}\x1b[0m", p);

            let command = match rl.readline(p) {
                Ok(c) => match c.as_str() {
                    "exit" => break,
                    "e" => break,
                    "quit" => break,
                    "" => continue,
                    _ => c,
                },
                Err(ReadlineError::Interrupted) => {
                    println!("CTRL-C");
                    break;
                }
                Err(ReadlineError::Eof) => {
                    println!("CTRL-D");
                    break;
                }
                Err(e) => bail!(e),
            };

            rl.add_history_entry(&command);

            let args = parse_args();
            // Don't return when an error occurs because it would break the loop and disrupt the workflow.
            // Errors returned from here mostly come from parsing the in-shell command line options.
            let matches = match args.try_get_matches_from(command.split_whitespace()) {
                Ok(m) => m,
                Err(e) => {
                    println!("{}", e);
                    continue;
                }
            };

            // Error raised here are probably parse errors from faulty user input
            let mode = match Mode::new(&matches) {
                Ok(m) => m,
                // Ok(m) => Rc::new(m),
                Err(e) => {
                    println!("{}", e);
                    continue;
                }
            };

            // Print errors instead of returning them
            if let Err(e) = dispatch(mode, &mut pdb, filename) {
                println! {"{}", e};
            }
        }
    } else if pdbman_match.is_present("File") {
        let inpfile = match pdbman_match.value_of("File") {
            Some(s) => s,
            None => {
                bail!(
                    "{}\n\n{}",
                    "No input file path for the '--file' option was given!".red(),
                    HELP_SHORT
                )
            }
        };
        let input = fs::read_to_string(inpfile)
            .with_context(|| format!("File '{}' could not be found", inpfile).red())?;
        let args = input.split('/');

        let mut pdb = read_pdb()?;

        for arg in args {
            let matches = parse_args().try_get_matches_from(arg.trim().split_whitespace())?;
            let mode = Mode::new(&matches)?;
            dispatch(mode, &mut pdb, filename)?
        }
    } else {
        let args = env::args().skip(2).join(" ");

        ensure!(
            !args.trim().is_empty(),
            "{}\n\n{}",
            "No actionable arguments were provided!".red(),
            HELP_SHORT
        );

        let mut pdb = read_pdb()?;
        let args_split = args.split('/');

        for arg in args_split {
            let app = parse_args();
            let matches = app.try_get_matches_from(arg.trim().split_whitespace())?;
            let mode = Mode::new(&matches)?;
            dispatch(mode, &mut pdb, filename)?
        }
    }

    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("{}", e);
        process::exit(1);
    }
}
