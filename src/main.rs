#![allow(clippy::float_cmp)]

#[macro_use]
extern crate clap;

// #[macro_use]
// extern crate prettytable;

#[macro_use]
extern crate lazy_static;

#[macro_use]
extern crate anyhow;

#[macro_use]
extern crate pest_derive;

mod dispatch;
mod functions;
mod help;
mod options;
mod residue_ascii;
mod revertable;
mod selection;
mod shell;

use std::env;
use std::fs;
use std::path::Path;
use std::process;

use anyhow::Context;
use anyhow::Result;
use clap::{AppSettings, Arg};
use colored::Colorize;
use pdbtbx::StrictnessLevel;
use rustyline::completion::FilenameCompleter;
use rustyline::config::OutputStreamType;
use rustyline::error::ReadlineError;
use rustyline::highlight::MatchingBracketHighlighter;
use rustyline::hint::HistoryHinter;
use rustyline::validate::MatchingBracketValidator;
use rustyline::{Cmd, ColorMode, CompletionType, Config, EditMode, Editor, KeyEvent};

use dispatch::dispatch;
use help::{HELP_INTER, HELP_LONG, HELP_SHORT};
use options::{clap_args, Mode};
use revertable::Revertable;
use shell::ShellHelper;

struct PDBCacher<T>
where
    T: Fn() -> Result<pdbtbx::PDB, anyhow::Error>,
{
    calculation: T,
    value: Option<Result<pdbtbx::PDB, anyhow::Error>>,
}

impl<T> PDBCacher<T>
where
    T: Fn() -> Result<pdbtbx::PDB, anyhow::Error>,
{
    fn new(calculation: T) -> PDBCacher<T> {
        PDBCacher {
            calculation,
            value: None,
        }
    }

    fn get_pdb(&mut self) -> &mut Result<pdbtbx::PDB, anyhow::Error> {
        if self.value.is_none() {
            self.value = Some((self.calculation)());
        }
        self.value.as_mut().unwrap()
    }
}

fn run() -> Result<(), anyhow::Error> {
    let pdbman_match = app_from_crate!()
        .setting(AppSettings::DisableVersionFlag)
        .setting(AppSettings::IgnoreErrors)
        .setting(AppSettings::NoAutoHelp)
        .arg(Arg::new("PDBFILE").help("Path to PDB file").required(true))
        .arg(
            Arg::new("Interactive")
                .help("Interactive Mode")
                .long("interactive")
                .short('i')
                .conflicts_with("File"),
        )
        .arg(
            Arg::new("File")
                .help("Read Commands from file")
                .long("file")
                .short('f')
                .takes_value(true)
                .conflicts_with("Interactive"),
        )
        .arg(
            Arg::new("Help")
                .help("Display help message")
                .long("help")
                .short('h'),
        )
        .get_matches();

    let given_args: Vec<String> = env::args().collect();

    // Print help if pdbman is called with no arguments
    if given_args.len() == 1 {
        println!("{}", HELP_LONG);
        return Ok(());
    }

    // Handle printing of help with and without input file given and potential subcommand help required.
    // The workaround with searching for the presence of "--help" is due to the presence of the "IgnoreErrors"
    // setting in clap which will prevent the flag from being detected if unknown options were given.
    if given_args.contains(&"-h".to_owned()) || given_args.contains(&"--help".to_owned()) {
        match given_args.len() {
            2 => {
                println!("{}", HELP_LONG);
            }
            3 => {
                if Path::new(pdbman_match.value_of("PDBFILE").unwrap()).exists() {
                    println!("{}", HELP_LONG);
                } else if let Err(e) = clap_args().try_get_matches_from(given_args.iter().skip(1)) {
                    println!("{}", e);
                }
            }
            _ => {
                let skip_val = if Path::new(pdbman_match.value_of("PDBFILE").unwrap()).exists() {
                    2
                } else {
                    1
                };

                if let Err(e) = clap_args().try_get_matches_from(given_args.iter().skip(skip_val)) {
                    println!("{}", e);
                }
            }
        }
        return Ok(());
    }

    // Must be handled explicitly because clap errors are disabled
    let filename = match pdbman_match.value_of("PDBFILE") {
        Some(s) => s,
        None => bail!("NO PDB FILE PATH WAS GIVEN!".red()),
    };

    // Define now so it can either be called or cached depending on use mode
    let read_pdb = || -> Result<pdbtbx::PDB, anyhow::Error> {
        match pdbtbx::open_pdb(filename, StrictnessLevel::Strict) {
            Ok((pdb_read, errors)) => {
                errors.iter().for_each(|x| println!("{}", x));
                Ok(pdb_read)
            }
            Err(errors) => {
                errors.iter().for_each(|x| println!("{}", x));
                bail!("EXITING...".red())
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

        let mut edit_ops: Vec<Box<dyn Revertable>> = Vec::new();
        let mut edit_ops_index = 0;
        let mut undone = false;

        // Be careful not to return any error unnecessarily because they would break the loop
        loop {
            let p = "\npdbman> ";
            rl.helper_mut().unwrap().colored_prompt = format!("\x1b[1;32m{}\x1b[0m", p);

            // Get command from read line and deal with special commands
            let command = match rl.readline(p) {
                Ok(c) => match c.as_str() {
                    "exit" | "e" | "quit" => break,
                    "hist" | "history" => {
                        rl.history().iter().for_each(|e| println!("{}", e));
                        continue;
                    }
                    "undo" => {
                        rl.add_history_entry("undo");
                        if edit_ops_index == 0 {
                            println!("Nothing to undo");
                        } else {
                            edit_ops[edit_ops_index - 1].undo(&mut pdb);
                            edit_ops_index -= 1;
                            undone = true;
                        }
                        continue;
                    }
                    "redo" => {
                        rl.add_history_entry("redo");
                        if edit_ops.len() == edit_ops_index {
                            println!("Nothing to redo");
                        } else {
                            edit_ops[edit_ops_index].redo(&mut pdb);
                            edit_ops_index += 1;
                        }
                        continue;
                    }
                    h @ ("-h" | "--help" | "help") => {
                        rl.add_history_entry(h);
                        println!("{}", HELP_INTER);
                        continue;
                    }
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

            // Don't return when an error occurs because it would break the loop and disrupt the workflow.
            // Errors returned from here mostly come from parsing the in-shell command line options.
            let matches = match clap_args().try_get_matches_from(command.split_whitespace()) {
                Ok(m) => m,
                Err(e) => {
                    println!("{}", e);
                    continue;
                }
            };

            // Error raised here are probably parse errors from faulty user input
            let mode = match Mode::new(&matches) {
                Ok(m) => m,
                Err(e) => {
                    println!("{}", e);
                    continue;
                }
            };

            match dispatch(&mode, &mut pdb, filename) {
                Ok(opt_op) => {
                    if let Some(edit_op) = opt_op {
                        // ensures the undo/redo tree remains via cutting stale branches
                        if undone {
                            edit_ops.truncate(edit_ops_index);
                            undone = false;
                        }
                        edit_ops.push(edit_op);
                        edit_ops_index += 1;
                    }
                }
                Err(e) => println!("{}", e),
            }
            // println!("Index: {}", edit_ops_index);
            // println!("Length: {}", edit_ops.len());
            // println!("{:#?}", edit_ops)
        }
    } else {
        let input;
        let args_env;

        let args = if pdbman_match.is_present("File") {
            let inpfile = match pdbman_match.value_of("File") {
                Some(s) => s,
                None => {
                    bail!(
                        "\n{}\n\n{}",
                        "NO INPUT FILE PATH FOR THE '--file' OPTION WAS GIVEN".red(),
                        HELP_SHORT
                    )
                }
            };

            input = fs::read_to_string(inpfile).with_context(|| {
                format!(
                    "\n{}: '{}'",
                    "FILE COULD NOT BE FOUND".red(),
                    inpfile.blue()
                )
            })?;

            input.trim().split('\n')
        } else {
            args_env = given_args[2..].join(" ");

            ensure!(
                !args_env.trim().is_empty(),
                "{}\n\n{}\n\n{}",
                "NO COMMAND ARGUMENTS WERE PROVIDED".red(),
                "If you want to enter interactive mode, provide the '-i' flag.",
                HELP_SHORT
            );

            args_env.split('/')
        };

        // More convenient so the args can be reused without cloning
        let args_vec: Vec<&str> = args.map(str::trim).collect();
        let mut pdb_cache = PDBCacher::new(read_pdb);

        // Test for input errors before actually processing anything
        for (i, arg) in args_vec.iter().enumerate() {
            let matches = match clap_args().try_get_matches_from(arg.split_whitespace()) {
                Ok(m) => m,
                Err(e) => bail!(
                    "\n{}{}: '{}'\n\n{}",
                    "FAILURE WHILE PARSING COMMAND #".red(),
                    (i + 1).to_string().red(),
                    arg.blue(),
                    e
                ),
            };

            if let Err(e) = Mode::new(&matches) {
                bail!(
                    "\n{}{}: '{}'\n\n{}",
                    "FAILURE WHILE PARSING COMMAND #".red(),
                    (i + 1).to_string().red(),
                    arg.blue(),
                    e
                )
            };
        }

        // Do the processing now that all inputs have been checked
        for arg in args_vec {
            let matches = clap_args().get_matches_from(arg.split_whitespace());
            let mode = Mode::new(&matches).unwrap();

            let pdb = match pdb_cache.get_pdb().as_mut() {
                Ok(p) => p,
                Err(e) => bail!(e.to_string()),
            };

            if let Err(e) = dispatch(&mode, pdb, filename) {
                bail!(
                    "\n{}: '{}'\n\n{}",
                    "ERROR DURING PROCESSING OF INPUT".red(),
                    arg.blue(),
                    e
                )
            };
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
