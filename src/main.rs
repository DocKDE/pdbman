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

use std::process;
// use std::rc::Rc;

use anyhow::Result;
use clap::Arg;
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
// use options::Mode2;

fn run() -> Result<(), anyhow::Error> {
    let pdbman_match = app_from_crate!()
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
            bail!("Exiting...")
        }
    };

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
            Ok(c) => match c.as_str() {
                "exit" => break,
                "e" => break,
                "quit" => break,
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
        let matches = match args.try_get_matches_from(command.split_ascii_whitespace()) {
            Ok(m) => m,
            Err(e) => {
                println!("{}", e);
                continue;
            }
        };

        // let mode = Mode::new(&matches)?;
        // Error raised here are probably parse errors from faulty user input
        let mode = match Mode::new(&matches) {
            Ok(m) => m,
            // Ok(m) => Rc::new(m),
            Err(e) => {
                println!("{}", e);
                continue;
            }
        };
        // let mode = Rc::new(Mode::new(&matches)?);

        // Print errors instead of returning them
        if let Err(e) = dispatch(mode, &mut pdb, filename) {
            println! {"{}", e};
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
