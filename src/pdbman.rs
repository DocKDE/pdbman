use std::process;
use pdbman::run;

fn main() {
    if let Err(e) = run() {
        eprintln!("{}", e);
        process::exit(1);
    }
}
