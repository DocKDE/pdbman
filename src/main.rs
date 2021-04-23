use pdbman::run;
use std::error::Error;
use std::process;

fn main() -> Result<(), Box<dyn Error>> {
    if let Err(e) = run() {
        eprintln!("{}", e);
        process::exit(1);
    }
    Ok(())
}
