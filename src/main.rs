use pdbman::run;
use std::process;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    if let Err(e) = run() {
        eprintln!("{}", e);
        process::exit(1);
    }
    Ok(())
}
