use pdbman::run;
use std::process;

fn main() {
    // let matches = match parse_args() {
    //     Ok(result) => result,
    //     Err(e) => {
    //         eprintln!("{}", e);
    //         process::exit(1)
    //     }
    // };
    // let matches = parse_args()?;

    if let Err(e) = run() {
        eprintln!("{}", e);
        process::exit(1);
    }
}
