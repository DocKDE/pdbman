* write tests for argparsing

[x] minimize compilation of regexes
[x] extend verbose error printing to string list input
[x] improve handling of list input when atoms or residues are not present in PDB file
[x]   - propagate error from expand_residue to parse_residue_list and print all of it at once
[x] solve issue with duplicate residue numbers 
[x]   - add insert codes to PDB file (or use PDBx/mmCIF)
[x]   - parse insert codes from list inputs for residues (regex?)
[x]   - update number expansion function to incorporate insertion codes
[x]   - pass parsed list input at a vec of tuples (isize, Option<String>)
[x]   - rewrite residue functions to take vec of tuples and check for insertion code
[x] rewrite sphere finding functions with RTree
[x] add validation for sphere input
[x] write struct containing origin atom, radius, include_self and atom/residues for sphere jobs
[x] write funtion for processing sphere command line input
[x] find more elegant solution to mode.clone() in run() (maybe smart pointers?)
[x] find better solution to pdb.clone() in run()
[x] find out how to return error from within closure
[x] Accept ranges in command line list input
[x]   - write function to replace use of values_of_t
[x] move validation of command line arguments (i.e. List) to argparse module (clap validator_regex); for some reason clap doesn't recognize this so regular validator was used with separate function
[x] Implement querying Atom or Residue Names
[x] Implement sphere-based query
[x] change parsing of List input in add and remove modes
[x] move run to lib.rs
[x] improve error handling of reading pdb file
[x] Maybe refactor with subcommands
[x]   - add enums for modes
[x]   - write parsing function that returns an enum variant (based on mode chosen) which contains a struct of command line flags given (which in turn contain potential data)
[x] Write tests