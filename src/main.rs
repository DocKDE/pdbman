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

mod functions;
mod options;
mod residue_ascii;

use clap::App;
use pdbtbx::StrictnessLevel;
use std::error::Error;
use std::fs;
use std::io::Write;
use std::process;

use functions::*;
use options::*;

fn dispatch() -> Result<(), Box<dyn Error>> {
    let args = parse_args();
    let matches = args.clone().get_matches();
    let mode = Mode::new(&matches)?;
    let filename = matches.value_of("INPUT").unwrap();
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

    if check_residue_overflow(&pdb) {
        eprintln!(
            "WARNING: More than 9999 residues present and not all of them have insertion codes.\n\
        All output generated with this file is untrustworthy! Please add appropriate insertion codes!\n",
        );
    }

    if matches.is_present("interactive") {
        let new_args = args
            .mut_arg("INPUT", |a| a.required(false).hidden(true))
            .mut_arg("interactive", |a| a.exclusive(true).hidden(true))
            .setting(clap::AppSettings::NoBinaryName);
        interactive(new_args, filename, &mut pdb)?;
    } else {
        run(matches.clone(), mode, filename, &mut pdb)?;
        // println!("{:?}", matches);
    }
    Ok(())
}

fn interactive(args: App, filename: &str, pdb:&mut pdbtbx::PDB) -> Result<(), Box<dyn Error>> {
    'outer: loop {
        print!("\npdbman> ");
        std::io::stdout().flush()?;
        let command: String = text_io::read!("{}\n");

        if command == "exit" || command == "quit" || command == "q" {
            break 'outer;
        }

        let matches = args
            .clone()
            .get_matches_from(command.split_ascii_whitespace());
        let mode = Mode::new(&matches)?;
        run(matches, mode, filename, pdb)?;
    }
    Ok(())
}

// Run function that handles the logic of when to call which function given an enum with the
// command line options. Hands all occurring errors to main.
fn run(
    matches: clap::ArgMatches,
    mode: Mode,
    filename: &str,
    mut pdb: &mut pdbtbx::PDB,
) -> Result<(), Box<dyn Error>> {
    // let matches = parse_args().get_matches();
    // let mode = Mode::new(&matches)?;
    // let filename = matches.value_of("INPUT").unwrap();
    // let mut pdb;

    // let match_string = "pdbman myfile.pdb Q -tl 1";
    // let str_matches = parse_args().get_matches_from(match_string.split_ascii_whitespace());

    // match pdbtbx::open_pdb(filename, StrictnessLevel::Strict) {
    //     Ok((pdb_read, errors)) => {
    //         pdb = pdb_read;
    //         errors.iter().for_each(|x| println!("{}", x))
    //     }
    //     Err(errors) => {
    //         errors.iter().for_each(|x| println!("{}", x));
    //         return Err("Exiting".into());
    //     }
    // }

    if check_residue_overflow(&pdb) {
        eprintln!(
            "WARNING: More than 9999 residues present and not all of them have insertion codes.\n\
        All output generated with this file is untrustworthy! Please add appropriate insertion codes!\n",
        );

        // let filename_insert = filename.to_string() + "_insert";

        // if !std::path::Path::new(&filename_insert).exists() {
        //     eprintln!(
        //         "WARNING: More than 9999 residues present and not all of them have insertion codes.\n\
        //     All output generated with this file is untrustworthy! Please add appropriate insertion codes!\n",
        //     );
        // add_insertion_codes(&mut pdb)?;
        // if let Err(e) = pdbtbx::save_pdb(pdb.clone(), &filename_insert, StrictnessLevel::Loose)
        // {
        //     e.iter().for_each(|x| println!("{}", x))
        // };
        // } else {
        //     println!(
        //         "WARNING: More than 9999 residues present and not all of them have insertion codes.\n\
        //     PDB file with appropriate insertion code '{}_insert' already exists!\n\
        //     Please use this file, otherwise the correctness of the results cannot be guaranteed.\n",
        //         filename
        //     );
        // }
    }

    match mode {
        Mode::Query { source, target } => match source {
            Source::List => {
                let list = matches
                    .subcommand_matches("Query")
                    .ok_or("Something wrong with subcommand 'Query'")?
                    .value_of("List")
                    .ok_or("Something wrong with option 'List'")?;
                match target {
                    Target::Atoms => query_atoms(&pdb, parse_atomic_list(list, &pdb)?)?,
                    Target::Residues => query_residues(&pdb, parse_residue_list(list, &pdb)?)?,
                    Target::None => unreachable!(),
                }
            }
            Source::Sphere => {
                let sphere = Sphere::new(
                    matches
                        .subcommand_matches("Query")
                        .ok_or("Something wrong with option 'Query'")?
                        .values_of("Sphere")
                        .ok_or("Something wrong with option 'Sphere'")?,
                    &pdb,
                )?;

                let list = match target {
                    Target::Atoms => calc_atom_sphere(&pdb, &sphere.origin, sphere.radius, false)?,
                    Target::Residues => {
                        calc_residue_sphere(&pdb, &sphere.origin, sphere.radius, false)?
                    }
                    Target::None => unreachable!(),
                };

                query_atoms(&pdb, list)?;
            }
            _ => return Err("Please specify another input for a query.".into()),
        },
        Mode::Analyze {
            region,
            target,
            distance,
        } => {
            analyze(&pdb, region, target)?;
            if distance == Distance::Clashes || distance == Distance::Contacts {
                find_contacts(&pdb, distance)?.printstd();
            }
            // match distance {
            //     Distance::Clashes | Distance::Contacts => find_contacts(&pdb, distance)?.printstd(),
            //     Distance::None => 0,
            // };
        }
        Mode::Add {
            region,
            source,
            target,
            partial,
            output: _,
        }
        | Mode::Remove {
            region,
            source,
            target,
            partial,
            output: _,
        } => {
            let edit_value = match mode {
                Mode::Remove { .. } => 0.00,
                Mode::Add { .. } => match region {
                    Region::Active => 1.00,
                    Region::QM1 => 1.00,
                    Region::QM2 => 2.00,
                    Region::None => unreachable!(),
                },
                _ => unreachable!(),
            };

            match source {
                Source::Infile => {
                    println!("This is not implemented yet.")
                }
                Source::List => {
                    let list = matches
                        .subcommand_matches(mode.to_string())
                        .ok_or("Something wrong with subcommand 'Add' or 'Remove'")?
                        .value_of("List")
                        .ok_or("Something wrong with option 'List'")?;

                    match target {
                        Target::Atoms => {
                            let atomic_list = parse_atomic_list(list, &pdb)?;

                            match region {
                                Region::QM1 | Region::QM2 => {
                                    edit_qm_atoms(&mut pdb, edit_value, atomic_list)?
                                }
                                Region::Active => {
                                    edit_active_atoms(&mut pdb, edit_value, atomic_list)?
                                }
                                Region::None => unreachable!(),
                            }
                        }
                        Target::Residues => {
                            let pdb_clone = pdb.clone();
                            let residue_list = parse_residue_list(list, &pdb_clone)?;

                            match region {
                                Region::QM1 | Region::QM2 => {
                                    edit_qm_residues(&mut pdb, edit_value, residue_list, partial)?
                                }
                                Region::Active => edit_active_residues(
                                    &mut pdb,
                                    edit_value,
                                    residue_list,
                                    partial,
                                )?,
                                Region::None => unreachable!(),
                            }
                        }
                        Target::None => unreachable!(),
                    }
                }
                Source::Sphere => {
                    let sphere = Sphere::new(
                        matches
                            .subcommand_matches(mode.to_string())
                            .ok_or("Something wrong with option 'Add' or 'Remove'")?
                            .values_of("Sphere")
                            .ok_or("Something wrong with option 'Sphere'")?,
                        &pdb,
                    )?;

                    let list = match target {
                        Target::Atoms => {
                            calc_atom_sphere(&pdb, &sphere.origin, sphere.radius, true)?
                        }
                        Target::Residues => {
                            calc_residue_sphere(&pdb, &sphere.origin, sphere.radius, true)?
                        }
                        Target::None => unreachable!(),
                    };

                    match region {
                        Region::QM1 | Region::QM2 => edit_qm_atoms(&mut pdb, edit_value, list)?,
                        Region::Active => edit_active_atoms(&mut pdb, edit_value, list)?,
                        Region::None => unreachable!(),
                    }
                }
                Source::None => {
                    if mode.to_string() == "Remove"
                        && region == Region::None
                        && target == Target::None
                    {
                        remove_all(&mut pdb)?
                    } else if { region == Region::QM1 || region == Region::QM2 }
                        && target == Target::None
                    {
                        remove_qm(&mut pdb)?
                    } else if region == Region::Active && target == Target::None {
                        remove_active(&mut pdb)?
                    } else {
                        return Err("Please provide the approprate options (see --help).".into());
                    }
                }
            }
        }
        Mode::None => {
            unreachable!()
        }
    }

    if mode.to_string() == "Add" || mode.to_string() == "Remove" {
        if matches
            .subcommand_matches(mode.to_string())
            .ok_or("Something wrong with mode 'Add' or 'Remove'")?
            .is_present("Overwrite")
        {
            let filename_new = filename.to_string() + "_new";

            if let Err(e) = pdbtbx::save_pdb(pdb.clone(), &filename_new, StrictnessLevel::Loose) {
                e.iter().for_each(|x| println!("{}", x))
            };

            fs::remove_file(filename)?;
            fs::rename(filename_new, filename)?;
        } else if matches
            .subcommand_matches(mode.to_string())
            .ok_or("Something wrong with mode 'Add' or 'Remove'")?
            .is_present("Outfile")
        {
            if let Err(e) = pdbtbx::save_pdb(
                pdb.clone(),
                matches
                    .subcommand_matches(mode.to_string())
                    .ok_or("Something wrong with mode 'Add' or 'Remove'")?
                    .value_of("Outfile")
                    .ok_or("Value for Outfile could not be parsed")?,
                StrictnessLevel::Loose,
            ) {
                e.iter().for_each(|x| println!("{}", x))
            }
        } else {
            print_to_stdout(&pdb)?;
        }
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    if let Err(e) = dispatch() {
        eprintln!("{}", e);
        process::exit(1);
    }
    Ok(())
}
