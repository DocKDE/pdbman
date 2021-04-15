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

// #![allow(dead_code)]
#[macro_use]
extern crate clap;

#[macro_use]
extern crate prettytable;

use pdbtbx::StrictnessLevel;
use std::error::Error;
use std::fs;
use std::rc::Rc;

use crate::argparse::*;
use crate::functions::*;

pub mod argparse;
pub mod functions;

// Run function that handles the logic of when to call which function given an enum with the
// command line options. Hands all occurring errors to main.
pub fn run() -> Result<(), Box<dyn Error>> {
    let matches = parse_args()?;
    let mode = Rc::new(Mode::new(&matches)?);
    let filename = matches.value_of("INPUT").ok_or("No input file given")?;
    let mut pdb;

    match pdbtbx::open_pdb(filename, StrictnessLevel::Strict) {
        Ok((pdb_read, errors)) => {
            pdb = pdb_read;
            errors.iter().for_each(|x| println!("{}", x))
        }
        Err(errors) => {
            errors.iter().for_each(|x| println!("{}", x));
            return Err("Breaking error detected!".into());
        }
    }

    match *mode {
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
                Target::None => {
                    return Err("Please provide either the 'atoms' or 'residues' flag.".into())
                }
            }},
            Source::Sphere => {
                let sphere = Sphere::from_str(&matches, &mode, &pdb)?;

                let list = match target {
                    Target::Atoms => calc_atom_sphere(&pdb, sphere.origin, sphere.radius, false)?,
                    Target::Residues => calc_residue_sphere(&pdb, sphere.origin, sphere.radius, false)?,
                    Target::None => return Err("No target was given!".into()),
                };

                query_atoms(&pdb, list)?;
            }
            _ => println!("Please specifiy another input for a query."),
        },
        Mode::Analyze {
            region,
            target,
            distance,
        } => {
            analyze(&pdb, region, target)?;
            match distance {
                Distance::Clashes => find_contacts(&pdb, 0)?.printstd(),
                Distance::Contacts => find_contacts(&pdb, 1)?.printstd(),
                Distance::None => 0,
            };
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
            let edit_value = match *mode {
                Mode::Remove { .. } => 0.00,
                Mode::Add { .. } => match region {
                    Region::Active => 1.00,
                    Region::QM1 => 1.00,
                    Region::QM2 => 2.00,
                    Region::None => return Err("Region input is needed".into()),
                },
                _ => return Err("A mode must be selected".into()),
            };

            match source {
                Source::Infile => {
                    println!("This is not implemented yet.")
                }
                Source::List =>  {
                    let list = matches
                        .subcommand_matches(mode.to_string())
                        .ok_or("Something wrong with subcommand 'Add' or 'Remove'")?
                        .value_of("List")
                        .ok_or("Something wrong with option 'List'")?;

                match target {
                    Target::Atoms => {
                        let atomic_list = parse_atomic_list(list, &pdb)?;

                        match region {
                            Region::QM1 | Region::QM2 => edit_qm_atoms(
                                &mut pdb,
                                edit_value,
                                atomic_list,
                            )?,
                            Region::Active => edit_active_atoms(
                                &mut pdb,
                                edit_value,
                                atomic_list,
                            )?,
                            Region::None => {
                                return Err("Please give a region to add atoms or residues to.".into());
                            }
                        }
                    },
                    Target::Residues => {
                        let residue_list = parse_residue_list(list, &pdb)?;

                        match region {
                            Region::QM1 | Region::QM2 => edit_qm_residues(
                                &mut pdb,
                                edit_value,
                                residue_list,
                                partial,
                            )?,
                            Region::Active => edit_active_residues(
                                &mut pdb,
                                edit_value,
                                residue_list,
                                partial,
                            )?,
                            Region::None => return Err("Please give a region to modify.".into()),
                        }
                    },
                    Target::None => {
                        return Err("Please give either an 'Atoms' or 'Residues' flag.".into())
                    }
                }},
                Source::Sphere => {
                    let sphere = Sphere::from_str(&matches, &mode, &pdb)?;

                    let list = match target {
                        Target::Atoms => calc_atom_sphere(&pdb, sphere.origin, sphere.radius, true)?,
                        Target::Residues => calc_residue_sphere(&pdb, sphere.origin, sphere.radius, true)?,
                        Target::None => return Err("No target was given!".into()),
                    };

                    match region {
                        Region::QM1 | Region::QM2 => edit_qm_atoms(&mut pdb, edit_value, list)?,
                        Region::Active => edit_active_atoms(&mut pdb, edit_value, list)?,
                        Region::None => return Err("Please give a region.".into()),
                    }
                }
                Source::None => {
                    if mode.to_string() == "Remove"
                        && region == Region::None
                        && target == Target::None
                    {
                        remove_all(&mut pdb)?
                    } else {
                        return Err("Please provide a source (input file, list or sphere) for the addition of atoms or residues.".into());
                    }
                }
            }
        }
        Mode::None => {
            return Err("Please choose either 'Add', 'Analyze', 'Query' or 'Remove' mode using the respective subcommand.".into())
        }
    }

    if mode.to_string() == "Add" || mode.to_string() == "Remove" {
        if matches
            .subcommand_matches(mode.to_string())
            .ok_or("Something wrong with mode 'Add' or 'Remove'")?
            .is_present("Overwrite")
        {
            let filename_new = &(filename.to_string() + "_new");

            if let Err(e) = pdbtbx::save_pdb(pdb, filename_new, StrictnessLevel::Loose) {
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
                pdb,
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
