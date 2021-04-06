// #![allow(dead_code)]
#[macro_use]
extern crate clap;

#[macro_use]
extern crate prettytable;

use std::error::Error;
use std::fs;
// use clap::{load_yaml, App};
use pdbtbx::StrictnessLevel;

use crate::argparse::*;
use crate::functions::*;

pub mod argparse;
pub mod functions;

type Result<T> = std::result::Result<T, Box<dyn Error>>;

pub fn run() -> Result<()> {
    let matches = parse_args()?;
    let mode = Mode::new(&matches)?;
    let filename = matches.value_of("INPUT").unwrap();
    let mut pdb;
    // println!("{:?}", matches.subcommand_matches("Add").unwrap().value_of("List").unwrap());

    match pdbtbx::open_pdb(filename, StrictnessLevel::Strict) {
        Ok(result) => pdb = result.0,
        Err(e) => {
            let err_string: String = e
                .iter()
                .map(|x| x.short_description())
                .collect::<Vec<&str>>()
                .join("\n");
            return Err(err_string.into());
        }
    }

    for i in pdbtbx::validate_pdb(&pdb).iter() {
        println!("{}", i)
    }

    match mode.clone() {
        Mode::Query { source, target } => match source {
            Source::List => match target {
                Target::Atoms => query_atoms(
                    &pdb,
                    parse_atomic_list(
                        matches
                            .subcommand_matches("Query")
                            .unwrap()
                            .value_of("List")
                            .unwrap(),
                        &pdb,
                    )?,
                    // matches
                    //     .subcommand_matches("Query")
                    //     .unwrap()
                    //     .values_of_t("List")?,
                )?,
                Target::Residues => query_residues(
                    &pdb,
                    parse_residue_list(
                        matches
                            .subcommand_matches("Query")
                            .unwrap()
                            .value_of("List")
                            .unwrap(),
                        &pdb,
                    )?,
                    // matches
                    //     .subcommand_matches("Query")
                    //     .unwrap()
                    //     .values_of_t("List")?,
                )?,
                Target::None => {
                    return Err("Please provide either the 'atoms' or 'residues' flag.".into())
                }
            },
            Source::Sphere => {
                let sphere: Vec<_> = matches
                    .subcommand_matches("Query")
                    .unwrap()
                    .values_of("Sphere")
                    .unwrap()
                    .collect();

                let (origin_id, radius): (usize, f64) = if let [o, r] = sphere[..] {
                    (o.parse()?, r.parse()?)
                } else {
                    return Err("Error parsing sphere values!".into());
                };

                let origin_atom = pdb
                    .atoms()
                    .find(|x| x.serial_number() == origin_id)
                    .ok_or("No atom corresponding to the given ID could be found.")?;

                let list = match target {
                    Target::Atoms => calc_atom_sphere(&pdb, origin_atom, radius, false)?,
                    Target::Residues => calc_residue_sphere(&pdb, origin_atom, radius, false)?,
                    _ => return Err("No target was given!".into()),
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
            let verbosity = match target {
                Target::Atoms => 2,
                Target::Residues => 1,
                _ => 0,
            };
            functions::analyze(&pdb, &region.to_string(), verbosity)?;
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
            let edit_value = match mode {
                Mode::Remove { .. } => 0.00,
                Mode::Add { .. } => match region {
                    Region::Active => 1.00,
                    Region::QM1 => 1.00,
                    Region::QM2 => 2.00,
                    Region::None => 0.00,
                },
                _ => 0.00,
            };
            match source {
                Source::Infile => {
                    println!("This is not implemented yet.")
                }
                Source::List => match target {
                    Target::Atoms => match region {
                        Region::QM1 | Region::QM2 => edit_qm_atoms(
                            &mut pdb,
                            edit_value,
                            matches
                                .subcommand_matches(mode.to_string())
                                .unwrap()
                                .values_of_t("List")?,
                        )?,
                        Region::Active => edit_active_atoms(
                            &mut pdb,
                            edit_value,
                            matches
                                .subcommand_matches(mode.to_string())
                                .unwrap()
                                .values_of_t("List")?,
                        )?,
                        Region::None => {
                            return Err("Please give a region to add atoms or residues to.".into());
                        }
                    },
                    Target::Residues => match region {
                        Region::QM1 | Region::QM2 => edit_qm_residues(
                            &mut pdb,
                            edit_value,
                            matches
                                .subcommand_matches(mode.to_string())
                                .unwrap()
                                .values_of_t("List")?,
                            &partial.to_string(),
                        )?,
                        Region::Active => edit_active_residues(
                            &mut pdb,
                            edit_value,
                            matches
                                .subcommand_matches(mode.to_string())
                                .unwrap()
                                .values_of_t("List")?,
                            &partial.to_string(),
                        )?,
                        Region::None => return Err("Please give a region to modify.".into()),
                    },
                    Target::None => {
                        return Err("Please give either an 'Atoms' or 'Residues' flag.".into())
                    }
                },
                Source::Sphere => {
                    let sphere: Vec<_> = matches
                        .subcommand_matches(mode.to_string())
                        .unwrap()
                        .values_of("Sphere")
                        .unwrap()
                        .collect();

                    let (origin_id, radius): (usize, f64) = if let [o, r] = sphere[..] {
                        (o.parse()?, r.parse()?)
                    } else {
                        return Err("Error parsing sphere values!".into());
                    };

                    let origin_atom = pdb
                        .atoms()
                        .find(|x| x.serial_number() == origin_id)
                        .ok_or("No atom corresponding to the given ID could be found.")?;

                    let list = match target {
                        Target::Atoms => calc_atom_sphere(&pdb, origin_atom, radius, true)?,
                        Target::Residues => calc_residue_sphere(&pdb, origin_atom, radius, true)?,
                        _ => return Err("No target was given!".into()),
                    };

                    match region {
                        Region::QM1 | Region::QM2 => edit_qm_atoms(&mut pdb, edit_value, list)?,
                        Region::Active => edit_active_atoms(&mut pdb, edit_value, list)?,
                        Region::None => return Err("Please give a target region.".into()),
                    }
                }
                Source::None => {
                    if mode.to_string() == "Remove".to_string()
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
            return Err("Please choose either 'Add', 'Analyze', 'Query' or 'Remove' mode.".into())
        }
    }

    if mode.to_string() == "Add".to_string() || mode.to_string() == "Remove".to_string() {
        if matches
            .subcommand_matches(mode.to_string())
            .unwrap()
            .is_present("Outfile")
        {
            pdbtbx::save_pdb(
                pdb,
                matches
                    .subcommand_matches(mode.to_string())
                    .unwrap()
                    .value_of("Outfile")
                    .ok_or("Value for Outfile could not be parsed")?,
                StrictnessLevel::Loose,
            )
            .unwrap();
        } else if matches
            .subcommand_matches(mode.to_string())
            .unwrap()
            .is_present("Overwrite")
        {
            let filename_new = &(filename.to_string() + "_new");
            pdbtbx::save_pdb(pdb, filename_new, StrictnessLevel::Loose).unwrap();
            fs::remove_file(filename)?;
            fs::rename(filename_new, filename)?;
        } else {
            print_to_stdout(&pdb)?;
        }
    }
    Ok(())
}
