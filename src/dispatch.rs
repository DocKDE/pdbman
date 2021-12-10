use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

use anyhow::Context;
use colored::Colorize;
use pdbtbx::{save_pdb, ContainsAtomConformer, ContainsAtomConformerResidue};
use rayon::iter::ParallelIterator;

use crate::functions::*;
use crate::options::{Mode, Output, Partial, Region, Source, Target};
use crate::revertable::{EditOp, Revertable};

// Run function that handles the logic of when to call which function given an enum with the
// command line options. Hands all occurring errors to caller.
pub fn dispatch(
    mode: Mode,
    pdb: &mut pdbtbx::PDB,
    infile: &str,
) -> Result<Option<Box<dyn Revertable>>, anyhow::Error> {
    let mut edit_op: Option<Box<dyn Revertable>> = None;

    match &mode {
        Mode::Query { source, target } => match source {
            Source::List(list) => match target {
                Target::Atoms => query_atoms(pdb, &parse_atomic_list(list, pdb)?)?,
                Target::Residues => query_residues(pdb, &parse_residue_list(list, pdb)?)?,
            },
            Source::Sphere(origin_id, radius) => {
                let sphere_origin = pdb
                    .atoms_with_hierarchy()
                    .find(|a| a.atom().serial_number() == *origin_id)
                    .ok_or_else::<_, _>(|| {
                        anyhow!(
                            "{}: '{}'",
                            "\nNO ATOM WITH FOUND WITH SERIAL NUMBER".red(),
                            origin_id.to_string().blue(),
                        )
                    })?;

                let list = match target {
                    Target::Atoms => calc_atom_sphere(pdb, sphere_origin.atom(), *radius, false)?,
                    Target::Residues => calc_residue_sphere(pdb, sphere_origin, *radius, false)?,
                };

                query_atoms(pdb, &list)?;
            }
            _ => bail!("Please don't do this to me..."),
        },
        Mode::Analyze {
            region,
            target,
            distance,
        } => {
            analyze(pdb, *region, *target)?;
            if let Some(d) = *distance {
                find_contacts(pdb, d)?.printstd();
            }
        }
        Mode::Add {
            region,
            source,
            target,
            partial,
        }
        | Mode::Remove {
            region,
            source,
            target,
            partial,
        } => {
            let mut input_list: Vec<usize> = Vec::new();
            match source {
                // If source is Some(_), clap requires target and region arguments as well, hence
                // calling unwrap on them is fine.
                Some(Source::List(_)) | Some(Source::Infile(_)) => {
                    let list = match *source.as_ref().unwrap() {
                        Source::List(l) => l.to_string(),
                        Source::Infile(f) => {
                            let file = BufReader::new(
                                File::open(f).context("\nNO SUCH FILE OR DIRECTORY".red())?,
                            );
                            file.lines()
                                .enumerate()
                                .map(|(i, l)| -> Result<String, anyhow::Error> {
                                    let s = l
                                        .context(format!(
                                            "{}: {}",
                                            "COULDN'T READ LINE FROM FILE".red(),
                                            i.to_string().blue()
                                        ))?
                                        .trim()
                                        .to_owned();
                                    Ok(s)
                                })
                                .collect::<Result<Vec<String>, anyhow::Error>>()?
                                .join(",")
                        }
                        _ => unreachable!(),
                    };

                    input_list.extend(match target.unwrap() {
                        Target::Atoms => parse_atomic_list(&list, pdb)?,
                        Target::Residues => {
                            let residue_list = parse_residue_list(&list, pdb)?;
                            match partial {
                                None => pdb
                                    .atoms_with_hierarchy()
                                    .filter(|a| residue_list.contains(&a.residue().serial_number()))
                                    .map(|a| a.atom().serial_number())
                                    .collect(),
                                Some(p) => pdb
                                    .atoms_with_hierarchy()
                                    .filter(|a| match p {
                                        Partial::Backbone => a.is_backbone(),
                                        Partial::Sidechain => a.is_sidechain(),
                                    } && residue_list.contains(&a.residue().serial_number()))
                                    // .filter(|a| match p {
                                    //     Partial::Backbone => a.is_backbone(),
                                    //     Partial::Sidechain => a.is_sidechain(),
                                    // })
                                    .map(|a| a.atom().serial_number())
                                    .collect(),
                            }
                        }
                    })
                }
                // Some(Source::Infile(f)) => {
                //     let file =
                //         BufReader::new(File::open(f).context("\nNO SUCH FILE OR DIRECTORY".red())?);
                //     let list = file
                //         .lines()
                //         .enumerate()
                //         .map(|(i, l)| -> Result<String, anyhow::Error> {
                //             let s = l
                //                 .context(format!(
                //                     "{}: {}",
                //                     "COULDN'T READ LINE FROM FILE".red(),
                //                     i.to_string().blue()
                //                 ))?
                //                 .trim()
                //                 .to_owned();
                //             Ok(s)
                //         })
                //         .collect::<Result<Vec<String>, anyhow::Error>>()?
                //         .join(",");
                //     input_list.extend(match target.unwrap() {
                //         Target::Atoms => parse_atomic_list(&list, pdb)?,
                //         Target::Residues => {
                //             let residue_list = parse_residue_list(&list, pdb)?;
                //             match partial {
                //                 None => pdb
                //                     .atoms_with_hierarchy()
                //                     .filter(|a| residue_list.contains(&a.residue().serial_number()))
                //                     .map(|a| a.atom().serial_number())
                //                     .collect(),
                //                 Some(p) => pdb
                //                     .atoms_with_hierarchy()
                //                     .filter(|a| residue_list.contains(&a.residue().serial_number()))
                //                     .filter(|a| match p {
                //                         Partial::Backbone => a.is_backbone(),
                //                         Partial::Sidechain => a.is_sidechain(),
                //                     })
                //                     .map(|a| a.atom().serial_number())
                //                     .collect(),
                //             }
                //         }
                //     })
                // }
                Some(Source::Sphere(origin_id, radius)) => {
                    let sphere_origin = pdb
                        .atoms_with_hierarchy()
                        .find(|a| a.atom().serial_number() == *origin_id)
                        .ok_or_else::<_, _>(|| {
                            anyhow!(
                                "{}: '{}'",
                                "NO ATOM FOUND WITH SERIAL NUMBER".red(),
                                origin_id.to_string().blue(),
                            )
                        })?;

                    input_list.extend(match target.unwrap() {
                        Target::Atoms => {
                            calc_atom_sphere(pdb, sphere_origin.atom(), *radius, true)?
                        }
                        Target::Residues => calc_residue_sphere(pdb, sphere_origin, *radius, true)?,
                    });
                }
                // Some(Source::List(_)) | Some(Source::Infile(_)) => {
                //     let list = match *source.as_ref().unwrap() {
                //         Source::List(l) => l.to_string(),
                //         Source::Infile(f) => {
                //             let file = BufReader::new(
                //                 File::open(f).context("\nNO SUCH FILE OR DIRECTORY".red())?,
                //             );
                //             let list = file
                //                 .lines()
                //                 .enumerate()
                //                 .map(|(i, l)| -> Result<String, anyhow::Error> {
                //                     let s = l
                //                         .context(format!(
                //                             "{}: {}",
                //                             "COULDN'T READ LINE FROM FILE".red(),
                //                             i.to_string().blue()
                //                         ))?
                //                         .trim()
                //                         .to_owned();
                //                     Ok(s)
                //                 })
                //                 .collect::<Result<Vec<String>, anyhow::Error>>()?;
                //             list.join(",")
                //         }
                //         _ => unreachable!(),
                //     };
                //     match target.unwrap() {
                //         Target::Atoms => {
                //             // let atomic_list = parse_atomic_list(&list, pdb)?;

                //             // Determine whether already existing are overwritten with the add/remove command
                //             // and factor those out. Necessary to be able to undo the action.
                //             let atom_set: HashSet<usize> =
                //                 HashSet::from_iter(parse_atomic_list(&list, pdb)?.iter().copied());
                //             let set_of_existing: HashSet<usize> = pdb
                //                 .par_atoms()
                //                 .filter(|a| match region.unwrap() {
                //                     Region::QM1 => a.occupancy() == 1.00,
                //                     Region::QM2 => a.occupancy() == 2.00,
                //                     Region::Active => a.b_factor() == 1.00,
                //                 })
                //                 .map(|a| a.serial_number())
                //                 .collect();
                //             let actual_op_list: Vec<usize> = match mode.to_string().as_str() {
                //                 "Add" => atom_set.difference(&set_of_existing).copied().collect(),
                //                 "Remove" => atom_set.intersection(&set_of_existing).copied().collect(),
                //                 _ => unreachable!(),
                //             };

                //             edit_op = Some(Box::new(match mode.to_string().as_str() {
                //                 "Add" => EditOp::ToAdd {
                //                     region: region.unwrap(),
                //                     atoms: actual_op_list.iter().copied().collect(),
                //                 },
                //                 "Remove" => EditOp::ToRemove {
                //                     region: region.unwrap(),
                //                     atoms: actual_op_list.iter().copied().collect(),
                //                 },
                //                 _ => unreachable!(),
                //             }));

                //             match region.unwrap() {
                //                 Region::QM1 | Region::QM2 => {
                //                     edit_atoms(pdb, &actual_op_list, &mode.to_string(), region.unwrap())
                //                 }
                //                 Region::Active => {
                //                     edit_atoms(pdb, &actual_op_list, &mode.to_string(), region.unwrap())
                //                 }
                //             }
                //         }
                //         Target::Residues => {
                //             let residue_list = parse_residue_list(&list, pdb)?;

                //             let atom_list: Vec<usize> = match partial {
                //                 None => pdb
                //                     .atoms_with_hierarchy()
                //                     .filter(|a| residue_list.contains(&a.residue().serial_number()))
                //                     .map(|a| a.atom().serial_number())
                //                     .collect(),
                //                 Some(p) => pdb
                //                     .atoms_with_hierarchy()
                //                     .filter(|a| residue_list.contains(&a.residue().serial_number()))
                //                     .filter(|a| match p {
                //                         Partial::Backbone => a.is_backbone(),
                //                         Partial::Sidechain => a.is_sidechain(),
                //                     })
                //                     .map(|a| a.atom().serial_number())
                //                     .collect(),
                //             };

                //             let atom_set: HashSet<usize> =
                //                 HashSet::from_iter(atom_list.iter().copied());
                //             let set_of_existing: HashSet<usize> = pdb
                //                 .par_atoms()
                //                 .filter(|a| match region.unwrap() {
                //                     Region::QM1 => a.occupancy() == 1.00,
                //                     Region::QM2 => a.occupancy() == 2.00,
                //                     Region::Active => a.b_factor() == 1.00,
                //                 })
                //                 .map(|a| a.serial_number())
                //                 .collect();
                //             let actual_op_list: Vec<usize> = match mode.to_string().as_str() {
                //                 "Add" => atom_set.difference(&set_of_existing).copied().collect(),
                //                 "Remove" => atom_set.intersection(&set_of_existing).copied().collect(),
                //                 _ => unreachable!(),
                //             };

                //             edit_op = Some(Box::new(match mode.to_string().as_str() {
                //                 "Add" => EditOp::ToAdd {
                //                     region: region.unwrap(),
                //                     atoms: actual_op_list.clone(),
                //                 },
                //                 "Remove" => EditOp::ToRemove {
                //                     region: region.unwrap(),
                //                     atoms: actual_op_list.clone(),
                //                 },
                //                 _ => unreachable!(),
                //             }));

                //             match region.unwrap() {
                //                 Region::QM1 | Region::QM2 => {
                //                     edit_atoms(pdb, &actual_op_list, &mode.to_string(), region.unwrap())
                //                 }
                //                 Region::Active => {
                //                     edit_atoms(pdb, &actual_op_list, &mode.to_string(), region.unwrap())
                //                 }
                //             }
                //         }
                //     }
                // }
                // Some(Source::Sphere(origin_id, radius)) => {
                //     let sphere_origin = pdb
                //         .atoms_with_hierarchy()
                //         .find(|a| a.atom().serial_number() == *origin_id)
                //         .ok_or_else::<_, _>(|| {
                //             // anyhow!("No Atom with serial number {} could be found", origin_id)
                //             anyhow!(
                //                 "{}: '{}'",
                //                 "NO ATOM FOUND WITH SERIAL NUMBER".red(),
                //                 origin_id.to_string().blue(),
                //             )
                //         })?;

                //     let list = match target.unwrap() {
                //         Target::Atoms => calc_atom_sphere(pdb, sphere_origin.atom(), *radius, true)?,
                //         Target::Residues => calc_residue_sphere(pdb, sphere_origin, *radius, true)?,
                //     };

                //     let atom_set: HashSet<usize> = HashSet::from_iter(list.iter().copied());
                //     let set_of_existing: HashSet<usize> = pdb
                //         .par_atoms()
                //         .filter(|a| match region.unwrap() {
                //             Region::QM1 => a.occupancy() == 1.00,
                //             Region::QM2 => a.occupancy() == 2.00,
                //             Region::Active => a.b_factor() == 1.00,
                //         })
                //         .map(|a| a.serial_number())
                //         .collect();

                //     let actual_op_list: Vec<usize> = match mode.to_string().as_str() {
                //         "Add" => atom_set.difference(&set_of_existing).copied().collect(),
                //         "Remove" => atom_set.intersection(&set_of_existing).copied().collect(),
                //         _ => unreachable!(),
                //     };

                //     edit_op = Some(Box::new(match mode.to_string().as_str() {
                //         "Add" => EditOp::ToAdd {
                //             region: region.unwrap(),
                //             atoms: actual_op_list.iter().copied().collect(),
                //         },
                //         "Remove" => EditOp::ToRemove {
                //             region: region.unwrap(),
                //             atoms: actual_op_list.iter().copied().collect(),
                //         },
                //         _ => unreachable!(),
                //     }));

                //     match region.unwrap() {
                //         Region::QM1 | Region::QM2 => {
                //             edit_atoms(pdb, &actual_op_list, &mode.to_string(), region.unwrap())
                //         }
                //         Region::Active => {
                //             edit_atoms(pdb, &actual_op_list, &mode.to_string(), region.unwrap())
                //         }
                //     }
                // }
                None => {
                    if mode.to_string() == "Remove" && *region == None && *target == None {
                        let mut qm1_atoms = Vec::new();
                        let mut qm2_atoms = Vec::new();
                        let mut active_atoms = Vec::new();

                        for atom in pdb.atoms() {
                            if atom.occupancy() == 1.00 {
                                qm1_atoms.push(atom.serial_number())
                            }
                            if atom.occupancy() == 2.00 {
                                qm2_atoms.push(atom.serial_number())
                            }
                            if atom.b_factor() == 1.00 {
                                active_atoms.push(atom.serial_number())
                            }
                        }

                        let mut remove_ops = Vec::with_capacity(3);
                        if !qm1_atoms.is_empty() {
                            remove_ops.push(EditOp::ToRemove {
                                region: Region::QM1,
                                atoms: qm1_atoms,
                            })
                        }
                        if !qm2_atoms.is_empty() {
                            remove_ops.push(EditOp::ToRemove {
                                region: Region::QM2,
                                atoms: qm2_atoms,
                            })
                        }
                        if !active_atoms.is_empty() {
                            remove_ops.push(EditOp::ToRemove {
                                region: Region::Active,
                                atoms: active_atoms,
                            })
                        }

                        if !remove_ops.is_empty() {
                            edit_op = Some(Box::new(remove_ops));
                        }

                        remove_region(pdb, None);
                    } else if region.is_some() && *target == None {
                        let region_atoms: Vec<usize> = pdb
                            .atoms()
                            .filter(|a| match region.unwrap() {
                                Region::QM1 => a.occupancy() == 1.00,
                                Region::QM2 => a.occupancy() == 2.00,
                                Region::Active => a.b_factor() == 1.00,
                            })
                            .map(|a| a.serial_number())
                            .collect();

                        if !region_atoms.is_empty() {
                            edit_op = Some(Box::new(EditOp::ToRemove {
                                region: region.unwrap(),
                                atoms: region_atoms,
                            }))
                        }

                        remove_region(pdb, Some(region.unwrap()))
                    } else {
                        bail!("Please provide the approprate options (see --help).".red())
                    }
                }
            }
            if source.is_some() {
                let atom_set: HashSet<usize> = HashSet::from_iter(input_list.iter().copied());
                let set_of_existing: HashSet<usize> = pdb
                    .par_atoms()
                    .filter(|a| match region.unwrap() {
                        Region::QM1 => a.occupancy() == 1.00,
                        Region::QM2 => a.occupancy() == 2.00,
                        Region::Active => a.b_factor() == 1.00,
                    })
                    .map(|a| a.serial_number())
                    .collect();

                let actual_op_list: Vec<usize> = match mode.to_string().as_str() {
                    "Add" => atom_set.difference(&set_of_existing).copied().collect(),
                    "Remove" => atom_set.intersection(&set_of_existing).copied().collect(),
                    _ => unreachable!(),
                };

                if !actual_op_list.is_empty() {
                    edit_op = Some(Box::new(match mode.to_string().as_str() {
                        "Add" => EditOp::ToAdd {
                            region: region.unwrap(),
                            atoms: actual_op_list.iter().copied().collect(),
                        },
                        "Remove" => EditOp::ToRemove {
                            region: region.unwrap(),
                            atoms: actual_op_list.iter().copied().collect(),
                        },
                        _ => unreachable!(),
                    }));

                    match region.unwrap() {
                        Region::QM1 | Region::QM2 => {
                            edit_atoms(pdb, &actual_op_list, &mode.to_string(), region.unwrap())
                        }
                        Region::Active => {
                            edit_atoms(pdb, &actual_op_list, &mode.to_string(), region.unwrap())
                        }
                    }
                }
            }
        }
        Mode::Write {
            output,
            region,
            target,
        } => match output {
            None => match region {
                None => print_pdb_to_stdout(pdb)?,
                _ => {
                    let stdout = io::stdout();
                    let mut handle = stdout.lock();
                    // target can be unwrapped because required to be Some by clap if region is Some
                    match target.unwrap() {
                        Target::Atoms => {
                            for num in get_atomlist(pdb, region.unwrap())? {
                                writeln!(handle, "{}", num)
                                    .context("FAILED TO WRITE LIST OF ATOMS TO STDOUT".red())?
                            }
                        }
                        Target::Residues => {
                            for num in get_residuelist(pdb, region.unwrap())? {
                                writeln!(handle, "{}", num)
                                    .context("FAILED TO WRITE LIST OF RESIDUES TO STDOUT".red())?
                            }
                        }
                    }
                }
            },
            Some(Output::Outfile(f)) => match region {
                None => {
                    if let Err(e) = save_pdb(pdb, f, pdbtbx::StrictnessLevel::Loose) {
                        e.into_iter().for_each(|e| println!("{}", e))
                    }
                }
                _ => {
                    let mut file = BufWriter::new(File::create(f)?);
                    // target can be unwrapped because required to be Some by clap if region is Some
                    match target.unwrap() {
                        Target::Atoms => {
                            for num in get_atomlist(pdb, region.unwrap())? {
                                writeln!(file, "{}", num)
                                    .context("FAILED TO WRITE LIST OF ATOMS TO FILE".red())?;
                            }
                        }
                        Target::Residues => {
                            for num in get_residuelist(pdb, region.unwrap())? {
                                writeln!(file, "{}", num)
                                    .context("FAILED TO WRITE LIST OF RESIDUES TO FILE".red())?;
                            }
                        }
                    };
                }
            },
            Some(Output::Overwrite) => {
                if let Err(e) = save_pdb(pdb, infile, pdbtbx::StrictnessLevel::Loose) {
                    e.into_iter().for_each(|e| println!("{}", e))
                }
            }
        },
    }
    Ok(edit_op)
}
