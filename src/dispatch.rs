use std::error::Error;
use std::rc::Rc;

use crate::functions::*;
use crate::options::*;

// Run function that handles the logic of when to call which function given an enum with the
// command line options. Hands all occurring errors to calöer.
pub fn dispatch(
    // matches: clap::ArgMatches,
    mode: Rc<Mode>,
    // mode: Mode,
    mut pdb: &mut pdbtbx::PDB,
) -> Result<(), Box<dyn Error>> {
    match &*mode {
        Mode::Query { source, target } => match source {
            // Source::List => {
            Source::List(list) => {
                // let list = matches
                //     .subcommand_matches("Query").unwrap()
                //     .value_of("List").unwrap();
                match target {
                    Target::Atoms => query_atoms(pdb, parse_atomic_list(list, pdb)?)?,
                    Target::Residues => query_residues(pdb, parse_residue_list(list, pdb)?)?,
                    Target::None => unreachable!(),
                }
            }
            Source::Sphere(origin_str, radius_str) => {
                let sphere = Sphere::new(origin_str, radius_str, pdb)?;
                // matches
                //     .subcommand_matches("Query")
                //     .unwrap()
                //     .values_of("Sphere")
                //     .unwrap(),
                //     pdb,
                // )?;

                let list = match target {
                    Target::Atoms => calc_atom_sphere(pdb, &sphere.origin, sphere.radius, false)?,
                    Target::Residues => {
                        calc_residue_sphere(pdb, &sphere.origin, sphere.radius, false)?
                    }
                    Target::None => unreachable!(),
                };

                query_atoms(pdb, list)?;
            }
            _ => return Err("Please specify another input for a query.".into()),
        },
        Mode::Analyze {
            region,
            target,
            distance,
        } => {
            analyze(pdb, *region, *target)?;
            if distance == &Distance::Clashes || distance == &Distance::Contacts {
                find_contacts(pdb, *distance)?.printstd();
            }
        }
        Mode::Add {
            region,
            source,
            target,
            partial,
            // output: _,
        }
        | Mode::Remove {
            region,
            source,
            target,
            partial,
            // output: _,
        } => {
            let edit_value = match *mode {
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
                    todo!()
                }
                // Source::List => {
                Source::List(list) => {
                    // let list = matches
                    //     .subcommand_matches(mode.to_string()).unwrap()
                    //     .value_of("List").unwrap();

                    match target {
                        Target::Atoms => {
                            let atomic_list = parse_atomic_list(list, pdb)?;

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
                            // let pdb = pdb.clone();
                            // let ref_pdb = std::cell::RefCell::new(pdb);
                            // let pdb_borrow = ref_pdb.borrow();
                            // let residue_list = parse_residue_list(list, &pdb_borrow)?;
                            // let mut pdb_mut = ref_pdb.borrow_mut();
                            let residue_list = parse_residue_list(list, &pdb_clone)?;

                            match region {
                                Region::QM1 | Region::QM2 => {
                                    // edit_qm_residues(&mut pdb_mut, edit_value, residue_list, *partial)?
                                    edit_qm_residues(pdb, edit_value, residue_list, *partial)?
                                }
                                Region::Active => edit_active_residues(
                                    // &mut pdb_mut,
                                    pdb,
                                    edit_value,
                                    residue_list,
                                    *partial,
                                )?,
                                Region::None => unreachable!(),
                            }
                        }
                        Target::None => unreachable!(),
                    }
                }
                Source::Sphere(origin_str, radius_str) => {
                    let sphere = Sphere::new(origin_str, radius_str, pdb)?;
                    //     matches
                    //         .subcommand_matches(mode.to_string())
                    //         .unwrap()
                    //         .values_of("Sphere")
                    //         .unwrap(),
                    //     pdb,
                    // )?;

                    let list = match target {
                        Target::Atoms => {
                            calc_atom_sphere(pdb, &sphere.origin, sphere.radius, true)?
                        }
                        Target::Residues => {
                            calc_residue_sphere(pdb, &sphere.origin, sphere.radius, true)?
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
                        && *region == Region::None
                        && *target == Target::None
                    {
                        remove_all(&mut pdb)?
                    } else if { *region == Region::QM1 || *region == Region::QM2 }
                        && target == &Target::None
                    {
                        remove_qm(&mut pdb)?
                    } else if *region == Region::Active && *target == Target::None {
                        remove_active(&mut pdb)?
                    } else {
                        return Err("Please provide the approprate options (see --help).".into());
                    }
                }
            }
        } // Mode::None => {
          //     unreachable!()
          // }
    }

    Ok(())
}
