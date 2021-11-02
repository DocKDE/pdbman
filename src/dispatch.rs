use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
// use std::path::Path;
// use std::rc::Rc;

// use itertools::Itertools;
use rayon::iter::ParallelIterator;

use crate::functions::*;
use crate::options::*;

// Run function that handles the logic of when to call which function given an enum with the
// command line options. Hands all occurring errors to caller.
pub fn dispatch(mode: Mode, mut pdb: &mut pdbtbx::PDB, infile: &str) -> Result<(), Box<dyn Error>> {
    match &mode {
        Mode::Query { source, target } => match source {
            Source::List(list) => match target {
                Target::Atoms => query_atoms(pdb, parse_atomic_list(list, pdb)?),
                Target::Residues => query_residues(pdb, parse_residue_list(list, pdb)?)?,
                // Target::None => unreachable!(),
            },
            Source::Sphere(origin_id, radius) => {
                let sphere_origin = pdb
                    .par_atoms_with_hierarchy()
                    .find_any(|a| a.atom.serial_number() == *origin_id)
                    .ok_or_else::<_, _>(|| {
                        format!("No Atom with serial number {} could be found", origin_id)
                    })?;

                let list = match target {
                    Target::Atoms => calc_atom_sphere(pdb, &sphere_origin, *radius, false)?,
                    Target::Residues => calc_residue_sphere(pdb, &sphere_origin, *radius, false)?,
                    // Target::None => unreachable!(),
                };

                query_atoms(pdb, list);
            }
            _ => return Err("Please specify another input for a query.".into()),
        },
        Mode::Analyze {
            region,
            target,
            distance,
        } => {
            analyze(pdb, *region, *target)?;
            if *distance == Some(Distance::Clashes) || *distance == Some(Distance::Contacts) {
                find_contacts(pdb, *distance)?.printstd();
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
        } => match source {
            Some(Source::List(_)) | Some(Source::Infile(_)) => {
                let list = match *source.as_ref().unwrap() {
                    Source::List(l) => l.to_string(),
                    Source::Infile(f) => {
                        let file = BufReader::new(File::open(f)?);
                        let list = file
                            .lines()
                            .map(|l| -> Result<String, anyhow::Error> {
                                let s = l?.trim().to_owned();
                                Ok(s)
                            })
                            .collect::<Result<Vec<String>, anyhow::Error>>()?;
                        list.join(",")
                    }
                    _ => unreachable!(),
                };
                match target.unwrap() {
                    Target::Atoms => {
                        let atomic_list = parse_atomic_list(&list, pdb)?;

                        match region.unwrap() {
                            Region::QM1 | Region::QM2 => {
                                edit_atoms(&mut pdb, atomic_list, *region)?
                            }
                            Region::Active => edit_atoms(&mut pdb, atomic_list, *region)?,
                            // None => unreachable!(),
                        }
                    }
                    Target::Residues => {
                        let residue_list = parse_residue_list(&list, pdb)?;

                        match region.unwrap() {
                            Region::QM1 | Region::QM2 => {
                                edit_residues(pdb, residue_list, *partial, *region)?
                            }
                            Region::Active => edit_residues(pdb, residue_list, *partial, *region)?,
                            // None => unreachable!(),
                        }
                    } // Target::None => unreachable!(),
                }
            }
            Some(Source::Sphere(origin_id, radius)) => {
                let sphere_origin = pdb
                    .par_atoms_with_hierarchy()
                    .find_any(|a| a.atom.serial_number() == *origin_id)
                    .ok_or_else::<_, _>(|| {
                        format!("No Atom with serial number {} could be found", origin_id)
                    })?;

                let list = match target.unwrap() {
                    Target::Atoms => calc_atom_sphere(pdb, &sphere_origin, *radius, true)?,
                    Target::Residues => calc_residue_sphere(pdb, &sphere_origin, *radius, true)?,
                    // Target::None => unreachable!(),
                };

                match region.unwrap() {
                    Region::QM1 | Region::QM2 => edit_atoms(&mut pdb, list, *region)?,
                    Region::Active => edit_atoms(&mut pdb, list, *region)?,
                    // None => unreachable!(),
                }
            }
            None => {
                if mode.to_string() == "Remove" && *region == None && *target == None {
                    remove_region(&mut pdb, None)
                } else if { *region == Some(Region::QM1) || *region == Some(Region::QM2) }
                    && *target == None
                {
                    remove_region(&mut pdb, Some(Region::QM1))
                } else if *region == Some(Region::Active) && *target == None {
                    remove_region(&mut pdb, Some(Region::Active))
                } else {
                    return Err("Please provide the approprate options (see --help).".into());
                }
            }
        },
        Mode::Write { output, region } => match output {
            None => match region {
                None => print_pdb_to_stdout(pdb)?,
                _ => {
                    let stdout = io::stdout();
                    let mut handle = stdout.lock();
                    for num in get_atomlist(pdb, *region)? {
                        writeln!(handle, "{}", num)?
                    }
                }
            },
            Some(Output::Outfile(f)) => match region {
                None => print_pdb_to_file(pdb, f)?,
                _ => {
                    let mut file = BufWriter::new(File::create(f)?);
                    for num in get_atomlist(pdb, *region)? {
                        writeln!(file, "{}", num)?;
                    }
                }
            },
            Some(Output::Overwrite) => print_pdb_to_file(pdb, infile)?,
        },
    }
    Ok(())
}
