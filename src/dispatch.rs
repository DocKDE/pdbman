// use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

use anyhow::Context;
use colored::Colorize;
use pdbtbx::save_pdb;
use rayon::iter::ParallelIterator;

use crate::functions::*;
use crate::options::*;

// Run function that handles the logic of when to call which function given an enum with the
// command line options. Hands all occurring errors to caller.
pub fn dispatch(mode: Mode, mut pdb: &mut pdbtbx::PDB, infile: &str) -> Result<(), anyhow::Error> {
    match &mode {
        Mode::Query { source, target } => match source {
            Source::List(list) => match target {
                Target::Atoms => query_atoms(pdb, parse_atomic_list(list, pdb)?)?,
                Target::Residues => query_residues(pdb, parse_residue_list(list, pdb)?)?,
            },
            Source::Sphere(origin_id, radius) => {
                let sphere_origin = pdb
                    .par_atoms_with_hierarchy()
                    .find_any(|a| a.atom.serial_number() == *origin_id)
                    .ok_or_else::<_, _>(|| {
                        // anyhow!("No Atom with serial number {} could be found", origin_id)
                        anyhow!(
                            "{} {}",
                            "No Atom with serial number found with serial number:".red(),
                            origin_id.to_string().red(),
                        )
                    })?;

                let list = match target {
                    Target::Atoms => calc_atom_sphere(pdb, &sphere_origin, *radius, false)?,
                    Target::Residues => calc_residue_sphere(pdb, &sphere_origin, *radius, false)?,
                };

                query_atoms(pdb, list)?;
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
        } => match source {
            // If source is Some(_), clap requires target and region arguments as well, hence
            // calling unwrap on them is fine.
            Some(Source::List(_)) | Some(Source::Infile(_)) => {
                let list = match *source.as_ref().unwrap() {
                    Source::List(l) => l.to_string(),
                    Source::Infile(f) => {
                        let file =
                            BufReader::new(File::open(f).context("No such file or directory")?);
                        let list = file
                            .lines()
                            .enumerate()
                            .map(|(i, l)| -> Result<String, anyhow::Error> {
                                let s = l
                                    .context(format!("Couldn't read line {} from file", i).red())?
                                    .trim()
                                    .to_owned();
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
                            Region::QM1 | Region::QM2 => edit_atoms(
                                &mut pdb,
                                atomic_list,
                                &mode.to_string(),
                                region.unwrap(),
                            ),
                            Region::Active => edit_atoms(
                                &mut pdb,
                                atomic_list,
                                &mode.to_string(),
                                region.unwrap(),
                            ),
                        }
                    }
                    Target::Residues => {
                        let residue_list = parse_residue_list(&list, pdb)?;

                        match region.unwrap() {
                            Region::QM1 | Region::QM2 => edit_residues(
                                pdb,
                                residue_list,
                                &mode.to_string(),
                                *partial,
                                region.unwrap(),
                            ),
                            Region::Active => edit_residues(
                                pdb,
                                residue_list,
                                &mode.to_string(),
                                *partial,
                                region.unwrap(),
                            ),
                        }
                    }
                }
            }
            Some(Source::Sphere(origin_id, radius)) => {
                let sphere_origin = pdb
                    .par_atoms_with_hierarchy()
                    .find_any(|a| a.atom.serial_number() == *origin_id)
                    .ok_or_else::<_, _>(|| {
                        // anyhow!("No Atom with serial number {} could be found", origin_id)
                        anyhow!(
                            "{} {}",
                            "No Atom with serial number found with serial number:".red(),
                            origin_id.to_string().red(),
                        )
                    })?;

                let list = match target.unwrap() {
                    Target::Atoms => calc_atom_sphere(pdb, &sphere_origin, *radius, true)?,
                    Target::Residues => calc_residue_sphere(pdb, &sphere_origin, *radius, true)?,
                };

                match region.unwrap() {
                    Region::QM1 | Region::QM2 => {
                        edit_atoms(&mut pdb, list, &mode.to_string(), region.unwrap())
                    }
                    Region::Active => {
                        edit_atoms(&mut pdb, list, &mode.to_string(), region.unwrap())
                    }
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
                    bail!("Please provide the approprate options (see --help).".green())
                }
            }
        },
        Mode::Write { output, region } => match output {
            None => match region {
                None => print_pdb_to_stdout(pdb)?,
                _ => {
                    let stdout = io::stdout();
                    let mut handle = stdout.lock();
                    for num in get_atomlist(pdb, region.unwrap())? {
                        writeln!(handle, "{}", num)
                            .context("Failed to write list of atoms to stdout".red())?
                    }
                }
            },
            Some(Output::Outfile(f)) => match region {
                None => {
                    if let Err(e) = save_pdb(pdb, f, pdbtbx::StrictnessLevel::Loose) {
                        e.into_iter().for_each(|e| println!("{}", e))
                    }
                }
                // None => print_pdb_to_file(pdb, f)?,
                _ => {
                    let mut file = BufWriter::new(File::create(f)?);
                    for num in get_atomlist(pdb, region.unwrap())? {
                        writeln!(file, "{}", num)
                            .context("Failed to write list of atoms to stdout".red())?;
                    }
                }
            },
            Some(Output::Overwrite) => {
                if let Err(e) = save_pdb(pdb, infile, pdbtbx::StrictnessLevel::Loose) {
                    e.into_iter().for_each(|e| println!("{}", e))
                }
            } // Some(Output::Overwrite) => print_pdb_to_file(pdb, infile)?,
        },
    }
    Ok(())
}
