use std::error::Error;
// use std::rc::Rc;

use rayon::iter::ParallelIterator;

use crate::functions::*;
use crate::options::*;

// Run function that handles the logic of when to call which function given an enum with the
// command line options. Hands all occurring errors to calöer.
pub fn dispatch(
    mode: Mode,
    mut pdb: &mut pdbtbx::PDB,
    infile: &str,
) -> Result<(), Box<dyn Error>> {
    match &mode {
        Mode::Query { source, target } => match source {
            Source::List(list) => match target {
                Target::Atoms => query_atoms(pdb, parse_atomic_list(list, pdb)?),
                Target::Residues => query_residues(pdb, parse_residue_list(list, pdb)?),
                Target::None => unreachable!(),
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
                    Target::None => unreachable!(),
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
            if distance == &Distance::Clashes || distance == &Distance::Contacts {
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
            Source::List(_) | Source::Infile(_) => {
                let list = match source {
                    Source::List(l) => l.to_string(),
                    Source::Infile(f) => std::fs::read_to_string(f)?,
                    _ => unreachable!(),
                };
                match target {
                    Target::Atoms => {
                        let atomic_list = parse_atomic_list(&list, pdb)?;

                        match region {
                            Region::QM1 | Region::QM2 => {
                                edit_atoms(&mut pdb, atomic_list, *region)?
                            }
                            Region::Active => edit_atoms(&mut pdb, atomic_list, *region)?,
                            Region::None => unreachable!(),
                        }
                    }
                    Target::Residues => {
                        let residue_list = parse_residue_list(&list, pdb)?;

                        match region {
                            Region::QM1 | Region::QM2 => {
                                edit_residues(pdb, residue_list, *partial, *region)?
                            }
                            Region::Active => edit_residues(pdb, residue_list, *partial, *region)?,
                            Region::None => unreachable!(),
                        }
                    }
                    Target::None => unreachable!(),
                }
            }
            Source::Sphere(origin_id, radius) => {
                let sphere_origin = pdb
                    .par_atoms_with_hierarchy()
                    .find_any(|a| a.atom.serial_number() == *origin_id)
                    .ok_or_else::<_, _>(|| {
                        format!("No Atom with serial number {} could be found", origin_id)
                    })?;

                let list = match target {
                    Target::Atoms => calc_atom_sphere(pdb, &sphere_origin, *radius, true)?,
                    Target::Residues => calc_residue_sphere(pdb, &sphere_origin, *radius, true)?,
                    Target::None => unreachable!(),
                };

                match region {
                    Region::QM1 | Region::QM2 => edit_atoms(&mut pdb, list, *region)?,
                    Region::Active => edit_atoms(&mut pdb, list, *region)?,
                    Region::None => unreachable!(),
                }
            }
            Source::None => {
                if mode.to_string() == "Remove"
                    && *region == Region::None
                    && *target == Target::None
                {
                    remove_region(&mut pdb, Region::None)
                } else if { *region == Region::QM1 || *region == Region::QM2 }
                    && target == &Target::None
                {
                    remove_region(&mut pdb, Region::QM1)
                } else if *region == Region::Active && *target == Target::None {
                    remove_region(&mut pdb, Region::Active)
                } else {
                    return Err("Please provide the approprate options (see --help).".into());
                }
            }
        },
        Mode::Write { output, region } => match output {
            Output::None => match region {
                Region::None => print_pdb_to_stdout(pdb),
                _ => println!("{}", get_atomlist(pdb, *region)?),
            },
            Output::Outfile(f) => match region {
                Region::None => print_pdb_to_file(pdb, f)?,
                _ => std::fs::write(f, get_atomlist(pdb, *region)?)?,
            },
            Output::Overwrite => print_pdb_to_file(pdb, infile)?,
        },
    }
    Ok(())
}
