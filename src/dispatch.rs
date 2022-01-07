use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

use anyhow::Context;
use colored::Colorize;
use nom::error::convert_error;
use nom::Finish;
// use nom::error::context;
// use nom_supreme::final_parser::final_parser;
use pdbtbx::save_pdb;
// use rayon::iter::ParallelIterator;

use crate::functions::*;
use crate::options::{Distance, Mode, Output, Region, Source, Target};
use crate::revertable::{EditOp, Revertable};
use crate::select::{full_list, Conjunction};

// Run function that handles the logic of when to call which function given an enum with the
// command line options. Hands all occurring errors to caller.
pub fn dispatch(
    mode: Mode,
    pdb: &mut pdbtbx::PDB,
    pdb_path: &str,
) -> Result<Option<Box<dyn Revertable>>, anyhow::Error> {
    let mut edit_op: Option<Box<dyn Revertable>> = None;

    match &mode {
        Mode::Query { input } => {
            // let (sel, conj) = final_parser(full_list)(input)?;
            // let (_, (sele_vec, conj_vec)) = context("test", full_list)(input).map_err(|e| e.to_owned()).finish()?;
            // let (_, (sele_vec, conj_vec)) = full_list(input).map_err(|e| e.to_owned())?;
            // let (_, (sele_vec, conj_vec)) = full_list(input)?;
            let (_, (sele_vec, conj_vec)) = match full_list(input).finish() {
                Ok(r) => r,
                Err(e) => {
                    println!("{:#?}", convert_error(input.as_ref(), e));
                    bail!("Something went wrong during parsing")
                }
            };

            let mut sele_iter = sele_vec.into_iter();
            let initial_sel = sele_iter.next().unwrap();
            let initial_atomvec = get_atoms_from_selection(initial_sel, pdb, None)?;

            match conj_vec {
                // if vec of conjunctions is not present, the vec of selections can only have one element
                None => {
                    query_atoms(pdb, &initial_atomvec)?.printstd();
                }
                Some(cv) => {
                    let mut prev_set: HashSet<usize> = HashSet::from_iter(initial_atomvec);

                    // let mut new_set: HashSet<usize>;

                    for (sele, conj) in sele_iter.zip(cv) {
                        let current_atomvec = get_atoms_from_selection(sele, pdb, None)?;
                        let current_set: HashSet<usize> = HashSet::from_iter(current_atomvec);

                        prev_set = match conj {
                            Conjunction::And => {
                                prev_set.intersection(&current_set).copied().collect()
                            }
                            Conjunction::Or => prev_set.union(&current_set).copied().collect(),
                        };
                    }

                    query_atoms(pdb, &prev_set.into_iter().collect::<Vec<usize>>())?.printstd();
                }
            }
            // println!("{:?}", sele_vec);
            // if let Some(c) = conj_vec {
            //     println!("{:?}", c)
            // }
        }
        // Mode::Query { source, target } => match source {
        //     Source::List(list) => match target {
        //         Target::Atoms => {
        //             query_atoms(pdb, &parse_atomic_list(list, pdb)?)?.printstd();
        //         }
        //         Target::Residues => {
        //             query_atoms(pdb, &parse_residue_list(list, pdb, None)?)?.printstd();
        //         }
        //     },
        //     Source::Sphere(origin_id, radius) => {
        //         let list = match target {
        //             Target::Atoms => get_atom_sphere(pdb, *origin_id, *radius, false)?,
        //             Target::Residues => get_residue_sphere(pdb, *origin_id, *radius, false)?,
        //         };

        //         query_atoms(pdb, &list)?.printstd();
        //     }
        //     _ => bail!("Please don't do this to me..."),
        // },
        Mode::Analyze {
            region,
            target,
            distance,
        } => {
            let (basic_table, detailed_table) = analyze(pdb, *region, *target)?;
            basic_table.printstd();

            if let Some(t) = detailed_table {
                // target must be present if detailed_table is Some
                let target_str = target.unwrap().to_string();
                writeln!(io::stdout(), "\n{} {}", region.unwrap(), target_str).with_context(
                    || format!(" Failed to print {} to stdout", target_str.to_lowercase()),
                )?;
                t.printstd();
            };

            if let Some(d) = *distance {
                let table = find_contacts(pdb, d)?;
                match d {
                    Distance::Clashes => {
                        writeln!(io::stdout(), "\nClash Analysis")
                            .context("Failed to print clash analysis to stdout.")?;
                    }
                    Distance::Contacts => {
                        writeln!(io::stdout(), "\nContact Analysis")
                            .context("Failed to print contact analysis to stdout.")?;
                    }
                }
                table.printstd();
            }
        }
        Mode::Add {
            region,
            source,
            // target,
            partial,
            // input,
        }
        | Mode::Remove {
            region,
            source,
            // target,
            partial,
            // input,
        } => {
            let mut input_list: Vec<usize> = Vec::new();
            match source {
                // If source is Some(_), clap requires target and region arguments as well, hence
                // calling unwrap on them is fine.
                Some(s) => {
                    let input = match s {
                        Source::List(l) => l.to_owned(),
                        Source::Infile(f) => {
                            let file = BufReader::new(
                                File::open(f).context("\nNo such file or directory".red())?,
                            );
                            file.lines()
                                .enumerate()
                                .map(|(i, l)| -> Result<String, anyhow::Error> {
                                    Ok(l.context(format!(
                                        "{}: {}",
                                        "COULDN'T READ LINE FROM FILE".red(),
                                        i.to_string().blue()
                                    ))?
                                    .trim()
                                    .to_owned())
                                    // Ok(s)
                                })
                                .collect::<Result<Vec<String>, anyhow::Error>>()?
                                .join(",")
                        }
                    };

                    let (_, (sele_vec, conj_vec)) = match full_list(&input).finish() {
                        Ok(r) => r,
                        Err(e) => {
                            println!("{:#?}", convert_error(input.as_ref(), e));
                            bail!("\nSomething went wrong during parsing".red())
                        }
                    };

                    let mut sele_iter = sele_vec.into_iter();
                    let initial_sel = sele_iter.next().unwrap();
                    let initial_atomvec = get_atoms_from_selection(initial_sel, pdb, *partial)?;

                    match conj_vec {
                        // if vec of conjunctions is not present, the vec of selections can only have one element
                        None => input_list.extend(initial_atomvec),
                        Some(cv) => {
                            let mut prev_set: HashSet<usize> = HashSet::from_iter(initial_atomvec);

                            // let mut new_set: HashSet<usize>;

                            for (sele, conj) in sele_iter.zip(cv) {
                                let current_atomvec =
                                    get_atoms_from_selection(sele, pdb, *partial)?;
                                let current_set: HashSet<usize> =
                                    HashSet::from_iter(current_atomvec);

                                prev_set = match conj {
                                    Conjunction::And => {
                                        prev_set.intersection(&current_set).copied().collect()
                                    }
                                    Conjunction::Or => {
                                        prev_set.union(&current_set).copied().collect()
                                    }
                                };
                            }

                            input_list.extend(prev_set)
                        }
                    }

                    // Some(Source::List(_)) | Some(Source::Infile(_)) => {
                    //     let input = match source.as_ref().unwrap() {
                    //         Source::List(l) => l.to_string(),
                    //         Source::Infile(f) => {
                    //             let file = BufReader::new(
                    //                 File::open(f).context("\nNO SUCH FILE OR DIRECTORY".red())?,
                    //             );
                    //             file.lines()
                    //                 .enumerate()
                    //                 .map(|(i, l)| -> Result<String, anyhow::Error> {
                    //                     Ok(l.context(format!(
                    //                         "{}: {}",
                    //                         "COULDN'T READ LINE FROM FILE".red(),
                    //                         i.to_string().blue()
                    //                     ))?
                    //                     .trim()
                    //                     .to_owned())
                    //                     // Ok(s)
                    //                 })
                    //                 .collect::<Result<Vec<String>, anyhow::Error>>()?
                    //                 .join(",")
                    //         }
                    //         _ => unreachable!(),
                    //     };

                    // input_list.extend(match target.unwrap() {
                    //     Target::Atoms => parse_atomic_list(&list, pdb)?,
                    //     Target::Residues => {
                    //         let res_list = parse_residue_list(&list, pdb, *partial)?;
                    //         get_atomlist_from_residuelist(&res_list, pdb)
                    //     }
                    // })
                }
                // Some(Source::Sphere(origin_id, radius)) => {
                //     input_list.extend(match target.unwrap() {
                //         Target::Atoms => get_atom_sphere(pdb, *origin_id, *radius, true)?,
                //         Target::Residues => get_residue_sphere(pdb, *origin_id, *radius, true)?,
                //     });
                // }
                None => {
                    if mode.to_string() == "Remove" && *region == None {
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
                    } else if region.is_some() {
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
                edit_op = match mode.to_string().as_str() {
                    // Necessary when adding QM1 atoms over QM2 atoms or vice versa
                    "Add" => match region.unwrap() {
                        r @ Region::QM1 | r @ Region::QM2 => {
                            let other_region = match r {
                                Region::QM1 => Region::QM2,
                                Region::QM2 => Region::QM1,
                                _ => unreachable!(),
                            };

                            // Try removing atoms from other QM region when adding any to see if anything would be
                            // overwritten.
                            let actual_other =
                                edit_atoms_checked(pdb, &input_list, "Remove", other_region);
                            let actual_self = edit_atoms_checked(pdb, &input_list, "Add", r)?;

                            if let Ok(actual) = actual_other {
                                Some(Box::new(vec![
                                    EditOp::ToRemove {
                                        region: other_region,
                                        atoms: actual,
                                    },
                                    EditOp::ToAdd {
                                        region: r,
                                        atoms: actual_self,
                                    },
                                ]))
                            } else {
                                Some(Box::new(EditOp::ToAdd {
                                    region: r,
                                    atoms: actual_self,
                                }))
                            }
                        }
                        Region::Active => Some(Box::new(EditOp::ToAdd {
                            region: Region::Active,
                            atoms: edit_atoms_checked(pdb, &input_list, "Add", Region::Active)?,
                        })),
                    },
                    "Remove" => Some(Box::new(EditOp::ToRemove {
                        region: region.unwrap(),
                        atoms: edit_atoms_checked(pdb, &input_list, "Remove", region.unwrap())?,
                    })),
                    _ => unreachable!(),
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
                if let Err(e) = save_pdb(pdb, pdb_path, pdbtbx::StrictnessLevel::Loose) {
                    e.into_iter().for_each(|e| println!("{}", e))
                }
            }
        },
    }
    Ok(edit_op)
}

#[cfg(test)]
mod tests {
    use crate::options::clap_args;

    use super::*;
    use clap::ArgMatches;
    use itertools::Itertools;
    use lazy_regex::regex;
    use pdbtbx::{StrictnessLevel, PDB};

    fn test_pdb(path: &str) -> PDB {
        let (pdb, _) = pdbtbx::open_pdb(path, StrictnessLevel::Strict).unwrap();
        pdb
    }

    fn get_matches<'a>(args: impl Iterator<Item = &'a str>) -> ArgMatches {
        clap_args().get_matches_from(args)
    }

    fn get_edit_action<'a>(pdb_path: &str, args: impl Iterator<Item = &'a str>) -> String {
        let mut pdb = test_pdb(pdb_path);
        let matches = clap_args().get_matches_from(args);
        let mode = Mode::new(&matches).unwrap();
        format!("{:?}", dispatch(mode, &mut pdb, pdb_path).unwrap().unwrap())
    }

    fn get_editop(text: &str) -> (&str, &str, &str) {
        let re = regex!(r"(\w+)\s\{\sregion:\s(\w+),\satoms:\s\[((\d(, )?)+)\]\s");
        let caps = re.captures(text).unwrap();
        let (edit_action, region, atoms) = (
            caps.get(1).unwrap().as_str(),
            caps.get(2).unwrap().as_str(),
            caps.get(3).unwrap().as_str(),
        );
        (edit_action, region, atoms)
    }

    fn get_atomvec(text: &str) -> Vec<usize> {
        let atom_vec: Vec<usize> = text
            .split(',')
            .map(|a| a.trim().parse::<usize>().unwrap())
            .sorted()
            .collect();
        atom_vec
    }

    #[test]
    fn return_nothing() {
        let pdb_path = "tests/test_blank.pdb";
        let mut pdb = test_pdb(pdb_path);
        let matches = get_matches(["Y"].into_iter());
        let mode = Mode::new(&matches).unwrap();
        assert!(dispatch(mode, &mut pdb, pdb_path).unwrap().is_none());

        // let matches = get_matches(["Q", "-rl", "12"].into_iter());
        // let mode = Mode::new(&matches).unwrap();
        // assert!(dispatch(mode, &mut pdb, pdb_path).unwrap().is_none())
    }

    #[test]
    fn add_qm1() {
        let edit_action = get_edit_action(
            "tests/test_blank.pdb",
            ["A", "-ql", "id", "4,1,9"].into_iter(),
        );
        let (edit_action, region, atoms) = get_editop(&edit_action);
        let atom_vec = get_atomvec(atoms);

        assert_eq!(edit_action, "ToAdd");
        assert_eq!(region, "QM1");
        assert_eq!(atom_vec, vec![1, 4, 9])
    }

    #[test]
    fn add_qm2() {
        let edit_action = get_edit_action(
            "tests/test_blank.pdb",
            ["A", "-ol", "id", "22,17,8"].into_iter(),
        );
        let (edit_action, region, atoms) = get_editop(&edit_action);
        let atom_vec = get_atomvec(atoms);

        assert_eq!(edit_action, "ToAdd");
        assert_eq!(region, "QM2");
        assert_eq!(atom_vec, vec![8, 17, 22])
    }

    #[test]
    fn add_active() {
        let edit_action = get_edit_action(
            "tests/test_blank.pdb",
            ["A", "-al", "id", "22,17,8"].into_iter(),
        );
        let (edit_action, region, atoms) = get_editop(&edit_action);
        let atom_vec = get_atomvec(atoms);

        assert_eq!(edit_action, "ToAdd");
        assert_eq!(region, "Active");
        assert_eq!(atom_vec, vec![8, 17, 22])
    }

    #[test]
    fn remove_qm1() {
        let edit_action = get_edit_action(
            "tests/test_full.pdb",
            ["R", "-ql", "id", "22,17,8"].into_iter(),
        );
        let (edit_action, region, atoms) = get_editop(&edit_action);
        let atom_vec = get_atomvec(atoms);

        assert_eq!(edit_action, "ToRemove");
        assert_eq!(region, "QM1");
        assert_eq!(atom_vec, vec![8, 17, 22])
    }

    #[test]
    fn remove_active() {
        let edit_action = get_edit_action(
            "tests/test_full.pdb",
            ["R", "-al", "id", "22,17,8"].into_iter(),
        );
        let (edit_action, region, atoms) = get_editop(&edit_action);
        let atom_vec = get_atomvec(atoms);

        assert_eq!(edit_action, "ToRemove");
        assert_eq!(region, "Active");
        assert_eq!(atom_vec, vec![8, 17, 22])
    }

    #[test]
    #[should_panic]
    fn remove_empty() {
        get_edit_action(
            "tests/test_blank.pdb",
            ["R", "-al", "id", "22,17,8"].into_iter(),
        );
    }

    #[test]
    #[should_panic]
    fn add_existing() {
        get_edit_action(
            "tests/test_full.pdb",
            ["A", "-ql", "id", "22,17,8"].into_iter(),
        );
    }

    #[test]
    fn checked_add() {
        let edit_action = get_edit_action(
            "tests/test_overwrite.pdb",
            ["A", "-ql", "id", "1,3,4"].into_iter(),
        );
        let (edit_action, region, atoms) = get_editop(&edit_action);
        let atom_vec = get_atomvec(atoms);

        assert_eq!(edit_action, "ToAdd");
        assert_eq!(region, "QM1");
        assert_eq!(atom_vec, vec![3, 4])
    }

    #[test]
    fn checked_remove() {
        let edit_action = get_edit_action(
            "tests/test_overwrite.pdb",
            ["R", "-al", "id", "1-3"].into_iter(),
        );
        let (edit_action, region, atoms) = get_editop(&edit_action);
        let atom_vec = get_atomvec(atoms);

        assert_eq!(edit_action, "ToRemove");
        assert_eq!(region, "Active");
        assert_eq!(atom_vec, vec![3])
    }

    #[test]
    fn overwrite_qm() {
        let edit_action = get_edit_action(
            "tests/test_overwrite.pdb",
            ["A", "-ol", "id", "9,1-3"].into_iter(),
        );
        let re = regex!(r"(\w+)\s\{\sregion:\s(\w+),\satoms:\s\[((\d(, )?)+)\]\s");
        let mut matches = re.find_iter(&edit_action);

        let mut qm1_vec = Vec::new();
        let match1 = matches.next().unwrap();
        let caps1 = re.captures(match1.as_str()).unwrap();

        for i in 1..=3 {
            qm1_vec.push(caps1.get(i).unwrap().as_str());
        }

        let mut qm2_vec = Vec::new();
        let match2 = matches.next().unwrap();
        let caps2 = re.captures(match2.as_str()).unwrap();

        for i in 1..=3 {
            qm2_vec.push(caps2.get(i).unwrap().as_str());
        }

        assert_eq!(qm1_vec[..2], vec!["ToRemove", "QM1"]);
        assert_eq!(qm2_vec[..2], vec!["ToAdd", "QM2"]);

        let qm1_atoms: Vec<usize> = qm1_vec[2]
            .split(',')
            .map(|a| a.trim().parse::<usize>().unwrap())
            .sorted()
            .collect();
        let qm2_atoms: Vec<usize> = qm2_vec[2]
            .split(',')
            .map(|a| a.trim().parse::<usize>().unwrap())
            .sorted()
            .collect();

        assert_eq!(qm1_atoms, vec![1]);
        assert_eq!(qm2_atoms, vec![1, 3, 9]);
    }
}
