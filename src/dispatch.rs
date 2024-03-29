use std::fs::File;
use std::io::{self, BufWriter, Write};

use anyhow::Context;
use colored::Colorize;
use comfy_table::modifiers::{UTF8_ROUND_CORNERS, UTF8_SOLID_INNER_BORDERS};
use comfy_table::presets::UTF8_FULL;
use comfy_table::{Row, Table};
use itertools::Itertools;
use pdbtbx::{save_pdb, Atom, ContainsAtomConformer, ContainsAtomConformerResidue};
use rayon::iter::ParallelIterator;

use crate::functions;
use crate::options::{Distance, MeasureTarget, Mode, Output, Region};
use crate::revertable::{EditOp, Revertable};

// Run function that handles the logic of when to call which function given an enum with the
// command line options. Hands all occurring errors to caller.
pub fn dispatch(
    mode: &Mode,
    pdb: &mut pdbtbx::PDB,
    pdb_path: &str,
) -> Result<Option<Revertable>, anyhow::Error> {
    let mut edit_op: Option<Revertable> = None;

    match mode {
        Mode::Query { input } => {
            let atomlist = functions::get_atomlist_from_input(input, pdb, None)?;
            let (table, res) = functions::query_atoms(pdb, &atomlist)?;
            if let Some(s) = res {
                writeln!(io::stdout(), "{}", s)
                    .context("Failed to print residue depiction to stdout")?
            }
            writeln!(io::stdout(), "{}", table).context("Failed to write table to stdout")?;
        }
        Mode::Analyze {
            region,
            target,
            distance,
        } => {
            let (basic_table, detailed_table) = functions::analyze(pdb, *region, *target)?;
            writeln!(io::stdout(), "{}", basic_table).context("Failed to write table to stdout")?;

            if let Some(t) = detailed_table {
                // target must be present if detailed_table is Some
                let target_str = target.unwrap().to_string();
                writeln!(io::stdout(), "\n{} {}", region.unwrap(), target_str).with_context(
                    || format!(" Failed to print {} to stdout", target_str.to_lowercase()),
                )?;
                writeln!(io::stdout(), "{}", t).context("Failed to write table to stdout")?;
            };

            if let Some(d) = *distance {
                let table = functions::find_contacts(pdb, d)?;
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
                writeln!(io::stdout(), "{}", table).context("Failed to print table to stdout")?;
            }
        }
        Mode::Measure { measure_target } => match measure_target {
            MeasureTarget::Atoms(atoms) => {
                let (table, geom) = functions::get_measurements(atoms, pdb)?;
                writeln!(io::stdout(), "{}\n{}", table, geom)
                    .context("Failed to print table to stdout")?;
            }
            MeasureTarget::Sphere(origin_id, radius) => {
                let origin_atom = pdb
                    .par_atoms()
                    .find_first(|a| a.serial_number() == *origin_id)
                    .ok_or_else(|| {
                        anyhow!(
                            "{}: '{}'",
                            "\nNO ATOM WITH FOUND WITH SERIAL NUMBER".red(),
                            origin_id.to_string().blue(),
                        )
                    })?;

                let tree = pdb.create_hierarchy_rtree();
                let sphere_iter = tree.nearest_neighbor_iter_with_distance_2(&origin_atom.pos());
                let mut table = Table::new();
                table
                    .load_preset(UTF8_FULL)
                    .apply_modifier(UTF8_ROUND_CORNERS)
                    .apply_modifier(UTF8_SOLID_INNER_BORDERS);
                table.set_header(Row::from(vec![
                    "Atom ID",
                    "Atom name",
                    "Residue ID",
                    "Residue Name",
                    "QM",
                    "Active",
                    "Distance",
                ]));

                for (atom_hier, mut dist) in sphere_iter {
                    dist = dist.sqrt();
                    if dist <= *radius {
                        table.add_row(Row::from(vec![
                            atom_hier.atom().serial_number().to_string(),
                            atom_hier.atom().name().to_string(),
                            atom_hier.residue().serial_number().to_string()
                                + atom_hier.residue().insertion_code().unwrap_or(""),
                            atom_hier.residue().name().unwrap_or("").to_owned(),
                            atom_hier.atom().occupancy().to_string(),
                            atom_hier.atom().b_factor().to_string(),
                            format!("{:.3}", dist),
                        ]));
                    } else {
                        break;
                    };
                }

                ensure!(
                    table.row_iter().peekable().peek().is_some(),
                    "No atoms within the given radius"
                );
                writeln!(
                    io::stdout(),
                    "\nAtoms within {:.3} of atom with ID {}\n{}",
                    radius,
                    origin_atom.serial_number(),
                    table
                )
                .context("Failed to print table to stdout")?;
            }
        },
        Mode::Add {
            region,
            selection,
            partial,
        }
        | Mode::Remove {
            region,
            selection,
            partial,
        } => {
            let mut input_list: Vec<usize> = Vec::new();
            match selection {
                Some(s) => {
                    input_list.extend(functions::get_atomlist_from_input(s, pdb, *partial)?);
                }
                None => {
                    if mode.to_string() == "Remove" && *region == None {
                        let mut qm1_atoms = Vec::new();
                        let mut qm2_atoms = Vec::new();
                        let mut active_atoms = Vec::new();

                        for atom in pdb.atoms() {
                            if atom.occupancy() == 1.00 {
                                qm1_atoms.push(atom.serial_number())
                            } else if atom.occupancy() == 2.00 {
                                qm2_atoms.push(atom.serial_number())
                            }
                            if atom.b_factor() == 1.00 {
                                active_atoms.push(atom.serial_number())
                            }
                        }

                        let mut remove_ops = Vec::with_capacity(3);

                        for (list, region) in [qm1_atoms, qm2_atoms, active_atoms]
                            .into_iter()
                            .zip([Region::QM1, Region::QM2, Region::Active].into_iter())
                        {
                            if !list.is_empty() {
                                remove_ops.push(EditOp::ToRemove {
                                    region,
                                    atoms: list,
                                })
                            }
                        }

                        if !remove_ops.is_empty() {
                            edit_op = Some(Revertable::Many(remove_ops));
                        }

                        functions::remove_region(pdb, None);
                    } else if let Some(r) = region {
                        let region_atoms: Vec<usize> = pdb
                            .atoms()
                            .filter(|a| match r {
                                Region::QM1 => a.occupancy() == 1.00,
                                Region::QM2 => a.occupancy() == 2.00,
                                Region::Active => a.b_factor() == 1.00,
                            })
                            .map(Atom::serial_number)
                            .collect();

                        if !region_atoms.is_empty() {
                            edit_op = Some(Revertable::One(EditOp::ToRemove {
                                region: region.unwrap(),
                                atoms: region_atoms,
                            }))
                        }

                        functions::remove_region(pdb, Some(region.unwrap()))
                    } else {
                        bail!("Please provide the approprate options (see --help).".red())
                    }
                }
            }
            if selection.is_some() {
                edit_op = match mode.to_string().as_str() {
                    // Necessary when adding QM1 atoms over QM2 atoms or vice versa
                    "Add" => match region.unwrap() {
                        r @ (Region::QM1 | Region::QM2) => {
                            let other_region = match r {
                                Region::QM1 => Region::QM2,
                                Region::QM2 => Region::QM1,
                                Region::Active => unreachable!(),
                            };

                            // Try removing atoms from other QM region when adding any to see if anything would be
                            // overwritten.
                            let actual_other = functions::edit_atoms_checked(
                                pdb,
                                &input_list,
                                "Remove",
                                other_region,
                            );
                            let actual_self =
                                functions::edit_atoms_checked(pdb, &input_list, "Add", r)?;

                            if let Ok(actual) = actual_other {
                                Some(Revertable::Many(vec![
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
                                Some(Revertable::One(EditOp::ToAdd {
                                    region: r,
                                    atoms: actual_self,
                                }))
                            }
                        }
                        Region::Active => Some(Revertable::One(EditOp::ToAdd {
                            region: Region::Active,
                            atoms: functions::edit_atoms_checked(
                                pdb,
                                &input_list,
                                "Add",
                                Region::Active,
                            )?,
                        })),
                    },
                    "Remove" => Some(Revertable::One(EditOp::ToRemove {
                        region: region.unwrap(),
                        atoms: functions::edit_atoms_checked(
                            pdb,
                            &input_list,
                            "Remove",
                            region.unwrap(),
                        )?,
                    })),
                    _ => unreachable!(),
                }
            }
        }
        Mode::Write { output, state } => match output {
            None => {
                if *state {
                    let stdout = io::stdout();
                    let mut handle = stdout.lock();

                    writeln!(handle, "R")?;
                    for (region, string) in [Region::QM1, Region::QM2, Region::Active]
                        .into_iter()
                        .zip(["-q", "-o", "-a"].into_iter())
                    {
                        if let Ok(l) = functions::get_atomlist(pdb, region) {
                            writeln!(
                                handle,
                                "A {} id {}",
                                string,
                                l.into_iter().map(|n| n.to_string()).join(",")
                            )?;
                        };
                    }
                    writeln!(handle, "W -w")?;
                } else {
                    functions::print_pdb_to_stdout(pdb)?;
                }
            }
            Some(Output::Outfile(f)) => {
                if *state {
                    let mut file = BufWriter::new(File::create(f)?);

                    writeln!(file, "R")?;
                    for (region, string) in [Region::QM1, Region::QM2, Region::Active]
                        .into_iter()
                        .zip(["-q", "-o", "-a"].into_iter())
                    {
                        if let Ok(l) = functions::get_atomlist(pdb, region) {
                            writeln!(
                                file,
                                "A {} id {}",
                                string,
                                l.into_iter().map(|n| n.to_string()).join(",")
                            )?;
                        };
                    }
                    writeln!(file, "W -w")?;
                } else if let Err(e) = save_pdb(pdb, f, pdbtbx::StrictnessLevel::Loose) {
                    e.into_iter()
                        .try_for_each(|e| writeln!(io::stdout(), "{}", e))?;
                }
            }
            Some(Output::Overwrite) => {
                if let Err(e) = save_pdb(pdb, pdb_path, pdbtbx::StrictnessLevel::Loose) {
                    e.into_iter()
                        .try_for_each(|e| writeln!(io::stdout(), "{}", e))?;
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
        format!(
            "{:?}",
            dispatch(&mode, &mut pdb, pdb_path).unwrap().unwrap()
        )
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
        assert!(dispatch(&mode, &mut pdb, pdb_path).unwrap().is_none());

        // let matches = get_matches(["Q", "-rl", "12"].into_iter());
        // let mode = Mode::new(&matches).unwrap();
        // assert!(dispatch(mode, &mut pdb, pdb_path).unwrap().is_none())
    }

    #[test]
    fn add_qm1() {
        let edit_action = get_edit_action(
            "tests/test_blank.pdb",
            ["A", "-q", "id", "4,1,9"].into_iter(),
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
            ["A", "-o", "id", "22,17,8"].into_iter(),
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
            ["A", "-a", "id", "22,17,8"].into_iter(),
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
            ["R", "-q", "id", "22,17,8"].into_iter(),
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
            ["R", "-a", "id", "22,17,8"].into_iter(),
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
            ["R", "-a", "id", "22,17,8"].into_iter(),
        );
    }

    #[test]
    #[should_panic]
    fn add_existing() {
        get_edit_action(
            "tests/test_full.pdb",
            ["A", "-q", "id", "22,17,8"].into_iter(),
        );
    }

    #[test]
    fn checked_add() {
        let edit_action = get_edit_action(
            "tests/test_overwrite.pdb",
            ["A", "-q", "id", "1,3,4"].into_iter(),
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
            ["R", "-a", "id", "1-3"].into_iter(),
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
            ["A", "-o", "id", "9,1-3"].into_iter(),
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
