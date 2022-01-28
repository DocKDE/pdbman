use crate::options::{Distance, Region, Target};

use anyhow::Result;
use comfy_table::modifiers::{UTF8_ROUND_CORNERS, UTF8_SOLID_INNER_BORDERS};
use comfy_table::presets::{UTF8_BORDERS_ONLY, UTF8_FULL};
use comfy_table::{Row, Table};
use pdbtbx::{ContainsAtomConformer, ContainsAtomConformerResidue, PDB};

/// Finds and prints all contacts present in the PDB file structure. Definition of
/// 'contact' is given by the 'level' arg which is 1.0A for Clashes and depends
/// on the atomic radius of the involved atoms for Contacts.
pub fn find_contacts(pdb: &PDB, level: Distance) -> Result<Table, anyhow::Error> {
    let mut table = Table::new();
    // table.set_format(*format::consts::FORMAT_NO_LINESEP_WITH_TITLE);
    table
        .load_preset(UTF8_BORDERS_ONLY)
        .apply_modifier(UTF8_ROUND_CORNERS);
    table.set_header(Row::from(vec![
        "Atom ID 1",
        "Atom Name 1",
        "Residue Name 1",
        "Atom ID 2",
        "Atom Name 2",
        "Residue Name 2",
        "Distance",
    ]));

    let tree = pdb.create_hierarchy_rtree();

    for atom_hier in pdb.atoms_with_hierarchy() {
        let radius: f64 = match level {
            Distance::Clashes => 1.0,
            Distance::Contacts => atom_hier
                .atom()
                .atomic_radius()
                .ok_or_else(|| {
                    anyhow!(
                        "{} {}",
                        "No radius found for given atom type:",
                        atom_hier.atom().element()
                    )
                })?
                .powf(2.0),
        };

        let contacts = tree.locate_within_distance(atom_hier.atom().pos(), radius);

        for other_atom_hier in contacts {
            // This eliminates duplicate entries
            if other_atom_hier.atom() < atom_hier.atom()
            // This eliminates atoms from same residue
                && other_atom_hier.residue() != atom_hier.residue()
                // This eliminates neighboring residues
                && !(other_atom_hier.atom().name() == "C"
                    && atom_hier.atom().name() == "N"
                    && other_atom_hier.residue().serial_number() + 1
                        == atom_hier.residue().serial_number())
            {
                let distance = other_atom_hier.atom().distance(atom_hier.atom());

                // if distance <= 0.75 && distance > 0.5 {
                //     table.add_row(Row::from(vec![
                //                    bByFd =>
                //                    other_atom_hier.atom().serial_number(),
                //                    other_atom_hier.atom().name(),
                //                    other_atom_hier.residue().name().unwrap_or(""),
                //                    atom_hier.atom().serial_number(),
                //                    atom_hier.atom().name(),
                //                    atom_hier.residue().name().unwrap_or(""),
                //                    format!("{:.2}", distance)
                //     ]));
                // } else if distance <= 0.5 {
                //     table.add_row(row![
                //         bBrFd =>
                //         other_atom_hier.atom().serial_number(),
                //         other_atom_hier.atom().name(),
                //         other_atom_hier.residue().name().unwrap_or(""),
                //         atom_hier.atom().serial_number(),
                //         atom_hier.atom().name(),
                //         atom_hier.residue().name().unwrap_or(""),
                //         format!("{:.2}", distance)
                //     ]);
                // } else {
                table.add_row(Row::from(vec![
                    other_atom_hier.atom().serial_number().to_string(),
                    other_atom_hier.atom().name().to_owned(),
                    other_atom_hier.residue().name().unwrap_or("").to_owned(),
                    atom_hier.atom().serial_number().to_string(),
                    atom_hier.atom().name().to_owned(),
                    atom_hier.residue().name().unwrap_or("").to_owned(),
                    format!("{:.2}", distance).to_string(),
                ]));
                // }
            }
        }
    }

    // Header and delimiter lines also count
    ensure!(table.lines().count() > 4, "No contacts found!");
    Ok(table)
}

pub fn analyze(
    pdb: &PDB,
    region: Option<Region>,
    target: Option<Target>,
) -> Result<(Table, Option<Table>), anyhow::Error> {
    let mut qm1_residue_list = Vec::new();
    let mut qm1_atom_list = Vec::new();
    let mut qm2_residue_list = Vec::new();
    let mut qm2_atom_list = Vec::new();
    let mut active_residue_list = Vec::new();
    let mut active_atom_list = Vec::new();
    let mut atom_num: u32 = 0;
    let mut res_num: u32 = 0;

    for residue in pdb.residues() {
        res_num += 1;
        for atom in residue.atoms() {
            atom_num += 1;
            if atom.occupancy() == 1.00 {
                qm1_residue_list.push(residue);
                qm1_atom_list.push(atom);
            } else if atom.occupancy() == 2.00 {
                qm2_residue_list.push(residue);
                qm2_atom_list.push(atom);
            }

            if atom.b_factor() == 1.00 {
                active_residue_list.push(residue);
                active_atom_list.push(atom);
            }
        }
    }

    // It's much faster to just dedup the vecs once than check whether
    // the respective item is already present every time. Also since
    // the residues are iterated serially, sorting should not be necessary.

    // qm1_residue_list.sort();
    // qm2_residue_list.sort();
    // active_residue_list.sort();

    qm1_residue_list.dedup();
    qm2_residue_list.dedup();
    active_residue_list.dedup();

    let mut basic_table = Table::new();
    basic_table
        .load_preset(UTF8_FULL)
        .apply_modifier(UTF8_ROUND_CORNERS)
        .apply_modifier(UTF8_SOLID_INNER_BORDERS);

    basic_table.set_header(Row::from(vec!["", "# of Atoms", "# of Residues"]));
    basic_table.add_row(Row::from(vec![
        "QM1".to_owned(),
        qm1_atom_list.len().to_string(),
        qm1_residue_list.len().to_string(),
    ]));
    basic_table.add_row(Row::from(vec![
        "QM2".to_owned(),
        qm2_atom_list.len().to_string(),
        qm2_residue_list.len().to_string(),
    ]));
    basic_table.add_row(Row::from(vec![
        "Active".to_owned(),
        active_atom_list.len().to_string(),
        active_residue_list.len().to_string(),
    ]));
    basic_table.add_row(Row::from(vec![
        "Total".to_owned(),
        atom_num.to_string(),
        res_num.to_string(),
    ]));
    // let basic_table = table!(
    //     ["", "# of Atoms", "# of Residues"],
    //     ["QM1", qm1_atom_list.len(), qm1_residue_list.len()],
    //     ["QM2", qm2_atom_list.len(), qm2_residue_list.len()],
    //     ["Active", active_atom_list.len(), active_residue_list.len()],
    //     ["Total", atom_num, res_num]
    // );

    let mut detailed_table = None;

    if target == Some(Target::Residues) {
        let residue_list = match region {
            Some(Region::QM1) => qm1_residue_list,
            Some(Region::QM2) => qm2_residue_list,
            Some(Region::Active) => active_residue_list,
            // Impossible because if target is Some(..), a region is required by clap
            None => unreachable!(),
        };

        ensure!(
            !residue_list.is_empty(),
            "No Residues found in given region!"
        );

        let mut residue_table = Table::new();
        residue_table
            .load_preset(UTF8_FULL)
            .apply_modifier(UTF8_ROUND_CORNERS)
            .apply_modifier(UTF8_SOLID_INNER_BORDERS);
        residue_table.set_header(Row::from(vec![
            "Residue ID",
            "Residue Name",
            "# of Atoms",
            match region {
                Some(Region::QM1) => "# of QM1 Atoms",
                Some(Region::QM2) => "# of QM2 Atoms",
                Some(Region::Active) => "# of Active Atoms",
                None => unreachable!(),
            },
        ]));

        for residue in residue_list {
            let mut resid_atoms: u32 = 0;
            let mut atom_counter: u32 = 0;
            for atom in residue.atoms() {
                atom_counter += 1;
                if (region == Some(Region::QM1) && atom.occupancy() == 1.00)
                    || (region == Some(Region::QM2) && atom.occupancy() == 2.00)
                    || (region == Some(Region::Active) && atom.b_factor() == 1.00)
                {
                    resid_atoms += 1;
                }
            }

            residue_table.add_row(Row::from(vec![
                residue.serial_number().to_string() + residue.insertion_code().unwrap_or(""),
                residue.name().unwrap_or("").to_owned(),
                atom_counter.to_string(),
                resid_atoms.to_string(),
            ]));
        }
        detailed_table = Some(residue_table);
    } else if target == Some(Target::Atoms) {
        let (atom_list, residue_list) = match region {
            Some(Region::QM1) => (qm1_atom_list, qm1_residue_list),
            Some(Region::QM2) => (qm2_atom_list, qm2_residue_list),
            Some(Region::Active) => (active_atom_list, active_residue_list),
            None => unreachable!(),
        };

        ensure!(!atom_list.is_empty(), "No Atoms found in given region!");

        let mut atom_table = Table::new();
        atom_table
            .load_preset(UTF8_FULL)
            .apply_modifier(UTF8_ROUND_CORNERS)
            .apply_modifier(UTF8_SOLID_INNER_BORDERS);
        atom_table.set_header(Row::from(vec![
            "Atom ID",
            "Atom name",
            "Residue ID",
            "Residue Name",
            "QM",
            "Active",
        ]));

        for residue in residue_list {
            for atom in residue.atoms() {
                if atom_list.contains(&atom) {
                    atom_table.add_row(Row::from(vec![
                        atom.serial_number().to_string(),
                        atom.name().to_owned(),
                        residue.serial_number().to_string()
                            + residue.insertion_code().unwrap_or(""),
                        residue.name().unwrap_or("").to_owned(),
                        atom.occupancy().to_string(),
                        atom.b_factor().to_string(),
                    ]));
                }
            }
        }
        detailed_table = Some(atom_table);
    }
    Ok((basic_table, detailed_table))
}

#[cfg(test)]
mod tests {
    use super::*;
    use pdbtbx::StrictnessLevel;

    fn test_pdb(path: &str) -> PDB {
        let (pdb, _) = pdbtbx::open_pdb(path, StrictnessLevel::Strict).unwrap();
        pdb
    }

    #[test]
    fn analyze_atoms_test() {
        let pdb = test_pdb("tests/test_overwrite.pdb");
        let (basic, qm1_atoms) = analyze(&pdb, Some(Region::QM1), Some(Target::Atoms)).unwrap();
        let mut basic_table = Table::new();
        basic_table
            .load_preset(UTF8_FULL)
            .apply_modifier(UTF8_ROUND_CORNERS)
            .apply_modifier(UTF8_SOLID_INNER_BORDERS);
        basic_table.set_header(Row::from(vec!["", "# of Atoms", "# of Residues"]));
        basic_table.add_row(Row::from(vec![
            "QM1".to_owned(),
            1_u8.to_string(),
            1_u8.to_string(),
        ]));
        basic_table.add_row(Row::from(vec![
            "QM2".to_owned(),
            1_u8.to_string(),
            1_u8.to_string(),
        ]));
        basic_table.add_row(Row::from(vec![
            "Active".to_owned(),
            1_u8.to_string(),
            1_u8.to_string(),
        ]));
        basic_table.add_row(Row::from(vec![
            "Total".to_owned(),
            83_u8.to_string(),
            7_u8.to_string(),
        ]));

        assert_eq!(format!("{}", basic), format!("{}", basic_table));

        let mut atom_table = Table::new();
        atom_table
            .load_preset(UTF8_FULL)
            .apply_modifier(UTF8_ROUND_CORNERS)
            .apply_modifier(UTF8_SOLID_INNER_BORDERS);
        atom_table.set_header(Row::from(vec![
            "Atom ID",
            "Atom name",
            "Residue ID",
            "Residue Name",
            "QM",
            "Active",
        ]));

        atom_table.add_row(Row::from(vec![
            1_u8.to_string(),
            "N".to_owned(),
            1_u8.to_string(),
            "HIE".to_owned(),
            1.00.to_string(),
            0.00.to_string(),
        ]));
        assert_eq!(format!("{}", qm1_atoms.unwrap()), format!("{}", atom_table));
    }

    #[test]
    fn analyze_residues_test() {
        let pdb = test_pdb("tests/test_get_residuelist.pdb");
        let (basic, qm2_residues) =
            analyze(&pdb, Some(Region::QM2), Some(Target::Residues)).unwrap();
        let mut basic_table = Table::new();
        basic_table
            .load_preset(UTF8_FULL)
            .apply_modifier(UTF8_ROUND_CORNERS)
            .apply_modifier(UTF8_SOLID_INNER_BORDERS);
        basic_table.set_header(Row::from(vec!["", "# of Atoms", "# of Residues"]));
        basic_table.add_row(Row::from(vec![
            "QM1".to_owned(),
            20_u8.to_string(),
            2_u8.to_string(),
        ]));
        basic_table.add_row(Row::from(vec![
            "QM2".to_owned(),
            2_u8.to_string(),
            2_u8.to_string(),
        ]));
        basic_table.add_row(Row::from(vec![
            "Active".to_owned(),
            14_u8.to_string(),
            2_u8.to_string(),
        ]));
        basic_table.add_row(Row::from(vec![
            "Total".to_owned(),
            83_u8.to_string(),
            7_u8.to_string(),
        ]));

        assert_eq!(format!("{}", basic), format!("{}", basic_table));

        let mut residue_table = Table::new();
        residue_table
            .load_preset(UTF8_FULL)
            .apply_modifier(UTF8_ROUND_CORNERS)
            .apply_modifier(UTF8_SOLID_INNER_BORDERS);
        residue_table.set_header(Row::from(vec![
            "Residue ID",
            "Residue Name",
            "# of Atoms",
            "# of QM2 Atoms",
        ]));
        residue_table.add_row(Row::from(vec![
            3_u8.to_string(),
            "TYR".to_owned(),
            21_u8.to_string(),
            1_u8.to_string(),
        ]));
        residue_table.add_row(Row::from(vec![
            4_u8.to_string(),
            "VAL".to_owned(),
            16_u8.to_string(),
            1_u8.to_string(),
        ]));
        assert_eq!(
            format!("{}", qm2_residues.unwrap()),
            format!("{}", residue_table)
        );
    }

    // #[test]
    // fn contacts_test() {
    //     let pdb = test_pdb("tests/test_clash.pdb");
    //     let clashes = find_contacts(&pdb, Distance::Clashes).unwrap();
    //     let contacts = find_contacts(&pdb, Distance::Contacts).unwrap();

    //     // Because prettytable table formatting differs between tables created from
    //     // scratch and tables read from csv the test is run by creating a table, saving
    //     // it as csv, re-reading it and checking against previously created csv tables.
    //     // This is alsy necessary because the color coding used in the binary is not saved
    //     // to csv.
    //     let clashes_out = File::create("tests/test_clashes.csv_tmp").unwrap();
    //     let contacts_out = File::create("tests/test_contacts.csv_tmp").unwrap();

    //     clashes.to_csv(clashes_out).unwrap();
    //     contacts.to_csv(contacts_out).unwrap();

    //     let clashes_in = Table::from_csv_file("tests/test_clashes.csv_tmp").unwrap();
    //     let contacts_in = Table::from_csv_file("tests/test_contacts.csv_tmp").unwrap();

    //     assert_eq!(
    //         clashes_in,
    //         Table::from_csv_file("tests/test_clashes.csv").unwrap()
    //     );
    //     assert_eq!(
    //         contacts_in,
    //         Table::from_csv_file("tests/test_contacts.csv").unwrap()
    //     );

    //     remove_file("tests/test_clashes.csv_tmp").unwrap();
    //     remove_file("tests/test_contacts.csv_tmp").unwrap();
    // }
}
