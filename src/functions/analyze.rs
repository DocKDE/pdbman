use crate::options::{Distance, Region, Target};

// use std::io;
// use std::io::prelude::Write;

use anyhow::Result;
use pdbtbx::{ContainsAtomConformer, ContainsAtomConformerResidue, PDB};
use prettytable::{format, Table};

/// Finds and prints all contacts present in the PDB file structure. Definition of
/// 'contact' is given by the 'level' arg which is 1.0A for Clashes and depends
/// on the atomic radius of the involved atoms for Contacts.
pub fn find_contacts(pdb: &PDB, level: Distance) -> Result<Table, anyhow::Error> {
    let mut table = Table::new();
    table.set_format(*format::consts::FORMAT_NO_LINESEP_WITH_TITLE);
    table.set_titles(row![
        "Atom ID 1",
        "Atom Name 1",
        "Residue Name 1",
        "Atom ID 2",
        "Atom Name 2",
        "Residue Name 2",
        "Distance"
    ]);

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

                if distance <= 0.75 && distance > 0.5 {
                    table.add_row(row![
                        bByFd =>
                        other_atom_hier.atom().serial_number(),
                        other_atom_hier.atom().name(),
                        other_atom_hier.residue().name().unwrap_or(""),
                        atom_hier.atom().serial_number(),
                        atom_hier.atom().name(),
                        atom_hier.residue().name().unwrap_or(""),
                        format!("{:.2}", distance)
                    ]);
                } else if distance <= 0.5 {
                    table.add_row(row![
                        bBrFd =>
                        other_atom_hier.atom().serial_number(),
                        other_atom_hier.atom().name(),
                        other_atom_hier.residue().name().unwrap_or(""),
                        atom_hier.atom().serial_number(),
                        atom_hier.atom().name(),
                        atom_hier.residue().name().unwrap_or(""),
                        format!("{:.2}", distance)
                    ]);
                } else {
                    table.add_row(row![
                        other_atom_hier.atom().serial_number(),
                        other_atom_hier.atom().name(),
                        other_atom_hier.residue().name().unwrap_or(""),
                        atom_hier.atom().serial_number(),
                        atom_hier.atom().name(),
                        atom_hier.residue().name().unwrap_or(""),
                        format!("{:.2}", distance)
                    ]);
                }
            }
        }
    }

    ensure!(!table.is_empty(), "No contacts found!");
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
    let mut atom_num = 0;
    let mut res_num = 0;

    for residue in pdb.residues() {
        res_num += 1;
        for atom in residue.atoms() {
            atom_num += 1;
            if atom.occupancy() == 1.00 {
                qm1_residue_list.push(residue);
                qm1_atom_list.push(atom)
            }

            if atom.occupancy() == 2.00 {
                qm2_residue_list.push(residue);
                qm2_atom_list.push(atom)
            }

            if atom.b_factor() == 1.00 {
                active_residue_list.push(residue);
                active_atom_list.push(atom)
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

    let basic_table = table!(
        ["", "# of Atoms", "# of Residues"],
        ["QM1", qm1_atom_list.len(), qm1_residue_list.len()],
        ["QM2", qm2_atom_list.len(), qm2_residue_list.len()],
        ["Active", active_atom_list.len(), active_residue_list.len()],
        ["Total", atom_num, res_num]
    );

    let mut detailed_table = None;

    if target == Some(Target::Residues) {
        let residue_list = match region {
            Some(Region::QM1) => qm1_residue_list,
            Some(Region::QM2) => qm2_residue_list,
            Some(Region::Active) => active_residue_list,
            // Impossible because if target is Some(..), a region is required by clap
            None => unreachable!(),
        };
        if !residue_list.is_empty() {
            let mut residue_table = Table::new();
            residue_table.add_row(row![
                "Residue ID",
                "Residue Name",
                "# of Atoms",
                match region {
                    Some(Region::QM1) => "# of QM1 Atoms",
                    Some(Region::QM2) => "# of QM2 Atoms",
                    Some(Region::Active) => "# of Active Atoms",
                    None => unreachable!(),
                }
            ]);

            for residue in residue_list {
                let mut resid_atoms = 0;
                let mut atom_counter = 0;
                for atom in residue.atoms() {
                    atom_counter += 1;
                    if (region == Some(Region::QM1) && atom.occupancy() == 1.00)
                        || (region == Some(Region::QM2) && atom.occupancy() == 2.00)
                        || (region == Some(Region::Active) && atom.b_factor() == 1.00)
                    {
                        resid_atoms += 1;
                    }
                }

                residue_table.add_row(row![
                    residue.serial_number().to_string() + residue.insertion_code().unwrap_or(""),
                    residue.name().unwrap_or(""),
                    atom_counter,
                    resid_atoms,
                ]);
            }
            detailed_table = Some(residue_table)
        } else {
            bail!("No Residues found in given region!");
        }
    } else if target == Some(Target::Atoms) {
        let (atom_list, residue_list) = match region {
            Some(Region::QM1) => (qm1_atom_list, qm1_residue_list),
            Some(Region::QM2) => (qm2_atom_list, qm2_residue_list),
            Some(Region::Active) => (active_atom_list, active_residue_list),
            None => unreachable!(),
        };

        if !atom_list.is_empty() {
            let mut atom_table = Table::new();
            atom_table.add_row(row![
                "Atom ID",
                "Atom name",
                "Residue ID",
                "Residue Name",
                "QM",
                "Active"
            ]);

            for residue in residue_list {
                for atom in residue.atoms() {
                    if atom_list.contains(&atom) {
                        atom_table.add_row(row![
                            atom.serial_number(),
                            atom.name(),
                            residue.serial_number().to_string()
                                + residue.insertion_code().unwrap_or(""),
                            residue.name().unwrap_or(""),
                            atom.occupancy(),
                            atom.b_factor(),
                        ]);
                    }
                }
            }
            detailed_table = Some(atom_table)
        } else {
            bail!("No Atoms found in given region!")
        }
    }
    Ok((basic_table, detailed_table))
}

#[cfg(test)]
mod tests {
    use super::*;
    use pdbtbx::StrictnessLevel;
    use std::fs::{remove_file, File};

    fn test_pdb(path: &str) -> PDB {
        let (pdb, _) = pdbtbx::open_pdb(path, StrictnessLevel::Strict).unwrap();
        pdb
    }

    #[test]
    fn analyze_atoms_test() {
        let pdb = test_pdb("tests/test_overwrite.pdb");
        let (basic, qm1_atoms) = analyze(&pdb, Some(Region::QM1), Some(Target::Atoms)).unwrap();
        let basic_table = table!(
            ["", "# of Atoms", "# of Residues"],
            ["QM1", 1, 1],
            ["QM2", 1, 1],
            ["Active", 1, 1],
            ["Total", 83, 7]
        );
        assert_eq!(basic, basic_table);

        let mut atom_table = Table::new();
        atom_table.add_row(row![
            "Atom ID",
            "Atom name",
            "Residue ID",
            "Residue Name",
            "QM",
            "Active"
        ]);
        atom_table.add_row(row!(1, "N", 1, "HIE", 1.00, 0.00));
        assert_eq!(qm1_atoms.unwrap(), atom_table);
    }

    #[test]
    fn analyze_residues_test() {
        let pdb = test_pdb("tests/test_overwrite.pdb");
        let (basic, qm2_residues) = analyze(&pdb, Some(Region::QM2), Some(Target::Residues)).unwrap();
        let basic_table = table!(
            ["", "# of Atoms", "# of Residues"],
            ["QM1", 1, 1],
            ["QM2", 1, 1],
            ["Active", 1, 1],
            ["Total", 83, 7]
        );
        assert_eq!(basic, basic_table);

        let mut residue_table = Table::new();
        residue_table.add_row(row![
            "Residue ID",
            "Residue Name",
            "# of Atoms",
            "# of QM2 Atoms",
        ]);

        residue_table.add_row(row!(1, "HIE", 18, 1));
        assert_eq!(qm2_residues.unwrap(), residue_table);
    }

    #[test]
    fn contacts_test() {
        let pdb = test_pdb("tests/test_clash.pdb");
        let clashes = find_contacts(&pdb, Distance::Clashes).unwrap();
        let contacts = find_contacts(&pdb, Distance::Contacts).unwrap();

        // Because prettytable table formatting differs between tables created from
        // scratch and tables read from csv the test is run by creating a table, saving
        // it as csv, re-reading it and checking against previously created csv tables.
        // This is alsy necessary because the color coding used in the binary is not saved
        // to csv.
        let clashes_out = File::create("tests/test_clashes.csv_tmp").unwrap();
        let contacts_out = File::create("tests/test_contacts.csv_tmp").unwrap();

        clashes.to_csv(clashes_out).unwrap();
        contacts.to_csv(contacts_out).unwrap();

        let clashes_in = Table::from_csv_file("tests/test_clashes.csv_tmp").unwrap();
        let contacts_in = Table::from_csv_file("tests/test_contacts.csv_tmp").unwrap();

        assert_eq!(
            clashes_in,
            Table::from_csv_file("tests/test_clashes.csv").unwrap()
        );
        assert_eq!(
            contacts_in,
            Table::from_csv_file("tests/test_contacts.csv").unwrap()
        );

        remove_file("tests/test_clashes.csv_tmp").unwrap();
        remove_file("tests/test_contacts.csv_tmp").unwrap();
    }
}
