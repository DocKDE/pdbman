use crate::residue_ascii::RESIDUE_ASCII;

use std::borrow::ToOwned;

use comfy_table::modifiers::{UTF8_ROUND_CORNERS, UTF8_SOLID_INNER_BORDERS};
use comfy_table::presets::UTF8_FULL;
use comfy_table::{Row, Table};
use itertools::Itertools;
use pdbtbx::{ContainsAtomConformer, ContainsAtomConformerResidue, PDB};

/// Query Molecule for information. Depending on the input this will print a table of
/// Residues and/or Atoms with available information that were asked for.
pub fn query_atoms<'a>(
    pdb: &'a PDB,
    atom_list: &'a [usize],
) -> Result<(Table, Option<&'a str>), anyhow::Error> {
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
    ]));

    let mut key: Option<&str> = None;
    let mut resname_vec = Vec::new();

    for atom_hier in pdb.atoms_with_hierarchy() {
        if atom_list.contains(&atom_hier.atom().serial_number()) {
            resname_vec.push(atom_hier.residue().name().map(ToOwned::to_owned));

            table.add_row(Row::from(vec![
                atom_hier.atom().serial_number().to_string(),
                atom_hier.atom().name().to_owned(),
                atom_hier.residue().serial_number().to_string()
                    + atom_hier.residue().insertion_code().unwrap_or(""),
                atom_hier.residue().name().unwrap_or("").to_owned(),
                atom_hier.atom().occupancy().to_string(),
                atom_hier.atom().b_factor().to_string(),
            ]));
        }
    }

    ensure!(!resname_vec.is_empty(), "No atoms found in given selection");

    if resname_vec.iter().all_equal() {
        key = resname_vec[0].as_deref();
    }

    if let Some(k) = key {
        Ok((
            table,
            RESIDUE_ASCII
                .get(&k.to_uppercase().as_ref())
                .map(ToOwned::to_owned),
        ))
    } else {
        Ok((table, None))
    }
}

// This cannot fail because if no residues can be queried, the
// parse_residue_list would have returned an error beforehand
// pub fn query_residues(pdb: &PDB, residue_list: &[isize]) -> Result<Table, anyhow::Error> {
//     let mut table = Table::new();

//     table.add_row(row![
//         "Atom ID",
//         "Atom name",
//         "Residue ID",
//         "Residue Name",
//         "QM",
//         "Active"
//     ]);

//     // let mut key = None;
//     let mut key = None;

//     for atom_hier in pdb.atoms_with_hierarchy() {
//         // let res = (
//         //     atom_hier.residue.serial_number(),
//         //     atom_hier.residue.insertion_code().map(ToOwned::to_owned),
//         // );

//         if residue_list.contains(&atom_hier.residue().serial_number()) {
//             if residue_list.len() == 1 {
//                 // key = atom_hier.residue().name();
//                 key = atom_hier.residue().name().map(|s| s.to_owned());
//             }
//             table.add_row(row![
//                 atom_hier.atom().serial_number(),
//                 atom_hier.atom().name(),
//                 atom_hier.residue().serial_number().to_string()
//                     + atom_hier.residue().insertion_code().unwrap_or(""),
//                 atom_hier.residue().name().unwrap_or(""),
//                 atom_hier.atom().occupancy(),
//                 atom_hier.atom().b_factor(),
//             ]);
//         }
//     }

//     if let Some(k) = key {
//         if let Some(res_ascii) = RESIDUE_ASCII.get(&k.to_uppercase().as_ref()) {
//             writeln!(io::stdout(), "{}", res_ascii)
//                 .context("Failed to print residue depiction to stdout")?
//         }
//     }

//     Ok(table)
// }

#[cfg(test)]
mod tests {
    use super::*;
    use pdbtbx::{StrictnessLevel, PDB};

    fn test_pdb(path: &str) -> PDB {
        let (pdb, _) = pdbtbx::open_pdb(path, StrictnessLevel::Strict).unwrap();
        pdb
    }

    #[test]
    fn query_atoms_test() {
        let pdb = test_pdb("tests/test_blank.pdb");
        let (table, _) = query_atoms(&pdb, &[1, 5, 19]).unwrap();

        let mut test_table = Table::new();
        test_table
            .load_preset(UTF8_FULL)
            .apply_modifier(UTF8_ROUND_CORNERS)
            .apply_modifier(UTF8_SOLID_INNER_BORDERS);
        test_table.set_header(Row::from(vec![
            "Atom ID",
            "Atom name",
            "Residue ID",
            "Residue Name",
            "QM",
            "Active",
        ]));

        test_table.add_row(Row::from(vec![
            1_u8.to_string(),
            "N".to_owned(),
            1_u8.to_string(),
            "HIE".to_owned(),
            0.00.to_string(),
            0.00.to_string(),
        ]));
        test_table.add_row(Row::from(vec![
            5_u8.to_string(),
            "HA".to_owned(),
            1_u8.to_string(),
            "HIE".to_owned(),
            0.00.to_string(),
            0.00.to_string(),
        ]));
        test_table.add_row(Row::from(vec![
            19_u8.to_string(),
            "N".to_owned(),
            2_u8.to_string(),
            "GLY".to_owned(),
            0.00.to_string(),
            0.00.to_string(),
        ]));

        assert_eq!(format!("{}", table), format!("{}", test_table))
    }

    // #[test]
    // fn query_residue_test() {
    //     let pdb = test_pdb("tests/test_blank.pdb");
    //     let table = query_residues(&pdb, &[2, 7]).unwrap();

    //     let mut test_table = Table::new();
    //     test_table.add_row(row![
    //         "Atom ID",
    //         "Atom name",
    //         "Residue ID",
    //         "Residue Name",
    //         "QM",
    //         "Active"
    //     ]);

    //     test_table.add_row(row![19, "N", 2, "GLY", 0.00, 0.00]);
    //     test_table.add_row(row![20, "H", 2, "GLY", 0.00, 0.00]);
    //     test_table.add_row(row![21, "CA", 2, "GLY", 0.00, 0.00]);
    //     test_table.add_row(row![22, "HA2", 2, "GLY", 0.00, 0.00]);
    //     test_table.add_row(row![23, "HA3", 2, "GLY", 0.00, 0.00]);
    //     test_table.add_row(row![24, "C", 2, "GLY", 0.00, 0.00]);
    //     test_table.add_row(row![25, "O", 2, "GLY", 0.00, 0.00]);
    //     test_table.add_row(row![81, "O", 7, "WAT", 0.00, 0.00]);
    //     test_table.add_row(row![82, "H1", 7, "WAT", 0.00, 0.00]);
    //     test_table.add_row(row![83, "H2", 7, "WAT", 0.00, 0.00]);

    //     assert_eq!(table, test_table)
    // }
}
