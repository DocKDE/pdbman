use crate::options::Region;

use pdbtbx::{Atom, PDB};
use rayon::prelude::ParallelIterator;

/// This function edits the q or b value of a given residue list (used by ORCA as input for active section
/// selection) by Residue.
// pub fn edit_residues(
//     pdb: &mut PDB,
//     list: &[isize],
//     mode: &str,
//     partial: Option<Partial>,
//     region: Region,
// ) {
//     let edit = |a: &mut Atom| match mode {
//         "Add" => match region {
//             // This cannot fail because new values are finite and positive
//             Region::QM1 => a.set_occupancy(1.00).unwrap(),
//             Region::QM2 => a.set_occupancy(2.00).unwrap(),
//             Region::Active => a.set_b_factor(1.00).unwrap(),
//         },
//         "Remove" => match region {
//             Region::QM1 | Region::QM2 => a.set_occupancy(0.00).unwrap(),
//             Region::Active => a.set_b_factor(0.00).unwrap(),
//         },
//         _ => unreachable!(),
//     };

//     pdb.atoms_with_hierarchy_mut().for_each(|mut res| {
//         if list.contains(&res.residue().serial_number()) {
//             match partial {
//                 None => edit(res.atom_mut()),
//                 Some(Partial::Sidechain) => {
//                     if !res.atom().is_backbone() {
//                         edit(res.atom_mut())
//                     }
//                 }
//                 Some(Partial::Backbone) => {
//                     if res.atom().is_backbone() {
//                         edit(res.atom_mut())
//                     }
//                 }
//             }
//         }
//     });
// }

/// This functions edits the q or b value of the PDB file (used by ORCA as input for QM region
/// selection) by Atom.
pub fn edit_atoms(pdb: &mut PDB, list: &[usize], mode: &str, region: Region) {
    let edit = |a: &mut Atom| match mode {
        "Add" => match region {
            // This cannot fail because new values are finite and positive
            Region::QM1 => a.set_occupancy(1.00).unwrap(),
            Region::QM2 => a.set_occupancy(2.00).unwrap(),
            Region::Active => a.set_b_factor(1.00).unwrap(),
        },
        "Remove" => match region {
            Region::QM1 | Region::QM2 => a.set_occupancy(0.00).unwrap(),
            Region::Active => a.set_b_factor(0.00).unwrap(),
        },
        _ => unreachable!(),
    };

    pdb.par_atoms_mut()
        .filter(|a| list.contains(&a.serial_number()))
        .for_each(|a| edit(a));
}

/// Removes a whole region from PDB file.
pub fn remove_region(pdb: &mut PDB, region: Option<Region>) {
    match region {
        Some(Region::QM1) => pdb
            .par_atoms_mut()
            .filter(|a| a.occupancy() == 1.00)
            .for_each(|a| a.set_occupancy(0.00).unwrap()),
        Some(Region::QM2) => pdb
            .par_atoms_mut()
            .filter(|a| a.occupancy() == 2.00)
            .for_each(|a| a.set_occupancy(0.00).unwrap()),
        Some(Region::Active) => pdb
            .par_atoms_mut()
            .for_each(|a| a.set_b_factor(0.00).unwrap()),
        None => pdb.par_atoms_mut().for_each(|a| {
            a.set_occupancy(0.00).unwrap();
            a.set_b_factor(0.00).unwrap()
        }),
    }
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
    fn edit_atoms_test() {
        let mut pdb = test_pdb("tests/test_blank.pdb");
        let atom_id_list = &[1, 5, 9];
        edit_atoms(&mut pdb, atom_id_list, "Add", Region::QM1);
        edit_atoms(&mut pdb, atom_id_list, "Add", Region::Active);

        let atom_list: Vec<&Atom> = pdb
            .atoms()
            .filter(|x| atom_id_list.contains(&x.serial_number()))
            .collect();
        for atom in atom_list.iter() {
            assert_eq!(atom.occupancy(), 1.00);
            assert_eq!(atom.b_factor(), 1.00);
        }

        edit_atoms(&mut pdb, atom_id_list, "Remove", Region::QM1);
        edit_atoms(&mut pdb, atom_id_list, "Remove", Region::Active);
        let atom_list: Vec<&Atom> = pdb
            .atoms()
            .filter(|x| atom_id_list.contains(&x.serial_number()))
            .collect();
        for atom in atom_list.iter() {
            assert_eq!(atom.occupancy(), 0.00);
            assert_eq!(atom.b_factor(), 0.00);
        }
    }

    #[test]
    // fn edit_residues_test_all() {
    //     let mut pdb = test_pdb("tests/test_blank.pdb");
    //     let res_id_list = &vec![2, 4];

    //     edit_residues(&mut pdb, res_id_list, "Add", None, Region::QM2);
    //     edit_residues(&mut pdb, res_id_list, "Add", None, Region::Active);

    //     let res_list = pdb
    //         .residues()
    //         .filter(|x| res_id_list.contains(&x.serial_number()));

    //     for residue in res_list {
    //         for atom in residue.atoms() {
    //             assert_eq!(atom.occupancy(), 2.00);
    //             assert_eq!(atom.b_factor(), 1.00);
    //         }
    //     }

    //     edit_residues(&mut pdb, res_id_list, "Remove", None, Region::QM2);
    //     edit_residues(&mut pdb, res_id_list, "Remove", None, Region::Active);

    //     let res_list = pdb
    //         .residues()
    //         .filter(|x| res_id_list.contains(&x.serial_number()));

    //     for residue in res_list {
    //         for atom in residue.atoms() {
    //             assert_eq!(atom.occupancy(), 0.00);
    //             assert_eq!(atom.b_factor(), 0.00);
    //         }
    //     }
    // }
    #[test]
    // fn edit_residues_test_side() {
    //     let mut pdb = test_pdb("tests/test_blank.pdb");
    //     let res_id_list = &vec![2, 4];

    //     edit_residues(
    //         &mut pdb,
    //         res_id_list,
    //         "Add",
    //         Some(Partial::Sidechain),
    //         Region::QM2,
    //     );
    //     edit_residues(
    //         &mut pdb,
    //         res_id_list,
    //         "Add",
    //         Some(Partial::Sidechain),
    //         Region::Active,
    //     );

    //     let res_list = pdb
    //         .residues()
    //         .filter(|x| res_id_list.contains(&x.serial_number()));

    //     let atom_list = res_list
    //         .flat_map(|x| x.atoms())
    //         .filter(|x| !x.is_backbone());

    //     for atom in atom_list {
    //         assert_eq!(atom.occupancy(), 2.00);
    //         assert_eq!(atom.b_factor(), 1.00);
    //     }
    // }
    #[test]
    // fn edit_residues_test_back() {
    //     let mut pdb = test_pdb("tests/test_blank.pdb");
    //     let res_id_list = &vec![2, 4];
    //     edit_residues(
    //         &mut pdb,
    //         res_id_list,
    //         "Add",
    //         Some(Partial::Backbone),
    //         Region::QM2,
    //     );
    //     edit_residues(
    //         &mut pdb,
    //         res_id_list,
    //         "Add",
    //         Some(Partial::Backbone),
    //         Region::Active,
    //     );

    //     let res_list = pdb
    //         .residues()
    //         .filter(|x| res_id_list.contains(&x.serial_number()));
    //     let atom_list = res_list.flat_map(|x| x.atoms()).filter(|x| x.is_backbone());

    //     for atom in atom_list {
    //         assert_eq!(atom.occupancy(), 2.00);
    //         assert_eq!(atom.b_factor(), 1.00);
    //     }
    // }
    #[test]
    fn remove_qm_region() {
        let mut pdb = test_pdb("tests/test_full.pdb");
        remove_region(&mut pdb, Some(Region::QM1));

        for atom in pdb.atoms() {
            assert_eq!(atom.occupancy(), 0.00);
            assert_eq!(atom.b_factor(), 1.00);
        }
    }

    #[test]
    fn remove_active_region() {
        let mut pdb = test_pdb("tests/test_full.pdb");
        remove_region(&mut pdb, Some(Region::Active));

        for atom in pdb.atoms() {
            assert_eq!(atom.occupancy(), 1.00);
            assert_eq!(atom.b_factor(), 0.00);
        }
    }

    #[test]
    fn remove_all_test() {
        let mut pdb = test_pdb("tests/test_full.pdb");
        remove_region(&mut pdb, None);

        for atom in pdb.atoms() {
            assert_eq!(atom.occupancy(), 0.00);
            assert_eq!(atom.b_factor(), 0.00);
        }
    }
}