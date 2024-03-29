use core::fmt;
use std::collections::HashSet;

use crate::{
    options::{Partial, Region},
    selection::{convert_result, parse_selection, Conjunction, Selection},
};

use super::{parse_atomic_list, parse_residue_list};
use anyhow::Result;
use colored::Colorize;
use comfy_table::{
    modifiers::{UTF8_ROUND_CORNERS, UTF8_SOLID_INNER_BORDERS},
    presets::UTF8_FULL,
    Row, Table,
};
use itertools::Itertools;
use pdbtbx::{
    Atom, AtomConformerResidueChainModel, ContainsAtomConformer, ContainsAtomConformerResidue,
    Residue, PDB,
};
use rayon::{iter::FromParallelIterator, prelude::ParallelIterator};

pub enum AtomMeasurement {
    Distance(f64),
    Angle(f64),
    Dihedral(f64),
}

impl fmt::Display for AtomMeasurement {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                AtomMeasurement::Distance(d) => format!("Distance: {:.3} \u{212B}", d),
                AtomMeasurement::Angle(a) => format!("Angle: {:.1}°", a),
                AtomMeasurement::Dihedral(d) => format!("Dihedral: {:.1}°", d),
            }
        )
    }
}

type AtomList = Vec<usize>;
// type ResidueList = Vec<isize>;

/// Takes an Atom struct as point of origin and a radius in A. Returns a Vector of Atom IDs
/// of Atoms within the given radius wrapped in a Result. Origin can be included or excluded.
fn get_atom_sphere(
    pdb: &PDB,
    origin_id: usize,
    // origin: &Atom,
    // origin: impl ContainsAtomConformerResidueChain,
    radius: f64,
    include_self: bool,
) -> Result<AtomList, anyhow::Error> {
    let origin_atom = pdb
        .par_atoms()
        .find_first(|a| a.serial_number() == origin_id)
        .ok_or_else::<_, _>(|| {
            anyhow!(
                "{}: '{}'",
                "\nNO ATOM WITH FOUND WITH SERIAL NUMBER".red(),
                origin_id.to_string().blue(),
            )
        })?;

    let tree = pdb.create_atom_rtree();
    let mut sphere_atoms: AtomList = tree
        .locate_within_distance(origin_atom.pos(), radius.powi(2))
        .map(|atom| atom.serial_number())
        .collect();
    // let test = tree.nearest_neighbor_iter_with_distance_2(&origin_atom.pos());
    // let mut sphere_atoms = Vec::new();
    // let mut distances = Vec::new();

    // for (atom, distance) in test {
    //     if distance <= radius {
    //         sphere_atoms.push(atom.serial_number());
    //         distances.push(distance);
    //     }
    // }

    if !include_self {
        sphere_atoms.retain(|&x| x != origin_atom.serial_number());
    }

    sphere_atoms.sort_unstable();

    ensure!(
        !sphere_atoms.is_empty(),
        "Calculated sphere doesn't contain any atoms."
    );
    Ok(sphere_atoms)
}

/// Takes an Atom struct as point of origin, a radius in A. Returns a Vector of Atom IDs
/// that belong to a Residue that had at least one Atom within the given radius wrapped in a Result.
/// Origin residue can be included or excluded.
fn get_residue_sphere(
    pdb: &PDB,
    origin_id: usize,
    // origin: impl ContainsAtomConformerResidueChain,
    radius: f64,
    include_self: bool,
) -> Result<AtomList, anyhow::Error> {
    let sphere_origin = pdb
        .atoms_with_hierarchy()
        .find(|a| a.atom().serial_number() == origin_id)
        .ok_or_else::<_, _>(|| {
            anyhow!(
                "{}: '{}'",
                "\nNO ATOM WITH FOUND WITH SERIAL NUMBER".red(),
                origin_id.to_string().blue(),
            )
        })?;
    let tree = pdb.create_hierarchy_rtree();

    let mut sphere_atoms: AtomList = tree
        .locate_within_distance(sphere_origin.atom().pos(), radius.powi(2))
        .flat_map(|atom_hier| atom_hier.residue().atoms().map(Atom::serial_number))
        .unique()
        .collect();

    let origin_res_atoms: AtomList = sphere_origin
        .residue()
        .atoms()
        .map(Atom::serial_number)
        .collect();

    if !include_self {
        sphere_atoms.retain(|x| !origin_res_atoms.contains(x));
    }

    ensure!(
        !sphere_atoms.is_empty(),
        "Calculated sphere doesn't contain any atoms."
    );
    Ok(sphere_atoms)
}

// Get list of atom IDs for printing to stdout or file
pub fn get_atomlist(pdb: &PDB, region: Region) -> Result<AtomList, anyhow::Error> {
    let filt_closure = match region {
        Region::QM1 => |a: &Atom| a.occupancy() == 1.00,
        Region::QM2 => |a: &Atom| a.occupancy() == 2.00,
        Region::Active => |a: &Atom| a.b_factor() == 1.00,
    };

    let num_vec = pdb
        .par_atoms()
        .filter(|&a| filt_closure(a))
        .map(Atom::serial_number)
        .collect::<Vec<usize>>();

    ensure!(!num_vec.is_empty(), "No atoms in the requested region!");
    Ok(num_vec)
}

// Get list of residue IDs for printing to stdout or file
// pub fn get_residuelist(pdb: &PDB, region: Region) -> Result<ResidueList, anyhow::Error> {
//     let filt_closure = match region {
//         Region::QM1 => |a: &AtomConformerResidueChainModel| a.atom().occupancy() == 1.00,
//         Region::QM2 => |a: &AtomConformerResidueChainModel| a.atom().occupancy() == 2.00,
//         Region::Active => |a: &AtomConformerResidueChainModel| a.atom().b_factor() == 1.00,
//     };

//     let num_vec = pdb
//         .atoms_with_hierarchy()
//         .filter(filt_closure)
//         .map(|a| a.residue().serial_number())
//         .dedup()
//         .collect::<Vec<isize>>();

//     ensure!(!num_vec.is_empty(), "No residues in the requested region!");
//     Ok(num_vec)
// }

fn get_atomlist_from_residuelist(
    list: &[isize],
    pdb: &PDB,
    partial: Option<Partial>,
) -> Vec<usize> {
    let residue_set: HashSet<&isize> = list.iter().collect();
    match partial {
        None => pdb
            .atoms_with_hierarchy()
            .filter(|a| residue_set.contains(&a.residue().serial_number()))
            .map(|a| a.atom().serial_number())
            .collect(),
        Some(p) => pdb
            .atoms_with_hierarchy()
            .filter(|a| match p {
                Partial::Backbone => a.is_backbone(),
                Partial::Sidechain => a.is_sidechain(),
            } && residue_set.contains(&a.residue().serial_number()))
            .map(|a| a.atom().serial_number())
            .collect(),
    }
}

fn get_inverted(atomlist: &[usize], pdb: &PDB) -> Vec<usize> {
    let atom_set: HashSet<&usize> = atomlist.iter().collect();
    pdb.par_atoms()
        .filter(|a| !atom_set.contains(&a.serial_number()))
        .map(Atom::serial_number)
        .collect()
}

fn verify_atomlist(list: &[usize], pdb: &PDB) -> Result<(), anyhow::Error> {
    let output_set: HashSet<usize> = list.iter().copied().collect();
    let pdb_set: HashSet<usize> = HashSet::from_par_iter(pdb.par_atoms().map(Atom::serial_number));
    let mut missing_atoms = output_set.difference(&pdb_set).peekable();

    ensure!(
        missing_atoms.peek().is_none(),
        "No atom(s) found with serial number(s): {}",
        missing_atoms.format(",")
    );

    Ok(())
}

fn verify_residuelist(list: &[isize], pdb: &PDB) -> Result<(), anyhow::Error> {
    let output_set: HashSet<isize> = list.iter().copied().collect();
    let pdb_set: HashSet<isize> =
        HashSet::from_par_iter(pdb.par_residues().map(Residue::serial_number));
    let mut missing_atoms = output_set.difference(&pdb_set).peekable();

    ensure!(
        missing_atoms.peek().is_none(),
        "No residue(s) found with serial number(s): {}",
        missing_atoms.format(",")
    );

    Ok(())
}

fn get_atoms_from_selection(
    s: Selection,
    pdb: &PDB,
    partial: Option<Partial>,
) -> Result<Vec<usize>, anyhow::Error> {
    let to_invert: bool;
    let mut atomvec = match s {
        Selection::ID { atomlist, invert } => {
            to_invert = invert;
            verify_atomlist(&atomlist, pdb)?;
            atomlist
        }
        Selection::Name { atomlist, invert } => {
            to_invert = invert;
            parse_atomic_list(&atomlist.join(","), pdb)?
        }
        Selection::Resid { reslist, invert } => {
            to_invert = invert;
            verify_residuelist(&reslist, pdb)?;
            get_atomlist_from_residuelist(&reslist, pdb, partial)
        }
        Selection::Resname { reslist, invert } => {
            to_invert = invert;
            let res_list = parse_residue_list(&reslist.join(","), pdb)?;
            get_atomlist_from_residuelist(&res_list, pdb, partial)
        }
        Selection::Sphere {
            origin,
            radius,
            invert,
        } => {
            to_invert = invert;
            get_atom_sphere(pdb, origin, radius, false)?
        }
        Selection::ResSphere {
            origin,
            radius,
            invert,
        } => {
            to_invert = invert;
            get_residue_sphere(pdb, origin, radius, false)?
        }
    };

    if to_invert {
        atomvec = get_inverted(&atomvec, pdb);
    }
    Ok(atomvec)
}

pub fn get_atomlist_from_input(
    input: &str,
    pdb: &PDB,
    partial: Option<Partial>,
) -> Result<Vec<usize>, anyhow::Error> {
    // Add a space to the given user input. This will make pest parse the
    // last character as a finished word resulting in more meaningful error messages.
    let input = input.to_owned() + " ";
    let (sele_vec, conj_vec) = convert_result(parse_selection(&input), &input)?;

    let mut sele_iter = sele_vec.into_iter();
    let initial_sel = sele_iter.next().unwrap();
    let initial_atomvec = get_atoms_from_selection(initial_sel, pdb, partial)?;

    if conj_vec.is_empty()
    // if vec of conjunctions is not present, the vec of selections can only have one element
    {
        Ok(initial_atomvec)
    } else {
        let atom_set = sele_iter.zip(conj_vec.iter()).try_fold(
            HashSet::from_iter(initial_atomvec),
            |acc, (sele, conj)| -> Result<HashSet<usize>, anyhow::Error> {
                match conj {
                    Conjunction::Or => Ok(acc
                        .union(&HashSet::from_iter(get_atoms_from_selection(
                            sele, pdb, partial,
                        )?))
                        .copied()
                        .collect()),
                    Conjunction::And => Ok(acc
                        .intersection(&HashSet::from_iter(get_atoms_from_selection(
                            sele, pdb, partial,
                        )?))
                        .copied()
                        .collect()),
                }
            },
        )?;
        Ok(atom_set.into_iter().collect())
    }
}

pub fn get_measurements(
    atomlist: &[usize],
    pdb: &PDB,
) -> Result<(Table, AtomMeasurement), anyhow::Error> {
    let atom_vec = atomlist
        .iter()
        .map(|i| {
            pdb.atoms_with_hierarchy()
                .find(|a| a.atom().serial_number() == *i)
                .ok_or_else(|| anyhow!("No atom found with ID: {}", i))
        })
        .collect::<Result<Vec<AtomConformerResidueChainModel>, anyhow::Error>>()?;

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

    for atom in &atom_vec {
        table.add_row(Row::from(vec![
            atom.atom().serial_number().to_string(),
            atom.atom().name().to_owned(),
            atom.residue().serial_number().to_string()
                + atom.residue().insertion_code().unwrap_or(""),
            atom.residue().name().unwrap_or("").to_owned(),
            atom.atom().occupancy().to_string(),
            atom.atom().b_factor().to_string(),
        ]));
    }

    match atom_vec.len() {
        2 => Ok((
            table,
            AtomMeasurement::Distance(atom_vec[0].atom().distance(atom_vec[1].atom())),
        )),
        3 => {
            let a = atom_vec[0].atom();
            let b = atom_vec[1].atom();
            let c = atom_vec[2].atom();
            Ok((table, AtomMeasurement::Angle(a.angle(b, c))))
        }
        4 => {
            let a = atom_vec[0].atom();
            let b = atom_vec[1].atom();
            let c = atom_vec[2].atom();
            let d = atom_vec[3].atom();
            Ok((table, AtomMeasurement::Dihedral(a.dihedral(b, c, d))))
        }
        _ => unreachable!(),
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
    fn atom_sphere_test() {
        let pdb = test_pdb("tests/test_blank.pdb");
        // let origin = pdb
        //     .atoms_with_hierarchy()
        //     .find(|x| x.atom().serial_number() == 26)
        //     .unwrap();
        let atom_list_incl = vec![21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 45, 46];
        let atom_list_excl = vec![21, 22, 23, 24, 25, 27, 28, 29, 30, 32, 45, 46];

        assert_eq!(get_atom_sphere(&pdb, 26, 4.0, true).unwrap().len(), 22);
        assert_eq!(
            get_atom_sphere(&pdb, 26, 3.0, true).unwrap(),
            atom_list_incl
        );

        assert_eq!(get_atom_sphere(&pdb, 26, 4.0, false).unwrap().len(), 21);
        assert_eq!(
            get_atom_sphere(&pdb, 26, 3.0, false).unwrap(),
            atom_list_excl
        );
    }

    #[test]
    fn residue_sphere_test() {
        let pdb = test_pdb("tests/test_blank.pdb");
        // let origin = pdb
        //     .atoms_with_hierarchy()
        //     .find(|x| x.atom().serial_number() == 26)
        //     .unwrap();
        let atom_list_excl = vec![
            19, 20, 21, 22, 23, 24, 25, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61,
            62,
        ];

        assert_eq!(get_residue_sphere(&pdb, 26, 4.0, true).unwrap().len(), 44);
        assert_eq!(get_residue_sphere(&pdb, 26, 4.0, false).unwrap().len(), 23);
        assert_eq!(
            get_residue_sphere(&pdb, 26, 4.0, false).unwrap(),
            atom_list_excl
        );
    }

    #[test]
    fn get_atomlist_test() {
        let pdb = test_pdb("tests/test_get_atomlist.pdb");
        let qm1_atoms = get_atomlist(&pdb, Region::QM1).unwrap();
        let qm2_atoms = get_atomlist(&pdb, Region::QM2).unwrap();
        let active_atoms = get_atomlist(&pdb, Region::Active).unwrap();

        assert_eq!(qm1_atoms, vec![1, 2, 4, 5, 6]);
        assert_eq!(qm2_atoms, vec![8, 9, 11, 12]);
        assert_eq!(active_atoms, vec![1, 2, 3, 5, 6, 8, 9]);
    }

    // #[test]
    // fn get_residuelist_test() {
    //     let pdb = test_pdb("tests/test_get_residuelist.pdb");
    //     let qm1_residues = get_residuelist(&pdb, Region::QM1).unwrap();
    //     let qm2_residues = get_residuelist(&pdb, Region::QM2).unwrap();
    //     let active_residues = get_residuelist(&pdb, Region::Active).unwrap();

    //     assert_eq!(qm1_residues, vec![1, 7]);
    //     assert_eq!(qm2_residues, vec![3, 4]);
    //     assert_eq!(active_residues, vec![2, 3]);
    // }
}
