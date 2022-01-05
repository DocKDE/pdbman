use std::collections::HashSet;

use crate::options::Region;

use anyhow::Result;
use colored::Colorize;
use itertools::Itertools;
use pdbtbx::{
    Atom, AtomConformerResidueChainModel, ContainsAtomConformer, ContainsAtomConformerResidue, PDB,
};
use rayon::prelude::ParallelIterator;

type AtomList = Vec<usize>;
type ResidueList = Vec<isize>;

/// Takes an Atom struct as point of origin and a radius in A. Returns a Vector of Atom IDs
/// of Atoms within the given radius wrapped in a Result. Origin can be included or excluded.
pub fn get_atom_sphere(
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
        .locate_within_distance(origin_atom.pos(), radius.powf(2.0))
        .map(|atom| atom.serial_number())
        .collect();

    if !include_self {
        sphere_atoms.retain(|&x| x != origin_atom.serial_number())
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
pub fn get_residue_sphere(
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
        .locate_within_distance(sphere_origin.atom().pos(), radius.powf(2.0))
        .flat_map(|atom_hier| atom_hier.residue().atoms().map(|atom| atom.serial_number()))
        .unique()
        .collect();

    let origin_res_atoms: AtomList = sphere_origin
        .residue()
        .atoms()
        .map(|atom| atom.serial_number())
        .collect();

    if !include_self {
        sphere_atoms.retain(|x| !origin_res_atoms.contains(x))
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
        .map(|a| a.serial_number())
        .collect::<Vec<usize>>();

    ensure!(!num_vec.is_empty(), "No atoms in the requested region!");
    Ok(num_vec)
}

// Get list of residue IDs for printing to stdout or file
pub fn get_residuelist(pdb: &PDB, region: Region) -> Result<ResidueList, anyhow::Error> {
    let filt_closure = match region {
        Region::QM1 => |a: &AtomConformerResidueChainModel| a.atom().occupancy() == 1.00,
        Region::QM2 => |a: &AtomConformerResidueChainModel| a.atom().occupancy() == 2.00,
        Region::Active => |a: &AtomConformerResidueChainModel| a.atom().b_factor() == 1.00,
    };

    let num_vec = pdb
        .atoms_with_hierarchy()
        .filter(filt_closure)
        .map(|a| a.residue().serial_number())
        .dedup()
        .collect::<Vec<isize>>();

    ensure!(!num_vec.is_empty(), "No residues in the requested region!");
    Ok(num_vec)
}

pub fn get_inverted(atomlist: &[usize], pdb: &PDB) -> Vec<usize> {
    let atom_set: HashSet<&usize> = HashSet::from_iter(atomlist);
    pdb.par_atoms()
        .filter(|a| !atom_set.contains(&a.serial_number()))
        .map(|a| a.serial_number())
        .collect()
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

    #[test]
    fn get_residuelist_test() {
        let pdb = test_pdb("tests/test_get_residuelist.pdb");
        let qm1_residues = get_residuelist(&pdb, Region::QM1).unwrap();
        let qm2_residues = get_residuelist(&pdb, Region::QM2).unwrap();
        let active_residues = get_residuelist(&pdb, Region::Active).unwrap();

        assert_eq!(qm1_residues, vec![1, 7]);
        assert_eq!(qm2_residues, vec![3, 4]);
        assert_eq!(active_residues, vec![2, 3]);
    }
}
