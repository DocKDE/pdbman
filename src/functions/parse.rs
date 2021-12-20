use std::collections::HashSet;

use anyhow::Result;
use itertools::Itertools;
use lazy_regex::regex;
use pdbtbx::{ContainsAtomConformer, ContainsAtomConformerResidue, PDB};
use rayon::iter::FromParallelIterator;
use rayon::prelude::ParallelIterator;

use crate::options::Partial;

type AtomList = Vec<usize>;
// type ResidueList = Vec<isize>;

/// Takes a comma-separated list (usually from command line input) as string and parses it into
/// a vector of Atom IDs. The Input may be atom IDs or Atom Names
pub fn parse_atomic_list(input: &str, pdb: &PDB) -> Result<AtomList, anyhow::Error> {
    let mut output_vec: AtomList = vec![];

    match input
        .split(&[',', '-', ':'][..])
        .all(|x| x.parse::<usize>().is_ok())
    {
        true => {
            // let input_vec = input.split(',');

            for i in input.split(',') {
                // All the unwraps work because the iterator
                // is guaranteed to be passed two usize values
                if i.contains(&['-', ':'][..]) {
                    let mut iter = i
                        .split(&['-', ':'][..])
                        .map(|x| x.parse::<usize>().unwrap());
                    let range = iter.next().unwrap()..=iter.next().unwrap();
                    output_vec.extend(range);
                } else {
                    output_vec.push(i.parse().unwrap())
                }
            }

            let output_set: HashSet<usize> = HashSet::from_iter(output_vec.iter().copied());
            let pdb_set: HashSet<usize> =
                HashSet::from_par_iter(pdb.par_atoms().map(|a| a.serial_number()));
            let mut missing_atoms = output_set.difference(&pdb_set).peekable();

            ensure!(
                missing_atoms.peek().is_none(),
                "No atom(s) found with serial number(s): {}",
                missing_atoms.format(",")
            );
        }
        false => {
            let input_vec: Vec<&str> = input.split(',').collect();

            // This might not be necessary
            ensure!(
                !regex!(r"[:-]").is_match(input),
                "Ranges are not allowed for atom name inputs."
            );

            let input_set: HashSet<String> = input_vec.iter().map(|s| s.to_lowercase()).collect();
            let pdb_set: HashSet<String> =
                HashSet::from_par_iter(pdb.par_atoms().map(|a| a.name().to_lowercase()));
            let mut missing_atoms = input_set.difference(&pdb_set).peekable();

            ensure!(
                missing_atoms.peek().is_none(),
                "No atom(s) found with identifier(s): {}",
                missing_atoms.format(",")
            );

            output_vec = pdb
                .atoms()
                .filter(|x| {
                    input_vec
                        .iter()
                        .any(|y| x.name().to_lowercase() == y.to_lowercase())
                })
                .map(|x| x.serial_number())
                .collect();
        }
    }

    Ok(output_vec)
}

/// Parses a string (usually taken from command line) and returns a list of residues given by a tuple
/// of serial numbers and insertion codes. The input can be either a comma-separated list of serial numbers
/// and insertion codes or residues names.
pub fn parse_residue_list(
    input: &str,
    pdb: &PDB,
    partial: Option<Partial>,
) -> Result<AtomList, anyhow::Error> {
    let mut residue_set: HashSet<isize> = HashSet::new();
    let atom_vec: AtomList;

    match input
        .split(&[',', '-', ':'][..])
        .all(|x| x.parse::<isize>().is_ok())
    {
        true => {
            for i in input.split(',') {
                // All the unwraps work because the iterator
                // is guaranteed to be passed two isize values
                if i.contains(&['-', ':'][..]) {
                    let mut iter = i
                        .split(&['-', ':'][..])
                        .map(|x| x.parse::<isize>().unwrap());
                    let range = iter.next().unwrap()..=iter.next().unwrap();
                    residue_set.extend(range);
                } else {
                    residue_set.insert(i.parse().unwrap());
                }
            }

            // let residue_set: HashSet<isize> = HashSet::from_iter(residue_vec.iter().copied());

            let pdb_set: HashSet<isize> =
                HashSet::from_par_iter(pdb.par_residues().map(|a| a.serial_number()));
            let mut missing_residues = residue_set.difference(&pdb_set).peekable();

            ensure!(
                missing_residues.peek().is_none(),
                "No residue(s) found with serial number(s): {}",
                missing_residues.format(",")
            );

            atom_vec = {
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
            };
        }
        false => {
            let residue_set: HashSet<String> = input.split(',').map(|s| s.to_owned()).collect();

            // This might not be necessary
            ensure!(
                !regex!(r"[:-]").is_match(input),
                "Ranges are not allowed for residue name inputs."
            );

            // let input_set: HashSet<String> = input_vec.iter().map(|s| s.to_lowercase()).collect();
            let pdb_set: HashSet<String> = HashSet::from_par_iter(
                pdb.par_residues().map(|a| a.name().unwrap().to_lowercase()),
            );
            let mut missing_residues = residue_set.difference(&pdb_set).peekable();

            ensure!(
                missing_residues.peek().is_none(),
                "No residue(s) found with identifier(s): {}",
                missing_residues.format(",")
            );

            // output_vec = pdb
            //     .residues()
            //     .filter(|x| {
            //         input_vec
            //             .iter()
            //             .any(|y| x.name().unwrap().to_lowercase() == y.to_lowercase())
            //     })
            //     .map(|r| r.serial_number())
            //     .collect();
            atom_vec = {
                match partial {
                    None => pdb
                        .atoms_with_hierarchy()
                        .filter(|a| residue_set.contains(&a.residue().name().unwrap().to_lowercase()))
                        .map(|a| a.atom().serial_number())
                        .collect(),
                    Some(p) => pdb
                        .atoms_with_hierarchy()
                        .filter(|a| match p {
                            Partial::Backbone => a.is_backbone(),
                            Partial::Sidechain => a.is_sidechain(),
                        } && residue_set.contains(a.residue().name().unwrap()))
                        .map(|a| a.atom().serial_number())
                        .collect(),
                }
            };
        }
    }

    Ok(atom_vec)
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
    fn parse_atomic_list_test() {
        let num_list = "1,2:5,7,9-11";
        let str_list = "OH,HH";
        let pdb = test_pdb("tests/test_blank.pdb");

        assert_eq!(
            parse_atomic_list(num_list, &pdb).unwrap(),
            vec!(1, 2, 3, 4, 5, 7, 9, 10, 11)
        );
        assert_eq!(parse_atomic_list(str_list, &pdb).unwrap(), vec![39, 40]);
    }

    #[test]
    fn parse_residue_list_test() {
        let num_list = "1,2:5,6-7";
        let str_list = "gly,wat";
        let pdb = test_pdb("tests/test_blank.pdb");

        assert_eq!(
            parse_atomic_list(num_list, &pdb).unwrap(),
            vec!(1, 2, 3, 4, 5, 6, 7)
        );
        assert_eq!(
            parse_residue_list(str_list, &pdb, None).unwrap(),
            vec![19, 20, 21, 22, 23, 24, 25, 78, 79, 80, 81, 82, 83]
        )
    }
}
