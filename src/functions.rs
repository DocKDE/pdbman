use crate::options::{Distance, Partial, Region, Target};
use crate::residue_ascii::RESIDUE_ASCII;

use std::collections::HashSet;
use std::io;
use std::io::prelude::Write;

use anyhow::{Context, Result};
use itertools::Itertools;
use lazy_regex::regex;
use pdbtbx::{
    Atom, AtomConformerResidueChainModel, ContainsAtomConformer, ContainsAtomConformerMut,
    ContainsAtomConformerResidue, ContainsAtomConformerResidueChain, PDB,
};
use prettytable::{format, Table};
use rayon::iter::FromParallelIterator;
use rayon::prelude::ParallelIterator;

type AtomList = Vec<usize>;
type ResidueList = Vec<isize>;

/// This function edits the q or b value of a given residue list (used by ORCA as input for active section
/// selection) by Residue.
pub fn edit_residues(
    pdb: &mut PDB,
    list: &[isize],
    mode: &str,
    partial: Option<Partial>,
    region: Region,
) {
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

    pdb.atoms_with_hierarchy_mut().for_each(|mut res| {
        if list.contains(&res.residue().serial_number()) {
            match partial {
                None => edit(res.atom_mut()),
                Some(Partial::Sidechain) => {
                    if !res.atom().is_backbone() {
                        edit(res.atom_mut())
                    }
                }
                Some(Partial::Backbone) => {
                    if res.atom().is_backbone() {
                        edit(res.atom_mut())
                    }
                }
            }
        }
    });
}

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

    pdb.par_atoms_mut().for_each(|a| {
        if list.contains(&a.serial_number()) {
            edit(a)
        }
    })
}

/// Removes a whole region from PDB file.
pub fn remove_region(pdb: &mut PDB, region: Option<Region>) {
    pdb.par_atoms_mut().for_each(|atom| match region {
        Some(Region::QM1) | Some(Region::QM2) => atom.set_occupancy(0.00).unwrap(),
        Some(Region::Active) => atom.set_b_factor(0.00).unwrap(),
        None => {
            atom.set_b_factor(0.00).unwrap();
            atom.set_occupancy(0.00).unwrap()
        }
    })
}

/// Prints all Atoms in Molecule to stdout in PDB file format
pub fn print_pdb_to_stdout(pdb: &PDB) -> Result<(), anyhow::Error> {
    let stdout = io::stdout();
    let mut handle = stdout.lock();
    for atom_hier in pdb.atoms_with_hierarchy() {
        writeln!(
            handle,
            "ATOM  {:>5} {:<4} {:>3} {:1}{:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}",
            atom_hier.atom().serial_number(),
            atom_hier.atom().name(),
            atom_hier.residue().name().unwrap_or(""),
            atom_hier.chain().id(),
            atom_hier.residue().serial_number(),
            atom_hier.atom().x(),
            atom_hier.atom().y(),
            atom_hier.atom().z(),
            atom_hier.atom().occupancy(),
            atom_hier.atom().b_factor(),
            atom_hier.atom().element()
        ).context("Failed to print PDB to stdout")?
    }
    Ok(())
}

/// Print all Atoms in PDB to a file given in path
// pub fn print_pdb_to_file(pdb: &PDB, path: &str) -> Result<(), std::io::Error> {
//     let mut file = BufWriter::new(File::create(path)?);

//     for atom_hier in pdb.atoms_with_hierarchy() {
//         // Deal with residue numbers above 9999
//         let mut resn = atom_hier.residue.serial_number();

//         if resn >= 10_000 {
//             resn -= (resn / 10_000) * 10_000;
//         };

//         // Deal with atom numbers above 99999
//         let mut atomn = atom_hier.atom.serial_number();
//         if atomn >= 100_000 {
//             atomn -= (atomn / 100_000) * 100_000
//         }

//         file.write_all(
//             format!(
//             "ATOM  {:>5} {:<4} {:>3} {:1}{:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}\n",
//             atomn,
//             atom_hier.atom.name(),
//             atom_hier.residue.name().unwrap_or(""),
//             atom_hier.chain.id(),
//             resn,
//             atom_hier.atom.x(),
//             atom_hier.atom.y(),
//             atom_hier.atom.z(),
//             atom_hier.atom.occupancy(),
//             atom_hier.atom.b_factor(),
//             atom_hier.atom.element(),
//         )
//             .as_bytes(),
//         )?;
//     }
//     Ok(())
// }

/// Takes an Atom struct as point of origin and a radius in A. Returns a Vector of Atom IDs
/// of Atoms within the given radius wrapped in a Result. Origin can be included or excluded.
pub fn calc_atom_sphere(
    pdb: &PDB,
    origin: &Atom,
    // origin: impl ContainsAtomConformerResidueChain,
    radius: f64,
    include_self: bool,
) -> Result<AtomList, anyhow::Error> {
    let tree = pdb.create_atom_rtree();
    let mut sphere_atoms: AtomList = tree
        .locate_within_distance(origin.pos(), radius.powf(2.0))
        .map(|atom| atom.serial_number())
        .collect();

    if !include_self {
        sphere_atoms.retain(|&x| x != origin.serial_number())
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
pub fn calc_residue_sphere(
    pdb: &PDB,
    origin: impl ContainsAtomConformerResidueChain,
    radius: f64,
    include_self: bool,
) -> Result<AtomList, anyhow::Error> {
    let tree = pdb.create_hierarchy_rtree();

    let mut sphere_atoms: AtomList = tree
        .locate_within_distance(origin.atom().pos(), radius.powf(2.0))
        .flat_map(|atom_hier| atom_hier.residue().atoms().map(|atom| atom.serial_number()))
        .unique()
        .collect();

    let origin_res_atoms: AtomList = origin
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

    if !table.is_empty() {
        if level == Distance::Clashes {
            writeln!(io::stdout(), "\nClash Analysis")
                .context("Failed to print clash analysis to stdout.")?;
        } else if level == Distance::Contacts {
            writeln!(io::stdout(), "\nContact Analysis")
                .context("Failed to print contact analysis to stdout.")?;
        }
        Ok(table)
    } else {
        bail!("No contacts found!")
    }
}

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

            let input_set: HashSet<String> = input_vec
                .iter()
                .map(|s| s.to_lowercase())
                .collect();
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
pub fn parse_residue_list(input: &str, pdb: &PDB) -> Result<ResidueList, anyhow::Error> {
    let mut output_vec: ResidueList = vec![];

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
                    output_vec.extend(range);
                } else {
                    output_vec.push(i.parse().unwrap())
                }
            }
            let output_set: HashSet<isize> = HashSet::from_iter(output_vec.iter().copied());
            let pdb_set: HashSet<isize> =
                HashSet::from_par_iter(pdb.par_residues().map(|a| a.serial_number()));
            let mut missing_residues = output_set.difference(&pdb_set).peekable();

            ensure!(
                missing_residues.peek().is_none(),
                "No residue(s) found with serial number(s): {}",
                missing_residues.format(",")
            );
        }
        false => {
            let input_vec: Vec<&str> = input.split(',').collect();

            // This might not be necessary
            ensure!(
                !regex!(r"[:-]").is_match(input),
                "Ranges are not allowed for residue name inputs."
            );

            let input_set: HashSet<String> = input_vec
                .iter()
                .map(|s| s.to_lowercase())
                .collect();
            let pdb_set: HashSet<String> = HashSet::from_par_iter(
                pdb.par_residues().map(|a| a.name().unwrap().to_lowercase()),
            );
            let mut missing_residues = input_set.difference(&pdb_set).peekable();

            ensure!(
                missing_residues.peek().is_none(),
                "No residue(s) found with identifier(s): {}",
                missing_residues.format(",")
            );

            output_vec = pdb
                .residues()
                .filter(|x| {
                    input_vec
                        .iter()
                        .any(|y| x.name().unwrap().to_lowercase() == y.to_lowercase())
                })
                .map(|r| r.serial_number())
                .collect();
        }
    }

    Ok(output_vec)
}

/// Query Molecule for information. Depending on the input this will print a table of
/// Residues and/or Atoms with available information that were asked for.
pub fn query_atoms(pdb: &PDB, atom_list: &[usize]) -> Result<(), anyhow::Error> {
    let mut table = Table::new();
    table.add_row(row![
        "Atom ID",
        "Atom name",
        "Residue ID",
        "Residue Name",
        "QM",
        "Active"
    ]);

    let mut key: Option<&str> = None;
    let mut resname_vec = Vec::new();

    for atom_hier in pdb.atoms_with_hierarchy() {
        if atom_list.contains(&atom_hier.atom().serial_number()) {
            resname_vec.push(atom_hier.residue().name().map(|s| s.to_owned()));

            table.add_row(row![
                atom_hier.atom().serial_number(),
                atom_hier.atom().name(),
                atom_hier.residue().serial_number().to_string()
                    + atom_hier.residue().insertion_code().unwrap_or(""),
                atom_hier.residue().name().unwrap_or(""),
                atom_hier.atom().occupancy(),
                atom_hier.atom().b_factor(),
            ]);
        }
    }

    if !resname_vec.is_empty() && resname_vec.iter().all_equal() {
        key = resname_vec[0].as_deref()
    }

    if let Some(k) = key {
        if let Some(res_ascii) = RESIDUE_ASCII.get(&k.to_uppercase().as_ref()) {
            writeln!(io::stdout(), "{}", res_ascii)
                .context("Failed to print residue depiction to stdout")?
        }
    }

    table.printstd();
    Ok(())
}

pub fn get_atomlist(pdb: &PDB, region: Region) -> Result<Vec<String>, anyhow::Error> {
    let filt_closure = match region {
        Region::QM1 => |a: &Atom| a.occupancy() == 1.00,
        Region::QM2 => |a: &Atom| a.occupancy() == 2.00,
        Region::Active => |a: &Atom| a.b_factor() == 1.00,
    };

    let str_vec = pdb
        .par_atoms()
        .filter(|&a| filt_closure(a))
        .map(|a| a.serial_number().to_string())
        .collect::<Vec<String>>();

    ensure!(!str_vec.is_empty(), "No atoms in the requested region!");
    Ok(str_vec)
}

// This cannot fail because if no residues can be queried, the
// parse_residue_list would have returned an error beforehand
pub fn query_residues(pdb: &PDB, residue_list: &[isize]) -> Result<(), anyhow::Error> {
    let mut table = Table::new();

    table.add_row(row![
        "Atom ID",
        "Atom name",
        "Residue ID",
        "Residue Name",
        "QM",
        "Active"
    ]);

    // let mut key = None;
    let mut key = None;

    for atom_hier in pdb.atoms_with_hierarchy() {
        // let res = (
        //     atom_hier.residue.serial_number(),
        //     atom_hier.residue.insertion_code().map(ToOwned::to_owned),
        // );

        if residue_list.contains(&atom_hier.residue().serial_number()) {
            if residue_list.len() == 1 {
                // key = atom_hier.residue().name();
                key = atom_hier.residue().name().map(|s| s.to_owned());
            }
            table.add_row(row![
                atom_hier.atom().serial_number(),
                atom_hier.atom().name(),
                atom_hier.residue().serial_number().to_string()
                    + atom_hier.residue().insertion_code().unwrap_or(""),
                atom_hier.residue().name().unwrap_or(""),
                atom_hier.atom().occupancy(),
                atom_hier.atom().b_factor(),
            ]);
        }
    }

    if let Some(k) = key {
        if let Some(res_ascii) = RESIDUE_ASCII.get(&k.to_uppercase().as_ref()) {
            writeln!(io::stdout(), "{}", res_ascii)
                .context("Failed to print residue depiction to stdout")?
        }
    }

    table.printstd();
    Ok(())
}

pub fn get_residuelist(pdb: &PDB, region: Region) -> Result<Vec<String>, anyhow::Error> {
    let filt_closure = match region {
        Region::QM1 => |a: &AtomConformerResidueChainModel| a.atom().occupancy() == 1.00,
        Region::QM2 => |a: &AtomConformerResidueChainModel| a.atom().occupancy() == 2.00,
        Region::Active => |a: &AtomConformerResidueChainModel| a.atom().b_factor() == 1.00,
    };

    let str_vec = pdb
        .atoms_with_hierarchy()
        .filter(|a| filt_closure(a))
        .map(|a| a.residue().serial_number().to_string())
        .dedup()
        .collect::<Vec<String>>();

    ensure!(!str_vec.is_empty(), "No residues in the requested region!");
    Ok(str_vec)
}

pub fn analyze(
    pdb: &PDB,
    region: Option<Region>,
    target: Option<Target>,
) -> Result<(), anyhow::Error> {
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

    basic_table.printstd();

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
            writeln!(io::stdout(), "\n{} Residues", region.unwrap().to_owned())
                .context("Failed to print residues to stdout")?;
            residue_table.printstd();
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
            writeln!(io::stdout(), "\n{} Atoms", region.unwrap().to_owned())
                .context("Failed to print atoms to stdout")?;
            atom_table.printstd();
        } else {
            bail!("No Atoms found in given region!")
        }
    }
    Ok(())
}

// pub fn check_residue_overflow(pdb: &PDB) -> bool {
//     pdb.par_residues().count() > 9999
//         && pdb
//             .par_residues()
//             .map(|r| r.insertion_code())
//             .any(|i| i.is_none())
// }

// pub fn add_insertion_codes(pdb: &mut PDB) {
//     let insertion_codes = ('A'..='Z').cycle();
//     for i in pdb
//         .residues_mut()
//         .chunks(9999)
//         .into_iter()
//         .zip(insertion_codes)
//     {
//         for res in i.0 {
//             res.set_insertion_code(&i.1.to_string());
//         }
//     }
// }

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
    fn parse_atomic_list_test() {
        let num_list = "1,2:5,7,9-11";
        let str_list = "OH,HH";
        let pdb = test_pdb("tests/test_blank.pdb");

        assert_eq!(
            parse_atomic_list(num_list, &pdb).unwrap(),
            vec!(1, 2, 3, 4, 5, 7, 9, 10, 11)
        );
        assert_eq!(parse_atomic_list(str_list, &pdb).unwrap(), vec!(39, 40));
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
        assert_eq!(parse_residue_list(str_list, &pdb).unwrap(), vec![2, 6, 7])
    }

    #[test]
    fn atom_sphere_test() {
        let pdb = test_pdb("tests/test_blank.pdb");
        let origin = pdb
            .atoms_with_hierarchy()
            .find(|x| x.atom().serial_number() == 26)
            .unwrap();
        let atom_list_incl = vec![21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 45, 46];
        let atom_list_excl = vec![21, 22, 23, 24, 25, 27, 28, 29, 30, 32, 45, 46];

        assert_eq!(
            calc_atom_sphere(&pdb, origin.atom(), 4.0, true)
                .unwrap()
                .len(),
            22
        );
        assert_eq!(
            calc_atom_sphere(&pdb, origin.atom(), 3.0, true).unwrap(),
            atom_list_incl
        );

        assert_eq!(
            calc_atom_sphere(&pdb, origin.atom(), 4.0, false)
                .unwrap()
                .len(),
            21
        );
        assert_eq!(
            calc_atom_sphere(&pdb, origin.atom(), 3.0, false).unwrap(),
            atom_list_excl
        );
    }

    #[test]
    fn residue_sphere_test() {
        let pdb = test_pdb("tests/test_blank.pdb");
        let origin = pdb
            .atoms_with_hierarchy()
            .find(|x| x.atom().serial_number() == 26)
            .unwrap();
        let atom_list_excl = vec![
            19, 20, 21, 22, 23, 24, 25, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61,
            62,
        ];

        assert_eq!(
            calc_residue_sphere(&pdb, origin.clone(), 4.0, true)
                .unwrap()
                .len(),
            44
        );
        assert_eq!(
            calc_residue_sphere(&pdb, origin.clone(), 4.0, false)
                .unwrap()
                .len(),
            23
        );
        assert_eq!(
            calc_residue_sphere(&pdb, origin.clone(), 4.0, false).unwrap(),
            atom_list_excl
        );
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
    fn edit_residues_test_all() {
        let mut pdb = test_pdb("tests/test_blank.pdb");
        let res_id_list = &vec![2, 4];

        edit_residues(&mut pdb, res_id_list, "Add", None, Region::QM2);
        edit_residues(&mut pdb, res_id_list, "Add", None, Region::Active);

        let res_list = pdb
            .residues()
            .filter(|x| res_id_list.contains(&x.serial_number()));

        for residue in res_list {
            for atom in residue.atoms() {
                assert_eq!(atom.occupancy(), 2.00);
                assert_eq!(atom.b_factor(), 1.00);
            }
        }

        edit_residues(&mut pdb, res_id_list, "Remove", None, Region::QM2);
        edit_residues(&mut pdb, res_id_list, "Remove", None, Region::Active);

        let res_list = pdb
            .residues()
            .filter(|x| res_id_list.contains(&x.serial_number()));

        for residue in res_list {
            for atom in residue.atoms() {
                assert_eq!(atom.occupancy(), 0.00);
                assert_eq!(atom.b_factor(), 0.00);
            }
        }
    }

    #[test]
    fn edit_residues_test_side() {
        let mut pdb = test_pdb("tests/test_blank.pdb");
        let res_id_list = &vec![2, 4];

        edit_residues(
            &mut pdb,
            res_id_list,
            "Add",
            Some(Partial::Sidechain),
            Region::QM2,
        );
        edit_residues(
            &mut pdb,
            res_id_list,
            "Add",
            Some(Partial::Sidechain),
            Region::Active,
        );

        let res_list = pdb
            .residues()
            .filter(|x| res_id_list.contains(&x.serial_number()));

        let atom_list = res_list
            .flat_map(|x| x.atoms())
            .filter(|x| !x.is_backbone());

        for atom in atom_list {
            assert_eq!(atom.occupancy(), 2.00);
            assert_eq!(atom.b_factor(), 1.00);
        }
    }

    #[test]
    fn edit_residues_test_back() {
        let mut pdb = test_pdb("tests/test_blank.pdb");
        let res_id_list = &vec![2, 4];
        edit_residues(
            &mut pdb,
            res_id_list,
            "Add",
            Some(Partial::Backbone),
            Region::QM2,
        );
        edit_residues(
            &mut pdb,
            res_id_list,
            "Add",
            Some(Partial::Backbone),
            Region::Active,
        );

        let res_list = pdb
            .residues()
            .filter(|x| res_id_list.contains(&x.serial_number()));
        let atom_list = res_list.flat_map(|x| x.atoms()).filter(|x| x.is_backbone());

        for atom in atom_list {
            assert_eq!(atom.occupancy(), 2.00);
            assert_eq!(atom.b_factor(), 1.00);
        }
    }

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

    #[test]
    fn get_atomlist_test() {
        let pdb = test_pdb("tests/test_get_atomlist.pdb");
        let qm1_atoms = get_atomlist(&pdb, Region::QM1).unwrap();
        let qm2_atoms = get_atomlist(&pdb, Region::QM2).unwrap();
        let active_atoms = get_atomlist(&pdb, Region::Active).unwrap();

        assert_eq!(qm1_atoms, vec!["1", "2", "4", "5", "6"]);
        assert_eq!(qm2_atoms, vec!["8", "9", "11", "12"]);
        assert_eq!(active_atoms, vec!["1", "2", "3", "5", "6", "8", "9"]);
    }
}
