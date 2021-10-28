use crate::options::{Distance, Partial, Region, Target};
use crate::residue_ascii::RESIDUE_ASCII;

use std::fs::File;
use std::io::prelude::Write;
use std::io::LineWriter;

use anyhow::Result;
use itertools::Itertools;
use lazy_regex::regex;
use pdbtbx::{Atom, AtomWithHierarchy, PDB};
use prettytable::{format, Table};
use rayon::prelude::ParallelIterator;

type AtomList = Vec<usize>;
type ResidueList = Vec<(isize, Option<String>)>;

#[derive(PartialEq, Debug, Clone)]
pub struct Sphere<'a> {
    pub origin: pdbtbx::AtomWithHierarchy<'a>,
    pub radius: f64,
}

/// This function parses the values of option 'Sphere' into a usize and an f64 and returns
/// a reference to the corresponding Atom and the radius. The return values are organized
/// in a Sphere struct to facilitate usage.
impl<'a> Sphere<'a> {
    pub fn new(
        origin_str: &str,
        radius_str: &str,
        pdb: &'a PDB,
    ) -> Result<Sphere<'a>, anyhow::Error> {
        // The sphere validator rejects anything containing something else than digits and dots.
        // Since it validates both arguments the same, it will not reject a decimal point in the Atom
        // ID which is taken care of here.
        let origin_id: usize = match origin_str.parse() {
            Ok(v) => v,
            Err(_) => bail!("Decimal points are not allowed in Atom IDs."),
        };

        // No more error handling should be necessary since the sphere validator already
        // takes care of everything that is not a float. Everything that can be constructed
        // with dots and digits is valid input for radius.
        let radius: f64 = radius_str.parse()?;

        let origin = pdb
            .par_atoms_with_hierarchy()
            .find_any(|x| x.atom.serial_number() == origin_id)
            .ok_or_else(|| anyhow!("No atom corresponding to the given ID could be found."))?;

        Ok(Sphere { origin, radius })
    }
}

/// This function edits the q value of the Molecule (used by ORCA as input for QM section
/// selection) by Residue. It takes a mutable reference of a PDB struct, the desired value q should
/// be set to, and a list of Residue serial numbers for the Atoms to edit.
/// The last option recognizes 'backbone', 'sidechain' and 'None' as options and will select the
/// correspdonding Atoms from each Residue.
pub fn edit_qm_residues(
    pdb: &mut PDB,
    qm_val: f64,
    list: ResidueList,
    partial: Partial,
) -> Result<(), String> {
    pdb.par_residues_mut()
        .try_for_each(|res| -> Result<(), String> {
            if list.contains(&(
                res.serial_number(),
                res.insertion_code().map(ToOwned::to_owned),
            )) {
                res.par_atoms_mut()
                    .try_for_each(|atom| -> Result<(), String> {
                        match partial {
                            Partial::None => atom.set_occupancy(qm_val)?,
                            Partial::Sidechain => {
                                if !atom.is_backbone() {
                                    atom.set_occupancy(qm_val)?
                                }
                            }
                            Partial::Backbone => {
                                if atom.is_backbone() {
                                    atom.set_occupancy(qm_val)?
                                }
                            }
                        }
                        Ok(())
                    })?;
            };
            Ok(())
        })?;
    Ok(())
}

/// This function edits the b value of the Molecule (used by ORCA as input for active section
/// selection) by Residue. It takes a mutable reference of a PDB struct, the desired value b should
/// be set to, and a list of Residue IDs for the Residues to edit.
/// The last option recognizes 'backbone', 'sidechain' and 'None' as options and will select the
/// correspdonding Atoms from each Residue.
pub fn edit_active_residues(
    pdb: &mut PDB,
    active_val: f64,
    list: ResidueList,
    partial: Partial,
) -> Result<(), String> {
    pdb.par_residues_mut()
        .try_for_each(|res| -> Result<(), String> {
            if list.contains(&(
                res.serial_number(),
                res.insertion_code().map(ToOwned::to_owned),
            )) {
                // if list.contains(&(res.id())) {
                res.par_atoms_mut()
                    .try_for_each(|atom| -> Result<(), String> {
                        match partial {
                            Partial::None => atom.set_b_factor(active_val)?,
                            Partial::Sidechain => {
                                if !atom.is_backbone() {
                                    atom.set_b_factor(active_val)?
                                }
                            }
                            Partial::Backbone => {
                                if atom.is_backbone() {
                                    atom.set_b_factor(active_val)?
                                }
                            }
                        }
                        Ok(())
                    })?;
            };
            Ok(())
        })?;

    Ok(())
}

/// This functions edits the q value of the PDB file (used by ORCA as input for QM region
/// selection) by Atom. It takes a mutable reference to a PDB struct, the desired value
/// q should be set to and a vector of Atom IDs.
pub fn edit_qm_atoms(pdb: &mut PDB, qm_val: f64, list: AtomList) -> Result<(), String> {
    pdb.par_atoms_mut()
        .try_for_each(|atom: &mut Atom| -> Result<(), String> {
            if list.contains(&atom.serial_number()) {
                atom.set_occupancy(qm_val)?;
            };
            Ok(())
        })?;
    Ok(())
}

/// This functions edits the b value of the PDB file (used by ORCA as input for QM region
/// selection) by Atom. It takes a mutable reference to a PDB struct, the desired value
/// b should be set to and a vector of Atom IDs.
pub fn edit_active_atoms(pdb: &mut PDB, active_val: f64, list: AtomList) -> Result<(), String> {
    pdb.par_atoms_mut()
        .try_for_each(|atom: &mut Atom| -> Result<(), String> {
            if list.contains(&atom.serial_number()) {
                atom.set_b_factor(active_val)?;
            };
            Ok(())
        })?;
    Ok(())
}

/// Sets all q and b values to zero which serves as a fresh start.
pub fn remove_all(pdb: &mut PDB) -> Result<(), String> {
    pdb.par_atoms_mut()
        .try_for_each(|atom| -> Result<(), String> {
            atom.set_occupancy(0.00)?;
            atom.set_b_factor(0.00)?;
            Ok(())
        })?;
    Ok(())
}

/// Sets all q values to zero which serves as a reset of the QM region.
pub fn remove_qm(pdb: &mut PDB) -> Result<(), String> {
    pdb.par_atoms_mut()
        .try_for_each(|atom| -> Result<(), String> {
            atom.set_occupancy(0.00)?;
            Ok(())
        })?;
    Ok(())
}

/// Sets all b values to zero which serves as a reset of the active region.
pub fn remove_active(pdb: &mut PDB) -> Result<(), String> {
    pdb.par_atoms_mut()
        .try_for_each(|atom| -> Result<(), String> {
            atom.set_b_factor(0.00)?;
            Ok(())
        })?;
    Ok(())
}

/// Prints all Atoms in Molecule to stdout in PDB file format
pub fn print_pdb_to_stdout(pdb: &PDB) -> Result<(), String> {
    for atom_hier in pdb.atoms_with_hierarchy() {
        println!(
            "ATOM  {:>5} {:<4} {:>3}  {:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}",
            atom_hier.atom.serial_number(),
            atom_hier.atom.name(),
            atom_hier.residue.name().ok_or("No Residue Name")?,
            atom_hier.residue.serial_number(),
            atom_hier.atom.x(),
            atom_hier.atom.y(),
            atom_hier.atom.z(),
            atom_hier.atom.occupancy(),
            atom_hier.atom.b_factor(),
            atom_hier.atom.element()
        );
    }
    Ok(())
}

/// Print all Atoms in PDB to a file given in path
pub fn print_pdb_to_file(pdb: &PDB, path: &str) -> Result<(), anyhow::Error> {
    let file = File::create(path)?;
    let mut file = LineWriter::new(file);

    for atom_hier in pdb.atoms_with_hierarchy() {
        // Deal with residue numbers above 9999
        let mut resn = atom_hier.residue.serial_number();

        if resn >= 10_000 {
            resn -= (resn / 10_000)  * 10_000;
        };

        // Deal with atom numbers above 99999
        let mut atomn = atom_hier.atom.serial_number();
        if atomn >= 100_000 {
            atomn -= (atomn / 100_000) * 100_000
        }

        file.write_all(
            format!(
            "ATOM  {:>5} {:<4} {:>3}  {:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}\n",
            // atom_hier.atom.serial_number(),
            atomn,
            atom_hier.atom.name(),
            atom_hier.residue.name().ok_or("No Residue Name").unwrap(),
            // atom_hier.residue.serial_number(),
            resn,
            atom_hier.atom.x(),
            atom_hier.atom.y(),
            atom_hier.atom.z(),
            atom_hier.atom.occupancy(),
            atom_hier.atom.b_factor(),
            atom_hier.atom.element(),
        )
            .as_bytes(),
        )?;
    }
    Ok(())
}

/// Takes an Atom struct as point of origin and a radius in A. Returns a Vector of Atom IDs
/// of Atoms within the given radius wrapped in a Result. Origin can be included or excluded.
pub fn calc_atom_sphere(
    pdb: &PDB,
    origin: &AtomWithHierarchy,
    radius: f64,
    include_self: bool,
) -> Result<AtomList, String> {
    let tree = pdb.create_atom_rtree();
    let mut sphere_atoms: AtomList = tree
        .locate_within_distance(origin.atom.pos_array(), radius.powf(2.0))
        .map(|atom| atom.serial_number())
        .collect();

    if !include_self {
        sphere_atoms.retain(|&x| x != origin.atom.serial_number())
    }

    sphere_atoms.sort_unstable();

    if !sphere_atoms.is_empty() {
        Ok(sphere_atoms)
    } else {
        Err("Calculated sphere doesn't contain any atoms.".into())
    }
}

/// Takes an Atom struct as point of origin, a radius in A. Returns a Vector of Atom IDs
/// that belong to a Residue that had at least one Atom within the given radius wrapped in a Result.
/// Origin residue can be included or excluded.
pub fn calc_residue_sphere(
    pdb: &PDB,
    origin: &AtomWithHierarchy,
    radius: f64,
    include_self: bool,
) -> Result<AtomList, String> {
    let tree = pdb.create_atom_with_hierarchy_rtree();

    let mut sphere_atoms: AtomList = tree
        .locate_within_distance(origin.atom.pos_array(), radius.powf(2.0))
        .flat_map(|atom_hier| atom_hier.residue.atoms().map(|atom| atom.serial_number()))
        .sorted()
        .dedup()
        .collect();

    let origin_res_atoms: AtomList = origin
        .residue
        .atoms()
        .map(|atom| atom.serial_number())
        .collect();

    if !include_self {
        sphere_atoms.retain(|x| !origin_res_atoms.contains(x))
    }

    if !sphere_atoms.is_empty() {
        Ok(sphere_atoms)
    } else {
        Err("Calculated sphere doesn't contain any atoms.".into())
    }
}

/// Finds and prints all contacts present in the PDB file structure. Definition of
/// 'contact' is given by the 'level' arg which is 1.0A for 'level' = 0 and the
/// respective atomic radius for 'level' = 1.
pub fn find_contacts(pdb: &PDB, level: Distance) -> Result<Table, String> {
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

    let tree = pdb.create_atom_with_hierarchy_rtree();

    for atom_hier in pdb.atoms_with_hierarchy() {
        let radius: f64 = match level {
            Distance::Clashes => 1.0,
            Distance::Contacts => atom_hier
                .atom
                .atomic_radius()
                .ok_or(format!(
                    "No radius found for given atom: {}",
                    atom_hier.atom
                ))?
                .powf(2.0),
            _ => return Err("Invalid".into()),
        };

        let contacts = tree.locate_within_distance(atom_hier.atom.pos_array(), radius);

        for other_atom_hier in contacts {
            // This eliminates duplicate entries
            if other_atom_hier.atom < atom_hier.atom
            // This eliminates atoms from same residue
                && other_atom_hier.residue != atom_hier.residue
                // This eliminates neighboring residues
                && !(other_atom_hier.atom.name() == "C"
                    && atom_hier.atom.name() == "N"
                    && other_atom_hier.residue.serial_number() + 1
                        == atom_hier.residue.serial_number())
            {
                let distance = other_atom_hier.atom.distance(atom_hier.atom);

                if distance <= 0.75 && distance > 0.5 {
                    table.add_row(row![
                        bByFd =>
                        other_atom_hier.atom.serial_number(),
                        other_atom_hier.atom.name(),
                        other_atom_hier.residue.name().ok_or("No Residue Name.")?,
                        atom_hier.atom.serial_number(),
                        atom_hier.atom.name(),
                        atom_hier.residue.name().ok_or("No Residue Name")?,
                        format!("{:.2}", distance)
                    ]);
                } else if distance <= 0.5 {
                    table.add_row(row![
                        bBrFd =>
                        other_atom_hier.atom.serial_number(),
                        other_atom_hier.atom.name(),
                        other_atom_hier.residue.name().ok_or("No Residue Name.")?,
                        atom_hier.atom.serial_number(),
                        atom_hier.atom.name(),
                        atom_hier.residue.name().ok_or("No Residue Name")?,
                        format!("{:.2}", distance)
                    ]);
                } else {
                    table.add_row(row![
                        other_atom_hier.atom.serial_number(),
                        other_atom_hier.atom.name(),
                        other_atom_hier.residue.name().ok_or("No Residue Name.")?,
                        atom_hier.atom.serial_number(),
                        atom_hier.atom.name(),
                        atom_hier.residue.name().ok_or("No Residue Name")?,
                        format!("{:.2}", distance)
                    ]);
                }
            }
        }
    }

    if !table.is_empty() {
        if level == Distance::Clashes {
            println!("\nClash Analysis");
        } else if level == Distance::Contacts {
            println!("\nContact Analysis");
        }
        Ok(table)
    } else {
        Err("No contacts found!".into())
    }
}

/// Takes a comma-separated list (usually from command line input) as string and parses it into
/// a vector of Atom IDs. The Input may be atom IDs or Atom Names
pub fn parse_atomic_list(input: &str, pdb: &PDB) -> Result<AtomList, anyhow::Error> {
    let mut output_vec: AtomList = vec![];
    let re = regex!(r"(?P<num>\d+)?(?P<str>[A-Za-z])?");

    // This test catches list inputs that contain digits followd by a letter.
    // This is valid input for residues but not atoms but since both get their
    // input from the list input validator this possibility is accounted for here.
    let mut invalid_chars = input
        .split(&[',', '-', ':'][..])
        .filter(|x| {
            let caps = re.captures(x).unwrap();
            caps.name("num").is_some() && caps.name("str").is_some()
        })
        .peekable();

    ensure!(
        invalid_chars.peek().is_none(),
        "Letters are not allowed in atom list input: {}",
        invalid_chars.join(",")
    );

    match input
        .split(&[',', '-', ':'][..])
        .all(|x| x.parse::<usize>().is_ok())
    {
        true => {
            let input_vec = input.split(',');

            for i in input_vec {
                if i.contains(&['-', ':'][..]) {
                    let vector: AtomList =
                        i.split(&['-', ':'][..]).map(|x| x.parse()).try_collect()?;

                    output_vec.extend(&mut (vector[0]..vector[1] + 1));
                } else {
                    output_vec.push(i.parse()?)
                }
            }

            let mut missing_atoms = output_vec
                .iter()
                .filter(|&x| !pdb.atoms().any(|atom| &atom.serial_number() == x))
                .sorted()
                .dedup()
                .peekable();

            ensure!(
                missing_atoms.peek().is_none(),
                "No atom(s) found with serial number(s): {}",
                missing_atoms.format(",")
            );
        }
        false => {
            let input_vec: Vec<&str> = input.split(',').collect();
            // let mut input_vec = input.split(',');

            let mut missing_atoms = input_vec
                .iter()
                .filter(|&x| {
                    !pdb.atoms()
                        .any(|atom| atom.name().to_lowercase() == x.to_lowercase())
                })
                .sorted()
                .dedup()
                .peekable();

            ensure!(
                missing_atoms.peek().is_none(),
                "No atom(s) found with serial number(s): {}",
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
pub fn parse_residue_list<'a>(input: &'a str, pdb: &'a PDB) -> Result<ResidueList> {
    // ) -> Result<ResidueList<'a>, anyhow::Error> {
    let re_num =
        regex!(r"^(?P<id1>\d+)(?P<insert1>[A-Za-z]?)([:-](?P<id2>\d+)(?P<insert2>[A-Za-z]?))?$");

    let mut output_vec: ResidueList = vec![];
    let mut err_vec: Vec<String> = vec![];

    match input.split(',').all(|x| re_num.is_match(x)) {
        true => {
            let input_vec: Vec<&str> = input.split(',').collect();

            for i in input_vec {
                if i.contains(&[':', '-'][..]) {
                    let re_num = regex!(
                        r"^(?P<id1>\d+)(?P<insert1>[A-Za-z]?)[:-](?P<id2>\d+)(?P<insert2>[A-Za-z]?)$"
                    );
                    let caps = re_num.captures(i).unwrap();

                    let id1: isize = caps.name("id1").unwrap().as_str().parse()?;
                    let insert1 = match caps.name("insert1").map(|x| x.as_str()) {
                        Some("") => None,
                        Some(x) => Some(x),
                        _ => unreachable!(),
                    };

                    let id2: isize = caps.name("id2").unwrap().as_str().parse()?;
                    let insert2 = match caps.name("insert2").map(|x| x.as_str()) {
                        Some("") => None,
                        Some(x) => Some(x),
                        _ => unreachable!(),
                    };

                    let mut start = None;
                    let mut end = None;
                    let mut parse_res = false;

                    // Determine start and end indices while simultaneously populating output_vec.
                    // This way the PDB is only traversed once.
                    for (index, residue) in pdb.residues().enumerate() {
                        if residue.id() == (id1, insert1) {
                            start = Some(index);
                            parse_res = true;
                        }

                        if parse_res {
                            output_vec.push((
                                residue.serial_number(),
                                residue.insertion_code().map(ToOwned::to_owned),
                            ))
                        }

                        if residue.id() == (id2, insert2) {
                            // Return prematurely if start is None which means that end < start which is invalid.
                            ensure!(start.is_some(), "Invalid range given: {}{}-{}{}. Left entry must preceed right one in PDB file!", 
                                id1, insert1.unwrap_or(""), id2, insert2.unwrap_or(""));
                            end = Some(index);
                            parse_res = false;
                        }
                    }

                    if start.is_none() {
                        err_vec.push(format!("{}{}", id1, insert1.unwrap_or("")));
                    }
                    if end.is_none() {
                        err_vec.push(format!("{}{}", id2, insert2.unwrap_or("")));
                    }
                } else {
                    let re = regex!(r"^(?P<resid>\d+)(?P<insert>[A-Za-z])?$");
                    let caps = re.captures(i).unwrap();
                    let resid: isize = caps.name("resid").unwrap().as_str().parse()?;
                    let insert = caps.name("insert").map(|x| x.as_str());

                    if !pdb.par_residues().any(|x| x.id() == (resid, insert)) {
                        err_vec.push(format!("{}{}", resid, insert.unwrap_or("")));
                    } else {
                        output_vec.push((resid, insert.map(ToOwned::to_owned)))
                    }
                }
            }
        }
        false => {
            let input_vec: Vec<String> = input.split(',').map(|x| x.to_lowercase()).collect();

            for i in &input_vec {
                if !pdb
                    .par_residues()
                    .any(|x| &x.name().unwrap_or("").to_lowercase() == i)
                {
                    err_vec.push(i.to_owned());
                }
            }

            output_vec = pdb
                .residues()
                .map(|x| {
                    Ok(input_vec
                        .contains(
                            &x.name()
                                .ok_or_else(|| anyhow!("No Residue Name"))?
                                .to_lowercase(),
                        )
                        .then(|| (x.serial_number(), x.insertion_code().map(ToOwned::to_owned))))
                })
                .filter_map(Result::transpose)
                .collect::<Result<Vec<(isize, Option<String>)>, anyhow::Error>>()?;
        }
    };

    err_vec.sort();
    err_vec.dedup();

    ensure!(
        err_vec.is_empty(),
        "No residue(s) found with identifier(s): {}",
        err_vec.join(",")
    );

    Ok(output_vec)
}

/// Query Molecule for information. Depending on the input this will print a table of
/// Residues and/or Atoms will available information that were asked for.
pub fn query_atoms(pdb: &PDB, atom_list: AtomList) -> Result<(), String> {
    let mut table = Table::new();
    table.add_row(row![
        "Atom ID",
        "Atom name",
        "Residue ID",
        "Residue Name",
        "QM",
        "Active"
    ]);

    for atom_hier in pdb.atoms_with_hierarchy() {
        if atom_list.contains(&atom_hier.atom.serial_number()) {
            table.add_row(row![
                atom_hier.atom.serial_number(),
                atom_hier.atom.name(),
                atom_hier.residue.serial_number().to_string()
                    + atom_hier.residue.insertion_code().unwrap_or(""),
                atom_hier.residue.name().ok_or("No Residue Name")?,
                atom_hier.atom.occupancy(),
                atom_hier.atom.b_factor(),
            ]);
        }
    }

    table.printstd();
    Ok(())
}

pub fn get_atomlist(pdb: &PDB, region: Region) -> Result<String, anyhow::Error> {
    let filt_value = match region {
        Region::QM1 => 1.00,
        Region::QM2 => 2.00,
        Region::Active => 1.00,
        _ => bail!("Invalid input for 'region' argument."),
    };

    let str_vec = pdb
        .par_atoms()
        .filter(|a| match region {
            Region::QM1 | Region::QM2 => a.occupancy() == filt_value,
            Region::Active => a.b_factor() == filt_value,
            Region::None => false,
        })
        .map(|a| a.serial_number().to_string())
        .collect::<Vec<String>>();

    if str_vec.is_empty() {
        bail!("No atoms in the requested region!")
    } else {
        Ok(str_vec.join(","))
    }
}

// pub fn get_active_atomlist(pdb: &PDB) -> Result<String, anyhow::Error> {
//     let str_vec = pdb
//         .par_atoms()
//         .filter(|a| a.b_factor() == 1.00)
//         .map(|a| a.serial_number().to_string()).collect::<Vec<String>>();

//     if str_vec.is_empty() {
//         bail!("No atoms in the requested region!")
//     } else {
//         Ok(str_vec.join(","))
//     }
// }

pub fn query_residues(pdb: &PDB, residue_list: ResidueList) -> Result<(), String> {
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

    for atom_hier in pdb.atoms_with_hierarchy() {
        let res = (
            atom_hier.residue.serial_number(),
            atom_hier.residue.insertion_code().map(ToOwned::to_owned),
        );

        if residue_list.contains(&res) {
            if residue_list.len() == 1 {
                key = atom_hier.residue.name();
            }
            table.add_row(row![
                atom_hier.atom.serial_number(),
                atom_hier.atom.name(),
                atom_hier.residue.serial_number().to_string()
                    + atom_hier.residue.insertion_code().unwrap_or(""),
                atom_hier.residue.name().ok_or("No Residue name")?,
                atom_hier.atom.occupancy(),
                atom_hier.atom.b_factor(),
            ]);
        }
    }

    if let Some(k) = key {
        println!(
            "{}",
            RESIDUE_ASCII.get(&k.to_uppercase().as_ref()).unwrap_or(&"")
        );
    }

    table.printstd();
    Ok(())
}

pub fn analyze(pdb: &PDB, region: Region, target: Target) -> Result<(), anyhow::Error> {
    let mut qm1_residue_list = Vec::new();
    let mut qm1_atom_list = Vec::new();
    let mut qm2_residue_list = Vec::new();
    let mut qm2_atom_list = Vec::new();
    let mut active_residue_list = Vec::new();
    let mut active_atom_list = Vec::new();

    for residue in pdb.residues() {
        for atom in residue.atoms() {
            if atom.occupancy() == 1.00 {
                if !qm1_residue_list.contains(&residue) {
                    qm1_residue_list.push(residue);
                }
                qm1_atom_list.push(atom)
            }

            if atom.occupancy() == 2.00 {
                if !qm2_residue_list.contains(&residue) {
                    qm2_residue_list.push(residue);
                }
                qm2_atom_list.push(atom)
            }

            if atom.b_factor() == 1.00 {
                if !active_residue_list.contains(&residue) {
                    active_residue_list.push(residue);
                }
                active_atom_list.push(atom)
            }
        }
    }

    let basic_table = table!(
        ["", "# of Atoms", "# of Residues"],
        ["QM1", qm1_atom_list.len(), qm1_residue_list.len()],
        ["QM2", qm2_atom_list.len(), qm2_residue_list.len()],
        ["Active", active_atom_list.len(), active_residue_list.len()]
    );

    basic_table.printstd();

    if target == Target::Residues {
        let residue_list = match region {
            Region::QM1 => qm1_residue_list,
            Region::QM2 => qm2_residue_list,
            Region::Active => active_residue_list,
            // Region::None => return Err("Invalid argument for 'region'".into()),
            Region::None => bail!("Invalid argument for 'region'"),
        };
        if !residue_list.is_empty() {
            let mut residue_table = Table::new();
            residue_table.add_row(row![
                "Residue ID",
                "Residue Name",
                "# of Atoms",
                match region {
                    Region::QM1 => "# of QM1 Atoms",
                    Region::QM2 => "# of QM2 Atoms",
                    Region::Active => "# of Active Atoms",
                    // Region::None => return Err("Invalid argument for 'region'".into()),
                    Region::None => bail!("Invalid argument for 'region'"),
                }
            ]);

            for residue in residue_list {
                let mut resid_atoms = 0;
                let mut atom_counter = 0;
                for atom in residue.atoms() {
                    atom_counter += 1;
                    if (region == Region::QM1 && atom.occupancy() == 1.00)
                        || (region == Region::QM2 && atom.occupancy() == 2.00)
                        || (region == Region::Active && atom.b_factor() == 1.00)
                    {
                        resid_atoms += 1;
                    }
                }

                residue_table.add_row(row![
                    residue.serial_number().to_string() + residue.insertion_code().unwrap_or(""),
                    residue.name().ok_or_else(|| anyhow!("No Residue Name"))?,
                    atom_counter,
                    resid_atoms,
                ]);
            }
            println!("\n{} Residues", region.to_owned());
            residue_table.printstd();
        } else {
            bail!("No Residues found in given region!");
        }
    } else if target == Target::Atoms {
        let (atom_list, residue_list) = match region {
            Region::QM1 => (qm1_atom_list, qm1_residue_list),
            Region::QM2 => (qm2_atom_list, qm2_residue_list),
            Region::Active => (active_atom_list, active_residue_list),
            Region::None => bail!("Invalid argument for 'region'"),
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
                            residue.name().ok_or_else(|| anyhow!("No Residue Name"))?,
                            atom.occupancy(),
                            atom.b_factor(),
                        ]);
                    }
                }
            }
            println!("\n{} Atoms", region.to_owned());
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

// pub fn add_insertion_codes(pdb: &mut PDB) -> Result<(), GenericErr> {
//     let insertion_codes = [
//         "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R",
//         "S", "T", "U", "V", "W", "X", "Y", "Z",
//     ];
//     for i in pdb
//         .residues_mut()
//         .chunks(9999)
//         .into_iter()
//         .zip_longest(insertion_codes)
//     {
//         match i {
//             Both(chunk, icode) => {
//                 for res in chunk {
//                     res.set_insertion_code(icode);
//                 }
//             }
//             Right(_) => (),
//             Left(_) => {
//                 return Err(
//                     "Too many residues for Latin alphabet. Please contact the developer.".into(),
//                 )
//             }
//         }
//     }
//     Ok(())
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

    // #[test]
    // fn parse_residue_list_test() {
    //     let num_list1 = "1-2,3:4,6,7";
    //     let num_list2 = "9999A:1B,6B,8B-10B";
    //     let str_list = "cu,NA+";

    //     let pdb1 = test_pdb("tests/test_blank.pdb");
    //     let pdb2 = test_pdb("tests/test_insert.pdb");

    //     assert_eq!(
    //         parse_residue_list(num_list1, &pdb1).unwrap(),
    //         vec!(
    //             (1, None),
    //             (2, None),
    //             (3, None),
    //             (4, None),
    //             (6, None),
    //             (7, None),
    //         )
    //     );

    //     assert_eq!(
    //         parse_residue_list(num_list2, &pdb2).unwrap(),
    //         vec!(
    //             (9999, Some("A")),
    //             (0, Some("B")),
    //             (1, Some("B")),
    //             (6, Some("B")),
    //             (8, Some("B")),
    //             (9, Some("B")),
    //             (10, Some("B")),
    //         )
    //     );

    //     assert_eq!(
    //         parse_residue_list(str_list, &pdb2).unwrap(),
    //         vec!((228, Some("A")), (229, Some("A")))
    //     );
    // }

    #[test]
    fn atom_sphere_test() {
        let pdb = test_pdb("tests/test_blank.pdb");
        let origin = pdb
            .atoms_with_hierarchy()
            .find(|x| x.atom.serial_number() == 26)
            .unwrap();
        let atom_list_incl = vec![21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 45, 46];
        let atom_list_excl = vec![21, 22, 23, 24, 25, 27, 28, 29, 30, 32, 45, 46];

        assert_eq!(
            calc_atom_sphere(&pdb, &origin, 4.0, true).unwrap().len(),
            22
        );
        assert_eq!(
            calc_atom_sphere(&pdb, &origin, 3.0, true).unwrap(),
            atom_list_incl
        );

        assert_eq!(
            calc_atom_sphere(&pdb, &origin, 4.0, false).unwrap().len(),
            21
        );
        assert_eq!(
            calc_atom_sphere(&pdb, &origin, 3.0, false).unwrap(),
            atom_list_excl
        );
    }

    #[test]
    fn residue_sphere_test() {
        let pdb = test_pdb("tests/test_blank.pdb");
        let origin = pdb
            .atoms_with_hierarchy()
            .find(|x| x.atom.serial_number() == 26)
            .unwrap();
        let atom_list_excl = vec![
            19, 20, 21, 22, 23, 24, 25, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61,
            62,
        ];

        assert_eq!(
            calc_residue_sphere(&pdb, &origin, 4.0, true).unwrap().len(),
            44
        );
        assert_eq!(
            calc_residue_sphere(&pdb, &origin, 4.0, false)
                .unwrap()
                .len(),
            23
        );
        assert_eq!(
            calc_residue_sphere(&pdb, &origin, 4.0, false).unwrap(),
            atom_list_excl
        );
    }

    #[test]
    fn edit_atoms_test() {
        let mut pdb = test_pdb("tests/test_blank.pdb");
        let atom_id_list = vec![1, 5, 9];
        edit_qm_atoms(&mut pdb, 1.00, atom_id_list.clone()).unwrap();
        edit_active_atoms(&mut pdb, 1.00, atom_id_list.clone()).unwrap();

        let atom_list = pdb
            .atoms()
            .filter(|x| atom_id_list.contains(&x.serial_number()));
        for atom in atom_list {
            assert_eq!(atom.occupancy(), 1.00);
            assert_eq!(atom.b_factor(), 1.00);
        }
    }

    #[test]
    fn edit_residues_test_all() {
        let mut pdb = test_pdb("tests/test_blank.pdb");
        let res_id_list = vec![(2, None), (4, None)];

        edit_qm_residues(&mut pdb, 2.00, res_id_list.clone(), Partial::None).unwrap();
        edit_active_residues(&mut pdb, 1.00, res_id_list.clone(), Partial::None).unwrap();

        let res_list = pdb.residues().filter(|x| {
            res_id_list.contains(&(x.serial_number(), x.insertion_code().map(ToOwned::to_owned)))
        });

        for residue in res_list {
            for atom in residue.atoms() {
                assert_eq!(atom.occupancy(), 2.00);
                assert_eq!(atom.b_factor(), 1.00);
            }
        }
    }

    #[test]
    fn edit_residues_test_side() {
        let mut pdb = test_pdb("tests/test_blank.pdb");
        let res_id_list = vec![(2, None), (4, None)];

        edit_qm_residues(&mut pdb, 2.00, res_id_list.clone(), Partial::Sidechain).unwrap();
        edit_active_residues(&mut pdb, 1.00, res_id_list.clone(), Partial::Sidechain).unwrap();

        let res_list = pdb.residues().filter(|x| {
            res_id_list.contains(&(x.serial_number(), x.insertion_code().map(ToOwned::to_owned)))
        });

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
        let res_id_list = vec![(2, None), (4, None)];
        edit_qm_residues(&mut pdb, 2.00, res_id_list.clone(), Partial::Backbone).unwrap();
        edit_active_residues(&mut pdb, 1.00, res_id_list.clone(), Partial::Backbone).unwrap();

        let res_list = pdb.residues().filter(|x| {
            res_id_list.contains(&(x.serial_number(), x.insertion_code().map(ToOwned::to_owned)))
        });
        let atom_list = res_list.flat_map(|x| x.atoms()).filter(|x| x.is_backbone());

        for atom in atom_list {
            assert_eq!(atom.occupancy(), 2.00);
            assert_eq!(atom.b_factor(), 1.00);
        }
    }

    #[test]
    fn remove_all_test() {
        let mut pdb = test_pdb("tests/test_full.pdb");
        remove_all(&mut pdb).unwrap();

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
}
