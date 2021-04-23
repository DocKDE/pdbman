use std::convert::TryFrom;
use std::error::Error;

use pdbtbx::{Atom, Residue, PDB};
use prettytable::{format, Table};
// use regex::Regex;
use crate::{Distance, Partial, Region, Target};
use rstar::{primitives::PointWithData, RTree};

// type Result<T> = std::result::Result<T, Box<dyn Error>>;

const AMINOS: [&str; 29] = [
    "ARG", "HIS", "HIP", "HID", "HIM", "HIE", "LYS", "LYN", "ASP", "ASH", "GLU", "GLH", "SER",
    "THR", "ASN", "GLN", "CYS", "CYX", "SEC", "GLY", "PRO", "ALA", "VAL", "ILE", "LEU", "MET",
    "PHE", "TYR", "TRP",
];
const BACKBONE_ATOMS: [&str; 8] = ["C", "O", "N", "H", "CA", "HA", "HA2", "HA3"];

#[derive(PartialEq, Debug, Clone)]
pub struct Sphere<'a> {
    pub origin: &'a pdbtbx::Atom,
    pub radius: f64,
}

/// This function parses the values of option 'Sphere' into a usize and an f64 and returns
/// a reference to the corresponding Atom and the radius. The return values are organized
/// in a Sphere struct to facilitate usage. 
impl<'a> Sphere<'a> {
    pub fn new(mut inp_str: clap::Values, pdb: &'a PDB) -> Result<Sphere<'a>, Box<dyn Error>> {
        let origin_str = inp_str.next().ok_or("No origin string present")?;
        let radius_str = inp_str.next().ok_or("No radius string present")?;

        // The sphere validator rejects anything containing something else than digits and dots.
        // Since it validates both values the same, it will not reject a decimal point in the Atom
        // ID which is taken care of here.
        let origin_id: usize = match origin_str.parse() {
            Ok(v) => v,
            Err(_) => return Err("Decimal points are not allowed in Atom IDs.".into())
        };

        // No more error handling should be necessary since the sphere validator already
        // takes care of everything that is not a float. Everything that can be constructed
        // with dots and digits is valid input for radius.
        let radius: f64 = radius_str.parse()?;

        let origin_atom = pdb
            .atoms()
            .find(|x| x.serial_number() == origin_id)
            .ok_or("No atom corresponding to the given ID could be found.")?;

        Ok(Sphere {
            origin: origin_atom,
            radius: radius,
        })
    }
}

/// This function edits the q value of the Molecule (used by ORCA as input for QM section
/// selection) by Residue. It takes a mutable reference of a PDB struct, the desired value q should
/// be set to and a list of Residue serial numbers for the Atoms to edit.
/// The last option recognizes 'backbone', 'sidechain' and 'None' as options and will select the
/// correspdonding Atoms from each Residue.
pub fn edit_qm_residues<'a>(
    pdb: &mut PDB,
    qm_val: f64,
    // list: impl Iterator<Item=Residue>,
    list: Vec<isize>,
    // region: &str,
    // partial: &str,
    partial: Partial,
) -> Result<(), Box<dyn Error>> {
    for residue in pdb.residues_mut() {
        let serial_number = residue.serial_number();
        let name = residue.name().ok_or("No Residue Name")?.to_string();
        for atom in residue.atoms_mut() {
            match partial {
                Partial::None => {
                    if list.contains(&serial_number) {
                        atom.set_occupancy(qm_val)?;
                    }
                    // if list.find(|x| x == residue).is_some() {
                    //     atom.set_occupancy(qm_val)?;
                    // }
                }
                Partial::Sidechain => {
                    if list.contains(&serial_number)
                    // if list.find(|x| x == residue).is_some()
                        && AMINOS.contains(&name.as_str())
                        && !BACKBONE_ATOMS.contains(&atom.name())
                    {
                        atom.set_occupancy(qm_val)?;
                    }
                }
                Partial::Backbone => {
                    if list.contains(&serial_number)
                    // if list.find(|x| x == residue).is_some()
                        && AMINOS.contains(&name.as_str())
                        && BACKBONE_ATOMS.contains(&atom.name())
                    {
                        atom.set_occupancy(qm_val)?;
                    }
                }
            }
        }
    }
    Ok(())
}

/// This function edits the b value of the Molecule (used by ORCA as input for active section
/// selection) by Residue. It takes a mutable reference of a PDB struct, the desired value b should
/// be set to and a list of Residue IDs for the Residues to edit.
/// The last option recognizes 'backbone', 'sidechain' and 'None' as options and will select the
/// correspdonding Atoms from each Residue.
pub fn edit_active_residues(
    pdb: &mut PDB,
    qm_val: f64,
    list: Vec<isize>,
    // region: &str,
    // partial: &str,
    partial: Partial,
) -> Result<(), Box<dyn Error>> {
    for residue in pdb.residues_mut() {
        let serial_number = residue.serial_number();
        let name = residue.name().ok_or("No Residue Name")?.to_string();

        for atom in residue.atoms_mut() {
            match partial {
                Partial::None => {
                    if list.contains(&serial_number) {
                        atom.set_b_factor(qm_val)?;
                    }
                }
                Partial::Sidechain => {
                    if list.contains(&serial_number)
                        && AMINOS.contains(&name.as_str())
                        && !BACKBONE_ATOMS.contains(&atom.name())
                    {
                        atom.set_b_factor(qm_val)?;
                    }
                }
                Partial::Backbone => {
                    if list.contains(&serial_number)
                        && AMINOS.contains(&name.as_str())
                        && BACKBONE_ATOMS.contains(&atom.name())
                    {
                        atom.set_b_factor(qm_val)?;
                    }
                }
            }
        }
    }
    Ok(())
}

/// This functions edits the q value of the PDB file (used by ORCA as input for QM region
/// selection) by Atom. It takes a mutable reference to a PDB struct, the desired value
/// q should be set to and a vector of Atom IDs.
pub fn edit_qm_atoms(pdb: &mut PDB, qm_val: f64, list: Vec<usize>) -> Result<(), Box<dyn Error>> {
    for atom in pdb.atoms_mut() {
        // if list.find(|x| x == atom).is_some() {
        if list.contains(&atom.serial_number()) {
            atom.set_occupancy(qm_val)?
        }
    }
    Ok(())
}

/// This functions edits the b value of the PDB file (used by ORCA as input for QM region
/// selection) by Atom. It takes a mutable reference to a PDB struct, the desired value
/// b should be set to and a vector of Atom IDs.
pub fn edit_active_atoms(
    pdb: &mut PDB,
    active_val: f64,
    list: Vec<usize>,
) -> Result<(), Box<dyn Error>> {
    for atom in pdb.atoms_mut() {
        if list.contains(&atom.serial_number()) {
            atom.set_b_factor(active_val)?
        }
    }
    Ok(())
}

/// Sets all q and b values to zero which serves as a fresh start.
pub fn remove_all(pdb: &mut PDB) -> Result<(), Box<dyn Error>> {
    for atom in pdb.atoms_mut() {
        atom.set_occupancy(0.00)?;
        atom.set_b_factor(0.00)?;
    }
    Ok(())
}

/// Prints all Atoms in Molecule to stdout in PDB file format This can be redirected
/// to a file if desired. This may be obsolete as the functionality of printing to a file
/// already exists with a separate flag.
pub fn print_to_stdout(pdb: &PDB) -> Result<(), Box<dyn Error>> {
    for residue in pdb.residues() {
        for atom in residue.atoms() {
            println!("ATOM  {:>5} {:<4} {:>3}  {:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}", 
                atom.serial_number(), atom.name(), residue.name().ok_or("No Residue Name")?, residue.serial_number(), atom.x(),
                atom.y(), atom.z(), atom.occupancy(), atom.b_factor(), atom.element());
        }
    }
    Ok(())
}

/// Takes an Atom struct as point of origin and a radius in A. Returns a Vector of Atom IDs
/// of Atoms within the given radius wrapped in a Result. Origin can be included or excluded.
pub fn calc_atom_sphere<'a>(
    pdb: &PDB,
    origin: &Atom,
    radius: f64,
    include_self: bool,
) -> Result<Vec<usize>, Box<dyn Error>> {
    // ) -> Result<impl Iterator<Item = &'a Atom>, Box<dyn Error>> {
    let mut sphere_atoms: Vec<usize> = Vec::new();

    // if include_self {
    //     return Ok(pdb.atoms().filter(|x| x.distance(origin) <= radius))
    // } else {
    //     return Ok(pdb.atoms().filter(|x| {x.distance(origin) <= radius && x != &origin}))
    // }

    for atom in pdb.atoms() {
        if include_self {
            if atom.distance(origin) <= radius {
                sphere_atoms.push(atom.serial_number())
            }
        } else if atom.distance(origin) <= radius && atom.serial_number() != origin.serial_number()
        {
            sphere_atoms.push(atom.serial_number())
        }
    }

    if !sphere_atoms.is_empty() {
        Ok(sphere_atoms)
    } else {
        Err("Calculated sphere doesn't contain any atoms.".into())
    }
}

/// Takes an Atom struct as point of origin, a radius in A. Returns a Vector of Atom IDs
/// that belong to a Residue that had at least one Atom within the given radius wrapped in a Result.
/// Origin can be included or excluded.
pub fn calc_residue_sphere(
    pdb: &PDB,
    origin: &Atom,
    radius: f64,
    include_self: bool,
) -> Result<Vec<usize>, Box<dyn Error>> {
    let mut sphere_atoms: Vec<usize> = Vec::new();
    let mut sphere_residues: Vec<&Residue> = Vec::new();

    let mut origin_res_id = 0;
    for residue in pdb.residues() {
        for atom in residue.atoms() {
            if atom.serial_number() == origin.serial_number() {
                origin_res_id = residue.serial_number();
            }
        }
    }

    for residue in pdb.residues() {
        for atom in residue.atoms() {
            if include_self {
                if atom.distance(origin) <= radius {
                    sphere_residues.push(residue);
                }
            } else if atom.distance(origin) <= radius && residue.serial_number() != origin_res_id {
                sphere_residues.push(residue)
            }
        }
    }

    for residue in sphere_residues {
        for atom in residue.atoms() {
            sphere_atoms.push(atom.serial_number());
        }
    }

    sphere_atoms.sort_unstable();
    sphere_atoms.dedup();

    if !sphere_atoms.is_empty() {
        Ok(sphere_atoms)
    } else {
        Err("Calculated sphere doesn't contain any atoms.sphere_atoms".into())
    }
}

/// Finds and prints all contacts present in the PDB file structure. Definition of
/// 'contact' is given by the 'level' arg which is 1.0A for 'level' = 0 and the
/// respective atomic radius for 'level' = 1.
pub fn find_contacts(pdb: &PDB, level: Distance) -> Result<Table, Box<dyn Error>> {
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

    let mut treevec = Vec::new();
    for residue in pdb.residues() {
        for atom in residue.atoms() {
            treevec.push(PointWithData::new(
                (atom, residue),
                // [atom.x(), atom.y(), atom.z()],
                atom.pos_array(),
            ))
        }
    }

    let tree = RTree::bulk_load(treevec);

    for residue in pdb.residues() {
        for atom in residue.atoms() {
            let radius: f64 = match level {
                Distance::Clashes => 1.0,
                Distance::Contacts => atom
                    .atomic_radius()
                    .ok_or("No atomic radius found for Atom type!")?
                    .powf(2.0),
                _ => return Err("No Valid 'level' argument".into()),
            };
            let contacts = tree.locate_within_distance(atom.pos_array(), radius);

            for PointWithData {
                data: (other_atom, other_residue),
                ..
            } in contacts
            {
                if other_atom < &atom
                    && !(other_residue.name() == residue.name()
                        && other_residue.serial_number() == residue.serial_number())
                    && !(other_atom.name() == "C"
                        && atom.name() == "N"
                        && other_residue.serial_number() + 1 == residue.serial_number())
                {
                    let distance = other_atom.distance(atom);
                    if distance <= 0.75 && distance > 0.5 {
                        table.add_row(row![
                            bByFd =>
                            other_atom.serial_number(),
                            other_atom.name(),
                            other_residue.name().ok_or("No Residue Name.")?,
                            atom.serial_number(),
                            atom.name(),
                            residue.name().ok_or("No Residue Name")?,
                            format!("{:.2}", distance)
                        ]);
                    } else if distance <= 0.5 {
                        table.add_row(row![
                            bBrFd =>
                            other_atom.serial_number(),
                            other_atom.name(),
                            other_residue.name().ok_or("No Residue Name.")?,
                            atom.serial_number(),
                            atom.name(),
                            residue.name().ok_or("No Residue Name")?,
                            format!("{:.2}", distance)
                        ]);
                    } else {
                        table.add_row(row![
                            other_atom.serial_number(),
                            other_atom.name(),
                            other_residue.name().ok_or("No Residue Name.")?,
                            atom.serial_number(),
                            atom.name(),
                            residue.name().ok_or("No Residue Name")?,
                            format!("{:.2}", distance)
                        ]);
                    }
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

// ///Input may contain only one instance of a range-indicating character
fn expand_range(input: &str) -> Result<Vec<usize>, std::num::ParseIntError> {
    let vector: Vec<usize> = input
        .split(&['-', ':'][..])
        .map(|x| x.parse())
        .collect::<Result<_, _>>()?;
    Ok((vector[0]..=vector[1]).collect())
}

/// Takes a comma-separated list (usually from command line input) as string and parses it into
/// a vector of Atom IDs. The Input may be atom IDs or Atom Names
pub fn parse_atomic_list(input: &str, pdb: &PDB) -> Result<Vec<usize>, Box<dyn Error>> {
    let mut output_vec: Vec<usize> = vec![];

    match input
        .split(&[',', '-', ':'][..])
        .all(|x| x.parse::<usize>().is_ok())
    {
        true => {
            let input_vec: Vec<&str> = input.split(',').collect();

            for i in input_vec {
                if i.contains(&['-', ':'][..]) {
                    output_vec.extend(expand_range(i)?)
                } else {
                    output_vec.push(i.parse()?)
                }
            }
        }
        false => match input
            .split(&[',', '-', ':'][..])
            .all(|x| x.parse::<usize>().is_err())
        {
            true => {
                let input_vec: Vec<String> = input.split(',').map(|x| x.to_lowercase()).collect();
                output_vec = pdb
                    .atoms()
                    .filter(|x| input_vec.contains(&x.name().to_lowercase()))
                    .map(|x| x.serial_number())
                    .collect();
            }
            false => return Err("Input List contains mixed types which is not supported.".into()),
        },
    }

    Ok(output_vec)
}

pub fn parse_residue_list(input: &str, pdb: &PDB) -> Result<Vec<isize>, Box<dyn Error>> {
    let mut output_vec: Vec<isize> = vec![];

    match input
        .split(&[',', '-', ':'][..])
        .all(|x| x.parse::<isize>().is_ok())
    {
        true => {
            let input_vec: Vec<&str> = input.split(',').collect();

            for i in input_vec {
                if i.contains(&['-', ':'][..]) {
                    output_vec.extend(
                        expand_range(i)?
                            .into_iter()
                            .map(|x| isize::try_from(x))
                            .collect::<Result<Vec<isize>, _>>()?,
                    )
                } else {
                    output_vec.push(i.parse()?)
                }
            }
        }
        false => match input
            .split(&[',', '-', ':'][..])
            .all(|x| x.parse::<isize>().is_err())
        {
            true => {
                let input_vec: Vec<String> = input.split(',').map(|x| x.to_lowercase()).collect();
                output_vec = pdb
                    .residues()
                    .map(|x| {
                        Ok(input_vec
                            .contains(&x.name().ok_or("No Residue Name")?.to_lowercase())
                            .then(|| x.serial_number()))
                    })
                    .filter_map(Result::transpose)
                    .collect::<Result<Vec<isize>, Box<dyn Error>>>()?;
                // .filter(|x| input_vec.contains(&x.name().unwrap().to_lowercase()))
                // .map(|x| x.serial_number())
                // .collect();
            }
            false => return Err("Input List contains mixed types which is not supported.".into()),
        },
    }

    Ok(output_vec)
}

/// Query Molecule for information. Depending on the input this will print a table of
/// Residues and/or Atoms will available information that were asked for.
pub fn query_atoms(pdb: &PDB, atom_list: Vec<usize>) -> Result<(), Box<dyn Error>> {
    let mut table = Table::new();
    table.add_row(row![
        "Atom ID",
        "Atom name",
        "Residue ID",
        "Residue Name",
        "QM",
        "Active"
    ]);

    for residue in pdb.residues() {
        for atom in residue.atoms() {
            if atom_list.contains(&atom.serial_number()) {
                table.add_row(row![
                    atom.serial_number(),
                    atom.name(),
                    residue.serial_number(),
                    residue.name().ok_or("No Residue Name")?,
                    atom.occupancy(),
                    atom.b_factor(),
                ]);
            }
        }
    }

    if !atom_list.is_empty() {
        table.printstd();
        Ok(())
    } else {
        Err("No Atoms found!".into())
    }
}

pub fn query_residues(pdb: &PDB, residue_list: Vec<isize>) -> Result<(), Box<dyn Error>> {
    let mut table = Table::new();
    table.add_row(row![
        "Atom ID",
        "Atom name",
        "Residue ID",
        "Residue Name",
        "QM",
        "Active"
    ]);

    for residue in pdb.residues() {
        for atom in residue.atoms() {
            if residue_list.contains(&residue.serial_number()) {
                table.add_row(row![
                    atom.serial_number(),
                    atom.name(),
                    residue.serial_number(),
                    residue.name().ok_or("No Residue name")?,
                    atom.occupancy(),
                    atom.b_factor(),
                ]);
            }
        }
    }

    if !residue_list.is_empty() {
        table.printstd();
        Ok(())
    } else {
        Err("No Residues found!".into())
    }
}

pub fn analyze(pdb: &PDB, region: Region, target: Target) -> Result<(), Box<dyn Error>> {
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
            Region::None => return Err("Invalid argument for 'region'".into()),
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
                    Region::None => return Err("Invalid argument for 'region'".into()),
                }
            ]);

            for residue in residue_list {
                let mut resid_atoms = 0;
                let mut atom_counter = 0;
                for atom in residue.atoms() {
                    atom_counter += 1;
                    if region == Region::QM1 && atom.occupancy() == 1.00 {
                        resid_atoms += 1;
                    } else if region == Region::QM2 && atom.occupancy() == 2.00 {
                        resid_atoms += 1;
                    } else if region == Region::Active && atom.b_factor() == 1.00 {
                        resid_atoms += 1;
                    }
                }

                residue_table.add_row(row![
                    residue.serial_number(),
                    residue.name().ok_or("No Residue Name")?,
                    atom_counter,
                    resid_atoms,
                ]);
            }
            println!("\n{} Residues", region.to_string());
            residue_table.printstd();
        } else {
            return Err("No Residues found in given region!".into());
        }
    } else if target == Target::Atoms {
        let (atom_list, residue_list) = match region {
            Region::QM1 => (qm1_atom_list, qm1_residue_list),
            Region::QM2 => (qm2_atom_list, qm2_residue_list),
            Region::Active => (active_atom_list, active_residue_list),
            Region::None => return Err("Invalid argument for 'region'".into()),
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
                            residue.serial_number(),
                            residue.name().ok_or("No Residue Name")?,
                            atom.occupancy(),
                            atom.b_factor(),
                        ]);
                    }
                }
            }
            println!("\n{} Atoms", region.to_string());
            atom_table.printstd();
        } else {
            return Err("No Atoms found in given region!".into());
        }
    }
    Ok(())
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
    fn atom_sphere_test() {
        let pdb = test_pdb("tests/test_blank.pdb");
        let origin = pdb.atoms().find(|x| x.serial_number() == 26).unwrap();
        let atom_list_incl = vec![21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 45, 46];
        let atom_list_excl = vec![21, 22, 23, 24, 25, 27, 28, 29, 30, 32, 45, 46];

        assert_eq!(calc_atom_sphere(&pdb, origin, 4.0, true).unwrap().len(), 22);
        assert_eq!(
            calc_atom_sphere(&pdb, origin, 3.0, true).unwrap(),
            atom_list_incl
        );

        assert_eq!(
            calc_atom_sphere(&pdb, origin, 4.0, false).unwrap().len(),
            21
        );
        assert_eq!(
            calc_atom_sphere(&pdb, origin, 3.0, false).unwrap(),
            atom_list_excl
        );
    }

    #[test]
    fn residue_sphere_test() {
        let pdb = test_pdb("tests/test_blank.pdb");
        let origin = pdb.atoms().find(|x| x.serial_number() == 26).unwrap();
        let atom_list_excl = vec![
            19, 20, 21, 22, 23, 24, 25, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61,
            62,
        ];

        assert_eq!(
            calc_residue_sphere(&pdb, origin, 4.0, true).unwrap().len(),
            44
        );
        assert_eq!(
            calc_residue_sphere(&pdb, origin, 4.0, false).unwrap().len(),
            23
        );
        assert_eq!(
            calc_residue_sphere(&pdb, origin, 4.0, false).unwrap(),
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
        let res_id_list = vec![2, 4];

        edit_qm_residues(&mut pdb, 2.00, res_id_list.clone(), Partial::None).unwrap();
        edit_active_residues(&mut pdb, 1.00, res_id_list.clone(), Partial::None).unwrap();

        let res_list = pdb
            .residues()
            .filter(|x| res_id_list.contains(&x.serial_number()));

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
        let res_id_list = vec![2, 4];
        // let backbone_atoms = ["C", "O", "N", "H", "CA", "HA", "HA2", "HA3"];

        edit_qm_residues(&mut pdb, 2.00, res_id_list.clone(), Partial::Sidechain).unwrap();
        edit_active_residues(&mut pdb, 1.00, res_id_list.clone(), Partial::Sidechain).unwrap();

        let res_list = pdb
            .residues()
            .filter(|x| res_id_list.contains(&x.serial_number()));

        let atom_list = res_list
            .map(|x| x.atoms())
            .flatten()
            .filter(|x| !BACKBONE_ATOMS.contains(&x.name()));

        for atom in atom_list {
            assert_eq!(atom.occupancy(), 2.00);
            assert_eq!(atom.b_factor(), 1.00);
        }
    }

    #[test]
    fn edit_residues_test_back() {
        let mut pdb = test_pdb("tests/test_blank.pdb");
        let res_id_list = vec![2, 4];
        // let backbone_atoms = ["C", "O", "N", "H", "CA", "HA", "HA2", "HA3"];
        edit_qm_residues(&mut pdb, 2.00, res_id_list.clone(), Partial::Backbone).unwrap();
        edit_active_residues(&mut pdb, 1.00, res_id_list.clone(), Partial::Backbone).unwrap();

        let res_list = pdb
            .residues()
            .filter(|x| res_id_list.contains(&x.serial_number()));
        let atom_list = res_list
            .map(|x| x.atoms())
            .flatten()
            .filter(|x| BACKBONE_ATOMS.contains(&x.name()));

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
