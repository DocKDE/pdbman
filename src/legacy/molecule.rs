use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::prelude::{BufRead, Write};
use std::io::{BufReader, LineWriter};

use prettytable::{format, Table};
use rstar::primitives::PointWithData;
use rstar::RTree;

type Result<T> = std::result::Result<T, Box<dyn Error>>;

#[derive(Debug)]
/// The Molecule struct follows a hierarchy of three levels:
/// * Molecule
///     * Residue
///         * Atom
///
/// Each level holds information about the one below it (may be changed somewhat in the
/// future). The constructor method meant for use by the user is parse_pdb which reads
/// a PDB file and extracts the information from there. However, it should be noted that
/// this program is most useful for PDB files that have been preprocessed in some way
/// (cleaned, missing atoms added, alternate conformations and chains removed etc.). It
/// should still work on raw PDB files but I can give no guarantee for it.
/// The Molecule strcut implements a couple of useful methods to obtain information about
/// the residues and atoms it contains (see documentation of methods).
pub struct Molecule {
    residues: Vec<Residue>,
}

impl Molecule {
    /// This reads all lines that begin with 'ATOM' or 'HETATM' and parses the contents
    /// into a Molecule struct. Since PDB files only contain ASCII characters, string
    /// slicing is used here which should be safe for this purpose. The parsing is conducted
    /// in compliance with this: http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html
    fn parse_pdb(path: &str) -> Result<Vec<Residue>> {
        let f = File::open(path)?;
        let f = BufReader::new(f);
        let mut molecule: Vec<Residue> = Vec::new();
        let mut res_counter: usize = 0;

        for line in f.lines() {
            let l = line?;

            // Parse all lines beginning with atom identifiers
            if l.starts_with("ATOM") || l.starts_with("HETATM") {
                let res_id: usize = l[22..26].trim().parse()?;
                let atom = Atom {
                    atom_id: l[6..11].trim().parse()?,
                    // Convert all of these to type "String" so the struct owns them.
                    // TODO: Maybe an owned string isn't necessary here and makes comparisons trickier
                    atom_name: l[12..16].trim().to_string(),
                    res_name: l[17..20].trim().to_string(),
                    res_id,
                    coords: [
                        l[30..38].trim().parse()?,
                        l[38..46].trim().parse()?,
                        l[46..54].trim().parse()?,
                    ],
                    qm: l[54..60].trim().parse()?,
                    active: l[60..66].trim().parse()?,
                    element: l[76..78].trim().to_string(),
                };

                if res_id != res_counter {
                    let residue = Residue {
                        res_name: l[17..20].trim().to_string(),
                        res_id,
                        atoms: Vec::new(),
                    };
                    molecule.push(residue);
                }
                molecule.last_mut().unwrap().atoms.push(atom);
                res_counter = res_id;
            }
        }

        if molecule.is_empty() {
            return Err(
                "No Atoms could be parsed! Please make sure you have a valid PDB file!".into(),
            );
        }
        Ok(molecule)
    }

    /// Constructor method for Molecule using the parse_pdb method
    pub fn from_pdb(filename: &str) -> Result<Molecule> {
        Ok(Molecule {
            residues: Molecule::parse_pdb(filename)?,
        })
    }

    /// Given an Atom ID returns the corresponding Atom struct wrapped in an Option
    pub fn atom_from_id(&self, id: usize) -> Option<&Atom> {
        for residue in &self.residues {
            for atom in &residue.atoms {
                if atom.atom_id == id {
                    return Some(&atom);
                }
            }
        }
        None
    }

    /// Given a vector of Atom IDs returns a vector of the corresponding Atoms wrapped in an Option.
    pub fn atoms_from_idlist(&self, list: Vec<usize>) -> Option<Vec<&Atom>> {
        let mut atom_list: Vec<&Atom> = Vec::new();
        for residue in &self.residues {
            for atom in &residue.atoms {
                if list.contains(&atom.atom_id) {
                    atom_list.push(atom)
                }
            }
        }
        if !atom_list.is_empty() {
            Some(atom_list)
        } else {
            None
        }
    }

    /// Given a Vector of Residue IDs returns a Vector of the contained Atom Structs wrapped in an
    /// Option
    pub fn atoms_byres_from_idlist(&self, list: Vec<usize>) -> Option<Vec<&Atom>> {
        let mut atom_list: Vec<&Atom> = Vec::new();
        for residue in &self.residues {
            if list.contains(&residue.res_id) {
                for atom in &residue.atoms {
                    atom_list.push(atom)
                }
            }
        }
        if !atom_list.is_empty() {
            Some(atom_list)
        } else {
            None
        }
    }

    /// Given a Residue ID returns the correspnding Residue struct wrapped in an Option
    pub fn residue_from_id(&self, id: usize) -> Option<&Residue> {
        for residue in &self.residues {
            if residue.res_id == id {
                return Some(residue);
            }
        }
        None
    }

    /// Given a list of Residue IDs returns a Vector of the corresponding Residue structs wrapped
    /// in an Option
    pub fn residues_from_idlist(&self, list: Vec<usize>) -> Option<Vec<&Residue>> {
        let mut residue_list: Vec<&Residue> = Vec::new();
        for residue in &self.residues {
            if list.contains(&residue.res_id) {
                residue_list.push(residue)
            }
        }
        if !residue_list.is_empty() {
            Some(residue_list)
        } else {
            None
        }
    }

    /// Given an Atom ID returns the corresponding Residue struct wrapped in an Option
    pub fn residue_from_atom_id(&self, id: usize) -> Option<&Residue> {
        for residue in &self.residues {
            for atom in &residue.atoms {
                if atom.atom_id == id {
                    return Some(residue);
                }
            }
        }
        None
    }

    /// Given a Vector of Atom IDs returns a Vector of the corresponding Residues wrapped in an Option
    pub fn residues_byatom_from_list(&self, list: Vec<usize>) -> Option<Vec<&Residue>> {
        let mut residue_list: Vec<&Residue> = Vec::new();
        for residue in &self.residues {
            for atom in &residue.atoms {
                if list.contains(&atom.atom_id) {
                    residue_list.push(residue)
                }
            }
        }
        if !residue_list.is_empty() {
            Some(residue_list)
        } else {
            None
        }
    }

    /// This function edits the q value of the Molecule (used by ORCA as input for QM section
    /// selection) by Residue. It takes a mutable reference of a Molecule struct, the desired value q should
    /// be set to, a list of Residue IDs for the Atoms to edit and an optional target region.
    /// The last option recognizes 'backbone', 'sidechain' and 'None' as options and will select the
    /// correspdonding Atoms from each Residue.
    pub fn edit_qm_residues(&mut self, qm_val: f64, list: Vec<usize>, partial: &str) -> Result<()> {
        let amino_list = [
            "ARG".to_string(),
            "HIS".to_string(),
            "HIP".to_string(),
            "HID".to_string(),
            "HIM".to_string(),
            "HIE".to_string(),
            "LYS".to_string(),
            "LYN".to_string(),
            "ASP".to_string(),
            "ASH".to_string(),
            "GLU".to_string(),
            "GLH".to_string(),
            "SER".to_string(),
            "THR".to_string(),
            "ASN".to_string(),
            "GLN".to_string(),
            "CYS".to_string(),
            "CYX".to_string(),
            "SEC".to_string(),
            "GLY".to_string(),
            "PRO".to_string(),
            "ALA".to_string(),
            "VAL".to_string(),
            "ILE".to_string(),
            "LEU".to_string(),
            "MET".to_string(),
            "PHE".to_string(),
            "TYR".to_string(),
            "TRP".to_string(),
        ];
        let backbone_atoms = [
            "C".to_string(),
            "O".to_string(),
            "N".to_string(),
            "H".to_string(),
            "CA".to_string(),
            "HA".to_string(),
            // In case of Glycin
            "HA2".to_string(),
            "HA3".to_string(),
        ];
        for residue in &mut self.residues {
            for atom in &mut residue.atoms {
                match partial {
                    "None" => {
                        if list.contains(&residue.res_id) {
                            atom.qm = qm_val
                        }
                    }
                    "Sidechain" => {
                        if list.contains(&residue.res_id)
                            && amino_list.contains(&residue.res_name)
                            && !backbone_atoms.contains(&atom.atom_name)
                        {
                            atom.qm = qm_val;
                        }
                    }
                    "Backbone" => {
                        if list.contains(&residue.res_id)
                            && amino_list.contains(&residue.res_name)
                            && backbone_atoms.contains(&atom.atom_name)
                        {
                            atom.qm = qm_val;
                        }
                    }
                    &_ => return Err("Passed argument 'partial' not recognized".into()),
                }
            }
        }
        Ok(())
    }

    /// This function edits the b value of the Molecule (used by ORCA as input for active section
    /// selection) by Residue. It takes a mutable reference of a Molecule struct, the desired value b should
    /// be set to, a list of Residue IDs for the Residues to edit and an optional target region.
    /// The last option recognizes 'backbone', 'sidechain' and 'None' as options and will select the
    /// correspdonding Atoms from each Residue.
    pub fn edit_active_residues(
        &mut self,
        active_val: f64,
        list: Vec<usize>,
        partial: &str,
    ) -> Result<()> {
        //This list contains the AMBER style nomenclature
        let amino_list = [
            "ARG".to_string(),
            "HIS".to_string(),
            "HIP".to_string(),
            "HID".to_string(),
            "HIM".to_string(),
            "HIE".to_string(),
            "LYS".to_string(),
            "LYN".to_string(),
            "ASP".to_string(),
            "ASH".to_string(),
            "GLU".to_string(),
            "GLH".to_string(),
            "SER".to_string(),
            "THR".to_string(),
            "ASN".to_string(),
            "GLN".to_string(),
            "CYS".to_string(),
            "CYX".to_string(),
            "SEC".to_string(),
            "GLY".to_string(),
            "PRO".to_string(),
            "ALA".to_string(),
            "VAL".to_string(),
            "ILE".to_string(),
            "LEU".to_string(),
            "MET".to_string(),
            "PHE".to_string(),
            "TYR".to_string(),
            "TRP".to_string(),
        ];
        let backbone_atoms = [
            "C".to_string(),
            "O".to_string(),
            "N".to_string(),
            "H".to_string(),
            "CA".to_string(),
            "HA".to_string(),
            // In case of Glycin
            "HA2".to_string(),
            "HA3".to_string(),
        ];
        for residue in &mut self.residues {
            for atom in &mut residue.atoms {
                match partial {
                    "None" => {
                        if list.contains(&residue.res_id) {
                            atom.active = active_val
                        }
                    }
                    "Sidechain" => {
                        if list.contains(&residue.res_id)
                            && amino_list.contains(&residue.res_name)
                            && !backbone_atoms.contains(&atom.atom_name)
                        {
                            atom.active = active_val;
                        }
                    }
                    "Backbone" => {
                        if list.contains(&residue.res_id)
                            && amino_list.contains(&residue.res_name)
                            && backbone_atoms.contains(&atom.atom_name)
                        {
                            atom.active = active_val;
                        }
                    }
                    &_ => return Err("Passed argument 'partial' not recognized".into()),
                }
            }
        }
        Ok(())
    }

    /// This functions edits the q value of the Molecule (used by ORCA as input for QM region
    /// selection) by Atom. It takes a mutable reference to a Molecule struct, the desired value
    /// q should be set to a list of Atom IDs wrapped in an Option.
    pub fn edit_qm_atoms(&mut self, qm_val: f64, list: Vec<usize>) {
        for residue in &mut self.residues {
            for atom in &mut residue.atoms {
                if list.contains(&atom.atom_id) {
                    atom.qm = qm_val;
                }
            }
        }
    }

    /// This functions edits the b value of the Molecule (used by ORCA as input for QM region
    /// selection) by Atom. It takes a mutable reference to a Molecule struct, the desired value
    /// b should be set to a list of Atom IDs wrapped in an Option.
    pub fn edit_active_atoms(&mut self, active_val: f64, list: Vec<usize>) {
        for residue in &mut self.residues {
            for atom in &mut residue.atoms {
                if list.contains(&atom.atom_id) {
                    atom.active = active_val;
                }
            }
        }
    }

    /// Sets all q and b values to zero which serves as a fresh start.
    pub fn remove_all(&mut self) {
        for residue in &mut self.residues {
            for atom in &mut residue.atoms {
                atom.qm = 0.00;
                atom.active = 0.00;
            }
        }
    }

    /// Prints all Atoms in Molecule to stdout in PDB file format This can be redirected
    /// to a file if desired.
    pub fn print_to_stdout(&self) {
        for residue in &self.residues {
            for atom in &residue.atoms {
                println!("ATOM  {:>5} {:<4} {:>3}  {:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}", 
                    atom.atom_id, atom.atom_name, residue.res_name, residue.res_id, atom.coords[0],
                    atom.coords[1], atom.coords[2], atom.qm, atom.active, atom.element);
            }
        }
    }

    /// Prints all Atoms in Molecule to a given filename in the current directory.
    pub fn print_to_file(&self, path: &str) -> Result<()> {
        let file = File::create(path)?;
        let mut file = LineWriter::new(file);

        for residue in &self.residues {
            for atom in &residue.atoms {
                file.write_all(format!("ATOM  {:>5} {:<4} {:>3}  {:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}\n", 
                    atom.atom_id, atom.atom_name, residue.res_name, residue.res_id, atom.coords[0],
                    atom.coords[1], atom.coords[2], atom.qm, atom.active, atom.element).as_bytes())?;
            }
        }
        Ok(())
    }

    /// Takes an Atom struct as point of origin, a radius in A. Returns a Vector of Atom IDs
    /// of Atoms within the given radius wrapped in an Option. Origin can be included or excluded.
    pub fn calc_atom_sphere(
        &self,
        origin: &Atom,
        radius: f64,
        include_self: bool,
    ) -> Result<Vec<usize>> {
        let mut sphere_atoms: Vec<usize> = Vec::new();

        for residue in &self.residues {
            for atom in &residue.atoms {
                if include_self {
                    if atom.calc_distance(origin) <= radius {
                        sphere_atoms.push(atom.atom_id)
                    }
                } else if atom.calc_distance(origin) <= radius && atom.atom_id != origin.atom_id {
                    sphere_atoms.push(atom.atom_id)
                }
            }
        }

        if !sphere_atoms.is_empty() {
            Ok(sphere_atoms)
        } else {
            Err("Calculated sphere doesn't contain any atoms!".into())
        }
    }

    /// Takes an Atom struct as point of origin, a radius in A. Returns a Vector of Atom IDs
    /// that belong to a Residue that had at least one Atom within the given radius wrapped in an Option.
    /// Origin can be included or excluded.
    pub fn calc_residue_sphere(
        &self,
        origin: &Atom,
        radius: f64,
        include_self: bool,
    ) -> Result<Vec<usize>> {
        let mut sphere_atoms: Vec<usize> = Vec::new();
        let mut sphere_residues: Vec<&Residue> = Vec::new();

        for residue in &self.residues {
            for atom in &residue.atoms {
                if include_self {
                    if atom.calc_distance(origin) <= radius {
                        sphere_residues.push(residue);
                    }
                } else if atom.calc_distance(origin) <= radius && atom.res_id != origin.res_id {
                    sphere_residues.push(residue);
                }
            }
        }

        for residue in sphere_residues {
            for atom in &residue.atoms {
                sphere_atoms.push(atom.atom_id)
            }
        }

        sphere_atoms.sort_unstable();
        sphere_atoms.dedup();

        if !sphere_atoms.is_empty() {
            Ok(sphere_atoms)
        } else {
            Err("Calculated sphere doesn't contain any atoms!".into())
        }
    }

    /// Finds and prints all contacts present in the PDB file structure. Definition of
    /// 'contact' is given by the 'level' arg which is 1.0A for 'level' = 0 and the
    /// respective atomic radius for 'level' = 1.
    pub fn find_contacts(&self, level: u8) {
        let mut vdw_radii = HashMap::new();
        vdw_radii.insert("H".to_string(), 1.54);
        vdw_radii.insert("C".to_string(), 1.90);
        vdw_radii.insert("N".to_string(), 1.79);
        vdw_radii.insert("O".to_string(), 1.71);
        vdw_radii.insert("P".to_string(), 2.23);
        vdw_radii.insert("S".to_string(), 2.14);
        vdw_radii.insert("CL".to_string(), 2.06);
        vdw_radii.insert("NA".to_string(), 2.25);
        vdw_radii.insert("CU".to_string(), 2.17);

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
        for residue in &self.residues {
            for atom in &residue.atoms {
                treevec.push(PointWithData::new(atom,  atom.coords))
            }
        }

        let tree = RTree::bulk_load(treevec);

        for residue in &self.residues {
            for atom in &residue.atoms {
                let radius: f64 = match level {
                    0 => 1.0,
                    1 => {
                        let rad = vdw_radii
                            .get(&atom.element.to_uppercase())
                            .expect("No Radius found for given element.");
                        rad * rad
                    }
                    _ => panic!("Too high level given"),
                };
                let contacts = tree.locate_within_distance(atom.coords, radius);
                for item in contacts {
                    if item.data.atom_id < atom.atom_id
                        && !(item.data.res_id == atom.res_id && item.data.res_name == atom.res_name)
                        && !(item.data.atom_name == "C"
                            && atom.atom_name == "N"
                            && item.data.res_id + 1 == atom.res_id)
                    {
                        let distance = item.data.calc_distance(atom);
                        if distance <= 0.75 && distance > 0.5 {
                            table.add_row(row![
                                bByFd =>
                                item.data.atom_id,
                                item.data.atom_name,
                                item.data.res_name,
                                atom.atom_id,
                                atom.atom_name,
                                atom.res_name,
                                format!("{:.2}", distance)
                            ]);
                        } else if distance <= 0.5 {
                            table.add_row(row![
                                bBrFd =>
                                item.data.atom_id,
                                item.data.atom_name,
                                item.data.res_name,
                                atom.atom_id,
                                atom.atom_name,
                                atom.res_name,
                                format!("{:.2}", distance)
                            ]);
                        } else {
                            table.add_row(row![
                                item.data.atom_id,
                                item.data.atom_name,
                                item.data.res_name,
                                atom.atom_id,
                                atom.atom_name,
                                atom.res_name,
                                format!("{:.2}", distance)
                            ]);
                        }
                    }
                }
            }
        }
        if table.len() > 1 {
            if level == 0 {
                println!("\nClash Analysis");
            } else if level == 1 {
                println!("\nContact Analysis");
            }
            table.printstd();
        }
    }

    /// Query Molecule for information. Depending on the input this will print a table of
    /// Residues and/or Atoms will available information that were asked for.
    pub fn query_atoms(&self, list: Vec<usize>) {
        let mut table = Table::new();
        table.add_row(row![
            "Atom ID",
            "Atom name",
            "Residue ID",
            "Residue Name",
            "QM",
            "Active"
        ]);

        let atom_list = self.atoms_from_idlist(list).unwrap();

        for atom in atom_list {
            table.add_row(row![
                atom.atom_id,
                atom.atom_name,
                atom.res_id,
                atom.res_name,
                atom.qm,
                atom.active
            ]);
        }

        table.printstd();
    }

    pub fn query_residues(&self, list: Vec<usize>) {
        let mut table = Table::new();
        table.add_row(row![
            "Atom ID",
            "Atom name",
            "Residue ID",
            "Residue Name",
            "QM",
            "Active"
        ]);

        let res_list = self.residues_from_idlist(list).unwrap();

        for residue in res_list {
            for atom in &residue.atoms {
                table.add_row(row![
                    atom.atom_id,
                    atom.atom_name,
                    atom.res_id,
                    atom.res_name,
                    atom.qm,
                    atom.active
                ]);
            }
        }
        table.printstd();
    }

    pub fn analyze(&self, region: &str, verbosity: u8) -> Result<()> {
        let mut qm1_residue_list = Vec::new();
        let mut qm1_atom_list = Vec::new();
        let mut qm2_residue_list = Vec::new();
        let mut qm2_atom_list = Vec::new();
        let mut active_residue_list = Vec::new();
        let mut active_atom_list = Vec::new();

        for residue in &self.residues {
            for atom in &residue.atoms {
                if atom.qm == 1.00 {
                    if !qm1_residue_list.contains(&residue) {
                        qm1_residue_list.push(residue);
                    }
                    qm1_atom_list.push(atom)
                }
                if atom.qm == 2.00 {
                    if !qm2_residue_list.contains(&residue) {
                        qm2_residue_list.push(residue);
                    }
                    qm2_atom_list.push(atom)
                }
                if atom.active == 1.00 {
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

        if verbosity == 1 {
            let residue_list = match region {
                "QM1" => qm1_residue_list,
                "QM2" => qm2_residue_list,
                "Active" => active_residue_list,
                &_ => return Err("No region given.".into()),
            };

            if !residue_list.is_empty() {
                let mut residue_table = Table::new();
                residue_table.add_row(row![
                    "Residue ID",
                    "Residue Name",
                    "# of Atoms",
                    "# of QM1 Atoms"
                ]);

                for residue in residue_list {
                    let mut resid_atoms = 0;
                    for atom in &residue.atoms {
                        if region == "QM1" && atom.qm == 1.00 {
                            resid_atoms += 1;
                        } else if region == "QM2" && atom.qm == 2.00 {
                            resid_atoms += 1;
                        } else if region == "Active" && atom.active == 1.00 {
                            resid_atoms += 1;
                        }
                    }

                    residue_table.add_row(row![
                        residue.res_id,
                        residue.res_name,
                        residue.atoms.len(),
                        resid_atoms,
                    ]);
                }
                println!("\n{} Residues", region);
                residue_table.printstd();
            }
        }

        if verbosity == 2 {
            let atom_list = match region {
                "QM1" => qm1_atom_list,
                "QM2" => qm2_atom_list,
                "Active" => active_atom_list,
                &_ => return Err("No region given.".into()),
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

                for atom in atom_list {
                    atom_table.add_row(row![
                        atom.atom_id,
                        atom.atom_name,
                        atom.res_id,
                        atom.res_name,
                        atom.qm,
                        atom.active
                    ]);
                }
                println!("\n{} Atoms", region);
                atom_table.printstd();
            }
        }
        Ok(())
    }
}

#[derive(Debug)]
pub struct Residue {
    res_name: String,
    res_id: usize,
    atoms: Vec<Atom>,
}

// Implemented in order to be able to perform .contains() on a vector of Residues
// which uses a == comparison under the hood (I think)
impl std::cmp::PartialEq for Residue {
    fn eq(&self, other: &Residue) -> bool {
        self.res_id == other.res_id && self.res_name == other.res_name && self.atoms == other.atoms
    }
}

#[derive(Debug)]
pub struct Atom {
    // Type chosen so IDs can be used for slicing.
    atom_id: usize,
    atom_name: String,
    res_name: String,
    res_id: usize,
    coords: [f64; 3],
    qm: f64,
    active: f64,
    element: String,
}

// Implemented in order to be able to perform .contains() on a vector of Atoms
// which uses a == comparison under the hood (I think)
impl std::cmp::PartialEq for Atom {
    fn eq(&self, other: &Atom) -> bool {
        self.atom_id == other.atom_id
            && self.atom_name == other.atom_name
            && self.res_id == other.res_id
            && self.res_name == other.res_name
            && self.coords == other.coords
            && self.element == other.element
    }
}

impl Clone for Atom {
    fn clone(&self) -> Atom {
        Atom {
            atom_id: self.atom_id,
            atom_name: self.atom_name.clone(),
            res_name: self.res_name.clone(),
            res_id: self.res_id,
            coords: self.coords,
            qm: self.qm,
            active: self.active,
            element: self.element.clone(),
        }
    }
}

impl Atom {
    /// Returns the distance between two instances of the Atom struct
    fn calc_distance(&self, other: &Atom) -> f64 {
        let pos_a = &self.coords;
        let pos_b = other.coords;
        let dist = [
            pos_a[0] - pos_b[0],
            pos_a[1] - pos_b[1],
            pos_a[2] - pos_b[2],
        ];
        let dist_sq = dist[0].powf(2.0) + dist[1].powf(2.0) + dist[2].powf(2.0);
        dist_sq.sqrt()
    }

    // Returns Some(distance) if the distance between the given atoms is smaller than
    // a given value. This value is controlled by the 'level' argument and is supposed
    // to find LJ clashes (i.e. really close atoms) close contacts (e.g. bonds that are not
    // within a residue) or overlap in atomic radii as defined by: https://doi.org/10.1002/chem.201602949
    // fn calc_contact(&self, other: &Atom, level: u8) -> Option<f64> {
    //     let mut vdw_radii = HashMap::new();
    //     vdw_radii.insert("H".to_string(), 1.54);
    //     vdw_radii.insert("C".to_string(), 1.90);
    //     vdw_radii.insert("N".to_string(), 1.79);
    //     vdw_radii.insert("O".to_string(), 1.71);
    //     vdw_radii.insert("P".to_string(), 2.23);
    //     vdw_radii.insert("S".to_string(), 2.14);
    //     vdw_radii.insert("Cl".to_string(), 2.06);
    //     vdw_radii.insert("Na".to_string(), 2.25);
    //     vdw_radii.insert("Cu".to_string(), 2.17);

    //     let dist = self.calc_distance(other);
    //     let radius = match level {
    //         0 => 1.0,
    //         1 => *vdw_radii
    //             .get(&self.element)
    //             .expect("No Radius found for given element."),
    //         2 => {
    //             vdw_radii.get(&self.element).unwrap()
    //                 + vdw_radii
    //                     .get(&other.element)
    //                     .expect("No Radius found for given element.")
    //         }
    //         _ => panic!("Too high level given"),
    //     };

    //     if dist < radius && self != other {
    //         return Some(dist);
    //     } else {
    //         return None;
    //     }
    // }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_pdb() -> Molecule {
        let pdb = Molecule::from_pdb("tests/test_blank.pdb");
        pdb
    }

    #[test]
    fn atom_sphere_inclusive() {
        let pdb = test_pdb();
        // let origin = pdb.atom(26).unwrap();
        let origin = pdb.atom_from_id(26);

        // assert_eq!(calc_atom_sphere(&pdb, origin, 4.0, true).unwrap().len(), 22);
        assert_eq!(pdb.calc_atom_sphere(origin, 4.0, true), 22);
    }
}