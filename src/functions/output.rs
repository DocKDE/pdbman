use std::io;
use std::io::prelude::Write;

use anyhow::{Context, Result};
use pdbtbx::{
    ContainsAtomConformer, ContainsAtomConformerResidue, ContainsAtomConformerResidueChain, PDB,
};

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
