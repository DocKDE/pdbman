use crate::functions;
use crate::options;

type AtomList = Vec<usize>;
type ResidueList = Vec<isize>;

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub enum OpTarget {
    Atoms(AtomList),
    Residues(ResidueList),
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub enum EditOp {
    ToAdd {
        target: OpTarget,
        region: options::Region,
    },
    ToRemove {
        target: OpTarget,
        region: options::Region,
    },
}

pub trait Revertable {
    fn undo(&self, pdb: &mut pdbtbx::PDB);
    fn redo(&self, pdb: &mut pdbtbx::PDB);
}

impl Revertable for EditOp {
    fn undo(&self, pdb: &mut pdbtbx::PDB) {
        match self {
            EditOp::ToAdd { target, region } => match target {
                OpTarget::Atoms(list) => {
                    functions::edit_atoms(pdb, list, "Remove", region.to_owned())
                }
                OpTarget::Residues(list) => {
                    functions::edit_residues(pdb, list, "Remove", None, region.to_owned())
                }
            },
            EditOp::ToRemove { target, region } => match target {
                OpTarget::Atoms(list) => functions::edit_atoms(pdb, list, "Add", region.to_owned()),
                OpTarget::Residues(list) => {
                    functions::edit_residues(pdb, list, "Add", None, region.to_owned())
                }
            },
        }
    }
    fn redo(&self, pdb: &mut pdbtbx::PDB) {
        match self {
            EditOp::ToAdd { target, region } => match target {
                OpTarget::Atoms(list) => functions::edit_atoms(pdb, list, "Add", region.to_owned()),
                OpTarget::Residues(list) => {
                    functions::edit_residues(pdb, list, "Add", None, region.to_owned())
                }
            },
            EditOp::ToRemove { target, region } => match target {
                OpTarget::Atoms(list) => {
                    functions::edit_atoms(pdb, list, "Remove", region.to_owned())
                }
                OpTarget::Residues(list) => {
                    functions::edit_residues(pdb, list, "Remove", None, region.to_owned())
                }
            },
        }
    }
}

impl Revertable for Vec<EditOp> {
    fn undo(&self, pdb: &mut pdbtbx::PDB) {
        for i in self {
            i.undo(pdb)
        }
    }
    fn redo(&self, pdb: &mut pdbtbx::PDB) {
        for i in self {
            i.redo(pdb)
        }
    }
}
