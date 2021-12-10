use strum_macros::Display;

use crate::functions;
use crate::options;

type AtomList = Vec<usize>;

#[derive(Debug, Display, Clone, PartialEq, PartialOrd)]
pub enum EditOp {
    ToAdd {
        region: options::Region,
        atoms: AtomList,
    },
    ToRemove {
        region: options::Region,
        atoms: AtomList,
    },
}

pub trait Revertable: std::fmt::Debug {
    fn undo(&self, pdb: &mut pdbtbx::PDB);
    fn redo(&self, pdb: &mut pdbtbx::PDB);
}

impl Revertable for EditOp {
    fn undo(&self, pdb: &mut pdbtbx::PDB) {
        match self {
            EditOp::ToAdd { region, atoms } => functions::edit_atoms(pdb, atoms, "Remove", *region),
            EditOp::ToRemove { region, atoms } => functions::edit_atoms(pdb, atoms, "Add", *region),
        }
    }

    fn redo(&self, pdb: &mut pdbtbx::PDB) {
        match self {
            EditOp::ToAdd { region, atoms } => functions::edit_atoms(pdb, atoms, "Add", *region),
            EditOp::ToRemove { region, atoms } => {
                functions::edit_atoms(pdb, atoms, "Remove", *region)
            }
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

// impl core::fmt::Debug for dyn Revertable {
//     fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
//         write!(f, "{:?}", self)
//     }
// }
