// use clap::{ArgEnum, Parser};

// #[derive(ArgEnum, Parser, Debug, Clone, Copy)]
// #[clap(name = "")]
// pub enum Mode2 {
//     Query {
//         source: Source,
//         target: Target,
//     },
//     Analyze {
//         region: Region,
//         target: Target,
//         distance: Distance,
//     },
//     Add {
//         region: Region,
//         source: Source,
//         target: Target,
//         partial: Partial,
//         // output: Output,
//     },
//     Remove {
//         region: Region,
//         source: Source,
//         target: Target,
//         partial: Partial,
//         // output: Output,
//     },
//     // Interactive,
// }

// #[derive(Parser, Debug, Copy, Clone)]
// enum Region {
//     QM1,
//     QM2,
//     Active,
//     None,
// }

// #[derive(Parser, Debug, Copy, Clone)]
// enum Source {
//     Infile,
//     List,
//     Sphere,
//     None,
// }

// #[derive(Parser, Debug, Copy, Clone)]
// enum Target {
//     Atoms,
//     Residues,
//     None,
// }

// #[derive(Parser, Debug, Copy, Clone)]
// enum Partial {
//     Sidechain,
//     Backbone,
//     None,
// }

// #[derive(Parser, Debug, Copy, Clone)]
// enum Distance {
//     Clashes,
//     Contacts,
//     None,
// }