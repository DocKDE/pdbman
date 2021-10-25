// use clap::Parser;

// #[derive(Parser, Debug)]
// #[clap(name = "pdbman")]
// enum Mode {
//     #[derive(Subcommand, Debug)]
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
//     None,
// }

// #[derive(Debug)]
// enum Region {
//     QM1,
//     QM2,
//     Active,
//     None,
// }

// #[derive(Debug)]
// enum Source {
//     Infile,
//     List,
//     Sphere,
//     None,
// }

// #[derive(Debug)]
// enum Target {
//     Atoms,
//     Residues,
//     None,
// }

// #[derive(Debug)]
// enum Partial {
//     Sidechain,
//     Backbone,
//     None,
// }

// #[derive(Debug)]
// enum Distance {
//     Clashes,
//     Contacts,
//     None,
// }