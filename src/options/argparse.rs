// use anyhow::Result;
use clap::{App, AppSettings, Arg, ArgGroup};
// use itertools::Itertools;
// use lazy_regex::{regex, regex_is_match};

// use crate::HELP_LONG;

// fn sphere_valid(v: &str) -> Result<(), anyhow::Error> {
//     let err_chars = v
//         .chars()
//         .filter(|x| !regex_is_match!(r"[\d.]", &x.to_string()))
//         .collect::<String>();

//     ensure!(
//         err_chars.is_empty(),
//         "\nInvalid characters: '{}'",
//         err_chars
//     );

//     Ok(())
// }

// fn list_valid(v: &str) -> Result<(), anyhow::Error> {
//     let re_num = regex!(r"^(?P<id1>\d+)([:-](?P<id2>\d+))?$");
//     let re_str = regex!(r"^[A-Za-z]+$");
//     let re_chars = regex!(r"[\dA-Z-a-z:,-]");

//     let mut err_chars = v
//         .chars()
//         .filter(|x| !re_chars.is_match(&x.to_string()))
//         .sorted()
//         .dedup()
//         .peekable();

//     ensure!(
//         err_chars.peek().is_none(),
//         "\nInvalid character(s) found: '{}'",
//         err_chars.join("")
//     );

//     let mut numerical_inp = false;
//     let mut string_inp = false;

//     for i in v.split(',') {
//         if re_num.is_match(i) {
//             numerical_inp = true;

//             let caps = re_num.captures(i).unwrap();
//             let id1 = caps.name("id1").unwrap().as_str().parse::<i32>()?;
//             let id2 = caps.name("id2").map(|m| m.as_str().parse::<i32>().unwrap());

//             if let Some(id2) = id2 {
//                 ensure!(id1 < id2, "\nLeft number must be lower: '{}-{}'", id1, id2)
//             }
//         } else if re_str.is_match(i) {
//             string_inp = true;
//         } else {
//             bail!("\nInvalid range definition found: '{}'", i)
//         }
//     }

//     ensure!(
//         !(numerical_inp && string_inp),
//         "\nInput list containts mixed types!"
//     );

//     Ok(())
// }

/// Defines all Args, their configuration and all ArgGroups as defined by clap.
pub fn clap_args() -> App<'static> {
    App::new("")
        // .setting(AppSettings::ArgRequiredElseHelp)
        // .setting(AppSettings::SubcommandRequiredElseHelp)
        // .setting(AppSettings::VersionlessSubcommands)
        .setting(AppSettings::InferSubcommands)
        .setting(AppSettings::NoBinaryName)
        // .setting(AppSettings::NoAutoHelp)
        .subcommand(App::new("Query")
            .about("Query for atoms and residues")
            .visible_aliases(&["query"])
            .arg(Arg::new("Input")
                .help("Input for atom selection")
                .required(true)
                .multiple_values(true))
            // .arg(
            //     Arg::new("Residues")
            //         .help("Residue Mode")
            //         .long("residues")
            //         .short('r')
            // )
            // .arg(
            //     Arg::new("Atoms")
            //         .help("Atom Mode")
            //         .long("atoms")
            //         .short('t')
            // )
            // .arg(
            //     Arg::new("Infile")
            //         .help("File for list input")
            //         .long("file")
            //         .short('f')
            //         .takes_value(true),
            // )
            // .arg(
            //     Arg::new("List")
            //         .help("Command line list")
            //         .long("list")
            //         .short('l')
            //         .takes_value(true)
            //         .validator(list_valid)
            // )
            // .arg(
            //     Arg::new("Sphere")
            //         .help("Calculate sphere around atom. Requires Atom ID and radius in Angstrom.")
            //         .long("sphere")
            //         .short('s')
            //         .takes_value(true)
            //         .number_of_values(2)
            //         .value_names(&["Atom ID", "Radius"])
            //         .conflicts_with_all(&["Sidechain", "Backbone"])
            //         .validator(sphere_valid),
            // )
            // .group(ArgGroup::new("target").args(&["Residues", "Atoms"]).required(true))
            // .group(ArgGroup::new("source").args(&["Infile", "List", "Sphere"]).required(true))
        )
        .subcommand(App::new("Analyze")
            .about("Analysis QM and active regions")
            .visible_aliases(&["analyze", "Y", "y"])
            .arg(
                Arg::new("Residues")
                    .help("Residue Mode")
                    .long("residues")
                    .short('r')
            )
            .arg(
                Arg::new("Atoms")
                    .help("Atom Mode")
                    .long("atoms")
                    .short('t')
            )
            .arg(
                Arg::new("QM1")
                    .help("QM1 region")
                    .long("qm1")
                    .short('q')
            )
            .arg(
                Arg::new("QM2")
                    .help("QM2 region")
                    .long("qm2")
                    .short('o')
            )
            .arg(
                Arg::new("Active")
                    .help("Active region")
                    .long("active")
                    .short('a')
            )
            .arg(
                Arg::new("Clashes")
                    .help("Find LJ clashes")
                    .long("clashes")
                    .short('c')
            )
            .arg(
                Arg::new("Contacts")
                    .help("Find vdW contacts")
                    .long("contacts")
                    .short('n')
            )
            .group(ArgGroup::new("target").args(&["Residues", "Atoms"]).requires("region"))
            .group(ArgGroup::new("region").args(&["QM1", "QM2", "Active"]).requires("target"))
            .group(ArgGroup::new("distances").args(&["Clashes", "Contacts"]))
        )
        .subcommand(App::new("Add")
            .about("Add atoms or residues to QM or active regions")
            .visible_aliases(&["add", "A", "a"])
            // .arg(
            //     Arg::new("Residues")
            //         .help("Residue Mode")
            //         .long("residues")
            //         .short('r')
            // )
            // .arg(
            //     Arg::new("Atoms")
            //         .help("Atom Mode")
            //         .long("atoms")
            //         .short('t')
            // )
            .arg(
                Arg::new("Infile")
                    .help("File for list input")
                    .long("file")
                    .short('f')
                    .takes_value(true),
            )
            .arg(
                Arg::new("List")
                    .help("Command line list")
                    .long("list")
                    .short('l')
                    .takes_value(true)
                    .multiple_values(true)
                    // .validator(list_valid)
            )
            // .arg(
            //     Arg::new("Sphere")
            //         .help("Calculate sphere around atom. Requires Atom ID and radius in Angstrom.")
            //         .long("sphere")
            //         .short('s')
            //         .takes_value(true)
            //         .number_of_values(2)
            //         .value_names(&["Atom ID", "Radius"])
            //         .conflicts_with_all(&["Sidechain", "Backbone"])
            //         .validator(sphere_valid),
            // )
            .arg(
                Arg::new("QM1")
                    .help("QM1 region")
                    .long("qm1")
                    .short('q')
            )
            .arg(
                Arg::new("QM2")
                    .help("QM2 region")
                    .long("qm2")
                    .short('o')
            )
            .arg(
                Arg::new("Active")
                    .help("Active region")
                    .long("active")
                    .short('a')
            )
            .arg(
                Arg::new("Sidechain")
                    .help("Use only sidechains of protein")
                    .long("sidechain")
                    .short('d')
            )
            .arg(
                Arg::new("Backbone")
                    .help("Use only backbones of protein")
                    .long("backbone")
                    .short('b')
            )
            // .group(
            //     ArgGroup::new("target")
            //         .args(&["Residues", "Atoms"])
            //         .required(true))
            .group(
                ArgGroup::new("partial")
                    .args(&["Sidechain", "Backbone"])
                    .requires_all(&["region", "source"]),
                    // .requires("Residues")
            )
            .group(
                ArgGroup::new("region")
                    .args(&["QM1", "QM2", "Active"])
                    .required(true))
            .group(
                ArgGroup::new("source")
                    .args(&["Infile", "List"])
                    .required(true))
        )
        .subcommand(App::new("Remove")
            .about("Remove atoms or residues to QM or active regions")
            .visible_aliases(&["remove"])
            // .arg(
            //     Arg::new("Residues")
            //         .help("Residue Mode")
            //         .long("residues")
            //         .short('r')
            // )
            // .arg(
            //     Arg::new("Atoms")
            //         .help("Atom Mode")
            //         .long("atoms")
            //         .short('t')
            // )
            .arg(
                Arg::new("Infile")
                    .help("File for list input")
                    .long("file")
                    .short('f')
                    .takes_value(true),
            )
            .arg(
                Arg::new("List")
                    .help("Command line list")
                    .long("list")
                    .short('l')
                    .takes_value(true)
                    .multiple_values(true)
                    // .validator(list_valid)
            )
            // .arg(
            //     Arg::new("Sphere")
            //         .help("Calculate sphere around atom. Requires Atom ID and radius in Angstrom.")
            //         .long("sphere")
            //         .short('s')
            //         .takes_value(true)
            //         .number_of_values(2)
            //         .value_names(&["Atom ID", "Radius"])
            //         .conflicts_with_all(&["Sidechain", "Backbone"])
            //         .validator(sphere_valid),
            // )
            .arg(
                Arg::new("QM1")
                    .help("QM1 region")
                    .long("qm1")
                    .short('q')
            )
            .arg(
                Arg::new("QM2")
                    .help("QM2 region")
                    .long("qm2")
                    .short('o')
            )
            .arg(
                Arg::new("Active")
                    .help("Active region")
                    .long("active")
                    .short('a')
            )
            .arg(
                Arg::new("Sidechain")
                    .help("Use only sidechains of protein")
                    .long("sidechain")
                    .short('d')
            )
            .arg(
                Arg::new("Backbone")
                    .help("Use only backbones of protein")
                    .long("backbone")
                    .short('b')
            )
            // .group(
            //     ArgGroup::new("target")
            //         .args(&["Residues", "Atoms"])
            //         .requires_all(&["region", "source"]))
            .group(
                ArgGroup::new("partial")
                    .args(&["Sidechain", "Backbone"])
                    // .requires("Residues")
                    .requires_all(&["region", "source"]),
            )
            .group(
                ArgGroup::new("region")
                    .args(&["QM1", "QM2", "Active"])
                    // .requires_all(&["target", "source"])
                )
            .group(
                ArgGroup::new("source")
                    .args(&["Infile", "List"])
                    .requires_all(&["region"]))
        )
        .subcommand(App::new("Write")
            .about("Write PDB information")
            .visible_aliases(&["write"])
            .arg(
                Arg::new("QM1")
                    .help("Write QM1 region list")
                    .long("qm1")
                    .short('q')
            )
            .arg(
                Arg::new("QM2")
                    .help("Write QM2 region list")
                    .long("qm2")
                    .short('o')
            )
            .arg(
                Arg::new("Active")
                    .help("Write active region list")
                    .long("active")
                    .short('a')
            )
            .arg(
                Arg::new("Atoms")
                .help("Write atom list")
                .long("atoms")
                .short('t')
            )
            .arg(
                Arg::new("Residues")
                .help("Write residue list")
                .long("residues")
                .short('r')
            )
            .arg(
                Arg::new("Outfile")
                    .help("File path for writing output")
                    .long("file")
                    .short('f')
                    .takes_value(true),
            )
            .arg(
                Arg::new("Overwrite")
                    .help("Overwrite PDB input file")
                    .long("overwrite")
                    .short('w')
                    .conflicts_with_all(&["QM1", "QM2", "Active", "Atoms", "Residues"])
            )
            .group(
                ArgGroup::new("output")
                    .args(&["Outfile", "Overwrite"])
            )
            .group(
                ArgGroup::new("region")
                    .args(&["QM1", "QM2", "Active"])
                    .requires("target")
            )
            .group(
                ArgGroup::new("target")
                    .args(&["Atoms", "Residues"])
                    .requires("region")
            )
        )
}
