use clap::{App, AppSettings, Arg, ArgGroup};

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
            //         .multiple_values(true)
            // )
            .arg(Arg::new("Input")
                .help("Input for atom selection")
                .required(true)
                .multiple_values(true))
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
            // .group(
            //     ArgGroup::new("source")
            //         .args(&["Infile", "List"])
            //         .required(true))
        )
        .subcommand(App::new("Remove")
            .about("Remove atoms or residues to QM or active regions")
            .visible_aliases(&["remove"])
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
            //         .multiple_values(true)
            // )
            .arg(Arg::new("Input")
                .help("Input for atom selection")
                // .required(true)
                .multiple_values(true))
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
            // .group(
            //     ArgGroup::new("source")
            //         .args(&["Infile", "List"])
            //         .requires_all(&["region"]))
        )
        .subcommand(App::new("Write")
            .about("Write PDB information")
            .visible_aliases(&["write"])
            .arg(
                Arg::new("State")
                    .help("Write state of system")
                    .long("state")
                    .short('s')
            )
            // .arg(
            //     Arg::new("QM1")
            //         .help("Write QM1 region list")
            //         .long("qm1")
            //         .short('q')
            // )
            // .arg(
            //     Arg::new("QM2")
            //         .help("Write QM2 region list")
            //         .long("qm2")
            //         .short('o')
            // )
            // .arg(
            //     Arg::new("Active")
            //         .help("Write active region list")
            //         .long("active")
            //         .short('a')
            // )
            // .arg(
            //     Arg::new("Atoms")
            //     .help("Write atom list")
            //     .long("atoms")
            //     .short('t')
            // )
            // .arg(
            //     Arg::new("Residues")
            //     .help("Write residue list")
            //     .long("residues")
            //     .short('r')
            // )
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
            // .group(
            //     ArgGroup::new("region")
            //         .args(&["QM1", "QM2", "Active"])
            //         .requires("target")
            // )
            // .group(
            //     ArgGroup::new("target")
            //         .args(&["Atoms", "Residues"])
            //         .requires("region")
            // )
        )
}
