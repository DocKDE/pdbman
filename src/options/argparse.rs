use anyhow::Result;
use clap::{App, AppSettings, Arg, ArgGroup};
use itertools::Itertools;
use lazy_regex::{regex, regex_is_match};

fn sphere_valid(v: &str) -> Result<(), anyhow::Error> {
    let err_chars = v
        .chars()
        .filter(|x| !regex_is_match!(r"[\d.]", &x.to_string()))
        .collect::<String>();

    ensure!(
        err_chars.is_empty(),
        "\nInvalid characters: '{}'",
        err_chars
    );

    Ok(())
}

fn list_valid(v: &str) -> Result<(), anyhow::Error> {
    let re_num = regex!(r"^(?P<id1>\d+)([:-](?P<id2>\d+))?$");
    let re_str = regex!(r"^[A-Za-z]+$");
    let re_chars = regex!(r"[\dA-Z-a-z:,-]");

    let mut err_chars = v
        .chars()
        .filter(|x| !re_chars.is_match(&x.to_string()))
        .sorted()
        .dedup()
        .peekable();

    ensure!(
        err_chars.peek().is_none(),
        "\nInvalid character(s) found: '{}'",
        err_chars.join("")
    );

    let mut numerical_inp = false;
    let mut string_inp = false;

    for i in v.split(',') {
        if re_num.is_match(i) {
            numerical_inp = true;

            let caps = re_num.captures(i).unwrap();
            let id1 = caps.name("id1").unwrap().as_str().parse::<i32>()?;
            let id2 = caps.name("id2").map(|m| m.as_str().parse::<i32>().unwrap());

            if let Some(id2) = id2 {
                ensure!(id1 < id2, "\nLeft number must be lower: '{}-{}'", id1, id2)
            }
        } else if re_str.is_match(i) {
            string_inp = true;
        } else {
            bail!("\nInvalid range definition found: '{}'", i)
        }
    }

    ensure!(
        !(numerical_inp && string_inp),
        "\nInput list containts mixed types!"
    );

    Ok(())
}

/// Defines all Args, their configuration and all ArgGroups as defined by clap.
pub fn parse_args() -> App<'static> {
    App::new("")
        .setting(AppSettings::ArgRequiredElseHelp)
        .setting(AppSettings::SubcommandRequiredElseHelp)
        // .setting(AppSettings::VersionlessSubcommands)
        .setting(AppSettings::NoBinaryName)
        .subcommand(App::new("Query")
            .about("Query mode")
            .visible_aliases(&["query", "que", "Q", "q"])
            .arg(
                Arg::new("Residues")
                    .about("Residue Mode")
                    .long("residues")
                    .short('r')
            )
            .arg(
                Arg::new("Atoms")
                    .about("Atom Mode")
                    .long("atoms")
                    .short('t')
            )
            .arg(
                Arg::new("Infile")
                    .about("File for list input")
                    .long("file")
                    .short('f')
                    .takes_value(true),
            )
            .arg(
                Arg::new("List")
                    .about("Command line list")
                    .long("list")
                    .short('l')
                    .takes_value(true)
                    .validator(list_valid)
            )
            .arg(
                Arg::new("Sphere")
                    .about("Calculate sphere around atom. Requires Atom ID and radius in Angstrom.")
                    .long("sphere")
                    .short('s')
                    .takes_value(true)
                    .number_of_values(2)
                    .value_names(&["Atom ID", "Radius"])
                    .conflicts_with_all(&["Sidechain", "Backbone"])
                    .validator(sphere_valid),
            )
            .group(ArgGroup::new("target").args(&["Residues", "Atoms"]).required(true))
            .group(ArgGroup::new("source").args(&["Infile", "List", "Sphere"]).required(true))
        )
        .subcommand(App::new("Analyze")
            .about("Analysis mode")
            .visible_aliases(&["analyze", "ana", "Y", "y"])
            .arg(
                Arg::new("Residues")
                    .about("Residue Mode")
                    .long("residues")
                    .short('r')
            )
            .arg(
                Arg::new("Atoms")
                    .about("Atom Mode")
                    .long("atoms")
                    .short('t')
            )
            .arg(
                Arg::new("QM1")
                    .about("QM1 region")
                    .long("qm1")
                    .short('q')
            )
            .arg(
                Arg::new("QM2")
                    .about("QM2 region")
                    .long("qm2")
                    .short('o')
            )
            .arg(
                Arg::new("Active")
                    .about("Active region")
                    .long("active")
                    .short('a')
            )
            .arg(
                Arg::new("Clashes")
                    .about("Find LJ clashes")
                    .long("clashes")
                    .short('c')
            )
            .arg(
                Arg::new("Contacts")
                    .about("Find vdW contacts")
                    .long("contacts")
                    .short('n')
            )
            .group(ArgGroup::new("target").args(&["Residues", "Atoms"]).requires("region"))
            .group(ArgGroup::new("region").args(&["QM1", "QM2", "Active"]).requires("target"))
            .group(ArgGroup::new("distances").args(&["Clashes", "Contacts"]))
        )
        .subcommand(App::new("Add")
            .about("Add mode")
            .visible_aliases(&["add", "A", "a"])
            .arg(
                Arg::new("Residues")
                    .about("Residue Mode")
                    .long("residues")
                    .short('r')
            )
            .arg(
                Arg::new("Atoms")
                    .about("Atom Mode")
                    .long("atoms")
                    .short('t')
            )
            .arg(
                Arg::new("Infile")
                    .about("File for list input")
                    .long("file")
                    .short('f')
                    .takes_value(true),
            )
            .arg(
                Arg::new("List")
                    .about("Command line list")
                    .long("list")
                    .short('l')
                    .takes_value(true)
                    .validator(list_valid)
            )
            .arg(
                Arg::new("Sphere")
                    .about("Calculate sphere around atom. Requires Atom ID and radius in Angstrom.")
                    .long("sphere")
                    .short('s')
                    .takes_value(true)
                    .number_of_values(2)
                    .value_names(&["Atom ID", "Radius"])
                    .conflicts_with_all(&["Sidechain", "Backbone"])
                    .validator(sphere_valid),
            )
            .arg(
                Arg::new("QM1")
                    .about("QM1 region")
                    .long("qm1")
                    .short('q')
            )
            .arg(
                Arg::new("QM2")
                    .about("QM2 region")
                    .long("qm2")
                    .short('o')
            )
            .arg(
                Arg::new("Active")
                    .about("Active region")
                    .long("active")
                    .short('a')
            )
            .arg(
                Arg::new("Sidechain")
                    .about("Use only sidechains of protein")
                    .long("sidechain")
                    .short('d')
            )
            .arg(
                Arg::new("Backbone")
                    .about("Use only backbones of protein")
                    .long("backbone")
                    .short('b')
            )
            .group(
                ArgGroup::new("target")
                    .args(&["Residues", "Atoms"])
                    .required(true))
            .group(
                ArgGroup::new("partial")
                    .args(&["Sidechain", "Backbone"])
                    .requires("Residues")
            )
            .group(
                ArgGroup::new("region")
                    .args(&["QM1", "QM2", "Active"])
                    .required(true))
            .group(
                ArgGroup::new("source")
                    .args(&["Infile", "List", "Sphere"])
                    .required(true))
        )
        .subcommand(App::new("Remove")
            .about("Remove mode")
            .visible_aliases(&["remove", "rem", "R", "r"])
            .arg(
                Arg::new("Residues")
                    .about("Residue Mode")
                    .long("residues")
                    .short('r')
            )
            .arg(
                Arg::new("Atoms")
                    .about("Atom Mode")
                    .long("atoms")
                    .short('t')
            )
            .arg(
                Arg::new("Infile")
                    .about("File for list input")
                    .long("file")
                    .short('f')
                    .takes_value(true),
            )
            .arg(
                Arg::new("List")
                    .about("Command line list")
                    .long("list")
                    .short('l')
                    .takes_value(true)
                    .validator(list_valid)
            )
            .arg(
                Arg::new("Sphere")
                    .about("Calculate sphere around atom. Requires Atom ID and radius in Angstrom.")
                    .long("sphere")
                    .short('s')
                    .takes_value(true)
                    .number_of_values(2)
                    .value_names(&["Atom ID", "Radius"])
                    .conflicts_with_all(&["Sidechain", "Backbone"])
                    .validator(sphere_valid),
            )
            .arg(
                Arg::new("QM1")
                    .about("QM1 region")
                    .long("qm1")
                    .short('q')
            )
            .arg(
                Arg::new("QM2")
                    .about("QM2 region")
                    .long("qm2")
                    .short('o')
            )
            .arg(
                Arg::new("Active")
                    .about("Active region")
                    .long("active")
                    .short('a')
            )
            .arg(
                Arg::new("Sidechain")
                    .about("Use only sidechains of protein")
                    .long("sidechain")
                    .short('d')
            )
            .arg(
                Arg::new("Backbone")
                    .about("Use only backbones of protein")
                    .long("backbone")
                    .short('b')
            )
            .group(
                ArgGroup::new("target")
                    .args(&["Residues", "Atoms"])
                    .requires_all(&["region", "source"]))
            .group(
                ArgGroup::new("partial")
                    .args(&["Sidechain", "Backbone"])
                    .requires("Residues")
                    .requires_all(&["target", "region", "source"]),
            )
            .group(
                ArgGroup::new("region")
                    .args(&["QM1", "QM2", "Active"])
                    // .requires_all(&["target", "source"])
                )
            .group(
                ArgGroup::new("source")
                    .args(&["Infile", "List", "Sphere"])
                    .requires_all(&["target", "region"]))
        )
        .subcommand(App::new("Write")
            .about("Write Mode")
            .visible_aliases(&["write", "wri", "W", "w"])
            .arg(
                Arg::new("QM1")
                    .about("Write QM1 atom list")
                    .long("qm1")
                    .short('q')
            )
            .arg(
                Arg::new("QM2")
                    .about("Write QM2 atom list")
                    .long("qm2")
                    .short('o')
            )
            .arg(
                Arg::new("Active")
                    .about("Write active atom list")
                    .long("active")
                    .short('a')
            )
            .arg(
                Arg::new("Outfile")
                    .about("File path for writing PDB")
                    .long("file")
                    .short('f')
                    .takes_value(true),
            )
            .arg(
                Arg::new("Overwrite")
                    .about("Overwrite PDB input file")
                    .long("overwrite")
                    .short('w')
                    .conflicts_with_all(&["QM1", "QM2", "Active"])
            )
            .group(
                ArgGroup::new("output")
                    .args(&["Outfile", "Overwrite"])
            )
            .group(
                ArgGroup::new("region")
                    .args(&["QM1", "QM2", "Active"])
            )
        )
}
