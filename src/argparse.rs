use clap::ArgMatches;
use clap::{App, AppSettings, Arg, ArgGroup};
use itertools::Itertools;
use std::error::Error;
use std::str::FromStr;

use lazy_static::lazy_static;
use regex::Regex;
use strum::VariantNames;
use strum_macros::{Display, EnumString, EnumVariantNames};

// type Result<T> = std::result::Result<T, Box<dyn Error>>;

#[derive(Display, PartialEq, Clone, Debug, EnumVariantNames)]
pub enum Mode {
    Query {
        source: Source,
        target: Target,
    },
    Analyze {
        region: Region,
        target: Target,
        distance: Distance,
    },
    Add {
        region: Region,
        source: Source,
        target: Target,
        partial: Partial,
        output: Output,
    },
    Remove {
        region: Region,
        source: Source,
        target: Target,
        partial: Partial,
        output: Output,
    },
    // None,
}

#[derive(Display, PartialEq, Debug, Clone, Copy, EnumString, EnumVariantNames)]
pub enum Region {
    QM1,
    QM2,
    Active,
    None,
}

#[derive(Display, PartialEq, Debug, Clone, Copy, EnumString, EnumVariantNames)]
pub enum Source {
    Infile,
    List,
    Sphere,
    None,
}

#[derive(Display, PartialEq, Debug, Clone, Copy, EnumString, EnumVariantNames)]
pub enum Target {
    Atoms,
    Residues,
    None,
}

#[derive(Display, PartialEq, Debug, Clone, Copy, EnumString, EnumVariantNames)]
pub enum Output {
    Outfile,
    Overwrite,
    None,
}

#[derive(Display, PartialEq, Debug, Clone, Copy, EnumString, EnumVariantNames)]
pub enum Partial {
    Sidechain,
    Backbone,
    None,
}

#[derive(Display, PartialEq, Debug, Clone, Copy, EnumString, EnumVariantNames)]
pub enum Distance {
    Clashes,
    Contacts,
    None,
}

impl Mode {
    /// Creates new Mode enum from clap::ArgMatches struct. This is
    /// where the given command line options are stored for later use.
    pub fn new(matches: &clap::ArgMatches) -> Result<Mode, Box<dyn Error>> {
        match matches.subcommand_name() {
            Some("Query") => {
                let source_str = Source::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Query").unwrap().is_present(x))
                    .unwrap_or(&"None");
                let source = Source::from_str(source_str)?;

                let target_str = Target::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Query").unwrap().is_present(x))
                    .unwrap_or(&"None");
                let target = Target::from_str(target_str)?;

                Ok(Mode::Query { source, target })
            }
            Some("Analyze") => {
                let region_str = Region::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Analyze").unwrap().is_present(x))
                    .unwrap_or(&"None");
                let region = Region::from_str(region_str)?;

                let target_str = Target::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Analyze").unwrap().is_present(x))
                    .unwrap_or(&"None");
                let target = Target::from_str(target_str)?;

                let distance_str = Distance::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Analyze").unwrap().is_present(x))
                    .unwrap_or(&"None");
                let distance = Distance::from_str(distance_str)?;

                Ok(Mode::Analyze {
                    region,
                    target,
                    distance,
                })
            }
            Some("Add") => {
                let region_str = Region::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Add").unwrap().is_present(x))
                    .unwrap_or(&"None");
                let region = Region::from_str(region_str)?;

                let source_str = Source::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Add").unwrap().is_present(x))
                    .unwrap_or(&"None");
                let source = Source::from_str(source_str)?;

                let target_str = Target::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Add").unwrap().is_present(x))
                    .unwrap_or(&"None");
                let target = Target::from_str(target_str)?;

                let partial_str = Partial::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Add").unwrap().is_present(x))
                    .unwrap_or(&"None");
                let partial = Partial::from_str(partial_str)?;

                let output_str = Output::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Add").unwrap().is_present(x))
                    .unwrap_or(&"None");
                let output = Output::from_str(output_str)?;

                Ok(Mode::Add {
                    region,
                    source,
                    target,
                    partial,
                    output,
                })
            }
            Some("Remove") => {
                let region_str = Region::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Remove").unwrap().is_present(x))
                    .unwrap_or(&"None");
                let region = Region::from_str(region_str)?;

                let source_str = Source::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Remove").unwrap().is_present(x))
                    .unwrap_or(&"None");
                let source = Source::from_str(source_str)?;

                let target_str = Target::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Remove").unwrap().is_present(x))
                    .unwrap_or(&"None");
                let target = Target::from_str(target_str)?;

                let partial_str = Partial::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Remove").unwrap().is_present(x))
                    .unwrap_or(&"None");
                let partial = Partial::from_str(partial_str)?;

                let output_str = Output::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Remove").unwrap().is_present(x))
                    .unwrap_or(&"None");
                let output = Output::from_str(output_str)?;

                Ok(Mode::Remove {
                    region,
                    source,
                    target,
                    partial,
                    output,
                })
            }
            Some(&_) => unreachable!(),
            None => unreachable!(),
        }
    }
}

fn sphere_valid(v: &str) -> Result<(), String> {
    lazy_static! {
        static ref RE: Regex = Regex::new(r"[\d.]").unwrap();
    }
    // let re = Regex::new(r"[\d\.]")?;

    let err_chars = v
        .chars()
        .filter(|x| !RE.is_match(&x.to_string()))
        .collect::<String>();

    if !err_chars.is_empty() {
        return Err(format!("Invalid characters: {}", err_chars).into());
    }

    Ok(())
}

fn list_valid(v: &str) -> Result<(), Box<dyn Error>> {
    lazy_static! {
        static ref RE_NUM: Regex = Regex::new(
            r"^(?P<id1>\d+)(?P<insert1>[A-Za-z]?)([:-](?P<id2>\d+)(?P<insert2>[A-Za-z]?))?$"
        )
        .unwrap();
    }
    lazy_static! {
        static ref RE_STR: Regex = Regex::new(r"^[A-Za-z]+$").unwrap();
    }
    lazy_static! {
        static ref RE_CHARS: Regex = Regex::new(r"[\dA-Z-a-z:,-]").unwrap();
    }

    let err_chars = v
        .chars()
        .filter(|x| !RE_CHARS.is_match(&x.to_string()))
        .sorted()
        .dedup()
        .collect::<String>();

    if !err_chars.is_empty() {
        return Err(format!("Invalid characters: {:?}", err_chars).into());
    }

    let mut numerical_inp = false;
    let mut string_inp = false;

    for i in v.split(',') {
        if RE_NUM.is_match(i) {
            numerical_inp = true;

            let caps = RE_NUM.captures(i).unwrap();
            if caps.name("id2").is_some()
                && caps.name("id1").unwrap().as_str().parse::<i32>()?
                    > caps.name("id2").unwrap().as_str().parse::<i32>()?
                && caps.name("insert1").unwrap().as_str() >= caps.name("insert2").unwrap().as_str()
            {
                return Err(format!(
                        "Invalid range given: {}{}-{}{}. Left entry must preceed right one in PDB file!",
                        caps.name("id1").unwrap().as_str(),
                        caps.name("insert1").unwrap().as_str(),
                        caps.name("id2").unwrap().as_str(),
                        caps.name("insert2").unwrap().as_str(),
                    )
                    .into());
            }
        } else if RE_STR.is_match(i) {
            string_inp = true;
        }
    }

    if numerical_inp && string_inp {
        return Err("Input List contains mixed types.".into());
    }

    Ok(())
}

/// Defines all Args, their configuration and all ArgGroups as defined by clap.
pub fn parse_args() -> ArgMatches {
    // let matches = App::new(crate_name!())
    App::new(crate_name!())
        .about(crate_description!())
        .version(crate_version!())
        .author(crate_authors!())
        .setting(AppSettings::ArgRequiredElseHelp)
        .setting(AppSettings::SubcommandRequiredElseHelp)
        .setting(AppSettings::VersionlessSubcommands)
        .arg(Arg::new("INPUT").about("Path to PDB file").value_name("PDB file").required(true))
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
                    // .validator_regex(&re, "Invalid characters for list input")
                    // .multiple(true)
                    // .use_delimiter(true),
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
                    // .multiple(true)
                    // .use_delimiter(true),
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
            .arg(
                Arg::new("Outfile")
                    .about("Name of PDB output file")
                    .long("outfile")
                    .short('e')
                    .takes_value(true),
            )
            .arg(
                Arg::new("Overwrite")
                    .about("Overwrite input PDB file")
                    .long("overwrite")
                    .short('w'),
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
            .group(
                ArgGroup::new("output")
                    .args(&["Outfile", "Overwrite"]))
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
                    // .multiple(true)
                    // .use_delimiter(true),
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
            .arg(
                Arg::new("Outfile")
                    .about("Name of PDB output file")
                    .long("outfile")
                    .short('e')
                    .takes_value(true),
            )
            .arg(
                Arg::new("Overwrite")
                    .about("Overwrite input PDB file")
                    .long("overwrite")
                    .short('w'),
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
                    .requires_all(&["target", "source"]))
            .group(
                ArgGroup::new("source")
                    .args(&["Infile", "List", "Sphere"])
                    .requires_all(&["target", "region"]))
            .group(
                ArgGroup::new("output")
                    .args(&["Outfile", "Overwrite"]))
        )
    .get_matches()
    // Ok(matches)
}
