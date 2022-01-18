use std::str::FromStr;

use anyhow::Result;
use itertools::Itertools;
use strum::VariantNames;
use strum_macros::{Display, EnumString, EnumVariantNames};

#[derive(Display, PartialEq, Debug, Clone, EnumVariantNames)]
pub enum Mode<'a> {
    Query {
        input: String,
    },
    Analyze {
        region: Option<Region>,
        target: Option<Target>,
        distance: Option<Distance>,
    },
    Add {
        region: Option<Region>,
        // source: Option<Source<'a>>,
        partial: Option<Partial>,
        selection: Option<String>,
    },
    Remove {
        region: Option<Region>,
        // source: Option<Source<'a>>,
        partial: Option<Partial>,
        selection: Option<String>,
    },
    Write {
        output: Option<Output<'a>>,
        state: bool,
        // region: Option<Region>,
        // target: Option<Target>,
    },
}

#[derive(Display, PartialEq, Debug, Clone, Copy, PartialOrd, EnumString, EnumVariantNames)]
pub enum Region {
    QM1,
    QM2,
    Active,
}

// #[derive(Display, PartialEq, Debug, Clone, EnumString, EnumVariantNames)]
// pub enum Source<'a> {
//     Infile(&'a str),
//     List(String),
// }

#[derive(Display, PartialEq, Debug, Clone, Copy, PartialOrd, EnumString, EnumVariantNames)]
pub enum Target {
    Atoms,
    Residues,
}

#[derive(Display, PartialEq, Debug, Clone, EnumString, EnumVariantNames)]
pub enum Output<'a> {
    Outfile(&'a str),
    Overwrite,
}

#[derive(Display, PartialEq, Debug, Clone, Copy, EnumString, EnumVariantNames)]
pub enum Partial {
    Sidechain,
    Backbone,
}

#[derive(Display, PartialEq, Debug, Clone, Copy, EnumString, EnumVariantNames)]
pub enum Distance {
    Clashes,
    Contacts,
}

impl<'a> Mode<'a> {
    /// Creates new Mode enum from clap::ArgMatches struct. This is
    /// where the given command line options are stored for later use.
    pub fn new(matches: &clap::ArgMatches) -> Result<Mode, anyhow::Error> {
        match matches.subcommand_name() {
            Some("Query") => Ok(Mode::Query {
                input: matches
                    .subcommand_matches("Query")
                    .unwrap()
                    .values_of("Input")
                    .unwrap()
                    .join(" "),
            }),
            Some("Analyze") => {
                // If no other argument to the subcommand is given, the find method yields None.
                // Thus the subsequent map method may use unwrap, since the presence and
                // validity of potential arguments are ensured.
                let region = Region::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Analyze").unwrap().is_present(x))
                    .map(|s| Region::from_str(s).unwrap());

                let target = Target::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Analyze").unwrap().is_present(x))
                    .map(|s| Target::from_str(s).unwrap());

                let distance = Distance::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Analyze").unwrap().is_present(x))
                    .map(|s| Distance::from_str(s).unwrap());

                Ok(Mode::Analyze {
                    region,
                    target,
                    distance,
                })
            }
            Some("Add") => {
                let region = Region::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Add").unwrap().is_present(x))
                    .map(|s| Region::from_str(s).unwrap());

                // This argument is required for the Add subcommand so unwrap is fine here.
                // let source_str = Source::VARIANTS
                //     .iter()
                //     .find(|x| matches.subcommand_matches("Add").unwrap().is_present(x))
                //     .unwrap();

                // let source = match *source_str {
                //     "Infile" => Some(Source::Infile(
                //         matches
                //             .subcommand_matches("Add")
                //             .unwrap()
                //             .value_of("Infile")
                //             .unwrap(),
                //     )),
                //     "List" => Some(Source::List(
                //         matches
                //             .subcommand_matches("Add")
                //             .unwrap()
                //             .values_of("List")
                //             .unwrap()
                //             .join(" "),
                //     )),
                //     _ => unreachable!(),
                // };

                let partial = Partial::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Add").unwrap().is_present(x))
                    .map(|s| Partial::from_str(s).unwrap());

                Ok(Mode::Add {
                    region,
                    selection: matches
                        .subcommand_matches("Add")
                        .unwrap()
                        .values_of("Input")
                        .map(|mut i| i.join(" ")),
                    partial,
                })
            }
            Some("Remove") => {
                let region = Region::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Remove").unwrap().is_present(x))
                    .map(|s| Region::from_str(s).unwrap());

                // let source_str = Source::VARIANTS
                //     .iter()
                //     .find(|x| matches.subcommand_matches("Remove").unwrap().is_present(x))
                //     .unwrap_or(&"None");

                // let source = match *source_str {
                //     "Infile" => Some(Source::Infile(
                //         matches
                //             .subcommand_matches("Remove")
                //             .unwrap()
                //             .value_of("Infile")
                //             .unwrap(),
                //     )),
                //     "List" => Some(Source::List(
                //         matches
                //             .subcommand_matches("Remove")
                //             .unwrap()
                //             .values_of("List")
                //             .unwrap()
                //             .join(" "),
                //     )),
                //     _ => None,
                // };

                let partial = Partial::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Remove").unwrap().is_present(x))
                    .map(|s| Partial::from_str(s).unwrap());

                Ok(Mode::Remove {
                    region,
                    selection: matches
                        .subcommand_matches("Remove")
                        .unwrap()
                        .values_of("Input")
                        .map(|mut i| i.join(" ")),
                    // target,
                    partial,
                    // input,
                })
            }
            Some("Write") => {
                let output_str = Output::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Write").unwrap().is_present(x))
                    .unwrap_or(&"None");

                let output = match *output_str {
                    "Outfile" => Some(Output::Outfile(
                        matches
                            .subcommand_matches("Write")
                            .unwrap()
                            .value_of("Outfile")
                            .unwrap(),
                    )),
                    "Overwrite" => Some(Output::Overwrite),
                    _ => None,
                };

                // let region = Region::VARIANTS
                //     .iter()
                //     .find(|x| matches.subcommand_matches("Write").unwrap().is_present(x))
                //     .map(|s| Region::from_str(s).unwrap());

                // let target = Target::VARIANTS
                //     .iter()
                //     .find(|x| matches.subcommand_matches("Write").unwrap().is_present(x))
                //     .map(|s| Target::from_str(s).unwrap());

                Ok(Mode::Write {
                    output,
                    state: matches.subcommand_matches("Write").unwrap().is_present("State")
                    // region,
                    // target,
                })
            }
            _ => unreachable!(),
        }
    }
}
