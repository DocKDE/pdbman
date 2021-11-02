use std::str::FromStr;

use anyhow::Result;
use itertools::Itertools;
use strum::VariantNames;
use strum_macros::{Display, EnumString, EnumVariantNames};

#[derive(Display, PartialEq, Debug, EnumVariantNames)]
pub enum Mode<'a> {
    Query {
        source: Source<'a>,
        target: Target,
    },
    Analyze {
        region: Option<Region>,
        target: Option<Target>,
        distance: Option<Distance>,
    },
    Add {
        region: Option<Region>,
        source: Option<Source<'a>>,
        target: Option<Target>,
        partial: Option<Partial>,
    },
    Remove {
        region: Option<Region>,
        source: Option<Source<'a>>,
        target: Option<Target>,
        partial: Option<Partial>,
    },
    Write {
        output: Option<Output<'a>>,
        region: Option<Region>,
    },
}

#[derive(Display, PartialEq, Debug, Clone, Copy, EnumString, EnumVariantNames)]
pub enum Region {
    QM1,
    QM2,
    Active,
    // None,
}

#[derive(Display, PartialEq, Debug, Clone, EnumString, EnumVariantNames)]
pub enum Source<'a> {
    Infile(&'a str),
    List(&'a str),
    Sphere(usize, f64),
    // None,
}

#[derive(Display, PartialEq, Debug, Clone, Copy, EnumString, EnumVariantNames)]
pub enum Target {
    Atoms,
    Residues,
    // None,
}

#[derive(Display, PartialEq, Debug, Clone, EnumString, EnumVariantNames)]
pub enum Output<'a> {
    Outfile(&'a str),
    Overwrite,
    // None,
}

#[derive(Display, PartialEq, Debug, Clone, Copy, EnumString, EnumVariantNames)]
pub enum Partial {
    Sidechain,
    Backbone,
    // None,
}

#[derive(Display, PartialEq, Debug, Clone, Copy, EnumString, EnumVariantNames)]
pub enum Distance {
    Clashes,
    Contacts,
    // None,
}

impl<'a> Mode<'a> {
    /// Creates new Mode enum from clap::ArgMatches struct. This is
    /// where the given command line options are stored for later use.
    pub fn new(matches: &clap::ArgMatches) -> Result<Mode, anyhow::Error> {
        match matches.subcommand_name() {
            Some("Query") => {
                let source_str = Source::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Query").unwrap().is_present(x))
                    .unwrap_or(&"None");

                let source = match *source_str {
                    "Infile" => Source::Infile(
                        matches
                            .subcommand_matches("Query")
                            .unwrap()
                            .value_of("Infile")
                            .unwrap(), // .to_owned(),
                    ),
                    "List" => Source::List(
                        matches
                            .subcommand_matches("Query")
                            .unwrap()
                            .value_of("List")
                            .unwrap(), // .to_owned(),
                    ),
                    "Sphere" => {
                        let (origin_str, radius_str) = matches
                            .subcommand_matches("Query")
                            .unwrap()
                            .values_of("Sphere")
                            .unwrap()
                            .next_tuple()
                            .unwrap();
                        Source::Sphere(origin_str.parse()?, radius_str.parse()?)
                    }
                    _ => Source::from_str(source_str)?,
                };

                let target_str = Target::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Query").unwrap().is_present(x))
                    .unwrap_or(&"None");
                let target = Target::from_str(target_str)?;

                Ok(Mode::Query { source, target })
            }
            Some("Analyze") => {
                // let region_str = Region::VARIANTS
                let region= Region::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Analyze").unwrap().is_present(x))
                    .and_then(|s| Some(Region::from_str(s).unwrap()));
                    // .unwrap_or(&"None");
                // let region = Region::from_str(region_str)?;
                // let region = match region_str {
                //     Some(s) => Some(Region::from_str(s)?),
                //     None => None,
                // };

                // let target_str = Target::VARIANTS
                let target= Target::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Analyze").unwrap().is_present(x))
                    .and_then(|s| Some(Target::from_str(s).unwrap()));
                    // .unwrap_or(&"None");
                // let target = Target::from_str(target_str)?;
                // let target = match target_str {
                //     Some(s) => Some(Target::from_str(s)?),
                //     None => None,
                // };

                // let distance_str = Distance::VARIANTS
                let distance= Distance::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Analyze").unwrap().is_present(x))
                    .and_then(|s| Some(Distance::from_str(s).unwrap()));
                    // .unwrap_or(&"None");
                // let distance = Distance::from_str(distance_str)?;
                // let distance = match distance_str {
                //     Some(s) => Some(Distance::from_str(s)?),
                //     None => None,
                // };

                Ok(Mode::Analyze {
                    region,
                    target,
                    distance,
                })
            }
            Some("Add") => {
                // let region_str = Region::VARIANTS
                //     .iter()
                //     .find(|x| matches.subcommand_matches("Add").unwrap().is_present(x))
                //     .unwrap_or(&"None");
                // let region = Region::from_str(region_str)?;
                let region= Region::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Add").unwrap().is_present(x))
                    .and_then(|s| Some(Region::from_str(s).unwrap()));

                let source_str = Source::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Add").unwrap().is_present(x))
                    .unwrap();

                let source = match *source_str {
                    "Infile" => Some(Source::Infile(
                        matches
                            .subcommand_matches("Add")
                            .unwrap()
                            .value_of("Infile")
                            .unwrap(), // .to_owned(),
                    )),
                    "List" => Some(Source::List(
                        matches
                            .subcommand_matches("Add")
                            .unwrap()
                            .value_of("List")
                            .unwrap(), // .to_string(),
                    )),
                    "Sphere" => {
                        let (origin_str, radius_str) = matches
                            .subcommand_matches("Add")
                            .unwrap()
                            .values_of("Sphere")
                            .unwrap()
                            .next_tuple()
                            .unwrap();
                        Some(Source::Sphere(origin_str.parse()?, radius_str.parse()?))
                    }
                    _ => None
                };

                // let target_str = Target::VARIANTS
                let target= Target::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Add").unwrap().is_present(x))
                    .and_then(|s| Some(Target::from_str(s).unwrap()));
                    // .unwrap_or(&"None");
                // let target = Target::from_str(target_str)?;

                // let partial_str = Partial::VARIANTS
                let partial= Partial::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Add").unwrap().is_present(x))
                    .and_then(|s| Some(Partial::from_str(s).unwrap()));
                    // .unwrap_or(&"None");
                // let partial = Partial::from_str(partial_str)?;

                Ok(Mode::Add {
                    region,
                    source,
                    target,
                    partial,
                })
            }
            Some("Remove") => {
                // let region_str = Region::VARIANTS
                //     .iter()
                //     .find(|x| matches.subcommand_matches("Remove").unwrap().is_present(x))
                //     .unwrap_or(&"None");
                // let region = Region::from_str(region_str)?;
                let region= Region::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Remove").unwrap().is_present(x))
                    .and_then(|s| Some(Region::from_str(s).unwrap()));

                let source_str = Source::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Remove").unwrap().is_present(x))
                    .unwrap_or(&"None");

                let source = match *source_str {
                    "Infile" => Some(Source::Infile(
                        matches
                            .subcommand_matches("Remove")
                            .unwrap()
                            .value_of("Infile")
                            .unwrap(), // .to_owned(),
                    )),
                    "List" => Some(Source::List(
                        matches
                            .subcommand_matches("Remove")
                            .unwrap()
                            .value_of("List")
                            .unwrap(), // .to_string(),
                    )),
                    "Sphere" => {
                        let (origin_str, radius_str) = matches
                            .subcommand_matches("Remove")
                            .unwrap()
                            .values_of("Sphere")
                            .unwrap()
                            .next_tuple()
                            .unwrap();
                        Some(Source::Sphere(origin_str.parse()?, radius_str.parse()?))
                    }
                    _ => None
                };

                // let target_str = Target::VARIANTS
                let target= Target::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Remove").unwrap().is_present(x))
                    .and_then(|s| Some(Target::from_str(s).unwrap()));
                    // .unwrap_or(&"None");
                // let target = Target::from_str(target_str)?;

                // let partial_str = Partial::VARIANTS
                let partial= Partial::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Remove").unwrap().is_present(x))
                    .and_then(|s| Some(Partial::from_str(s).unwrap()));
                    // .unwrap_or(&"None");
                // let partial = Partial::from_str(partial_str)?;

                Ok(Mode::Remove {
                    region,
                    source,
                    target,
                    partial,
                })
            }
            Some("Write") => {
                // let region_str = Region::VARIANTS
                //     .iter()
                //     .find(|x| matches.subcommand_matches("Write").unwrap().is_present(x))
                //     .unwrap_or(&"None");
                // let region = Region::from_str(region_str)?;
                let region= Region::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Write").unwrap().is_present(x))
                    .and_then(|s| Some(Region::from_str(s).unwrap()));

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
                            .unwrap(), // .to_owned(),
                    )),
                    "Overwrite" => Some(Output::from_str(output_str)?),
                    _ => None
                };

                Ok(Mode::Write { output, region })
            }
            _ => unreachable!(),
        }
    }
}
