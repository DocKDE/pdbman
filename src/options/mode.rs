use std::str::FromStr;

use anyhow::{Context, Result};
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
        partial: Option<Partial>,
        selection: Option<String>,
    },
    Remove {
        region: Option<Region>,
        partial: Option<Partial>,
        selection: Option<String>,
    },
    Write {
        output: Option<Output<'a>>,
        state: bool,
    },
    Measure {
        measure_target: MeasureTarget,
    },
}

#[derive(Display, PartialEq, Debug, Clone, Copy, PartialOrd, EnumString, EnumVariantNames)]
pub enum Region {
    QM1,
    QM2,
    Active,
}

#[derive(Display, PartialEq, Debug, Clone, Copy, PartialOrd, EnumString, EnumVariantNames)]
pub enum Target {
    Atoms,
    Residues,
}

#[derive(Display, PartialEq, Debug, Clone, PartialOrd, EnumString, EnumVariantNames)]
pub enum MeasureTarget {
    Sphere(usize, f64),
    Atoms(Vec<usize>),
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
                    partial,
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

                Ok(Mode::Write {
                    output,
                    state: matches
                        .subcommand_matches("Write")
                        .unwrap()
                        .is_present("State"),
                })
            }
            Some("Measure") => {
                let measure_str = MeasureTarget::VARIANTS
                    .iter()
                    .find(|x| matches.subcommand_matches("Measure").unwrap().is_present(x))
                    .unwrap();

                let measure = match *measure_str {
                    "Sphere" => {
                        let mut sphere_values = matches
                            .subcommand_matches("Measure")
                            .unwrap()
                            .values_of("Sphere")
                            .unwrap();
                        let origin_str = sphere_values.next().unwrap();
                        let origin_id = origin_str
                            .parse::<usize>()
                            .context(format!("Invalid input for Atom ID: {}", origin_str))?;
                        let radius_str = sphere_values.next().unwrap();
                        let radius_float = radius_str
                            .parse::<f64>()
                            .context(format!("Invalid input for radius: {}", radius_str))?;
                        Ok(Mode::Measure {
                            measure_target: MeasureTarget::Sphere(origin_id, radius_float),
                        })
                    }
                    "Atoms" => {
                        let atom_str = matches
                            .subcommand_matches("Measure")
                            .unwrap()
                            .values_of("Atoms")
                            .unwrap();

                        let mut atom_ids = Vec::new();
                        for id in atom_str {
                            if let Ok(i) = id.parse::<usize>() {
                                atom_ids.push(i)
                            } else {
                                bail!("Invalid input for Atom ID: {}", id)
                            }
                        }

                        // Error out if duplicates are present because the input will make no sense
                        let id_count = atom_ids.iter().copied().sorted().dedup().count();
                        ensure!(
                            atom_ids.len() == id_count,
                            "Atom list for measurements contains duplicate elements"
                        );

                        Ok(Mode::Measure {
                            measure_target: MeasureTarget::Atoms(atom_ids),
                        })
                    }
                    _ => unreachable!(),
                };

                measure
            }
            _ => unreachable!(),
        }
    }
}
