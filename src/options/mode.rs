use std::error::Error;
use std::str::FromStr;

use strum::VariantNames;
use strum_macros::{Display, EnumString, EnumVariantNames};

#[derive(Display, PartialEq, Clone, Copy, Debug, EnumVariantNames)]
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
        // output: Output,
    },
    Remove {
        region: Region,
        source: Source,
        target: Target,
        partial: Partial,
        // output: Output,
    },
    // Interactive,
    None,
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

// #[derive(Display, PartialEq, Debug, Clone, Copy, EnumString, EnumVariantNames)]
// pub enum Output {
//     Outfile,
//     Overwrite,
//     None,
// }

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

                // let output_str = Output::VARIANTS
                //     .iter()
                //     .find(|x| matches.subcommand_matches("Add").unwrap().is_present(x))
                //     .unwrap_or(&"None");
                // let output = Output::from_str(output_str)?;

                Ok(Mode::Add {
                    region,
                    source,
                    target,
                    partial,
                    // output,
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

                // let output_str = Output::VARIANTS
                //     .iter()
                //     .find(|x| matches.subcommand_matches("Remove").unwrap().is_present(x))
                //     .unwrap_or(&"None");
                // let output = Output::from_str(output_str)?;

                Ok(Mode::Remove {
                    region,
                    source,
                    target,
                    partial,
                    // output,
                })
            }
            // Some("Interactive") => Ok(Mode::Interactive),
            Some(&_) => unreachable!(),
            None => Ok(Mode::None),
        }
    }
}
