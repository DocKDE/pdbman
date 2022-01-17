use colored::Colorize;
use itertools::Itertools;
use pest::error::{ErrorVariant, InputLocation};
use pest::Parser;

#[derive(Parser)]
#[grammar = "selection/selection.pest"]
pub struct SelectionParser;

#[derive(Debug, PartialEq)]
pub enum Conjunction {
    And,
    Or,
}

#[derive(Debug, PartialEq)]
pub enum Selection<'a> {
    ID {
        atomlist: Vec<usize>,
        invert: bool,
    },
    Name {
        atomlist: Vec<&'a str>,
        invert: bool,
    },
    Resid {
        reslist: Vec<isize>,
        invert: bool,
    },
    Resname {
        reslist: Vec<&'a str>,
        invert: bool,
    },
    Sphere {
        origin: usize,
        radius: f64,
        invert: bool,
    },
    ResSphere {
        origin: usize,
        radius: f64,
        invert: bool,
    },
}

pub fn parse_selection(
    input: &str,
) -> Result<(Vec<Selection>, Vec<Conjunction>), pest::error::Error<Rule>> {
    let success = SelectionParser::parse(Rule::selection_chain, input)?
        .next()
        .unwrap()
        .into_inner();

    let mut selections = Vec::new();
    let mut conjunctions = Vec::new();

    // println!("{:#?}", success);
    let mut invert = false;
    for i in success {
        match i.as_rule() {
            Rule::selection => {
                for j in i.into_inner() {
                    match j.as_rule() {
                        Rule::negate => invert = !j.as_str().is_empty(),
                        Rule::namesel => {
                            let mut pairs = j.into_inner();
                            let name_str = pairs.next().unwrap().as_str();
                            let name_list: Vec<&str> = pairs
                                .next()
                                .unwrap()
                                .into_inner()
                                .map(|pair| pair.as_str())
                                .collect();
                            selections.push(match name_str.to_lowercase().as_str() {
                                "name" => Selection::Name {
                                    atomlist: name_list,
                                    invert,
                                },
                                "resn" | "resname" => Selection::Resname {
                                    reslist: name_list,
                                    invert,
                                },
                                _ => unreachable!(),
                            })
                        }
                        Rule::idsel => {
                            let mut pairs = j.into_inner();
                            let num_str = pairs.next().unwrap().as_str();

                            let num_elements = pairs.next().unwrap().into_inner();
                            let mut num_vec = Vec::new();
                            for num_pair in num_elements {
                                if let Ok(r) = num_pair.as_str().parse::<usize>() {
                                    num_vec.push(r)
                                } else {
                                    let mut range =
                                        num_pair.into_inner().next().unwrap().into_inner();
                                    let start =
                                        range.next().unwrap().as_str().parse::<usize>().unwrap();

                                    let end =
                                        range.next().unwrap().as_str().parse::<usize>().unwrap();
                                    if start < end {
                                        num_vec.extend(start..=end)
                                    } else {
                                        num_vec.extend(end..=start)
                                    }
                                }
                            }
                            selections.push(match num_str.to_lowercase().as_str() {
                                "id" => Selection::ID {
                                    atomlist: num_vec,
                                    invert,
                                },
                                "resid" => Selection::Resid {
                                    reslist: num_vec.into_iter().map(|n| n as isize).collect(),
                                    invert,
                                },
                                _ => unreachable!(),
                            })
                        }
                        Rule::sphere => {
                            let mut pairs = j.into_inner();
                            let sphere_str = pairs.next().unwrap().as_str();
                            let mut values = pairs.next().unwrap().into_inner();
                            let origin: usize = values.next().unwrap().as_str().parse().unwrap();
                            let radius: f64 = values.next().unwrap().as_str().parse().unwrap();

                            selections.push(match sphere_str.to_lowercase().as_str() {
                                "s" | "sphere" => Selection::Sphere {
                                    origin,
                                    radius,
                                    invert,
                                },
                                "rs" | "ressphere" => Selection::ResSphere {
                                    origin,
                                    radius,
                                    invert,
                                },
                                _ => unreachable!(),
                            })
                        }
                        _ => unreachable!(),
                    }
                }
            }
            Rule::conjunction => match i.as_str().to_lowercase().as_str() {
                "and" | "&" => conjunctions.push(Conjunction::And),
                "or" | "|" => conjunctions.push(Conjunction::Or),
                _ => unreachable!(),
            },
            Rule::EOI => {}
            _ => unreachable!(),
        }
    }

    Ok((selections, conjunctions))
}

pub fn convert_result<'a>(
    input: Result<(Vec<Selection<'a>>, Vec<Conjunction>), pest::error::Error<Rule>>,
    text: &'a str,
) -> Result<(Vec<Selection<'a>>, Vec<Conjunction>), anyhow::Error> {
    let mut err_string = String::new();

    match input {
        Ok(o) => Ok(o),
        Err(e) => {
            match e.location {
                InputLocation::Pos(p) => {
                    err_string.push_str("Error while parsing selection input:\n\n");
                    err_string.push_str(text);
                    err_string.push('\n');
                    err_string.push_str(&(format!("{}{}", " ".repeat(p), "^---".green())));
                }
                _ => unreachable!(),
            };

            match e.variant {
                ErrorVariant::CustomError { message } => err_string.push_str(&message),
                ErrorVariant::ParsingError {
                    positives,
                    negatives: _,
                } => {
                    err_string.push_str("\nExpected: ");
                    err_string.push_str(
                        &positives
                            .into_iter()
                            .map(|rule| match rule {
                                Rule::EOI => "the end of input",
                                Rule::selection => {
                                    "selection keyword: 'id/name/resid/resname/sphere/ressphere'"
                                }
                                Rule::conjunction => "conjunction: 'and'/'or'",
                                Rule::range => "range of numbers, e.g. '1-4'",
                                Rule::sphere_values => "sphere values: origin atom ID and radius",
                                Rule::origin => "origin atom ID",
                                Rule::radius => "sphere radius",
                                Rule::namelist_element => "name of atom or residue",
                                Rule::range_start => {
                                    "number or a range of numbers, e.g. '23' or '4-8'"
                                }
                                Rule::range_end => "number for end of range",
                                _ => "something else",
                            })
                            .join(&format!("{}", " OR\n".red())),
                    )
                }
            }
            bail!(err_string)
        }
    }
}
