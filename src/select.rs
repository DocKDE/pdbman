use std::vec;

use nom::branch::alt;
use nom::character::complete::{alphanumeric1, char, space0, space1, u32};
use nom::combinator::{opt, recognize};
use nom::error::context;
use nom::multi::{separated_list0, separated_list1};
use nom::number::complete::double;
use nom::sequence::{pair, separated_pair, terminated};
use nom::{IResult, Parser};
use nom_supreme::error::ErrorTree;
use nom_supreme::final_parser::{final_parser, Location};
use nom_supreme::parser_ext::ParserExt;
use nom_supreme::tag::complete::{tag, tag_no_case};

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

type Res<T, U> = IResult<T, U, ErrorTree<T>>;

// fn identifier(input: &str) -> IResult<&str, &str> {
//     alt((
//         tag_no_case("ID"),
//         tag_no_case("Name"),
//         tag_no_case("Resid"),
//         tag_no_case("Resname"),
//         tag_no_case("Sphere"),
//         tag_no_case("ResSphere"),
//     ))(input)
// }

fn range_sep(input: &str) -> Res<&str, &str> {
    recognize(alt((char('-'), char(':'))))(input)
}

fn range(input: &str) -> Res<&str, Vec<u32>> {
    let (rem, (a, b)) = separated_pair(
        u32.context("Start of range"),
        range_sep.context("Separator"),
        u32.context("End of range"),
    )(input)?;
    if a < b {
        Ok((rem, (a..=b).collect()))
    } else {
        Ok((rem, (b..=a).collect()))
    }
}

fn digit_to_vec(input: &str) -> Res<&str, Vec<u32>> {
    let (rem, b) = u32(input)?;
    Ok((rem, vec![b]))
}

fn id_with_opt_range(input: &str) -> Res<&str, Vec<u32>> {
    context(
        "id",
        alt((
            range.context("Range"),
            digit_to_vec.context("Single number"),
        )),
    )(input)
}

fn num_list(input: &str) -> Res<&str, Vec<u32>> {
    let (rem, vec) =
        separated_list1(char(','), id_with_opt_range.context("Number list element"))(input)?;
    let mut flat_vec: Vec<u32> = vec.into_iter().flatten().collect();
    flat_vec.sort_unstable();
    flat_vec.dedup();
    Ok((rem, flat_vec))
}

fn str_list(input: &str) -> Res<&str, Vec<&str>> {
    let (rem, mut str_vec) = separated_list1(char(','), alphanumeric1.context("Name"))(input)?;
    str_vec.sort_unstable();
    str_vec.dedup();
    Ok((rem, str_vec))
}

fn sphere_values(input: &str) -> Res<&str, (usize, f64)> {
    let (rem, (origin, radius)) =
        separated_pair(u32.context("Origin"), space1, double.context("Radius"))(input)?;
    Ok((rem, (origin as usize, radius)))
}

fn sphere(input: &str) -> Res<&str, Selection> {
    let (rem, (neg, (sphere_str, (origin, radius)))) = pair(
        opt(negate.context("negate")),
        separated_pair(
            alt((
                alt((tag_no_case("Sphere"), tag_no_case("S"))),
                alt((tag_no_case("ResSphere"), tag_no_case("RS"))),
            )),
            space0,
            sphere_values,
        ),
    )(input)?;
    match sphere_str.to_lowercase().as_str() {
        "sphere" | "s" => Ok((
            rem,
            Selection::Sphere {
                origin,
                radius,
                invert: neg.is_some(),
            },
        )),
        "ressphere" | "rs" => Ok((
            rem,
            Selection::ResSphere {
                origin,
                radius,
                invert: neg.is_some(),
            },
        )),
        _ => unreachable!(),
    }
}

fn id_and_list(input: &str) -> Res<&str, Selection> {
    let (rem, (neg, (id_str, atom_vec))) = pair(
        opt(negate.context("Negate")),
        separated_pair(
            alt((tag_no_case("ID"), tag_no_case("Resid"))),
            space0,
            num_list.context("Number list"),
        ),
    )(input)?;
    match id_str.to_lowercase().as_str() {
        "id" => Ok((
            rem,
            Selection::ID {
                atomlist: atom_vec.into_iter().map(|n| n as usize).collect(),
                invert: neg.is_some(),
            },
        )),
        "resid" => Ok((
            rem,
            Selection::Resid {
                reslist: atom_vec.into_iter().map(|n| n as isize).collect(),
                invert: neg.is_some(),
            },
        )),
        _ => unreachable!(),
    }
}

fn str_and_list(input: &str) -> Res<&str, Selection> {
    let (rem, (neg, (name_str, atom_vec))) = pair(
        opt(negate.context("Negate")),
        separated_pair(
            alt((
                tag_no_case("Name"),
                alt((tag_no_case("Resname"), tag_no_case("Resn"))),
            )),
            space0,
            str_list.context("List of names"),
        ),
    )(input)?;
    match name_str.to_lowercase().as_str() {
        "name" => Ok((
            rem,
            Selection::Name {
                atomlist: atom_vec,
                invert: neg.is_some(),
            },
        )),
        "resname" | "resn" => Ok((
            rem,
            Selection::Resname {
                reslist: atom_vec,
                invert: neg.is_some(),
            },
        )),
        _ => unreachable!(),
    }
}

fn negate(input: &str) -> Res<&str, &str> {
    terminated(alt((tag("!"), tag_no_case("Not"))), space0)(input)
}

fn selection(input: &str) -> Res<&str, Selection> {
    alt((sphere, id_and_list, str_and_list))(input)
}

fn conjunction(input: &str) -> Res<&str, Conjunction> {
    let (rem, conj) = alt((
        alt((tag("&"), tag_no_case("and"))),
        alt((tag("|"), tag_no_case("or"))),
    ))(input)?;
    Ok((
        rem,
        match conj.to_lowercase().as_str() {
            "and" | "&" => Conjunction::And,
            "or" | "|" => Conjunction::Or,
            _ => unreachable!(),
        },
    ))
}

fn conj_then_sel(input: &str) -> Res<&str, (Conjunction, Selection)> {
    separated_pair(conjunction.context("Conjunction"), space1, selection)(input)
}

#[allow(clippy::type_complexity)]
fn optional_sels(input: &str) -> Res<&str, Option<(Vec<Conjunction>, Vec<Selection>)>> {
    // let list_parser = separated_list0(
    //     space1,
    //     separated_pair(conjunction.context("Conjunction"), space1, selection),
    // );
    let (rem, vec) = separated_list0(space1, conj_then_sel)(input)?;
    let mut conjunctions = Vec::new();
    let mut selections = Vec::new();
    if vec.is_empty() {
        Ok((rem, (None)))
    } else {
        for (conj, sel) in vec {
            conjunctions.push(conj);
            selections.push(sel);
        }
        Ok((rem, Some((conjunctions, selections))))
    }
}

fn full_list(input: &str) -> Res<&str, (Vec<Selection>, Option<Vec<Conjunction>>)> {
    let (rem, (base_sele, opt_sele)) =
        separated_pair(selection.cut(), space0, optional_sels.cut())(input)?;
    let mut final_sele_vec = vec![base_sele];

    if let Some((conj_vec, mut sele_vec)) = opt_sele {
        final_sele_vec.append(&mut sele_vec);
        Ok((rem, (final_sele_vec, Some(conj_vec))))
    } else {
        Ok((rem, (final_sele_vec, None)))
    }
}

fn final_parse<'a>(
    list_parse: impl Parser<&'a str, (Vec<Selection<'a>>, Option<Vec<Conjunction>>), ErrorTree<&'a str>>,
) -> impl FnMut(&'a str) -> Result<(Vec<Selection>, Option<Vec<Conjunction>>), ErrorTree<Location>>
{
    final_parser(list_parse)
}

#[allow(clippy::type_complexity)]
pub fn parse_expr(
    input: &str,
) -> Result<(Vec<Selection>, Option<Vec<Conjunction>>), ErrorTree<Location>> {
    final_parse(full_list.cut())(input)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn num_list_test() {
        let str = "1,2,4-5,6:9asdjlaks";
        assert_eq!(num_list(str).unwrap().0, "asdjlaks");
        assert_eq!(num_list(str).unwrap().1, vec![1, 2, 4, 5, 6, 7, 8, 9])
    }

    #[test]
    fn str_list_test() {
        let str = "his,wat,ala,cu";
        assert_eq!(str_list(str).unwrap().0, "");
        assert_eq!(str_list(str).unwrap().1, vec!["ala", "cu", "his", "wat"])
    }

    #[test]
    fn selection_test() {
        let id = "id 1,2,3-8,4:12";
        let resn = "resname his,wat,ala,cu";
        let sphere_str = "sphere  2589 5.9";

        let (id_rem, sele) = id_and_list(id).unwrap();
        assert_eq!(id_rem, "");
        assert_eq!(
            sele,
            Selection::ID {
                atomlist: vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                invert: false
            }
        );

        let (resn_rem, sele) = str_and_list(resn).unwrap();
        assert_eq!(resn_rem, "");
        assert_eq!(
            sele,
            Selection::Resname {
                reslist: vec!["ala", "cu", "his", "wat"],
                invert: false
            }
        );

        let (rem, sele) = sphere(sphere_str).unwrap();
        assert_eq!(rem, "");
        assert_eq!(
            sele,
            Selection::Sphere {
                origin: 2589,
                radius: 5.9,
                invert: false
            }
        )
    }

    #[test]
    fn negate_test() {
        let strings = ["!name cu", "! name cu", "not name cu", "notname cu"];

        for str in strings {
            let (rem, sele) = selection(str).unwrap();
            assert_eq!(rem, "");
            assert_eq!(
                sele,
                Selection::Name {
                    atomlist: vec!["cu"],
                    invert: true
                }
            )
        }
    }

    #[test]
    fn full_list_test() {
        let str = "resname his and id 1,2,3 or id 123 and !sphere 23 4.6";
        let (rem, sele) = full_list(str).unwrap();

        let selections_vec = vec![
            Selection::Resname {
                reslist: vec!["his"],
                invert: false,
            },
            Selection::ID {
                atomlist: vec![1, 2, 3],
                invert: false,
            },
            Selection::ID {
                atomlist: vec![123],
                invert: false,
            },
            Selection::Sphere {
                origin: 23,
                radius: 4.6,
                invert: true,
            },
        ];

        let conj_vec = vec![Conjunction::And, Conjunction::Or, Conjunction::And];

        assert_eq!(rem, "");
        assert_eq!(sele.0, selections_vec);
        assert_eq!(sele.1, Some(conj_vec));
    }

    #[test]
    fn alias_test() {
        let strings = ["resn cu", "rs 2589 6", "s 2589 6"];

        for str in strings {
            let (rem, _) = selection(str).unwrap();
            assert_eq!(rem, "");
        }
    }
}
