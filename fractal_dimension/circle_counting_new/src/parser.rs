#![allow(dead_code)]
use std::fs::read_to_string;

use nalgebra::DMatrix;
use nom::{
    bytes::complete::is_a,
    character::complete::{char, digit0},
    combinator::map,
    error::{ErrorKind, ParseError},
    multi::separated_list0,
    number::complete::double,
    sequence::delimited,
    IResult,
};

/// Combinator to skip "whitespace"
fn whitespace<'a, E: ParseError<&'a str>>(i: &'a str) -> IResult<&'a str, &'a str, E> {
    is_a(" \t\r\n`\\,")(i)
}

/// Combinator to parse a positive integer and return the integer it represents as a value
fn int<'a, E: ParseError<&'a str>>(i: &'a str) -> IResult<&'a str, usize, E> {
    map(digit0, |s: &str| {
        if let Ok(val) = s.parse::<usize>() {
            val
        } else {
            unreachable!()
        }
    })(i)
}

/// Function that, given a combinator representing the kind of element, returns a combinator that
/// parses a list delimited by braces
fn list<'a, F: 'a, O, E: ParseError<&'a str>>(
    element: F,
) -> impl FnMut(&'a str) -> IResult<&'a str, Vec<O>, E>
where
    F: FnMut(&'a str) -> IResult<&'a str, O, E>,
{
    delimited(char('{'), separated_list0(whitespace, element), char('}'))
}

/// Combinator that parses a matrix: i.e. a nested list of floats, and returns a `DMatrix`
fn matrix<'a, E: 'a + ParseError<&'a str>>(i: &'a str) -> IResult<&'a str, DMatrix<f64>, E> {
    let (tail, output) = list(list(double))(i)?;
    let n = output.len();
    let m = output[0].len();
    for row in &output {
        if row.len() != m {
            return Err(nom::Err::Failure(E::from_error_kind(i, ErrorKind::Fail)));
        }
    }
    Ok((tail, DMatrix::from_row_slice(n, m, &output.concat())))
}

/// Reads the file `filename`, which should first have the gram matrix, and then the face scheme,
/// parses it, and returns the gram matrix and face scheme
pub fn read_file(filename: &str) -> (DMatrix<f64>, Vec<Vec<usize>>) {
    let contents = read_to_string(filename)
        .unwrap_or_else(|e| panic!("Something went wrong reading {}: {}", filename, e));
    let (tail, gram_matrix) =
        matrix::<(&str, ErrorKind)>(contents.as_str()).expect("ERROR PARSING");
    let (tail, _) = whitespace::<(&str, ErrorKind)>(tail).expect("ERROR PARSING");
    let (_, faces) = list(list(int::<(&str, ErrorKind)>))(tail).expect("ERROR PARSING");
    (gram_matrix, faces)
}

#[cfg(test)]
mod tests {
    use nalgebra::DMatrix;
    use nom::error::{Error, ErrorKind};

    use crate::parser::{double, matrix, whitespace};

    use super::{int, list};

    #[test]
    fn test_whitespace() {
        let input = "   \n`\\  asdf";
        assert_eq!(
            whitespace::<(&str, ErrorKind)>(input).unwrap(),
            ("asdf", "   \n`\\  ")
        );
    }

    #[test]
    fn test_number() {
        let input = "-123.456 asdf";
        assert_eq!(
            double::<&str, Error<&str>>(input).unwrap(),
            (" asdf", -123.456)
        );

        let input = "-3 asdf";
        assert_eq!(double::<&str, Error<&str>>(input).unwrap(), (" asdf", -3.0));
    }

    #[test]
    fn test_list() {
        let mut int_list = list(int::<Error<&str>>);
        let mut int_list_list = list(list(int::<Error<&str>>));
        assert_eq!(
            int_list("{3, 2, 1, 42} asdf").unwrap(),
            (" asdf", vec![3, 2, 1, 42])
        );
        assert_eq!(
            int_list_list("{{3, 2, 1, 42}, {1, 2}} asdf").unwrap(),
            (" asdf", vec![vec![3, 2, 1, 42], vec![1, 2]])
        );
    }

    #[test]
    fn test_matrix() {
        let input = "{{1, 0}, {0, 1}} asdf";
        assert_eq!(
            matrix::<Error<&str>>(input).unwrap(),
            (
                " asdf",
                DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0])
            )
        );
        let input = "{{1, -1, -3, -1, -1, -3, -5, -3}, {-1,
    1, -1, -3, -3, -1, -3, -5}, {-3, -1,
    1, -1, -5, -3, -1, -3}, {-1, -3, -1,
    1, -3, -5, -3, -1}, {-1, -3, -5, -3,
    1, -1, -3, -1}, {-3, -1, -3, -5, -1,
    1, -1, -3}, {-5, -3, -1, -3, -3, -1,
    1, -1}, {-3, -5, -3, -1, -1, -3, -1, 1}}";
        assert_eq!(
            matrix::<Error<&str>>(input).unwrap(),
            (
                "",
                DMatrix::from_row_slice(
                    8,
                    8,
                    &[
                        1.0, -1.0, -3.0, -1.0, -1.0, -3.0, -5.0, -3.0, -1.0, 1.0, -1.0, -3.0, -3.0,
                        -1.0, -3.0, -5.0, -3.0, -1.0, 1.0, -1.0, -5.0, -3.0, -1.0, -3.0, -1.0,
                        -3.0, -1.0, 1.0, -3.0, -5.0, -3.0, -1.0, -1.0, -3.0, -5.0, -3.0, 1.0, -1.0,
                        -3.0, -1.0, -3.0, -1.0, -3.0, -5.0, -1.0, 1.0, -1.0, -3.0, -5.0, -3.0,
                        -1.0, -3.0, -3.0, -1.0, 1.0, -1.0, -3.0, -5.0, -3.0, -1.0, -1.0, -3.0,
                        -1.0, 1.0
                    ]
                )
            )
        );
    }
}
