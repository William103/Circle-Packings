#![allow(dead_code)]

use crate::command_args::DataFormats;
use nalgebra::DMatrix;
use std::fmt::{Debug, Display};

/// Formats a matrix according to the format specified by `format`.
pub fn format_matrix<T: 'static>(matrix: &DMatrix<T>, format: &DataFormats) -> String
where
    T: Display + Copy + Debug + PartialEq,
{
    match format {
        DataFormats::HumanReadable => format!("{}", matrix),
        DataFormats::Mathematica | DataFormats::CStyle => {
            let mut res = "{".to_string();
            for (i, row) in matrix.row_iter().enumerate() {
                res += "{";
                for (j, elem) in row.iter().enumerate() {
                    res += format!("{}", elem).as_str();
                    if j < row.len() - 1 {
                        res += ", ";
                    }
                }
                res += "}";
                if i < matrix.nrows() - 1 {
                    res += ",\n";
                }
            }

            res + "}"
        }
        DataFormats::Python => {
            let mut res = "numpy.array([".to_string();
            for (i, row) in matrix.row_iter().enumerate() {
                res += "[";
                for (j, elem) in row.iter().enumerate() {
                    res += format!("{}", elem).as_str();
                    if j < row.len() - 1 {
                        res += ", ";
                    }
                }
                res += "]";
                if i < matrix.nrows() - 1 {
                    res += ",\n";
                }
            }

            res + "])"
        }
    }
}

/// Formats a vector of matrices according to the format specified by `format`.
pub fn format_vec_of_matrices<T: 'static>(vec: &Vec<DMatrix<T>>, format: &DataFormats) -> String
where
    T: Display + Copy + Debug + PartialEq,
{
    match format {
        DataFormats::HumanReadable | DataFormats::Python => {
            let mut res = "[".to_string();
            for (i, elem) in vec.iter().enumerate() {
                res += format!("{}", format_matrix(elem, format)).as_str();
                if i < vec.len() - 1 {
                    res += ", ";
                }
            }
            res + "]"
        }
        DataFormats::CStyle | DataFormats::Mathematica => {
            let mut res = "{".to_string();
            for (i, elem) in vec.iter().enumerate() {
                res += format!("{}", format_matrix(elem, format)).as_str();
                if i < vec.len() - 1 {
                    res += ", ";
                }
            }
            res + "}"
        }
    }
}
