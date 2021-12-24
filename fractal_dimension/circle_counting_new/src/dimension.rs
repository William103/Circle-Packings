#![allow(dead_code)]

use crate::{search::{Searcher, Task}, time::Timer};
use ansi_term::Color::Yellow;
use linregress::{FormulaRegressionBuilder, RegressionDataBuilder};
use nalgebra::{DMatrix, DVector};

/// struct representing the circle counting task to be used by `crate::search::Searcher` to compute
/// the fractal dimension via the circle counting method
struct CircleCountingTask {
    counts: Vec<usize>,
    maxes: Vec<f64>,
}

impl Task for CircleCountingTask {
    /// Implementation of `Task` for `CircleCountingTask` that counts circles, unsurprisingly.
    fn recurse(&mut self, circle: &DVector<f64>) -> bool {
        let mut too_big = true;
        for (i, max) in self.maxes.iter().enumerate() {
            if circle[1] < *max {
                too_big = false;
                self.counts[i] += 1;
            }
        }
        !too_big
    }
}

impl CircleCountingTask {
    /// Create a new `CircleCountingTask`. `n` is the number of data points (skewed to larger
    /// maximum curvatures) while `max_curvature` is the overall maximum curvature.
    fn new(n: usize, max_curvature: f64) -> Self {
        Self {
            counts: vec![0; n],
            maxes: (1..=n)
                .map(|x| (x as f64 * max_curvature / (2.0 * n as f64) + max_curvature / 2.0))
                .collect(),
        }
    }
}

/// Compute the fractal dimension of a polyhedral circle packing with generators given by
/// `generators` and root tuple given by `root` by counting circles with curvature less than
/// `max_curvature`, using `n` data points in the regression. The value of `debug` determines
/// whether or not debug information is printed, and `max_depth` is an optional maximum recursion
/// depth.
pub fn fractal_dimension(
    generators: &[DMatrix<f64>],
    root: &DMatrix<f64>,
    max_curvature: f64,
    n: usize,
    debug: bool,
    max_depth: usize,
    timer: &Timer,
) -> Result<f64, linregress::Error> {
    if debug {
        println!("{}", Yellow.paint("Counting circles"));
    }
    let mut circle_counter = CircleCountingTask::new(n, max_curvature);
    let mut searcher = Searcher::new(&mut circle_counter, generators, max_depth);
    searcher.search(root, debug, timer);

    if debug {
        println!("{}", Yellow.paint("Done counting circles!"));
        println!("sample points:\n{:?}", circle_counter.maxes);
        println!(
            "\nnumber of circles fewer than each of those sample points:\n{:?}",
            circle_counter.counts
        );
    }

    let xs: Vec<f64> = circle_counter.maxes.iter().map(|x| x.ln()).collect();
    let ys: Vec<f64> = circle_counter
        .counts
        .iter()
        .map(|x| (*x as f64).ln())
        .collect();
    let data = vec![("Y", ys), ("X", xs)];
    let data = RegressionDataBuilder::new().build_from(data)?;

    let formula = "Y ~ X";
    let model = FormulaRegressionBuilder::new()
        .data(&data)
        .formula(formula)
        .fit()?;

    if debug {
        println!("\n{}", Yellow.paint("Built regression model!"));
        println!("Model info:");
        println!("\tr^2:\t\t{}", model.rsquared);
        println!("\tp-value:\t{}", model.pvalues.pairs()[0].1);
        println!("\tintercept:\t{}", model.parameters.intercept_value);
        println!("\tslope:\t{}", model.parameters.regressor_values[0]);
    }

    Ok(model.parameters.regressor_values[0])
}
