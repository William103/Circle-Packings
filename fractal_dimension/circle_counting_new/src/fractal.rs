use ansi_term::Color::Yellow;
use linregress::{FormulaRegressionBuilder, RegressionDataBuilder};
use nalgebra::{DMatrix, DVector};

use crate::search::Searcher;

pub fn fractal_dimension(
    generators: Vec<DMatrix<f64>>,
    root: Vec<DVector<f64>>,
    upper_bound: f64,
    n: usize,
    debug: bool,
    generations: usize,
    orthogonal_generators: Vec<Vec<usize>>,
) -> Result<f64, linregress::Error> {
    let xs: Vec<f64> = (1..=n)
        .map(|x| (x as f64 * upper_bound / (2.0 * n as f64) + upper_bound / 2.0))
        .collect();

    let mut searcher = Searcher::new(&xs, &generators);
    for circle in root {
        searcher.search(&circle, std::usize::MAX, 0, generations);
    }

    if debug {
        println!("{}", Yellow.paint("Done counting circles!"));
        println!("sample points:\n{:?}", xs);
        println!(
            "\nnumber of circles fewer than each of those sample points:\n{:?}",
            searcher.counts
        );
    }

    let xs: Vec<f64> = xs.iter().map(|x| x.ln()).collect();
    let ys: Vec<f64> = searcher.counts.iter().map(|x| (*x as f64).ln()).collect();
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

    let parameters = model.parameters;

    Ok(parameters.regressor_values[0])
}
