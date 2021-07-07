use linregress::{FormulaRegressionBuilder, RegressionDataBuilder};
use nalgebra::{DMatrix, DVector};
use ansi_term::Color::Yellow;

pub fn fractal_dimension(
    generators: Vec<DMatrix<f64>>,
    root: DVector<f64>,
    faces: Vec<Vec<usize>>,
    upper_bound: f64,
    n: usize,
    debug: bool,
) -> Result<f64, linregress::Error> {
    let mut totals = vec![root.len(); n];
    let mut current = vec![(root, std::usize::MAX)];
    let mut next = vec![];

    let xs: Vec<f64> = (1..=n)
        .map(|x| (x as f64 * upper_bound / (n as f64)))
        .collect();

    let mut i = 0;
    loop {
        next.clear();
        for (tuple, previous_generator) in &current {
            for (i, generator) in generators.iter().enumerate() {
                if i != *previous_generator {
                    let new_tuple = generator * tuple;
                    let mut add = false;
                    for (j, curvature) in new_tuple.iter().enumerate() {
                        let mut skip = false;
                        for face in &faces[i] {
                            if j == *face {
                                skip = true;
                                break;
                            }
                        }
                        if skip {
                            continue;
                        }
                        for (k, max) in xs.iter().enumerate() {
                            if curvature <= max {
                                if k == n - 1 {
                                    add = true;
                                }
                                totals[k] += 1;
                            }
                        }
                    }
                    if add {
                        next.push((new_tuple, i));
                    }
                }
            }
        }
        std::mem::swap(&mut current, &mut next);

        i += 1;
        if debug {
            println!("Generation {}:", i);
            println!("\tnumber of leaves:\t{}", current.len());
            if !current.is_empty() {
                println!("\trandom tuple:\t\t{}", current[current.len() / 2].0);
            }
            println!();
        }

        if current.is_empty() {
            break;
        }
    }

    if debug {
        println!("{}", Yellow.paint("Done counting circles!"));
        println!("sample points:\n{:?}", xs);
        println!("\nnumber of circles fewer than each of those sample points:\n{:?}", totals);
    }

    let xs: Vec<f64> = xs.iter().map(|x| x.ln()).collect();
    let ys: Vec<f64> = totals.iter().map(|x| (*x as f64).ln()).collect();
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
