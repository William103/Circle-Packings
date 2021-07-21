use ansi_term::Color::Yellow;
use linregress::{FormulaRegressionBuilder, RegressionDataBuilder};
use nalgebra::{DMatrix, DVector};

pub fn fractal_dimension(
    generators: Vec<DMatrix<f64>>,
    root: DVector<f64>,
    faces: Vec<Vec<usize>>,
    upper_bound: f64,
    n: usize,
    debug: bool,
    generations: usize,
    orthogonal_generators: Vec<Vec<usize>>,
) -> Result<f64, linregress::Error> {
    let mut totals = vec![root.len(); n];
    let mut current = vec![(root, std::usize::MAX, false)];
    let mut next = vec![];
    let mut nodes: u64 = 1;

    let xs: Vec<f64> = (1..=n)
        .map(|x| (x as f64 * upper_bound / (2.0 * n as f64) + upper_bound / 2.0))
        .collect();

    let mut i = 0;

    loop {
        next.clear();
        for (tuple, previous_generator, tuple_too_big) in &current {
            let mut add_children = false;
            let mut children = vec![];
            for (i, generator) in generators.iter().enumerate() {
                let mut skip = false;
                for orthogonal_pairs in &orthogonal_generators {
                    if i == orthogonal_pairs[1] && *previous_generator == orthogonal_pairs[0] {
                        skip = true;
                        break;
                    }
                }
                if skip {
                    continue;
                }
                if i != *previous_generator {
                    let new_tuple = generator * tuple;
                    if new_tuple.iter().sum::<f64>() < tuple.iter().sum() {
                        continue;
                    }
                    let mut too_big = true;
                    for (j, curvature) in new_tuple.iter().enumerate() {
                        let mut skip = false;
                        for vertex in &faces[i] {
                            if j == *vertex {
                                skip = true;
                                break;
                            }
                        }
                        if skip {
                            continue;
                        }
                        for (k, max) in xs.iter().enumerate() {
                            if curvature <= max {
                                add_children = true;
                                too_big = false;
                                totals[k] += 1;
                            }
                        }
                    }
                    children.push((new_tuple, i, too_big));
                }
            }
            if add_children || !tuple_too_big {
                nodes += 1;
                for child in children {
                    next.push(child);
                }
            }
        }
        std::mem::swap(&mut current, &mut next);

        i += 1;
        if generations != 0 && i > generations {
            break;
        }
        if debug {
            println!("Generation {}:", i);
            println!("\tnumber of leaves:\t{}", current.len());
            if !current.is_empty() {
                let tuple = &current[current.len() / 2];
                println!("\trandom tuple:\t\t{}", tuple.0);
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
        println!(
            "\nnumber of circles fewer than each of those sample points:\n{:?}",
            totals
            );
        println!("\nTotal number of nodes:\t{}", nodes);
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
