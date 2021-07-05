use nalgebra::{SMatrix, SVector};
use linregress::{FormulaRegressionBuilder, RegressionDataBuilder};

type Matrix4 = SMatrix<f64, 4, 4>;
type Vector4 = SVector<f64, 4>;

fn search(generators: &Vec<Matrix4>, root: Vector4, upper_bound: f64) -> usize {
    let mut current = vec![(root, 5)];
    let mut next = vec![];
    let mut total = 0;
    loop {
        next.clear();
        for (tuple, previous_generator) in &current {
            for (i, generator) in generators.iter().enumerate() {
                let new_tuple = generator * tuple;
                if !(i == 0 && *previous_generator == 3) && !(i == 1 && *previous_generator == 2) && i != *previous_generator && new_tuple[i] < upper_bound && new_tuple.iter().sum::<f64>() > tuple.iter().sum() {
                    next.push((new_tuple, i));
                }
            }
        }

        if current.is_empty() {
            break;
        }
        total += next.len();
        std::mem::swap(&mut current, &mut next);
    }
    2 * (total + root.len())
}

fn fractal_dimension(generators: Vec<Matrix4>, root: Vector4, upper_bound: f64, n: isize) -> Result<f64, linregress::Error> {
    let xs: Vec<f64> = (1..=n).map(|x| (x as f64 * upper_bound / (n as f64)).ln()).collect();
    let ys: Vec<f64> = xs.iter().map(|upper_bound| (search(&generators, root, upper_bound.exp()) as f64).ln()).collect();

    let data = vec![("Y", ys), ("X", xs)];
    let data = RegressionDataBuilder::new().build_from(data)?;

    let formula = "Y ~ X";
    let model = FormulaRegressionBuilder::new()
        .data(&data)
        .formula(formula)
        .fit()?;

    let parameters = model.parameters;
    Ok(parameters.regressor_values[0])
}

fn main() {
    let generators = vec![
        Matrix4::new(
            -1.0, 2.0, 2.0, 0.0,
            0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 1.0,
        ),
        Matrix4::new(
            1.0, 0.0, 0.0, 0.0,
            2.0, -1.0, 0.0, 2.0,
            0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 1.0,
        ),
        Matrix4::new(
            1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0,
            2.0, 0.0, -1.0, 2.0,
            0.0, 0.0, 0.0, 1.0,
        ),
        Matrix4::new(
            1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0,
            0.0, 2.0, 2.0, -1.0,
        ),
    ];

    let root = Vector4::new(-2., 4., 5., 9.);

    println!("{}", fractal_dimension(generators, root, 10_000_000.0, 50).unwrap());
}
