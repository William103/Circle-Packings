use gad::prelude::*;
use nalgebra::{SMatrix, SVector};
use num::Complex;

type F = f64;
type Matrix2x2 = SMatrix<Complex<F>, 2, 2>;

fn sum(
    g: &mut GraphN,
    ar: &Value<F>,
    ai: &Value<F>,
    br: &Value<F>,
    bi: &Value<F>,
) -> Result<(Value<F>, Value<F>)> {
    Ok((g.add(ar, br)?, g.add(ai, bi)?))
}

fn product(
    g: &mut GraphN,
    ar: &Value<F>,
    ai: &Value<F>,
    br: &Value<F>,
    bi: &Value<F>,
) -> Result<(Value<F>, Value<F>)> {
    let t1 = g.mul(ar, br)?;
    let t2 = g.mul(ai, bi)?;
    let t3 = g.sub(&t1, &t2)?;

    let s1 = g.mul(ar, bi)?;
    let s2 = g.mul(ai, br)?;
    let s3 = g.add(&s1, &s2)?;

    Ok((t3, s3))
}

fn division(
    g: &mut GraphN,
    ar: &Value<F>,
    ai: &Value<F>,
    br: &Value<F>,
    bi: &Value<F>,
) -> Result<(Value<F>, Value<F>)> {
    let numerator1 = {
        let t1 = g.mul(ar, br)?;
        let t2 = g.mul(ai, bi)?;
        g.add(&t1, &t2)?
    };
    let denominator1 = {
        let t1 = g.mul(br, br)?;
        let t2 = g.mul(bi, bi)?;
        g.add(&t1, &t2)?
    };

    let numerator2 = {
        let t1 = g.mul(ar, bi)?;
        let t2 = g.mul(ai, br)?;
        g.sub(&t1, &t2)?
    };
    let denominator2 = {
        let t1 = g.mul(br, br)?;
        let t2 = g.mul(bi, bi)?;
        g.add(&t1, &t2)?
    };

    let f1 = g.div(&numerator1, &denominator1)?;
    let f2 = g.div(&numerator2, &denominator2)?;

    Ok((f1, f2))
}

fn mobius_derivative(mat: Matrix2x2) -> Result<(Value<F>, Value<F>)> {
    let mut g = GraphN::new();
    let x = g.variable(0.0);
    let y = g.variable(0.0);
    let a11r = g.constant(mat[(0, 0)].re);
    let a12r = g.constant(mat[(0, 1)].re);
    let a21r = g.constant(mat[(1, 0)].re);
    let a22r = g.constant(mat[(1, 1)].re);

    let a11i = g.constant(mat[(0, 0)].im);
    let a12i = g.constant(mat[(0, 1)].im);
    let a21i = g.constant(mat[(1, 0)].im);
    let a22i = g.constant(mat[(1, 1)].im);

    let (numeratorr, numeratori) = {
        let (prodr, prodi) = product(&mut g, &x, &y, &a11r, &a11i)?;
        sum(&mut g, &a12r, &a12i, &prodr, &prodi)?
    };

    let (denominatorr, denominatori) = {
        let (prodr, prodi) = product(&mut g, &x, &y, &a21r, &a21i)?;
        sum(&mut g, &a22r, &a22i, &prodr, &prodi)?
    };

    let (resultr, resulti) = division(
        &mut g,
        &numeratorr,
        &numeratori,
        &denominatorr,
        &denominatori,
    )?;

    let x = x.gid()?;

    let one = g.constant(1.0);
    let one2 = g.constant(1.0);
    let gradients1 = g.compute_gradients(resultr.gid()?, one)?;
    let gradients2 = g.compute_gradients(resulti.gid()?, one2)?;

    let du_dx = gradients1.get(x).unwrap();
    let dv_dx = gradients2.get(x).unwrap();

    Ok(Complex::new(*du_dx.data(), -*dv_dx.data()))
}

fn next_order_derivative(
    g: &mut GraphN,
    expr: &GradientId<F>,
    var: &GradientId<F>,
) -> Result<Value<F>> {
    let dz = g.constant(1.0);
    let dz_d = g.compute_gradients(*expr, dz)?;
    Ok(dz_d.get(*var).unwrap().clone())
}

fn power_method<const N: usize>(
    vec: SVector<Complex<F>, N>,
    mat: SMatrix<Complex<F>, N, N>,
    iterations: usize,
) -> Complex<F> {
    let mut current = vec;
    let mut previous = vec;
    for _ in 0..iterations {
        previous = current;
        current = mat * current;
        println!("current guess: {}", current[0] / previous[0]);
    }

    current[0] / previous[0]
}

fn secant_method(f: fn(F) -> F, x0: F, x1: F, accuracy: F, iterations: usize) -> F {
    let mut x0 = x0;
    let mut x1 = x1;
    let mut y0 = f(x0);
    let mut y1 = f(x1);
    let mut count = 0;

    while y1.abs() >= accuracy && count < iterations {
        let new_x = x0 - y0 * (x1 - x0) / (y1 - y0);
        x0 = x1;
        x1 = new_x;
        y0 = y1;
        y1 = f(x1);
        println!("f({}) =\t{}", x1, y1);
        count += 1;
    }
    x0
}

fn main() {
    // let g = Matrix2x2::new(
    //     Complex::new(2.0 / 12.0, -2.0 / 12.0),
    //     Complex::new(-1.0 / 12.0, -5.0 / 12.0),
    //     Complex::new(4.0 / 12.0, -4.0 / 12.0),
    //     Complex::new(-2.0 / 12.0, -10.0 / 12.0),
    // );

    // let r = Matrix2x2::new(
    //     Complex::new(-6.0 / 12.0, 8.0 / 12.0),
    //     Complex::new(0.0, 11.0 / 12.0),
    //     Complex::new(0.0, 4.0 / 12.0),
    //     Complex::new(-6.0 / 12.0, -8.0 / 12.0),
    // );
}
