#![allow(dead_code)]

mod diff;
mod integration;

use diff::*;
use num_complex::Complex64;
use std::f64::consts::PI;
// use nalgebra::{SMatrix, SVector};

type F = f64;
type Matrix2x2 = nalgebra::SMatrix<Complex64, 2, 2>;
type MatrixBig = nalgebra::SMatrix<Complex64, NC2, NC2>;
type MatrixBigr = nalgebra::SMatrix<f64, NC2, NC2>;

const G: Matrix2x2 = Matrix2x2::new(
    Complex64::new(-2.0 / 12.0, -2.0 / 12.0),
    Complex64::new(-1.0 / 12.0, -5.0 / 12.0),
    Complex64::new(4.0 / 12.0, -4.0 / 12.0),
    Complex64::new(-2.0 / 12.0, -10.0 / 12.0),
);

const R: Matrix2x2 = Matrix2x2::new(
    Complex64::new(-6.0 / 12.0, 8.0 / 12.0),
    Complex64::new(0.0, 11.0 / 12.0),
    Complex64::new(0.0, 4.0 / 12.0),
    Complex64::new(-6.0 / 12.0, -8.0 / 12.0),
);

const NC: usize = 10;
const K0: i32 = 100;
const LC: usize = 15;
const NC2: usize = NC * NC;

fn fancy_l(q: f64) -> MatrixBigr {
    MatrixBigr::from_diagonal_element(2.0) * fancy_m(q).map(|z| z.re)
}

fn integral_choose(n: i32, k: i32) -> f64 {
    match k.cmp(&0) {
	std::cmp::Ordering::Less => 0.0,
	std::cmp::Ordering::Equal => 1.0,
	std::cmp::Ordering::Greater => {
	    (0..k).map(|i|
		     (n - i) as f64 / (k - i) as f64
	    ).product()
	}
    }
}

fn non_integral_choose(q: f64, k: i32) -> f64 {
    match k.cmp(&0) {
	std::cmp::Ordering::Less => 0.0,
	std::cmp::Ordering::Equal => 1.0,
	std::cmp::Ordering::Greater => {
	    (0..k).map(|i|
		     (q - i as f64) / (k - i) as f64
	    ).product()
	}
    }
}

fn normal_m(k: i32, q: f64, n: i32, s: i32) -> Complex64 {
    let ak = R + Matrix2x2::from_diagonal_element(Complex64::new(k as f64, 0.0)) * G;
    let a11 = ak[(0, 0)];
    let a12 = ak[(0, 1)];
    let a21 = ak[(1, 0)];
    let a22 = ak[(1, 1)];

    (0..=s)
	.map(|j| {
	    non_integral_choose(-n as f64 - q, j)
		* integral_choose(n, s - j)
		* a21.powi(k)
		* a22.powf(-n as f64 - q - j as f64)
		* a11.powi(s - j)
		* a12.powi(n - s + k)
	})
	.sum()
}

fn fancy_m(q: f64) -> MatrixBig {
    let mut fancy_m = MatrixBig::zeros();
    for m in 0..NC {
	for n in 0..NC {
	    for r in 0..NC {
		for s in 0..NC {
		    if m <= n {
			let mut sum = Complex64::new(0.0, 0.0);
			for k in 1..K0 {
			    sum += normal_m(k, q, n as i32, s as i32)
				* normal_m(k as i32, q, m as i32, r as i32).conj();
			}
			fancy_m[(m * NC + n, r * NC + s)] = sum;
		    } else {
			fancy_m[(m * NC + n, r * NC + s)] =
			    fancy_m[(n * NC + m, s * NC + r)].conj();
		    }
		}
	    }
	}
    }
    fancy_m
}

fn secant_method(f: fn(F) -> F, target: F, x0: F, x1: F, accuracy: F, iterations: usize) -> F {
    let mut x0 = x0;
    let mut x1 = x1;
    let mut y0 = f(x0);
    let mut y1 = f(x1);
    let mut count = 0;

    println!("f({}) =\t{}", x1, y1);
    while (y1 - target).abs() >= accuracy && count < iterations {
	let new_x = x0 - (y0 - target) * (x1 - x0) / (y1 - y0);
	x0 = x1;
	x1 = new_x;
	y0 = y1;
	y1 = f(x1);
	println!("f({}) =\t{}", x1, y1);
	count += 1;
    }
    x0
}

fn power_method<const N: usize>(
    vec: nalgebra::SVector<f64, N>,
    mat: nalgebra::SMatrix<f64, N, N>,
    iterations: usize,
) -> f64 {
    let mut current = vec;
    let mut previous = vec;
    for _ in 0..iterations {
	previous = current;
	current = mat * current;
	println!("current eigenvalue: {}", current[0] / previous[0]);
    }

    current[0] / previous[0]
}

fn dzdt(t: f64) -> Complex64 {
    Complex64::new(
	-2.0 * PI * (2.0 * PI * t).sin(),
	2.0 * PI * (2.0 * PI * t).cos(),
    )
}

fn lambda(q: f64) -> f64 {
    let lq = fancy_l(q);
    let mut phi0 = nalgebra::SVector::zeros();
    phi0[0] = 1.0;
    power_method(phi0, lq, 10)
}

fn main() {
    println!("{}", non_integral_choose(15.3, 3));
    println!("{}", normal_m(3, 1.3, 2, 1));
    // secant_method(lambda, 1.0, 1.3, 1.31, f64::EPSILON, 10);
}
