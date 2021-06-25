#![allow(dead_code)]

mod diff;

use diff::*;
use num_complex::Complex64;
use std::f64::consts::PI;
// use nalgebra::{SMatrix, SVector};

type F = f64;
// type NComplex = num::Complex<F>;
// type Matrix2x2 = SMatrix<num::Complex<F>, 2, 2>;
// type Complex = Complex64;

// fn power_method<const N: usize>(
//     vec: SVector<NComplex, N>,
//     mat: SMatrix<NComplex, N, N>,
//     iterations: usize,
// ) -> NComplex {
//     let mut current = vec;
//     let mut previous = vec;
//     for _ in 0..iterations {
//         previous = current;
//         current = mat * current;
//         println!("current guess: {}", current[0] / previous[0]);
//     }

//     current[0] / previous[0]
// }

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

fn dzdt(t: f64) -> Complex64 {
    Complex64::new(
	-2.0 * PI * (-2.0 * PI * t).sin(),
	 2.0 * PI * ( 2.0 * PI * t).cos(),
    )
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

    println!(
        "{}",
        diff(&|z: Complex64| z * z, Complex64::new(0.3, 0.3), 1).unwrap()
    );

    // println!(
    //     "{}",
    //     bacon_sci::integrate::integrate(
    //         0.0,
    //         1.0,
    //         &|t: f64| dzdt(t)
    //             / Complex64::new(
    //                 (std::f64::consts::PI * 2.0 * t).cos(),
    //                 (std::f64::consts::PI * 2.0 * t).sin()
    //             ),
    //         0.001
    //     )
    //     .unwrap()
    // );

    println!(
	"{}",
	integrate(&|z: Complex64| 1.0 / z, &unit_circle, &trapezoid)
    );
}
