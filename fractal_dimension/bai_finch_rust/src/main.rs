#![allow(dead_code)]

use num_complex::Complex64;

type F = f64;
type Matrix2x2 = nalgebra::SMatrix<Complex64, 2, 2>;
type MatrixBig = nalgebra::DMatrix<Complex64>;
type MatrixBigr = nalgebra::DMatrix<f64>;

const G: Matrix2x2 = Matrix2x2::new(
    Complex64::new(2.0 / 12.0, -2.0 / 12.0),
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

const NC: usize = 5;
const K0: i32 = 100;
const LC: usize = 3;
const NC2: usize = NC * NC;
const UPPER_BOUND: i32 = 10_000;

fn fancy_l(q: f64) -> MatrixBigr {
    let mut fancy_l = MatrixBigr::zeros(NC2, NC2);
    for m in 0..NC {
        for n in 0..NC {
            for r in 0..NC {
                for s in 0..NC {
                    if m <= n {
                        let mut sum_m = Complex64::new(0.0, 0.0);
                        for k in 1..K0 {
                            sum_m += normal_m(k, q, n as i32, s as i32)
                                * normal_m(k as i32, q, m as i32, r as i32).conj();
                        }
                        let sum_f = (0..=LC)
                            .map(|l: usize| -> Complex64 {
                                zeta(l as f64 + 2.0 * q, K0)
                                    * (0..=l)
                                        .map(|lp| {
                                            normal_f(q, m as i32, r as i32, (l - lp) as i32)
                                                * normal_f(q, n as i32, s as i32, lp as i32).conj()
                                        })
                                        .sum::<Complex64>()
                            })
                            .sum::<Complex64>();
                        fancy_l[(m * NC + n, r * NC + s)] = 2.0 * (sum_m.re + sum_f.re);
                    } else {
                        fancy_l[(m * NC + n, r * NC + s)] = fancy_l[(n * NC + m, s * NC + r)];
                    }
                }
            }
        }
    }
    fancy_l
    // MatrixBigr::from_diagonal_element(NC2, NC2, 2.0) * (fancy_m(q) + fancy_f(q)).map(|z| z.re)
}

fn integral_choose(n: i32, k: i32) -> f64 {
    match k.cmp(&0) {
        std::cmp::Ordering::Less => 0.0,
        std::cmp::Ordering::Equal => 1.0,
        std::cmp::Ordering::Greater => (0..k).map(|i| (n - i) as f64 / (k - i) as f64).product(),
    }
}

fn non_integral_choose(q: f64, k: i32) -> f64 {
    match k.cmp(&0) {
        std::cmp::Ordering::Less => 0.0,
        std::cmp::Ordering::Equal => 1.0,
        std::cmp::Ordering::Greater => (0..k).map(|i| (q - i as f64) / (k - i) as f64).product(),
    }
}

fn normal_f(q: f64, n: i32, s: i32, l: i32) -> Complex64 {
    let r11 = R[(0, 0)];
    let r12 = R[(0, 1)];
    let r21 = R[(1, 0)];
    let r22 = R[(1, 1)];

    let g11 = G[(0, 0)];
    let g12 = G[(0, 1)];
    let g21 = G[(1, 0)];
    let g22 = G[(1, 1)];
    (0..=s)
        .into_iter()
        .map(|j: i32| -> Complex64 {
            non_integral_choose(-n as f64 - q, j)
                * integral_choose(n, s - j)
                * (0..=j)
                    .map(|l1: i32| -> Complex64 {
                        (0..=(s - j))
                            .map(|l3: i32| -> Complex64 {
                                (0..=(n - s + j))
                                    .map(|l4: i32| -> Complex64 {
                                        integral_choose(j, l1)
                                            * non_integral_choose(
                                                -n as f64 - q - j as f64,
                                                l - l1 - l3 - l4,
                                            )
                                            * integral_choose(s - j, l3)
                                            * integral_choose(n - s + j, l4)
                                            * r21.powi(l1)
                                            * r22.powi(l - l1 - l3 - l4)
                                            * r11.powi(l3)
                                            * r12.powi(l4)
                                            * g21.powi(j - l1)
                                            * g22.powf(
                                                -n as f64 - q - j as f64 - l as f64
                                                    + l1 as f64
                                                    + l3 as f64
                                                    + l4 as f64,
                                            )
                                            * g11.powi(s - j - l3)
                                            * g12.powi(n - s + j - l4)
                                    })
                                    .sum()
                            })
                            .sum()
                    })
                    .sum::<Complex64>()
        })
        .sum()
}

fn zeta(s: f64, k0: i32) -> f64 {
    (k0..=UPPER_BOUND)
        .into_iter()
        .map(|j| (j as f64).powf(-s))
        .sum::<f64>()
        - (UPPER_BOUND as f64).powf(1.0 - s) / (1.0 - s)
}

fn fancy_f(q: f64) -> MatrixBig {
    let mut fancy_f = MatrixBig::zeros(NC2, NC2);
    for m in 0..NC {
        for n in 0..NC {
            for r in 0..NC {
                for s in 0..NC {
                    if m <= n {
                        fancy_f[(m * NC + n, r * NC + s)] = (0..=LC)
                            .map(|l: usize| -> Complex64 {
                                zeta(l as f64 + 2.0 * q, K0)
                                    * (0..=l)
                                        .map(|lp| {
                                            normal_f(q, m as i32, r as i32, (l - lp) as i32)
                                                * normal_f(q, n as i32, s as i32, lp as i32).conj()
                                        })
                                        .sum::<Complex64>()
                            })
                            .sum::<Complex64>();
                    } else {
                        fancy_f[(m * NC + n, r * NC + s)] =
                            fancy_f[(n * NC + m, s * NC + r)].conj();
                    }
                }
            }
        }
    }
    fancy_f
}

fn normal_m(k: i32, q: f64, n: i32, s: i32) -> Complex64 {
    let ak = R + Matrix2x2::from_diagonal_element(Complex64::new(k as f64, 0.0)) * G;
    let a11 = ak[(0, 0)];
    let a12 = ak[(0, 1)];
    let a21 = ak[(1, 0)];
    let a22 = ak[(1, 1)];

    (0..=s)
        .into_iter()
        .map(|j| {
            non_integral_choose(-n as f64 - q, j)
                * integral_choose(n, s - j)
                * a21.powi(j)
                * a22.powf(-n as f64 - q - j as f64)
                * a11.powi(s - j)
                * a12.powi(n - s + j)
        })
        .sum()
}

fn fancy_m(q: f64) -> MatrixBig {
    let mut fancy_m = MatrixBig::zeros(NC2, NC2);
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

fn power_method(
    vec: nalgebra::DVector<f64>,
    mat: nalgebra::DMatrix<f64>,
    iterations: usize,
    tolerance: f64,
) -> f64 {
    let mut previous_entry = vec[0];
    let mut previous_val = 0f64;
    let mut current = mat.clone() * vec;
    let mut current_val = current[0] / previous_entry;
    let mut count = 0;
    while count < iterations && (current_val - previous_val).abs() > tolerance {
        previous_val = current_val;
        previous_entry = current[0];
        current = mat.clone() * current;
        current_val = current[0] / previous_entry;
        count += 1;
    }

    println!("found eigenvalue {} in {} iterations", current_val, count);
    current_val
}

fn test(x: f64) -> f64 {
    x.cos() - x
}

fn lambda(q: f64) -> f64 {
    let lq = fancy_l(q);
    let mut phi0 = nalgebra::DVector::zeros(NC2);
    phi0[0] = 1.0;
    power_method(phi0, lq, 50, std::f64::EPSILON)
}

fn main() {
    // secant_method(lambda, 1.0, 1.3, 1.31, f64::EPSILON, 100);
    // secant_method(test, 0.0, 0.0, 1.0, f64::EPSILON, 100);
    println!("{}", normal_m(5, 1.3, 3, 3));
}
