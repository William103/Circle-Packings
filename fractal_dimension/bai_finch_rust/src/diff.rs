use num_complex::Complex64;
use simba::scalar::ComplexField;
use std::f64::consts::PI;

pub type CFunction<'a> = &'a dyn Fn(Complex64) -> Complex64;
pub type Path<'a> = &'a dyn Fn(f64) -> Complex64;
pub type Rule<'a> = &'a dyn Fn(CFunction, Complex64, Complex64) -> Complex64;

pub fn unit_circle(t: f64) -> Complex64 {
    num::Complex::new((2.0 * PI * t).cos(), (2.0 * PI * t).sin())
}

pub fn trapezoid(f: CFunction, a: Complex64, b: Complex64) -> Complex64 {
    0.5 * (b - a) * (f(a) + f(b))
}

fn factorial(n: i32) -> i32 {
    (1..=n).product()
}

fn test(t: f64) -> Complex64 {
    1.0 / unit_circle(t)
}

// pub fn diff(fun: CFunction, a: Complex64, n: i32) -> Result<Complex64, String> {
//     let f = |t: f64| -> Complex64 {
//         let z = unit_circle(t);
//         fun(z) / (z - a).powi(n + 1)
//     };
//     let integral = bacon_sci::integrate::integrate_fixed(0f64, 1f64, f, 10)?;
//     Ok(integral * Complex64::new(factorial(n) as f64, 0.0) / Complex64::new(0.0, 2.0 * PI))
// }

pub fn diff(fun: CFunction, a: Complex64, n: i32) -> Result<Complex64, String> {
    let f = |z: Complex64| -> Complex64 {
        fun(z) / (z - a).powi(n + 1)
    };
    let integral = integrate(&f, &unit_circle, &trapezoid);
    Ok(integral * Complex64::new(factorial(n) as f64, 0.0) / Complex64::new(0.0, 2.0 * PI))
}

pub fn integrate(f: CFunction, path: Path, rule: Rule) -> Complex64 {
    let increment = 0.00001;
    let mut integral = Complex64::new(0.0, 0.0);

    let mut i = increment;
    let mut previous_point = path(0.0);
    while i <= 1.0 {
        let point = path(i);
        integral += rule(f, previous_point, point);
        previous_point = point;

        i += increment;
    }

    integral
}
