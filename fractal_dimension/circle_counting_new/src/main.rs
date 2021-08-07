#![allow(dead_code)]

use crate::{
    fractal::fractal_dimension,
    gram_matrix::{algebraic_generators, geometric_generators, root_tuple},
    parser::read_file,
};
use nalgebra::DVector;
use std::process;
use structopt::StructOpt;

use ansi_term::Color::Yellow;

pub mod constants;
pub mod fractal;
pub mod gram_matrix;
pub mod parser;
pub mod search;

/// Compute fractal dimension of crystallographic packings via the circle counting method
#[derive(StructOpt)]
#[structopt(name = "Circle Counter")]
struct Opt {
    /// File containing data in the order generators, root, faces
    #[structopt(name = "data file")]
    data_file: String,

    /// Activate debug mode
    #[structopt(short, long)]
    debug: bool,

    /// Number of sample points in linear regression
    #[structopt(short, long, default_value = "50")]
    n: usize,

    /// Maximum curvature of circles
    #[structopt(short, long, default_value = "1000000")]
    max: f64,

    /// Whether or not to time excecution (automatically enabled by --debug)
    #[structopt(short, long)]
    time: bool,

    /// Cap on the recursion depth (0 means no cap)
    #[structopt(short, long, default_value = "0")]
    depth: usize,

    /// Whether the generators given are geometric or algebraic generators (defaults to algebraic)
    #[structopt(short, long)]
    geometric: bool,
}

fn main() {
    let opt = Opt::from_args();
    let debug = opt.debug;
    let time = debug || opt.time;
    let depth = opt.depth;

    let beginning = std::time::Instant::now();

    let (gram_matrix, faces) = read_file(&opt.data_file).unwrap_or_else(|err| {
        eprintln!("{}", err);
        process::exit(-1);
    });

    let after_parsing = std::time::Instant::now();
    if time {
        let duration = after_parsing.duration_since(beginning);
        println!(
            "Took {}s to parse {}",
            duration.as_secs_f64(),
            opt.data_file
        );
    }

    if debug {
        println!(
            "{} (parsed from command args):\t{}",
            Yellow.paint("max"),
            opt.max
        );
        println!(
            "{} (parsed from command args):\t{}",
            Yellow.paint("n"),
            opt.n
        );
    }

    if debug {
        println!(
            "{} (parsed from file {}):",
            Yellow.paint("Gram Matrix"),
            opt.data_file
        );
        println!("{}\n", gram_matrix);
        println!(
            "{} (parsed from file {}):",
            Yellow.paint("Faces"),
            opt.data_file
        );
        println!("{:?}\n", faces);
        /* println!(
            "{} (parsed from file {}):",
            Yellow.paint("Orthogonal Generators"),
            opt.data_file
        );
        println!("{:?}\n", orthogonal_generators); */
    }

    // finding a suitable root tuple
    let mut root = None;
    for face in faces.iter().skip(1) {
        for vertex in face {
            let mut valid = true;
            for vertex2 in &faces[0] {
                if *vertex == *vertex2 {
                    valid = false;
                    break;
                }
            }
            if valid {
                let temp = root_tuple(
                    &gram_matrix,
                    (faces[0][0], faces[0][1], faces[0][2]),
                    *vertex,
                );

                let mut valid = true;

                let generators = geometric_generators(&gram_matrix, &faces, &temp);
                for gen in generators {
                    if gen[(1, 0)].abs() <= 1e-8 || !gen[(1, 0)].is_finite() {
                        valid = false;
                        break;
                    }
                }

                if valid {
                    root = Some(temp);
                    println!("{}\t{:?}", vertex, faces[0]);
                    break;
                }
            }
        }
        if root.is_some() {
            break;
        }
    }

    if root.is_none() {
        panic!("Invalid face scheme!");
    }

    let root = root.unwrap();

    let generators = geometric_generators(&gram_matrix, &faces, &root);

    if debug {
        println!(
            "{} (computed using gram_matrix):",
            Yellow.paint("Root Matrix")
        );
        println!("{}", root);

        println!(
            "{} (computed using gram_matrix and root):",
            Yellow.paint("Geometric generators")
        );
        for generator in &generators {
            println!("{}", generator);
        }

        let generators = algebraic_generators(&gram_matrix, &faces);
        println!(
            "{} (computed using gram_matrix and root):",
            Yellow.paint("Algebraic generators")
        );
        for generator in &generators {
            println!("{}", generator);
        }
    }

    let mut temp = vec![];
    for circle in root.column_iter() {
        temp.push(DVector::from_iterator(
            4,
            circle.iter().map(|val| val.clone()),
        ));
    }
    let root = temp;

    let delta = fractal_dimension(
        generators,
        root,
        opt.max,
        opt.n,
        debug,
        depth,
        faces,
        vec![],
    )
    .unwrap();
    let after_computing = std::time::Instant::now();
    if time {
        let duration1 = after_computing.duration_since(after_parsing);
        let duration2 = after_computing.duration_since(beginning);

        let time1 = duration1.as_secs_f64();
        let time2 = duration2.as_secs_f64();

        let seconds1 = time1 % 60.0;
        let seconds2 = time2 % 60.0;

        let minutes1 = (time1 as isize) / 60;
        let minutes2 = (time2 as isize) / 60;

        if minutes1 > 0 {
            println!(
                "\nTook {}m {}s to compute fractal dimension; {}m {}s total",
                minutes1, seconds1, minutes2, seconds2
            );
        } else {
            println!(
                "\nTook {}s to compute fractal dimension; {}s total",
                seconds1, seconds2,
            );
        }
    }
    println!("\n{}", delta);
}
