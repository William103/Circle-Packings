#![allow(dead_code)]

use crate::{lib::fractal_dimension, parser::read_file};
use std::process;
use structopt::StructOpt;

use ansi_term::Color::Yellow;

mod lib;
mod parser;

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

    /// Cap on the number of generations (0 means no cap)
    #[structopt(short, long, default_value = "0")]
    generations: usize,
}

fn main() {
    let opt = Opt::from_args();
    let debug = opt.debug;
    let time = debug || opt.time;
    let generations = opt.generations;

    let beginning = std::time::Instant::now();

    let (generators, root, faces) = read_file(&opt.data_file).unwrap_or_else(|err| {
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
            "{} (parsed from file {})",
            Yellow.paint("Generators"),
            opt.data_file
        );
        for generator in &generators {
            println!("{}", generator);
        }
        println!(
            "{} (parsed from file {}):",
            Yellow.paint("Root Tuple"),
            opt.data_file
        );
        println!("{}", root);
        println!(
            "{} (parsed from file {}):",
            Yellow.paint("Faces"),
            opt.data_file
        );
        println!("{:?}\n", faces);
    }

    let delta = fractal_dimension(generators, root, faces, opt.max, opt.n, debug, generations).unwrap();
    let after_computing = std::time::Instant::now();
    if time {
        let duration1 = after_computing.duration_since(after_parsing);
        let duration2 = after_computing.duration_since(beginning);
        println!(
            "\nTook {}s to compute fractal dimension; {}s total",
            duration1.as_secs_f64(),
            duration2.as_secs_f64()
        );
    }
    println!("\n{}", delta);
}
