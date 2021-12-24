#![allow(unused_variables)]

mod command_args;
mod dimension;
mod gram_matrix;
mod parser;
mod picture;
mod search;
mod time;

use command_args::Polyhedral;
use dimension::fractal_dimension;
use gram_matrix::{bounded_root_tuple, geometric_generators};
use parser::read_file;
use picture::render_packing;
use structopt::StructOpt;
use time::Timer;

fn main() {
    let opt = Polyhedral::from_args();
    match opt {
        Polyhedral::Dimension {
            data_file,
            debug,
            n,
            max,
            time,
            recursion_depth,
        } => {
            let timer = Timer::new(time);
            let (gram_matrix, faces) = read_file(&data_file);
            let root = bounded_root_tuple(&gram_matrix, &faces);
            let dimension = fractal_dimension(
                &geometric_generators(&gram_matrix, &faces, &root),
                &root,
                max,
                n,
                debug,
                0,
                &timer,
            )
            .unwrap_or_else(|e| panic!("Regression went wrong: {}", e));
            if debug {
                println!("Fractal dimension: {}", dimension);
            } else {
                println!("{}", dimension);
            }
        }

        Polyhedral::Picture {
            data_file,
            debug,
            max,
            time,
            recursion_depth,
            output_file,
            width,
            height,
        } => {
            let timer = Timer::new(time);
            let (gram_matrix, faces) = read_file(&data_file);
            render_packing(
                width,
                height,
                max,
                &gram_matrix,
                &faces,
                &output_file,
                recursion_depth,
                debug,
                &timer,
            );
        }
    }
}
