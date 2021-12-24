#![allow(dead_code)]

use ansi_term::Colour::Yellow;
use std::{fs::File, io::Write, path::Path};

use crate::{
    gram_matrix::{geometric_generators, bounded_root_tuple},
    search::{Searcher, Task}, time::Timer,
};
use nalgebra::{DMatrix, DVector};

/// Struct keeping track of information needed to create a picture of a circle packing
struct PictureTask {
    svg: String,
    max_curvature: f64,
    scale: f64,
    width: f64,
    height: f64,
}

impl Task for PictureTask {
    /// Implementation of `Task` for `PictureTask` that draws each circle it's given in an svg file.
    fn recurse(&mut self, circle: &DVector<f64>) -> bool {
        if circle[1] >= self.max_curvature {
            return false;
        }
        // rescale and center the circle in pixel-space on the canvas
        let x = circle[2] / circle[1] * self.scale;
        let y = circle[3] / circle[1] * self.scale;
        let cx = x + self.width / 2.0;
        let cy = y + self.height / 2.0;
        let r = (self.scale / circle[1]).abs();

        self.svg.push_str(
            format!(
                "<circle cx=\"{}\" cy=\"{}\" r=\"{}\" fill=\"none\" stroke=\"black\"/>\n",
                cx, cy, r
            )
            .as_str(),
        );

        true
    }
}

impl PictureTask {
    /// Constructor for `PictureTask` that sets up all the data it's going to need.
    fn new(max_curvature: f64, width: usize, height: usize) -> PictureTask {
        PictureTask {
            svg: format!(
                "<?xml version=\"1.0\" standalone=\"no\"?>\n<svg width=\"{}\" height=\"{}\">\n",
                width, height
            ),
            max_curvature,
            scale: std::cmp::min(height, width) as f64 / 2.0,
            width: width as f64,
            height: height as f64,
        }
    }
}

/// Normalizes the root tuple such that the outer circle has radius 1 and is centered at the origin
fn normalize_root_tuple(root: &mut DMatrix<f64>) {
    // find the outer circle; must exist if root tuple is bounded
    let mut outer_circle = None;
    for circle in root.column_iter() {
        if circle[1] < 0.0 {
            outer_circle = Some(circle);
            break;
        }
    }
    let outer_circle = outer_circle.expect("something went very very wrong");

    // modify each of the circles in the root tuple to be normalized
    let outer_x = outer_circle[2] / outer_circle[1];
    let outer_y = outer_circle[3] / outer_circle[1];
    let outer_r = 1.0 / outer_circle[1];
    for mut circle in root.column_iter_mut() {
        let x = circle[2] / circle[1];
        let y = circle[3] / circle[1];
        let r = 1.0 / circle[1];
        let new_x = (x - outer_x) / outer_r.abs();
        let new_y = (y - outer_y) / outer_r.abs();
        let new_r = r / outer_r.abs();
        circle[1] = 1.0 / new_r;
        circle[2] = circle[1] * new_x;
        circle[3] = circle[1] * new_y;
        circle[0] = (circle[2] * circle[2] + circle[3] * circle[3] - 1.0) / circle[1];
    }
}

/// Renders an svg of the circle packing; `width` and `height` are the width and height of the svg
/// respectively. All the circles with curvature less than `max_curvature` will be rendered.
/// `gram_matrix` is, unsurprisingly, the gram matrix while `faces` is the list of faces
/// corresponding to the gram matrix. The svg will be written to `output_file`. `max_depth` is an
/// optional cap on the recursion depth (i.e. word length) that shouldn't ever really be necessary,
/// but can be useful if you don't want very much detail near the points of parabolicity for
/// whatever reason.
pub fn render_packing(
    width: usize,
    height: usize,
    max_curvature: f64,
    gram_matrix: &DMatrix<f64>,
    faces: &[Vec<usize>],
    output_file: &str,
    max_depth: usize,
    debug: bool,
    timer: &Timer,
) {
    // find a bounded root tuple and normalize it
    let mut root = bounded_root_tuple(gram_matrix, faces);
    normalize_root_tuple(&mut root);
    if debug {
        println!("{}: {}", Yellow.paint("Normalized root tuple: "), root);
    }

    // find generators based on this new root tuple
    let generators = geometric_generators(gram_matrix, faces, &root);

    // run picture_task on said normalized root tuple
    if debug {
        println!("{}", Yellow.paint("Drawing circles"));
    }
    let mut picture_task = PictureTask::new(max_curvature, width, height);
    let mut searcher = Searcher::new(&mut picture_task, &generators, max_depth);
    searcher.search(&root, debug, timer);
    picture_task.svg.push_str("</svg>");

    // write result to disk
    let mut file = File::create(Path::new(&output_file)).expect("file creation went wrong");
    file.write_all(picture_task.svg.as_bytes())
        .expect("write went wrong");
}
