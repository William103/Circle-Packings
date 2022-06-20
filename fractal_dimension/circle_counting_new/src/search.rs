#![allow(dead_code)]
use crate::time::Timer;
use nalgebra::{DMatrix, DVector};

/// If (a - b).abs() < TOLERANCE, then a is close enough to b to be considered equal, to get around
/// floating point errors.
const TOLERANCE: f64 = 1e-3;

/// A trait representing some kind of task to be done on each circle in a packing.
pub trait Task {
    /// Will be given `circle` from each circle in the packing by recursively applying generators.
    /// Should return `true` if we should recurse on this particular circle and false otherwise.
    fn recurse(&mut self, circle: &DVector<f64>) -> bool;
}

/// Struct that will actually do the searching, using the task.
pub struct Searcher<'a, T: Task> {
    task: &'a mut T,
    generators: &'a [DMatrix<f64>],
    max_depth: usize,
}

impl<'a, T: Task> Searcher<'a, T> {
    /// Create a new `Searcher`. `task` is the particular task that this `Searcher` will run.
    /// `generators` is the list of generators that will be used to find all the circles.
    /// `max_depth` is an optional maximum recursion depth; a value of `0` indicates no maximum
    /// depth.
    pub fn new(task: &'a mut T, generators: &'a [DMatrix<f64>], max_depth: usize) -> Self {
        Self {
            task,
            generators,
            max_depth,
        }
    }

    /// Actually runs the search on the packing with the particular root tuple given by `root`.
    pub fn search(&mut self, root: &DMatrix<f64>, debug: bool, timer: &Timer) {
        let n = root.ncols();
        for (i, circle) in root.column_iter().enumerate() {
            let circle = circle.clone_owned();
            if self.task.recurse(&circle) {
                self.search_helper(&circle, usize::max_value(), 0);
            }
            if debug && i < n - 1 {
                print!("{}%", ((i + 1) as f64 / n as f64 * 100.0).round());
            }
            if i < n - 1 {
                let (n, i) = (n as f64, i as f64);
                if debug {
                    timer.time_remaining(
                        (n - i - 1.0) / (i + 1.0),
                        if debug { "\t\t" } else { "" },
                        " remaining",
                    );
                }
            }
        }
        if debug {
            println!();
            timer.time_stamp("Took ", " to count circles\n");
        } else {
            timer.time_stamp("", "");
        }
    }

    /// Auxiliary function with additional arguments to aid in recursion and make the front-facing
    /// function less cluttered and easier to run.
    fn search_helper(&mut self, circle: &DVector<f64>, previous_generator: usize, depth: usize) {
        if self.max_depth > 0 && depth > self.max_depth {
            return;
        }
        for (i, generator) in self.generators.iter().enumerate() {
            if i != previous_generator {
                let new_circle = generator * circle;
                if new_circle[1] - circle[1] <= TOLERANCE {
                    continue;
                }
                if self.task.recurse(&new_circle) {
                    self.search_helper(&new_circle, i, depth + 1);
                }
            }
        }
    }
}
