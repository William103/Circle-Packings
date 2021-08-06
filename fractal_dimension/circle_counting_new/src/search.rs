use nalgebra::{DMatrix, DVector};

pub struct Searcher<'a> {
    pub counts: Vec<u64>,
    maxes: &'a Vec<f64>,
    generators: &'a Vec<DMatrix<f64>>,
    max_depth: usize,
    tolerance: f64,
    faces: &'a Vec<Vec<usize>>,
}

impl<'a> Searcher<'a> {
    pub fn new(
        maxes: &'a Vec<f64>,
        generators: &'a Vec<DMatrix<f64>>,
        max_depth: usize,
        n: u64,
        faces: &'a Vec<Vec<usize>>,
    ) -> Self {
        Self {
            counts: vec![n; maxes.len()],
            maxes,
            generators,
            max_depth,
            tolerance: 1e-8,
            faces,
        }
    }

    pub fn search(
        &mut self,
        circle: &DVector<f64>,
        seed: usize,
        previous_generator: usize,
        depth: usize,
    ) {
        if self.max_depth > 0 && depth > self.max_depth {
            return;
        }
        for (i, generator) in self.generators.iter().enumerate() {
            if depth == 0 {
                for vertex in &self.faces[i] {
                    if *vertex == seed {
                        continue;
                    }
                }
            }
            if i != previous_generator {
                let new_circle = generator * circle;
                if new_circle[1] - circle[1] <= self.tolerance {
                    continue;
                } else {
                    let mut seen = false;
                    for (j, max) in self.maxes.iter().enumerate() {
                        if seen || new_circle[1] <= *max {
                            seen = true;
                            self.counts[j] += 1;
                        }
                    }
                    if !seen {
                        continue;
                    }
                } // ensure seen is dropped before recursive call to minimize stack frame size
                self.search(&new_circle, seed, i, depth + 1);
            }
        }
    }
}
