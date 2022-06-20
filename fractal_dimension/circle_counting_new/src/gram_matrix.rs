#![allow(dead_code)]
use nalgebra::DMatrix;

/// Computes the extended gram matrix from G using Nooria's formula.
pub fn extended_gram_matrix(g: &DMatrix<f64>, faces: &[Vec<usize>]) -> DMatrix<f64> {
    let n = g.nrows();
    let m = faces.len();

    let mut top_right = DMatrix::zeros(n, m);

    for i in 0..n {
        for j in 0..m {
            let n = faces[j][0];
            let x = i;
            let one = faces[j][1];
            let two = faces[j][2];

            top_right[(i, j)] = (1.0
                - g[(n, x)].powi(2)
                - (g[(one, x)] + g[(n, x)]) / (1.0 - g[(two, n)])
                    * (g[(two, n)] * g[(n, x)]
                        - g[(two, n)] * g[(one, x)]
                        - g[(one, x)]
                        - g[(n, x)]
                        - 2.0 * g[(two, x)]))
                .abs()
                .sqrt();
        }
    }
    let bottom_right = top_right.transpose()
        * g.clone()
            .pseudo_inverse(1e-8)
            .expect("Invalid gram matrix!")
        * &top_right;
    let mut gext = DMatrix::zeros(n + m, n + m);

    for i in 0..n + m {
        for j in 0..n + m {
            if i < n && j < n {
                gext[(i, j)] = g[(i, j)];
            } else if i < n && j >= n {
                gext[(i, j)] = top_right[(i, j - n)];
            } else if i >= n && j < n {
                gext[(i, j)] = top_right[(j, i - n)];
            } else {
                gext[(i, j)] = bottom_right[(i - n, j - n)];
            }
        }
    }

    gext
}

/// Finds the matrix corresponding to inversion through the circle given by (bt, b, h1, h2).
fn invert_matrix(bt: f64, b: f64, h1: f64, h2: f64) -> DMatrix<f64> {
    DMatrix::from_row_slice(
        4,
        4,
        &[
            h1 * h1 + h2 * h2,
            bt * bt,
            -2.0 * h1 * bt,
            -2.0 * h2 * bt,
            b * b,
            h1 * h1 + h2 * h2,
            -2.0 * b * h1,
            -2.0 * b * h2,
            b * h1,
            h1 * bt,
            1.0 - 2.0 * h1 * h1,
            -2.0 * h1 * h2,
            b * h2,
            h2 * bt,
            -2.0 * h1 * h2,
            1.0 - 2.0 * h2 * h2,
        ],
    )
}

/// Given a gram matrix and face scheme, finds a bounded root tuple using `root_tuple`. Makes sure
/// to avoid extra symmetry in the packing, as that can confuse it later on. I don't actually think
/// that is a problem any more, but I'm too scared to delete it. Picking nice tuples has other
/// benefits too anyway.
pub fn bounded_root_tuple(gram_matrix: &DMatrix<f64>, faces: &[Vec<usize>]) -> DMatrix<f64> {
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
                    gram_matrix,
                    (faces[0][0], faces[0][1], faces[0][2]),
                    *vertex,
                );

                let mut valid = true;

                let generators = geometric_generators(gram_matrix, faces, &temp);
                for gen in generators {
                    if gen[(1, 0)].abs() <= 1e-8 || !gen[(1, 0)].is_finite() {
                        valid = false;
                        break;
                    }
                }

                if valid {
                    root = Some(temp);
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
    root.unwrap()
}

/// Given G, three vertices on a face, and a circle to invert through, it finds the corresponding
/// root tuple. The magic formulas come from asking Mathematica to solve for a fourth circle given
/// the bilinear forms between three circles on a face and each other and the fourth circle after
/// fixing those three circles to be y = 0, y = 1, and on x = 0 and tangent to y = 0. Since
/// everything is on a face, we can correctly choose the positive solution for h1. NOTE: the
/// columns of `root_tuple` are the circles, NOT the rows.
pub fn root_tuple(g: &DMatrix<f64>, face: (usize, usize, usize), invert: usize) -> DMatrix<f64> {
    let (c1, c2, c3) = face;
    let n = g.nrows();
    let mut c = DMatrix::zeros(4, n);

    c[(0, c1)] = 2.0;
    c[(1, c1)] = 0.0;
    c[(2, c1)] = 0.0;
    c[(3, c1)] = 1.0;

    c[(0, c2)] = 0.0;
    c[(1, c2)] = 0.0;
    c[(2, c2)] = 0.0;
    c[(3, c2)] = -1.0;

    c[(0, c3)] = (1.0 - g[(c2, c3)].powi(2)) / (g[(c1, c3)] + g[(c2, c3)]);
    c[(1, c3)] = -g[(c1, c3)] - g[(c2, c3)];
    c[(2, c3)] = 0.0;
    c[(3, c3)] = -g[(c2, c3)];

    let a13 = g[(c1, c3)];
    let a23 = g[(c2, c3)];

    for i in 0..n {
        if i != c1 && i != c2 && i != c3 {
            let a1 = g[(i, c1)];
            let a2 = g[(i, c2)];
            let a3 = g[(i, c3)];

            // b-tilde
            c[(0, i)] = (a1 * (-1.0 + a23 * a23) - a2 * (1.0 + 2.0 * a13 * a23 + a23 * a23)
                + 2.0 * (a13 + a23) * a3)
                / (a13 + a23).powi(2);

            // b
            c[(1, i)] = -a1 - a2;

            // h1 (positive solution)
            c[(2, i)] = (a2 * a2 - a13 * a13 * (-1.0 + a2 * a2) + a23 * a23
                - a1 * a1 * (-1.0 + a23 * a23)
                - 2.0 * a2 * a23 * a3
                + 2.0 * a13 * (a23 - a2 * a3)
                + 2.0 * a1 * (a2 + a13 * a2 * a23 - (a13 + a23) * a3))
                .sqrt()
                / (a13 + a23).abs();

            // h2
            c[(3, i)] = -a2;
        }
    }

    invert_matrix(
        c[(0, invert)],
        c[(1, invert)],
        c[(2, invert)],
        c[(3, invert)],
    ) * &c
}

/// Finds the algebraic generators given g and the faces.
pub fn algebraic_generators(g: &DMatrix<f64>, faces: &[Vec<usize>]) -> Vec<DMatrix<f64>> {
    let gext = extended_gram_matrix(g, faces);
    let n = g.nrows();
    let gt = g
        .clone()
        .pseudo_inverse(1e-8)
        .expect("Something wrong with G");

    let gens: Vec<DMatrix<f64>> = faces
        .iter()
        .enumerate()
        .map(|(i, _)| {
            let column = gext.column(i + n);
            let alpha_i = column.rows(0, n);
            DMatrix::identity(n, n) - 2.0 * alpha_i * alpha_i.transpose() * &gt
        })
        .collect();
    gens
}

/// Finds the geometric generators given g, the faces, and a root tuple
pub fn geometric_generators(
    g: &DMatrix<f64>,
    faces: &[Vec<usize>],
    root: &DMatrix<f64>,
) -> Vec<DMatrix<f64>> {
    let roott = root
        .clone()
        .pseudo_inverse(1e-8)
        .expect("Invalid root tuple!");

    algebraic_generators(g, faces)
        .iter()
        .map(|sigma| root * sigma.transpose() * &roott)
        .collect()
}

#[cfg(test)]
#[rustfmt::skip]
mod tests {
    use super::{root_tuple, extended_gram_matrix, algebraic_generators};
    use nalgebra::DMatrix;

    fn matrices_approx_equal(tolerance: f64, m1: &DMatrix<f64>, m2: &DMatrix<f64>) -> bool {
        if m1.nrows() != m2.nrows() || m1.ncols() != m2.ncols() {
            return false;
        }
        for (i, row) in m1.row_iter().enumerate() {
            for (j, val) in row.column_iter().enumerate() {
                if (val[0] - m2[(i, j)]).abs() >= tolerance {
                    return false;
                }
            }
        }
        return true;
    }

    #[test]
    fn approx_matrices() {
        assert!(matrices_approx_equal(
            1e-3,
            &DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]),
            &DMatrix::from_row_slice(2, 2, &[1.0, -0.0001, 0.0001, 1.0])
        ));
    }

    #[test]
    fn extended_gram_matrix_test() {
        let g = DMatrix::from_row_slice(4, 4,
            &[
                1.0, -1.0, -1.0, -1.0,
                -1.0, 1.0, -1.0, -1.0,
                -1.0, -1.0, 1.0, -1.0,
                -1.0, -1.0, -1.0, 1.0,
            ],
        );

        assert!(matrices_approx_equal(0.001,
            &extended_gram_matrix(
                &g,
                &vec![vec![0, 1, 2], vec![0, 1, 3], vec![0, 2, 3], vec![1, 2, 3]]
            ),
            &DMatrix::from_row_slice(8, 8,
            &[
                1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 2.0,
                -1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 2.0, 0.0,
                -1.0, -1.0, 1.0, -1.0, 0.0, 2.0, 0.0, 0.0,
                -1.0, -1.0, -1.0, 1.0, 2.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 2.0, 1.0, -1.0, -1.0, -1.0,
                0.0, 0.0, 2.0, 0.0, -1.0, 1.0, -1.0, -1.0,
                0.0, 2.0, 0.0, 0.0, -1.0, -1.0, 1.0, -1.0,
                2.0, 0.0, 0.0, 0.0, -1.0, -1.0, -1.0, 1.0,
            ]))
        );
    }

    #[test]
    fn root_tuple_test() {
        let octahedron = DMatrix::from_row_slice(6, 6, 
            &[
                1.0, -1.0, -1.0, -1.0, -3.0, -1.0,
                -1.0, 1.0, -1.0, -1.0, -1.0, -3.0,
                -1.0, -1.0, 1.0, -3.0, -1.0, -1.0,
                -1.0, -1.0, -3.0, 1.0, -1.0, -1.0,
                -3.0, -1.0, -1.0, -1.0, 1.0, -1.0,
                -1.0, -3.0, -1.0, -1.0, -1.0, 1.0,
            ],
        );
        let root = root_tuple(&octahedron, (0, 1, 3), 2);
        let p = DMatrix::from_row_slice(4, 4,
            &[
                0.0, -0.5, 0.0, 0.0,
                -0.5, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0,
                0.0, 0.0, 0.0, 1.0,
            ]);
        assert_eq!(root.row(1).iter().filter(|&x| *x < 0.0).collect::<Vec<&f64>>().len(), 1);
        assert!(matrices_approx_equal(0.001, &octahedron, &(root.transpose() * p * root)));
    }

    #[test]
    fn algebraic_generators_test() {
        let g = DMatrix::from_row_slice(4, 4,
            &[
                1.0, -1.0, -1.0, -1.0,
                -1.0, 1.0, -1.0, -1.0,
                -1.0, -1.0, 1.0, -1.0,
                -1.0, -1.0, -1.0, 1.0,
            ],
        );
        let gens = algebraic_generators(&g, &vec![vec![1, 2, 3], vec![0, 2, 3], vec![0, 1, 3], vec![0, 1, 2]]);
        assert!(matrices_approx_equal(0.001,
                &DMatrix::from_row_slice(4, 4,
                &[
                    -1.0, 2.0, 2.0, 2.0,
                    0.0, 1.0, 0.0, 0.0,
                    0.0, 0.0, 1.0, 0.0,
                    0.0, 0.0, 0.0, 1.0,
                ]),
                &gens[0]));
        assert!(matrices_approx_equal(0.001,
                &DMatrix::from_row_slice(4, 4,
                &[
                    1.0, 0.0, 0.0, 0.0,
                    2.0, -1.0, 2.0, 2.0,
                    0.0, 0.0, 1.0, 0.0,
                    0.0, 0.0, 0.0, 1.0,
                ]),
                &gens[1]));
        assert!(matrices_approx_equal(0.001,
                &DMatrix::from_row_slice(4, 4,
                &[
                    1.0, 0.0, 0.0, 0.0,
                    0.0, 1.0, 0.0, 0.0,
                    2.0, 2.0, -1.0, 2.0,
                    0.0, 0.0, 0.0, 1.0,
                ]),
                &gens[2]));
        assert!(matrices_approx_equal(0.001,
                &DMatrix::from_row_slice(4, 4,
                &[
                    1.0, 0.0, 0.0, 0.0,
                    0.0, 1.0, 0.0, 0.0,
                    0.0, 0.0, 1.0, 0.0,
                    2.0, 2.0, 2.0, -1.0,
                ]),
                &gens[3]));
    }
}
