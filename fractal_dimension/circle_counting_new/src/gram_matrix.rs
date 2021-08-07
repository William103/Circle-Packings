use nalgebra::DMatrix;

/// Computes the extended gram matrix from G using Nooria's formula.
pub fn extended_gram_matrix(g: &DMatrix<f64>, faces: &Vec<Vec<usize>>) -> DMatrix<f64> {
    let n = g.row_iter().len();
    let m = faces.len();

    let mut b = DMatrix::zeros(n, m);

    for i in 0..n {
        for j in 0..m {
            let n = faces[j][0];
            let x = i;
            let one = faces[j][1];
            let two = faces[j][2];

            b[(i, j)] = (1.0
                - g[(n, x)].powi(2)
                - (g[(one, x)] + g[(n, x)]) / (1.0 - g[(two, n)])
                    * (g[(two, n)] * g[(n, x)]
                        - g[(two, n)] * g[(one, x)]
                        - g[(one, x)]
                        - g[(n, x)]
                        - 2.0 * g[(two, x)]))
                .sqrt();
        }
    }
    let d = b.transpose()
        * g.clone()
            .pseudo_inverse(1e-8)
            .expect("Invalid gram matrix!")
        * &b;
    let mut gext = DMatrix::zeros(n + m, n + m);

    for i in 0..n + m {
        for j in 0..n + m {
            if i < n && j < n {
                gext[(i, j)] = g[(i, j)];
            } else if i < n && j >= n {
                gext[(i, j)] = b[(i, j - n)];
            } else if i >= n && j < n {
                gext[(i, j)] = b[(j, i - n)];
            } else {
                gext[(i, j)] = d[(i - n, j - n)];
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

/// Given G, three vertices on a face, and a circle to invert through, it finds the corresponding
/// root tuple. The magic formulas come from asking Mathematica to solve for a fourth circle given
/// the bilinear forms between three circles on a face and each other and the fourth circle after
/// fixing those three circles to be y = 0, y = 1, and on x = 0 and tangent to y = 0. Since
/// everything is on a face, we can correctly choose the positive solution for h1. NOTE: the
/// columns of `root_tuple` are the circles, NOT the rows.
pub fn root_tuple(g: &DMatrix<f64>, face: (usize, usize, usize), invert: usize) -> DMatrix<f64> {
    let (c1, c2, c3) = face;
    let n = g.row_iter().len();
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

    let newc = invert_matrix(
        c[(0, invert)],
        c[(1, invert)],
        c[(2, invert)],
        c[(3, invert)],
    ) * &c;

    newc
}

/// Finds the algebraic generators given g and the faces.
pub fn algebraic_generators(g: &DMatrix<f64>, faces: &Vec<Vec<usize>>) -> Vec<DMatrix<f64>> {
    let gext = extended_gram_matrix(g, &faces);
    let n = g.row_iter().len();
    let gt = g
        .clone()
        .pseudo_inverse(1e-8)
        .expect("Something wrong with G");

    faces
        .iter()
        .enumerate()
        .map(|(i, _)| {
            let column = gext.column(i + n);
            let alpha_i = column.rows(0, n);
            DMatrix::identity(n, n) - 2.0 * alpha_i * alpha_i.transpose() * &gt
        })
        .collect()
}

/// Finds the geometric generators given g, the faces, and a root tuple
pub fn geometric_generators(
    g: &DMatrix<f64>,
    faces: &Vec<Vec<usize>>,
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

        assert_ne!(
            extended_gram_matrix(
                &g,
                &vec![vec![0, 1, 2], vec![0, 1, 3], vec![0, 2, 3], vec![1, 2, 3]]
            ),
            DMatrix::from_row_slice(8, 8,
            &[
                1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 2.0,
                -1.0, 1.0, -1.0, -1.0, 0.0, 2.0, 0.0, 0.0,
                -1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 2.0, 0.0,
                -1.0, -1.0, -1.0, 1.0, 2.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 2.0, 1.0, -1.0, -1.0, -1.0,
                0.0, 2.0, 0.0, 0.0, -1.0, 1.0, -1.0, -1.0,
                0.0, 0.0, 2.0, 0.0, -1.0, -1.0, 1.0, -1.0,
                2.0, 0.0, 0.0, 0.0, -1.0, -1.0, -1.0, 1.0,
            ])
        );

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

        assert_ne!(
            extended_gram_matrix(
                &octahedron,
                &vec![
                    vec![0, 1, 3],
                    vec![0, 2, 1],
                    vec![0, 3, 5],
                    vec![0, 5, 2],
                    vec![1, 2, 4],
                    vec![1, 4, 3],
                    vec![2, 5, 4],
                    vec![3, 4, 5],
                ]
            ),
            g
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
        assert_ne!(
            root_tuple(&octahedron, (0, 1, 3), 2),
            octahedron,
        );
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
        assert_ne!(
            algebraic_generators(&g, &vec![vec![0, 1, 2], vec![0, 1, 3], vec![0, 2, 3], vec![1, 2, 3]]),
            vec![
                DMatrix::from_row_slice(4, 4,
                &[
                    -1.0, 2.0, 2.0, 2.0,
                    0.0, 1.0, 0.0, 0.0,
                    0.0, 0.0, 1.0, 0.0,
                    0.0, 0.0, 0.0, 1.0,
                ]),
                DMatrix::from_row_slice(4, 4,
                &[
                    1.0, 0.0, 0.0, 0.0,
                    2.0, -1.0, 2.0, 2.0,
                    0.0, 0.0, 1.0, 0.0,
                    0.0, 0.0, 0.0, 1.0,
                ]),
                DMatrix::from_row_slice(4, 4,
                &[
                    1.0, 0.0, 0.0, 0.0,
                    0.0, 1.0, 0.0, 0.0,
                    2.0, 2.0, -1.0, 2.0,
                    0.0, 0.0, 0.0, 1.0,
                ]),
                DMatrix::from_row_slice(4, 4,
                &[
                    1.0, 0.0, 0.0, 0.0,
                    0.0, 1.0, 0.0, 0.0,
                    0.0, 0.0, 1.0, 0.0,
                    2.0, 2.0, 2.0, -1.0,
                ]),
            ]
        );

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

        assert_ne!(
            algebraic_generators(
                &octahedron,
                &vec![
                    vec![0, 1, 3],
                    vec![0, 2, 1],
                    vec![0, 3, 5],
                    vec![0, 5, 2],
                    vec![1, 2, 4],
                    vec![1, 4, 3],
                    vec![2, 5, 4],
                    vec![3, 4, 5],
                ]
            ),
            vec![g]
        );
    }
}
