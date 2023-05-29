pub trait Interpolation {
    fn interpolate(&self, x: f64) -> f64;
}

pub struct Lagrange<'a> {
    // func: fn(f64) -> f64,
    nodes: &'a [f64],
    values: Vec<f64>,
    omega_prime: Vec<f64>,
}

impl Lagrange<'_> {
    pub fn new(func: fn(f64) -> f64, nodes: &'_ [f64]) -> Lagrange<'_> {
        let values = nodes.iter().copied().map(func).collect();
        dbg!(&values);

        let omega_prime = nodes.iter().enumerate().map(|(i, x_i)| {
            nodes[..i].iter().chain(&nodes[i + 1..])
                .fold(1.0, |acc, x_j| { acc * (x_i - x_j) })
        }).collect::<Vec<f64>>();
        dbg!(&omega_prime);

        Lagrange {
            // func,
            nodes,
            values,
            omega_prime,
        }
    }
}

impl Interpolation for Lagrange<'_> {
    fn interpolate(&self, x: f64) -> f64 {
        dbg!(x);

        let omega = self.nodes.iter()
            .map(|x_i| x - x_i)
            .fold(1.0, |acc, x_i| acc * x_i);
        dbg!(omega);

        let result = self.nodes.iter().enumerate()
            .map(|(k, x_k)| {
                let y_k = self.values[k];
                if omega.is_normal() && self.omega_prime[k].is_normal() {
                    y_k * omega / ((x - x_k) * self.omega_prime[k])
                } else {
                    // Can't use cached omega_prime
                    let x_k = self.nodes[k];
                    let l_k = self.nodes[..k].iter().chain(&self.nodes[k + 1..])
                        .fold(1.0, |acc, x_i|
                            acc * (x - x_i) / (x_k - x_i));
                    y_k * l_k
                }
            }).sum::<f64>();
        dbg!(result);

        result
    }
}

pub fn chebyshev_zeros(n: usize) -> Vec<f64> {
    let denominator = (2 * n) as f64;
    (0..n).map(|k| {
        let numerator = (2 * k + 1) as f64;
        -(numerator * std::f64::consts::PI / denominator).cos()
    }).collect()
}

pub struct Linear<'a> {
    // func: fn(f64) -> f64,
    nodes: &'a [f64],
    values: Vec<f64>,
}

impl Linear<'_> {
    pub fn new(func: fn(f64) -> f64, nodes: &'_ [f64]) -> Linear<'_> {
        assert!(nodes.len() > 1);

        // Validate that nodes are sorted
        // let all_sorted = nodes.array_windows().all(|[a, b]| a < b);
        assert!(nodes.is_sorted(), "Nodes must be sorted!");

        let values = nodes.iter().copied().map(func).collect();
        dbg!(&values);

        Linear {
            // func,
            nodes,
            values,
        }
    }
}

impl Interpolation for Linear<'_> {
    fn interpolate(&self, x: f64) -> f64 {
        let (a, b) = find_between(self.nodes, x);

        let x_a = self.nodes[a];
        let x_b = self.nodes[b];
        let y_a = self.values[a];
        let y_b = self.values[b];

        let result = y_a + (y_b - y_a) * (x - x_a) / (x_b - x_a);
        result
    }
}

use nalgebra::{DMatrix, DVector};

pub struct CubicSpline<'a> {
    // func: fn(f64) -> f64,
    nodes: &'a [f64],
    values: Vec<f64>,
    b: DVector<f64>,
}

pub enum SplineType {
    /// Natural spline
    ///
    /// S''(x_0) = S''(x_n) = 0
    Natural,
    /// Clamped spline
    ///
    /// S'(x_0) = f'(x_0) and S'(x_n) = f'(x_n)
    Clamped(f64, f64),
}

impl CubicSpline<'_> {
    pub fn new(func: fn(f64) -> f64, nodes: &'_ [f64], spline_type: SplineType) -> CubicSpline<'_> {
        assert!(nodes.len() > 1);

        // Validate that nodes are sorted
        // let all_sorted = nodes.array_windows().all(|[a, b]| a < b);
        assert!(nodes.is_sorted(), "Nodes must be sorted!");

        let values: Vec<f64> = nodes.iter().copied().map(func).collect();
        dbg!(&values);

        let n = nodes.len();
        let mut matrix = DMatrix::zeros(n, n);
        let mut b = DVector::zeros(n);

        // Fill matrix
        for i in 1..n - 1 {
            matrix[(i, i - 1)] = (nodes[i] - nodes[i - 1]) / 6.0;
            matrix[(i, i)] = (nodes[i + 1] - nodes[i - 1]) / 3.0;
            matrix[(i, i + 1)] = (nodes[i + 1] - nodes[i]) / 6.0;
        }
        for i in 1..n - 1 {
            b[i] = (values[i + 1] - values[i]) / (nodes[i + 1] - nodes[i]) - (values[i] - values[i - 1]) / (nodes[i] - nodes[i - 1]);
        }

        match spline_type {
            SplineType::Natural => {
                matrix[(0, 0)] = 1.0;
                matrix[(n - 1, n - 1)] = 1.0;
            }
            SplineType::Clamped(f_a, f_b) => {
                matrix[(0, 0)] = 2.0 * (nodes[1] - nodes[0]);
                matrix[(0, 1)] = nodes[1] - nodes[0];
                matrix[(n - 1, n - 2)] = nodes[n - 1] - nodes[n - 2];
                matrix[(n - 1, n - 1)] = 2.0 * (nodes[n - 1] - nodes[n - 2]);
                b[0] = (values[1] - values[0]) / (nodes[1] - nodes[0]) - f_a;
                b[n - 1] = f_b - (values[n - 1] - values[n - 2]) / (nodes[n - 1] - nodes[n - 2]);
            }
        }
        let solve_result = matrix.lu().solve_mut(&mut b);
        assert!(solve_result);
        // let b = crate::task3::steepest_descent(matrix.into(), b.into());

        CubicSpline {
            // func,
            nodes,
            values,
            b,
        }
    }
}

impl Interpolation for CubicSpline<'_> {
    fn interpolate(&self, x: f64) -> f64 {
        let (a, b) = find_between(self.nodes, x);
        let x_a = self.nodes[a];
        let x_b = self.nodes[b];
        let y_a = self.values[a];
        let y_b = self.values[b];
        let h = x_b - x_a;
        let b_a = self.b[a];
        let b_b = self.b[b];

        let result =
            (1.0 / 6.0 / h) * (b_a * (x_b - x).powi(3) + b_b * (x - x_a).powi(3))
                + (y_a / h - b_a * h / 6.0) * (x_b - x)
                + (y_b / h - b_b * h / 6.0) * (x - x_a);
        result
    }
}

fn find_between(slice: &[f64], x: f64) -> (usize, usize) {
    match slice.iter().enumerate().find_map(|(i, value)| (*value > x).then_some(i)) {
        Some(0) => (0, 1),
        Some(non_zero) => (non_zero - 1, non_zero),
        None => (slice.len() - 2, slice.len() - 1),
    }
}
