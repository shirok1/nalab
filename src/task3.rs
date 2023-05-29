use nalgebra::{Matrix3, SMatrix, SVector, UniformNorm, Vector3};

const ITERATION_THRESHOLD: f64 = 1e-8;

type SqMat<const S: usize> = SMatrix<f64, S, S>;
type Vector<const S: usize> = SVector<f64, S>;

fn dlu<const S: usize>(a: SqMat<S>) -> (SqMat<S>, SqMat<S>, SqMat<S>) {
    (
        SqMat::<S>::from_diagonal(&a.diagonal()),
        {
            let mut res = -a;
            res.fill_upper_triangle(0., 0);
            res
        },
        {
            let mut res = -a;
            res.fill_lower_triangle(0., 0);
            res
        },
    )
}

pub fn jacobi<const S: usize>(a: SqMat<S>, b: Vector<S>) -> Vector<S> {
    let (d, l, u) = dlu(a);

    let d_inv = d.try_inverse().unwrap();
    let j = d_inv * (l + u);
    let f = d_inv * b;

    abstract_iteration(|x| j * x + f)
}

pub fn gauss_seidel<const S: usize>(a: SqMat<S>, b: Vector<S>) -> Vector<S> {
    let (d, l, u) = dlu(a);

    let d_l_inv = (d - l).try_inverse().unwrap();
    let g = d_l_inv * u;
    let f = d_l_inv * b;

    abstract_iteration(|x| g * x + f)
}

pub fn omega_opt(spectral_radius: f64) -> f64 {
    2.0 / (1.0 + (1.0 - spectral_radius.powi(2)).sqrt())
}

pub fn successive_over_relaxation<const S: usize>(a: SqMat<S>, b: Vector<S>, omega: f64) -> Vector<S> {
    dbg!(omega);
    let (d, l, u) = dlu(a);

    let d_omega_l_inv = (d - omega * l).try_inverse().unwrap();
    let l_omega = d_omega_l_inv * ((1.0 - omega) * d + omega * u);
    let f = d_omega_l_inv * omega * b;

    abstract_iteration(|x| l_omega * x + f)
}

pub fn steepest_descent<const S: usize>(a: SqMat<S>, b: Vector<S>) -> Vector<S> {
    abstract_iteration(|x| {
        let r = b - a * x;
        let alpha = r.dot(&r) / r.dot(&(a * r));
        x + alpha * r
    })
}

fn abstract_iteration<const S: usize>(iter: impl Fn(Vector<S>) -> Vector<S>) -> Vector<S> {
    let mut x = Vector::<S>::zeros();
    let mut iteration = 0;
    loop {
        let x_new = iter(x);
        iteration += 1;
        let epsilon = x_new - x;
        let epsilon_norm = epsilon.apply_norm(&UniformNorm);
        println!("Iteration: {}, x: {:?}, epsilon: {:e}", iteration, x_new, epsilon_norm);
        assert!(epsilon_norm.is_finite());
        if epsilon_norm < ITERATION_THRESHOLD {
            return x_new;
        }
        x = x_new;
    }
}