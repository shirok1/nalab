#![feature(core_intrinsics)]

use std::env;
use std::intrinsics::type_name;
use nalgebra::{Matrix3, UniformNorm, Vector3};
use plotters::prelude::IntoLinspace;
use nalab::{task1, task2, task3};
use nalab::task1::{chebyshev_zeros, CubicSpline, Interpolation, Lagrange, Linear};

use plotters::prelude::*;

fn task1() -> Result<(), Box<dyn std::error::Error>> {
    let f = |x| -> f64 { 1.0 / (1.0 + x * x) };
    // let zone = -5.0..5.0;
    let naive_nodes = [-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0];

    let lagrange = Lagrange::new(f, &naive_nodes);
    let linear = Linear::new(f, &naive_nodes);

    let mut chebyshev_nodes = chebyshev_zeros(11);
    chebyshev_nodes.iter_mut().for_each(|x| *x *= 5.0);
    dbg!(&chebyshev_nodes);
    let chebyshev_lagrange = Lagrange::new(f, chebyshev_nodes.as_slice());

    // let chebyshev_lagrange = task1::CubicSpline::new(f,  chebyshev_nodes.as_slice(), task1::SplineType::Clamped(-5.0, 5.0));
    let cubic_spline = CubicSpline::new(f, &naive_nodes, task1::SplineType::Clamped(-5.0, 5.0));

    // let root = BitMapBackend::new("lagrange.png", (640, 1024)).into_drawing_area();
    let root = SVGBackend::new("interpolate.svg", (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("插值", ("Source Han Serif SC", 24).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(0.0..5.0, -0.5..2.0)?;

    chart.configure_mesh().draw()?;

    let lagrange_series = LineSeries::new(
        (0..=100).map(|x| x as f64 / 20.0)
            .map(|x| (x, lagrange.interpolate(x))),
        &RED,
    );
    let chebyshev_lagrange_series = LineSeries::new(
        (0..=100).map(|x| x as f64 / 20.0)
            .map(|x| (x, chebyshev_lagrange.interpolate(x))),
        &MAGENTA,
    );
    let linear_series = LineSeries::new(
        (0..=100).map(|x| x as f64 / 20.0)
            .map(|x| (x, linear.interpolate(x))),
        &GREEN,
    );
    let cubic_spline_series = LineSeries::new(
        (0..=100).map(|x| x as f64 / 20.0)
            .map(|x| (x, cubic_spline.interpolate(x))),
        &CYAN,
    );
    let truth_series = LineSeries::new(
        (0..=100).map(|x| x as f64 / 20.0)
            .map(|x| (x, f(x))),
        &BLUE,
    );
    chart
        .draw_series(lagrange_series)?
        .label(type_name::<Lagrange>())
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
    chart
        .draw_series(chebyshev_lagrange_series)?
        .label(format!("{} + chebyshev_nodes", type_name::<Lagrange>()))
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &MAGENTA));
    chart
        .draw_series(linear_series)?
        .label(type_name::<Linear>())
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));
    chart
        .draw_series(cubic_spline_series)?
        .label(type_name::<CubicSpline>())
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &CYAN));
    chart
        .draw_series(truth_series)?
        .label("truth")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .label_font(("Source Han Serif SC", 12).into_font())
        .draw()?;

    root.present()?;

    let nodes = (0..=100).map(|x| -0.5 + x as f64 / 10.0).collect::<Vec<_>>();
    let truth = nodes.iter().copied().map(f).collect::<Vec<_>>();
    let lagrange_max_error = nodes.iter().copied().map(|x| lagrange.interpolate(x)).zip(truth.iter()).map(|(x, y)| (x - y).abs())
        .max_by(|x, y| x.total_cmp(y)).unwrap();
    let chebyshev_lagrange_max_error = nodes.iter().copied().map(|x| chebyshev_lagrange.interpolate(x)).zip(truth.iter()).map(|(x, y)| (x - y).abs())
        .max_by(|x, y| x.total_cmp(y)).unwrap();
    let linear_max_error = nodes.iter().copied().map(|x| linear.interpolate(x)).zip(truth.iter()).map(|(x, y)| (x - y).abs())
        .max_by(|x, y| x.total_cmp(y)).unwrap();
    let cubic_spline_max_error = nodes.iter().copied().map(|x| cubic_spline.interpolate(x)).zip(truth.iter()).map(|(x, y)| (x - y).abs())
        .max_by(|x, y| x.total_cmp(y)).unwrap();
    println!("lagrange_max_error: {:e}", lagrange_max_error);
    println!("chebyshev_lagrange_max_error: {:e}", chebyshev_lagrange_max_error);
    println!("linear_max_error: {:e}", linear_max_error);
    println!("cubic_spline_max_error: {:e}", cubic_spline_max_error);


    Ok(())
}

fn task2() {
    println!("TASK 2 Numerical Integration");
    let f = |x: f64| (2.0 * x).exp();
    let zone = 3.0..5.0;
    let ground_truth = (f(zone.end) - f(zone.start)) / 2.0;
    // println!("ground truth: {}", ground_truth);
    dbg!(ground_truth);

    let result_1 = task2::composite_trapezoid(f, zone.clone(), 40);
    dbg!(result_1);
    dbg!((result_1 - ground_truth).abs());

    let result_2 = task2::composite_simpson(f, zone.clone(), 20);
    dbg!(result_2);
    dbg!((result_2 - ground_truth).abs());

    let result_3 = task2::simple_romberg(f, zone.clone(), 20, 1);
    let result_4 = task2::simple_romberg(f, zone.clone(), 10, 2);
    let result_5 = task2::simple_romberg(f, zone.clone(), 5, 3);

    dbg!(result_3);
    dbg!((result_3 - ground_truth).abs());
    dbg!(result_4);
    dbg!((result_4 - ground_truth).abs());
    dbg!(result_5);
    dbg!((result_5 - ground_truth).abs());

    let result_3 = task2::simple_romberg_2(f, zone.clone(), 20, 1);
    let result_4 = task2::simple_romberg_2(f, zone.clone(), 10, 2);
    let result_5 = task2::simple_romberg_2(f, zone.clone(), 5, 3);

    dbg!(result_3);
    dbg!((result_3 - ground_truth).abs());
    dbg!(result_4);
    dbg!((result_4 - ground_truth).abs());
    dbg!(result_5);
    dbg!((result_5 - ground_truth).abs());

    let result_6 = task2::composite_gauss_2(f, zone.clone(), 20);
    dbg!(result_6);
    dbg!((result_6 - ground_truth).abs());
    let result_7 = task2::composite_gauss_4(f, zone.clone(), 40);
    dbg!(result_7);
    dbg!((result_7 - ground_truth).abs());


    // let result_auto = task2::auto_romberg(f, zone.clone(), 20, 0.00000000001);
    // dbg!(result_auto);
    // dbg!((result_auto - ground_truth).abs());
}

fn task3() {
    let a = Matrix3::new(3., 1., 1.,
                         1., 3., 1.,
                         1., 1., 3.);
    let b = Vector3::new(3., 0., 2.);
    let truth = Vector3::new(1., -0.5, 0.5);
    let jacobi_result = task3::jacobi(a, b);
    dbg!(jacobi_result);
    dbg!((jacobi_result - truth).apply_norm(&UniformNorm));
    let gauss_seidel_result = task3::gauss_seidel(a, b);
    dbg!(gauss_seidel_result);
    dbg!((gauss_seidel_result - truth).apply_norm(&UniformNorm));
    let sor_result = task3::successive_over_relaxation(a, b, task3::omega_opt(1. / 3.));
    dbg!(sor_result);
    dbg!((sor_result - truth).apply_norm(&UniformNorm));
    let steepest_descent_result = task3::steepest_descent(a, b);
    dbg!(steepest_descent_result);
    dbg!((steepest_descent_result - truth).apply_norm(&UniformNorm));
}

fn main() {
    let args: Vec<String> = env::args().collect();
    assert_eq!(args.len(), 2);
    match args[1].as_str() {
        "task1" => task1().unwrap(),
        "task2" => task2(),
        "task3" => task3(),
        _ => panic!("unknown task, enter one of task1, task2, task3"),
    }
}