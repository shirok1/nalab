use std::cell::Cell;
use std::ops::Range;
use const_soft_float::soft_f64::SoftF64;

pub fn composite_trapezoid(f: fn(f64) -> f64, zone: Range<f64>, partition: usize) -> f64 {
    let h = (zone.end - zone.start) / partition as f64;
    h * ((f(zone.start) + f(zone.end)) / 2.0
        + (1..partition).map(|i| zone.start + i as f64 * h).map(f).sum::<f64>())
}

pub fn composite_simpson(f: fn(f64) -> f64, zone: Range<f64>, partition: usize) -> f64 {
    let h = (zone.end - zone.start) / partition as f64;
    let half_h = h / 2.0;
    h / 6.0 * ((f(zone.start) + f(zone.end))
        + 4.0 * (0..partition).map(|i| zone.start + i as f64 * h + half_h).map(f).sum::<f64>()
        + 2.0 * (1..partition).map(|i| zone.start + i as f64 * h).map(f).sum::<f64>()
    )
}

fn romberg_rec(f: fn(f64) -> f64, zone: Range<f64>, partition: usize, partition_exp: usize, acceleration: usize, buffer: &mut [Option<f64>]) -> f64 {
    let line = partition_exp + acceleration;
    assert_eq!(buffer.len(), (line + 1) * (line + 2) / 2 - partition_exp);
    let (current, rest) = buffer.split_last_mut().unwrap();
    *(current.get_or_insert_with(|| {
        if acceleration == 0 {
            composite_trapezoid(f, zone, partition * 2usize.pow(partition_exp as u32))
        } else {
            let numerator = 4u32.pow(acceleration as u32);
            let denominator = numerator - 1;
            let rest_len = rest.len();
            ((numerator as f64) * romberg_rec(f, zone.clone(), partition, partition_exp + 1, acceleration - 1, rest)
                - romberg_rec(f, zone.clone(), partition, partition_exp, acceleration - 1, rest[0..rest_len - line].as_mut()))
                / (denominator as f64)
        }
    }))
}

fn romberg_2(
    f: fn(f64) -> f64,
    zone: Range<f64>,
    partition: usize,
    partition_exp: usize,
    last_line: &[f64],
    mut this_line: &mut [f64],
) {
    assert_eq!(last_line.len() + 1, this_line.len());
    this_line[0] = composite_trapezoid(f, zone, partition * 2usize.pow(partition_exp as u32));
    for (i, (t_m_minus_1_k, [t_m_minus_1_k_plus_1, t_m_k])) in
    last_line.iter().zip(Cell::from_mut(this_line).as_slice_of_cells().array_windows()).enumerate() {
        let acceleration = i + 1;
        let numerator = 4u32.pow(acceleration as u32);
        let denominator = numerator - 1;
        t_m_k.set(
            ((numerator as f64) * t_m_minus_1_k_plus_1.get() - t_m_minus_1_k) / (denominator as f64)
        );
    }
}

pub fn simple_romberg_2(f: fn(f64) -> f64, zone: Range<f64>, init_partition: usize, acceleration: usize) -> f64 {
    let mut buffer_a = vec![0.; acceleration + 1];
    let mut buffer_b = vec![0.; acceleration + 1];

    let mut last_line = &mut buffer_a[..];
    let mut this_line = &mut buffer_b[..];

    for i in 0..=acceleration {
        romberg_2(f, zone.clone(), init_partition, i, &last_line[0..i], &mut this_line[0..i + 1]);
        std::mem::swap(&mut last_line, &mut this_line);
    }
    last_line[acceleration]
}

pub fn simple_romberg(f: fn(f64) -> f64, zone: Range<f64>, init_partition: usize, acceleration: usize) -> f64 {
    let mut buffer = vec![None; (acceleration + 1) * (acceleration + 2) / 2];
    romberg_rec(f, zone, init_partition, 0, acceleration, buffer.as_mut_slice())
}

pub fn auto_romberg(f: fn(f64) -> f64, zone: Range<f64>, init_partition: usize, eps: f64) -> f64 {
    let mut buffer = vec![None; 3 * 4 / 2];
    let mut acceleration = 0;
    let mut last = romberg_rec(f, zone.clone(), init_partition, 0, acceleration, buffer[0..1].as_mut());
    loop {
        acceleration += 1;
        dbg!(acceleration);
        buffer.resize((acceleration + 1) * (acceleration + 2) / 2, None);
        let current = romberg_rec(f, zone.clone(), init_partition, 0, acceleration, buffer.as_mut_slice());
        if (current - last).abs() < eps {
            return current;
        }
        last = current;
    }
}

fn composite_gauss(f: fn(f64) -> f64, zone: Range<f64>, partition: usize, weights_and_nodes: &[(f64, f64)]) -> f64 {
    let per_zone_h = (zone.end - zone.start) / partition as f64;
    let per_zone_h_half = per_zone_h / 2.0;

    // let zone_nodes = (0..=partition).map(|i| zone.start + i as f64 * per_zone_h).collect::<Vec<_>>();
    // let zone_nodes = (0..partition).map(|i| zone.start + i as f64 * per_zone_h);
    let zone_nodes = (0..partition).map(|i| zone.start + i as f64 * per_zone_h).collect::<Vec<_>>();

    // zone_nodes.array_windows()
    //     .map(|&[start, _]| per_zone_h_half * weights_and_nodes.iter()
    //         .map(|(w, t)| w * f(per_zone_h_half * (t + 1.0) + start)).sum::<f64>()).sum::<f64>()
    // per_zone_h_half * zone_nodes.array_windows()
    //     .map(|&[start, _]| weights_and_nodes.iter()
    //         .map(|(w, t)| w * f(per_zone_h_half * (t + 1.0) + start)).sum::<f64>()).sum::<f64>()
    // per_zone_h_half * zone_nodes.map(|start| weights_and_nodes.iter()
    //     .map(|(w, t)| w * f(per_zone_h_half * (t + 1.0) + start)).sum::<f64>()).sum::<f64>()

    per_zone_h_half * weights_and_nodes.iter()
        .map(|(w, t)| zone_nodes.iter()
            .map(|start| w * f(per_zone_h_half * (t + 1.0) + start)).sum::<f64>()).sum::<f64>()
}

const GAUSS_2_WEIGHTS_AND_NODES: [(f64, f64); 2] = [
    (1.0, -1.0 / SoftF64(3.0).sqrt().to_f64()),
    (1.0, 1.0 / SoftF64(3.0).sqrt().to_f64())
];

pub fn composite_gauss_2(f: fn(f64) -> f64, zone: Range<f64>, partition: usize) -> f64 {
    composite_gauss(f, zone, partition, &GAUSS_2_WEIGHTS_AND_NODES)
}

const GAUSS_4_WEIGHTS_AND_NODES: [(f64, f64); 4] = {
    const LEFT: f64 = 3.0 / 7.0;
    const RIGHT_SQUARE: f64 = (2.0 * 2.0 * 6.0) / (5.0 * 7.0 * 7.0);
    const RIGHT: f64 = SoftF64(RIGHT_SQUARE).sqrt().to_f64();
    [
        ((18.0 - SoftF64(30.0).sqrt().to_f64()) / 36.0, -SoftF64(LEFT + RIGHT).sqrt().to_f64()),
        ((18.0 + SoftF64(30.0).sqrt().to_f64()) / 36.0, -SoftF64(LEFT - RIGHT).sqrt().to_f64()),
        ((18.0 + SoftF64(30.0).sqrt().to_f64()) / 36.0, SoftF64(LEFT - RIGHT).sqrt().to_f64()),
        ((18.0 - SoftF64(30.0).sqrt().to_f64()) / 36.0, SoftF64(LEFT + RIGHT).sqrt().to_f64()),
    ]
};

pub fn composite_gauss_4(f: fn(f64) -> f64, zone: Range<f64>, partition: usize) -> f64 {
    composite_gauss(f, zone, partition, &GAUSS_4_WEIGHTS_AND_NODES)
}
