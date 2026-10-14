//! Orders of convergence, against functions known in closed form
//! (docs/design.ja.md section 8.4.2).
//!
//! The data of a smooth function is interpolated on finer and finer points,
//! and the largest error is measured. Halving the spacing `h` divides the
//! error of a method of order `p` by `2^p`, so `p` is read off two errors:
//! `log2(e(h) / e(h/2))`. A method that reproduces its data and its
//! polynomials can still have the wrong order (a spline with the wrong end
//! condition goes through every point and loses two orders), and no test of
//! one grid sees it.
//!
//! The orders asserted are those of the design. The widths around them are
//! fixed here, from the orders measured when the tests were written (given at
//! each case): wide enough for the spacings used, which are not yet in the
//! limit, and narrower than the distance to the next order.

mod common;

use common::functions::genz;
use interpx::{
    Akima1D, AkimaMethod, BSpline1D, BoundaryCondition, CubicHermiteSpline, CubicSpline,
    GridMethod, Interpolator, Linear1D, Nearest1D, NearestKind, PchipInterpolator,
    RegularGridInterpolator,
};

/// The right end of the interval the functions of one variable are on.
const END: f64 = 3.0;

/// A smooth function with no symmetry and with turning points: a sine and an
/// exponential.
fn f(x: f64) -> f64 {
    (1.3 * x + 0.4).sin() + 0.5 * (0.7 * x).exp()
}

/// Its first derivative.
fn df(x: f64) -> f64 {
    1.3 * (1.3 * x + 0.4).cos() + 0.35 * (0.7 * x).exp()
}

/// Its second derivative.
fn d2f(x: f64) -> f64 {
    -1.69 * (1.3 * x + 0.4).sin() + 0.245 * (0.7 * x).exp()
}

/// Its integral from 0 to `x`.
fn integral(x: f64) -> f64 {
    ((0.4_f64).cos() - (1.3 * x + 0.4).cos()) / 1.3 + ((0.7 * x).exp() - 1.0) / 1.4
}

/// `n + 1` points from 0 to [`END`]: evenly spaced, or moved smoothly so that
/// the widest interval is about three times the narrowest (the orders hold
/// for points that are not evenly spaced, with the widest interval for `h`).
fn points(n: usize, even: bool) -> Vec<f64> {
    (0..=n)
        .map(|i| {
            let t = END * i as f64 / n as f64;
            if even {
                t
            } else {
                t + 0.25 * (2.0 * std::f64::consts::PI * t / END).sin()
            }
        })
        .collect()
}

/// The largest distance between `got` and `want` on 3001 points of the
/// interval.
fn worst(got: impl Fn(f64) -> f64, want: impl Fn(f64) -> f64) -> f64 {
    (0..=3000)
        .map(|i| {
            let x = END * f64::from(i) / 3000.0;
            (got(x) - want(x)).abs()
        })
        .fold(0.0, f64::max)
}

/// The orders between each two of `errors`, which belong to spacings each
/// half the one before.
fn orders(errors: &[f64]) -> Vec<f64> {
    errors
        .windows(2)
        .map(|pair| (pair[0] / pair[1]).log2())
        .collect()
}

/// Holds the last two orders between `errors` against `expected`, give or
/// take `width`.
#[cfg(test)]
fn hold(what: &str, errors: &[f64], expected: f64, width: f64) {
    let orders = orders(errors);
    assert!(
        errors.iter().all(|e| e.is_finite() && *e > 0.0),
        "{what}: {errors:?}"
    );
    for order in &orders[orders.len() - 2..] {
        assert!(
            (order - expected).abs() <= width,
            "{what}: orders {orders:.2?} where {expected} is expected (errors {errors:?})"
        );
    }
}

/// The errors of an interpolant of `f` built by `build` on `sizes` intervals,
/// each twice the one before, measured by `error`.
fn errors<I>(
    sizes: &[usize],
    even: bool,
    build: impl Fn(&[f64], &[f64]) -> I,
    error: impl Fn(&I) -> f64,
) -> Vec<f64> {
    sizes
        .iter()
        .map(|n| {
            let x = points(*n, even);
            let y: Vec<f64> = x.iter().map(|x| f(*x)).collect();
            error(&build(&x, &y))
        })
        .collect()
}

/// The spacings of the cases of one variable.
const SIZES: [usize; 4] = [32, 64, 128, 256];

#[test]
fn the_types_of_one_variable_have_the_orders_of_the_design() {
    for even in [true, false] {
        let spacing = if even { "even" } else { "uneven" };
        let value = |eval: &dyn Fn(f64) -> Option<f64>| worst(|x| eval(x).unwrap_or(f64::NAN), f);

        // Measured: 1.98 to 1.99.
        let linear = errors(
            &SIZES,
            even,
            |x, y| Linear1D::new(x, y).expect("valid"),
            |i| value(&|x| i.eval(&x)),
        );
        hold(&format!("Linear1D, {spacing}"), &linear, 2.0, 0.1);

        // The nearest value is off by the slope times half an interval.
        // Measured: 0.97 to 1.02 (the largest error moves from interval to
        // interval, so the order is not as steady as the others).
        let nearest = errors(
            &SIZES,
            even,
            |x, y| Nearest1D::new(x, y, NearestKind::Nearest).expect("valid"),
            |i| value(&|x| i.eval(&x)),
        );
        hold(&format!("Nearest1D, {spacing}"), &nearest, 1.0, 0.15);

        // Measured: 3.96 to 3.99.
        let spline = |x: &[f64], y: &[f64]| {
            CubicSpline::new(x, y, BoundaryCondition::NotAKnot).expect("valid")
        };
        let cubic = errors(&SIZES, even, spline, |i| value(&|x| i.eval(&x)));
        hold(&format!("CubicSpline, {spacing}"), &cubic, 4.0, 0.1);

        // A derivative loses an order each time. Measured: 2.96 to 2.99, and
        // 1.96 to 2.0.
        let slope = errors(&SIZES, even, spline, |i| {
            worst(
                |x| i.as_ppoly().eval_derivative(x, 1).unwrap_or(f64::NAN),
                df,
            )
        });
        hold(&format!("CubicSpline, slope, {spacing}"), &slope, 3.0, 0.1);
        let bend = errors(&SIZES, even, spline, |i| {
            worst(
                |x| i.as_ppoly().eval_derivative(x, 2).unwrap_or(f64::NAN),
                d2f,
            )
        });
        hold(
            &format!("CubicSpline, second derivative, {spacing}"),
            &bend,
            2.0,
            0.1,
        );

        // The integral over the whole interval has the order of the values.
        // Measured: 3.97 to 4.04.
        let area = errors(&SIZES, even, spline, |i| {
            (i.as_ppoly().integrate(0.0, END).unwrap_or(f64::NAN) - integral(END)).abs()
        });
        hold(
            &format!("CubicSpline, integral, {spacing}"),
            &area,
            4.0,
            0.15,
        );

        // The natural end condition sets a second derivative to zero that
        // the function does not have: two orders are lost at the ends. (The
        // not-a-knot spline above must not be this one.) Measured: 2.00.
        let natural = errors(
            &SIZES,
            even,
            |x, y| CubicSpline::new(x, y, BoundaryCondition::natural()).expect("valid"),
            |i| value(&|x| i.eval(&x)),
        );
        hold(
            &format!("CubicSpline, natural, {spacing}"),
            &natural,
            2.0,
            0.1,
        );

        // Measured: 2.97 to 2.99.
        let pchip = errors(
            &SIZES,
            even,
            |x, y| PchipInterpolator::new(x, y).expect("valid"),
            |i| value(&|x| i.eval(&x)),
        );
        hold(&format!("PchipInterpolator, {spacing}"), &pchip, 3.0, 0.1);

        // Akima's slopes are weighted means of secants: second order on
        // points of any spacing, and so third order for the values. The
        // weights change from point to point, and the order with them.
        // Measured: 2.82 to 3.14 (Akima), 2.91 to 2.96 (makima).
        for (method, width) in [(AkimaMethod::Akima, 0.3), (AkimaMethod::Makima, 0.15)] {
            let akima = errors(
                &SIZES,
                even,
                |x, y| Akima1D::new(x, y, method).expect("valid"),
                |i| value(&|x| i.eval(&x)),
            );
            hold(
                &format!("Akima1D, {method:?}, {spacing}"),
                &akima,
                3.0,
                width,
            );
        }

        // With the slopes of the function itself, a cubic on each interval
        // is as good as a cubic gets. Measured: 4.00.
        let hermite = errors(
            &SIZES,
            even,
            |x, y| {
                let slopes: Vec<f64> = x.iter().map(|x| df(*x)).collect();
                CubicHermiteSpline::new(x, y, &slopes).expect("valid")
            },
            |i| value(&|x| i.eval(&x)),
        );
        hold(
            &format!("CubicHermiteSpline, {spacing}"),
            &hermite,
            4.0,
            0.05,
        );

        // Degree k: order k + 1. The spacings stop where the error of the
        // highest degree is still above the rounding of the solve (1.5e-11
        // at 128 intervals). Measured: within 0.1 of k + 1 for degrees 1 to
        // 4, 5.90 and 5.96 for degree 5 (5.72 and 5.90 on uneven points).
        for degree in 1..=5 {
            let spline = errors(
                &[16, 32, 64, 128],
                even,
                |x, y| BSpline1D::new(x, y, degree, BoundaryCondition::NotAKnot).expect("valid"),
                |i| value(&|x| i.eval(&x)),
            );
            let width = if degree == 5 { 0.3 } else { 0.2 };
            hold(
                &format!("BSpline1D, degree {degree}, {spacing}"),
                &spline,
                degree as f64 + 1.0,
                width,
            );
        }
    }
}

/// The errors of a grid of `f(x) f(y)` with `n` intervals an axis, for each
/// `n` of `sizes`.
#[cfg(test)]
fn grid_errors(sizes: &[usize], even: bool, method: GridMethod) -> Vec<f64> {
    sizes
        .iter()
        .map(|n| {
            let x = points(*n, even);
            let values: Vec<f64> = x
                .iter()
                .flat_map(|a| x.iter().map(move |b| f(*a) * f(*b)))
                .collect();
            let grid = RegularGridInterpolator::new([&x[..], &x[..]], &values, method);
            let grid = grid.expect("valid");
            let mut worst = 0.0_f64;
            for i in 0..=150 {
                for j in 0..=150 {
                    let (a, b) = (END * f64::from(i) / 150.0, END * f64::from(j) / 150.0);
                    let got = grid.eval(&[a, b]).unwrap_or(f64::NAN);
                    worst = worst.max((got - f(a) * f(b)).abs());
                }
            }
            worst
        })
        .collect()
}

#[test]
fn the_grid_methods_have_the_orders_of_the_design() {
    // In two dimensions the spacings stay coarser than in one, and the
    // orders are farther from their limit. Measured at the last two
    // halvings, on even and on uneven axes: nearest 0.94 to 1.09, linear
    // 1.87 to 1.96, cubic 3.60 to 3.93, quintic 5.05 to 5.87, pchip 2.85 to
    // 3.03.
    let cases = [
        (GridMethod::Nearest, 1.0, 0.2),
        (GridMethod::Linear, 2.0, 0.2),
        // The lowest order measured for the cubic is 3.60, on uneven axes:
        // the width leaves 0.1 under it.
        (GridMethod::Cubic, 4.0, 0.5),
        (GridMethod::Pchip, 3.0, 0.2),
    ];
    for even in [true, false] {
        for (method, expected, width) in cases {
            let errors = grid_errors(&[8, 16, 32, 64], even, method);
            hold(
                &format!("grid, {method:?}, even {even}"),
                &errors,
                expected,
                width,
            );
        }
        // The quintic is still coming up to its order at these spacings
        // (5.42, 5.77, 5.87 on even axes; 4.55, 5.05, 5.70 on uneven ones):
        // the last order is held, and that the orders rise toward it.
        let errors = grid_errors(&[8, 16, 32, 64], even, GridMethod::Quintic);
        let orders = orders(&errors);
        assert!(
            orders.windows(2).all(|pair| pair[0] < pair[1]),
            "quintic: {orders:.2?}"
        );
        assert!((orders[2] - 6.0).abs() <= 0.35, "quintic: {orders:.2?}");
    }
}

#[test]
fn a_grid_of_three_dimensions_has_the_orders_of_its_method() {
    // Genz's oscillatory function on the unit cube (design section 8.4.4).
    // Three samples to the narrowest cell, and not a multiple of the cells:
    // with one sample a cell or fewer the largest error (in the end cells
    // of a not-a-knot spline) falls between the samples, and the last order
    // read 5.9 and 7.0 where it is 4.
    const SAMPLES: u32 = 49;
    let genz = |p: [f64; 3]| genz("oscillatory", &p).expect("a family");
    let errors = |method: GridMethod, sizes: &[usize]| -> Vec<f64> {
        sizes
            .iter()
            .map(|n| {
                let axis: Vec<f64> = (0..=*n).map(|i| i as f64 / *n as f64).collect();
                let mut values = Vec::with_capacity(axis.len().pow(3));
                for a in &axis {
                    for b in &axis {
                        for c in &axis {
                            values.push(genz([*a, *b, *c]));
                        }
                    }
                }
                let grid = RegularGridInterpolator::new([&axis[..]; 3], &values, method);
                let grid = grid.expect("valid");
                let mut worst = 0.0_f64;
                for i in 0..=SAMPLES {
                    for j in 0..=SAMPLES {
                        for k in 0..=SAMPLES {
                            let p = [i, j, k].map(|m| f64::from(m) / f64::from(SAMPLES));
                            let got = grid.eval(&p).unwrap_or(f64::NAN);
                            worst = worst.max((got - genz(p)).abs());
                        }
                    }
                }
                worst
            })
            .collect()
    };
    hold(
        "grid of three dimensions, linear",
        &errors(GridMethod::Linear, &[4, 8, 16]),
        2.0,
        0.2,
    );
    hold(
        "grid of three dimensions, cubic",
        &errors(GridMethod::Cubic, &[4, 8, 16]),
        4.0,
        0.5,
    );
}

/// Point `i` (from 1) of the Halton sequence of the bases 2 and 3: points in
/// the unit square that fill it evenly and lie on no grid (design section
/// 8.4.5).
#[cfg(feature = "scattered")]
fn halton(i: usize) -> [f64; 2] {
    [2, 3].map(|base| {
        let (mut share, mut at, mut i) = (1.0, 0.0, i);
        while i > 0 {
            share /= base as f64;
            at += share * (i % base) as f64;
            i /= base;
        }
        at
    })
}

#[cfg(feature = "scattered")]
#[test]
fn scattered_interpolants_come_nearer_the_function_with_more_points() {
    use common::functions::franke;
    use interpx::{NearestNDInterpolator, RBFInterpolator, RbfKernel, RbfOptions};

    // Franke's first function, on Halton points; the error on a lattice
    // inside the square (away from its edges, where there is data on one
    // side only).
    let data = |n: usize| -> (Vec<[f64; 2]>, Vec<f64>) {
        let points: Vec<[f64; 2]> = (1..=n).map(halton).collect();
        let values = points
            .iter()
            .map(|p| franke(1, p[0], p[1]).expect("F1"))
            .collect();
        (points, values)
    };
    let worst = |eval: &dyn Fn([f64; 2]) -> Option<f64>| {
        let mut worst = 0.0_f64;
        for i in 0..=40 {
            for j in 0..=40 {
                let p = [i, j].map(|m| 0.1 + 0.8 * f64::from(m) / 40.0);
                let got = eval(p).unwrap_or(f64::NAN);
                worst = worst.max((got - franke(1, p[0], p[1]).expect("F1")).abs());
            }
        }
        worst
    };

    // Scattered points have no spacing to halve, and an RBF interpolant has
    // no one order: it depends on the kernel, and for the kernels with a
    // shape parameter on that too. What is held is that eight times the
    // points bring the error down by the factor given, and under the bound
    // given. Measured, from 100 points to 800 (errors at 800 in brackets):
    //
    //   kernel                 all points          30 neighbors
    //   linear                 12x  (4.0e-3)       9x   (5.2e-3)
    //   thin plate spline      26x  (7.3e-4)       22x  (8.6e-4)
    //   cubic                  95x  (2.0e-4)       44x  (4.5e-4)
    //   quintic                610x (3.1e-5)       110x (1.9e-4)
    //   multiquadric           2e4x (8.2e-7)       290x (7.0e-5)
    //   inverse multiquadric   4e3x (4.2e-6)       110x (1.7e-4)
    //   inverse quadratic      2e3x (1.2e-5)       67x  (2.7e-4)
    //   Gaussian               180x (9.4e-5)       350x (5.4e-5)
    //
    // The factors asserted are half of these or less, and the bounds three
    // times the errors. The Gaussian over all points is held to less: its
    // system is ill-conditioned at these sizes (coefficients of 1e9), and
    // its error moves with the order of the points: at 800 points between
    // 7.3e-6 and 9.4e-5 over 12 orders (SciPy's between 7.0e-6 and 6.2e-4
    // over 24; #159).
    let shape = 5.0;
    let cases = [
        (RbfKernel::Linear, 5.0, 1.2e-2, 4.0, 1.6e-2),
        (RbfKernel::ThinPlateSpline, 10.0, 2.2e-3, 10.0, 2.6e-3),
        (RbfKernel::Cubic, 40.0, 6.0e-4, 20.0, 1.4e-3),
        (RbfKernel::Quintic, 200.0, 1.0e-4, 50.0, 6.0e-4),
        (
            RbfKernel::Multiquadric { epsilon: shape },
            1000.0,
            2.5e-6,
            100.0,
            2.1e-4,
        ),
        (
            RbfKernel::InverseMultiquadric { epsilon: shape },
            1000.0,
            1.3e-5,
            50.0,
            5.1e-4,
        ),
        (
            RbfKernel::InverseQuadratic { epsilon: shape },
            500.0,
            3.6e-5,
            30.0,
            8.1e-4,
        ),
        (
            RbfKernel::Gaussian { epsilon: shape },
            20.0,
            1.0e-3,
            100.0,
            1.7e-4,
        ),
    ];
    for (kernel, all_factor, all_bound, near_factor, near_bound) in cases {
        for (options, factor, bound) in [
            (RbfOptions::new(), all_factor, all_bound),
            (
                RbfOptions::new().with_neighbors(30),
                near_factor,
                near_bound,
            ),
        ] {
            let error = |n: usize| {
                let (points, values) = data(n);
                let rbf = RBFInterpolator::with_options(&points, &values, kernel, &options);
                let rbf = rbf.expect("valid");
                worst(&|p| rbf.eval(&p))
            };
            let (few, many) = (error(100), error(800));
            assert!(
                many <= bound && few / many >= factor,
                "{kernel:?}, {options:?}: {few:.1e} at 100 points, {many:.1e} at 800"
            );
        }
    }

    // The nearest value is off by the slope times the distance to the
    // nearest point, which goes down as the square root of the number of
    // points: 16 times the points, a quarter of the error at best. But
    // the largest error is at the steepest place and moves with the points,
    // so the factor held is 2.
    let error = |n: usize| {
        let (points, values) = data(n);
        let nearest = NearestNDInterpolator::new(&points, &values).expect("valid");
        worst(&|p| nearest.eval(&p))
    };
    let (few, many) = (error(100), error(1600));
    assert!(
        few / many >= 2.0 && many <= 0.15,
        "nearest: {few:.1e}, {many:.1e}"
    );
}
