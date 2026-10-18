//! The interpx side of the comparison: one case, the phases of
//! `perf-harness`. Every evaluation is the serial batch, `eval_batch_into`,
//! into an array that is there already; or, when the process has `par` as
//! its second argument, the parallel batch, `eval_batch_into_par`.

use interpx::{
    Akima1D, AkimaMethod, BSpline, BSpline1D, BSplineCurve, BoundaryCondition,
    CloughTocher2DInterpolator, CubicSpline, GridMethod, Interpolator, Linear1D,
    LinearNDInterpolator, Nearest1D, NearestKind, NearestNDInterpolator, PPoly, ParInterpolator,
    PchipInterpolator, RBFInterpolator, RbfKernel, RbfOptions, RegularGridInterpolator,
};
use perf_harness::{run, run_given, unsupported, Case};
use std::cell::RefCell;

const LIBRARY: &str = "interpx";
/// The same library, built with `from_owned`: the data is handed over and
/// not copied. Run when the process has `owned` as its second argument, for
/// the types that have such a constructor.
const OWNED: &str = "interpx_owned";

fn owned() -> bool {
    std::env::args().nth(2).as_deref() == Some("owned")
}

/// The same library, evaluated with the parallel batch on the threads the
/// driver allows (`RAYON_NUM_THREADS`).
const PARALLEL: &str = "interpx_par";

/// Read once: this is asked at every batch, and reading the arguments of the
/// process takes a microsecond, which is ten times a small batch.
fn parallel() -> bool {
    static MODE: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *MODE.get_or_init(|| std::env::args().nth(2).as_deref() == Some("par"))
}

/// The name the results of this process go under.
fn library() -> &'static str {
    if parallel() {
        PARALLEL
    } else {
        LIBRARY
    }
}

/// The batch of the mode this process runs in.
fn batch<I>(built: &I, q: &[I::Point], out: &mut [I::Output], fill: I::Output)
where
    I: Interpolator + Sync,
    I::Point: Sync,
    I::Output: Send + Sync + Clone,
{
    let done = if parallel() {
        built.eval_batch_into_par(q, out, fill)
    } else {
        built.eval_batch_into(q, out, fill)
    };
    done.expect("same length");
}

/// Runs a scalar interpolator of one variable over the queries `q`.
fn scalar<I>(case: &Case, build: impl Fn() -> I)
where
    I: Interpolator<Point = f64, Output = f64> + Sync,
{
    let q = case.array("q");
    run(case, library(), q.len(), build, |built, out| {
        batch(built, &q, out, f64::NAN);
    });
}

fn oned(case: &Case) {
    let (x, y) = (case.array("x"), case.array("y"));
    let spline = |degree| {
        scalar(case, || {
            BSpline1D::new(&x, &y, degree, BoundaryCondition::NotAKnot).expect("valid")
        })
    };
    let nearest = |kind| scalar(case, || Nearest1D::new(&x, &y, kind).expect("valid"));
    let akima = |method| scalar(case, || Akima1D::new(&x, &y, method).expect("valid"));
    match case.text("method") {
        "linear" => scalar(case, || Linear1D::new(&x, &y).expect("valid")),
        "nearest" => nearest(NearestKind::Nearest),
        "previous" => nearest(NearestKind::Previous),
        "bspline1" => spline(1),
        "bspline2" => spline(2),
        "bspline3" => spline(3),
        "cubic_spline" => {
            scalar(case, || {
                CubicSpline::new(&x, &y, BoundaryCondition::NotAKnot).expect("valid")
            });
        }
        "pchip" => scalar(case, || PchipInterpolator::new(&x, &y).expect("valid")),
        "akima" => akima(AkimaMethod::Akima),
        "makima" => akima(AkimaMethod::Makima),
        other => unsupported(&format!("no method {other}")),
    }
}

fn curve(case: &Case) {
    let (t, q) = (case.array("t"), case.array("q"));
    let points = case.rows::<2>("points");
    let degree = case.count("degree");
    let build = || {
        let curve = BSplineCurve::with_parameters(&t, &points, degree, BoundaryCondition::NotAKnot);
        // The values of a curve are points: a place for them, kept with the
        // curve so that the timed evaluation does not allocate it.
        (curve.expect("valid"), RefCell::new(vec![[0.0; 2]; q.len()]))
    };
    run(
        case,
        library(),
        2 * q.len(),
        build,
        |(curve, space), out| {
            let mut space = space.borrow_mut();
            batch(curve, &q, &mut space, [f64::NAN; 2]);
            for (pair, point) in out.as_chunks_mut::<2>().0.iter_mut().zip(space.iter()) {
                *pair = *point;
            }
        },
    );
}

fn grid<const D: usize>(case: &Case) {
    let axes: [Vec<f64>; D] = core::array::from_fn(|d| case.array(&format!("axis{d}")));
    let values = case.array("values");
    let q = case.rows::<D>("q");
    let method = match case.text("method") {
        "linear" => GridMethod::Linear,
        "nearest" => GridMethod::Nearest,
        "cubic" => GridMethod::Cubic,
        "quintic" => GridMethod::Quintic,
        "pchip" => GridMethod::Pchip,
        other => unsupported(&format!("no method {other}")),
    };
    let eval = |built: &RegularGridInterpolator<f64, D>, out: &mut [f64]| {
        batch(built, &q, out, f64::NAN);
    };
    if owned() {
        let make = || (axes.clone(), values.clone());
        let build = |(axes, values)| {
            RegularGridInterpolator::from_owned(axes, values, method).expect("valid")
        };
        return run_given(case, OWNED, q.len(), make, build, eval);
    }
    let build = || {
        let slices: [&[f64]; D] = core::array::from_fn(|d| axes[d].as_slice());
        RegularGridInterpolator::new(slices, &values, method).expect("valid")
    };
    run(case, library(), q.len(), build, eval);
}

fn nearest_nd(case: &Case) {
    let (points, values, q) = (
        case.rows::<2>("points"),
        case.array("values"),
        case.rows::<2>("q"),
    );
    let eval = |built: &NearestNDInterpolator<2>, out: &mut [f64]| {
        batch(built, &q, out, f64::NAN);
    };
    if owned() {
        let build = |values| NearestNDInterpolator::from_owned(&points, values).expect("valid");
        return run_given(case, OWNED, q.len(), || values.clone(), build, eval);
    }
    let build = || NearestNDInterpolator::new(&points, &values).expect("valid");
    run(case, library(), q.len(), build, eval);
}

fn rbf(case: &Case) {
    let (points, values, q) = (
        case.rows::<2>("points"),
        case.array("values"),
        case.rows::<2>("q"),
    );
    let epsilon = case.number("epsilon");
    let kernel = match case.text("kernel") {
        "linear" => RbfKernel::Linear,
        "thin_plate_spline" => RbfKernel::ThinPlateSpline,
        "cubic" => RbfKernel::Cubic,
        "quintic" => RbfKernel::Quintic,
        "multiquadric" => RbfKernel::Multiquadric { epsilon },
        "inverse_multiquadric" => RbfKernel::InverseMultiquadric { epsilon },
        "inverse_quadratic" => RbfKernel::InverseQuadratic { epsilon },
        "gaussian" => RbfKernel::Gaussian { epsilon },
        other => unsupported(&format!("no kernel {other}")),
    };
    let options = match case.count_or_none("neighbors") {
        Some(k) => RbfOptions::new().with_neighbors(k),
        None => RbfOptions::new(),
    };
    let eval = |built: &RBFInterpolator<2>, out: &mut [f64]| {
        batch(built, &q, out, f64::NAN);
    };
    if owned() {
        let make = || (points.clone(), values.clone());
        let build = |(points, values)| {
            RBFInterpolator::from_owned(points, values, kernel, &options).expect("valid")
        };
        return run_given(case, OWNED, q.len(), make, build, eval);
    }
    let build =
        || RBFInterpolator::with_options(&points, &values, kernel, &options).expect("valid");
    run(case, library(), q.len(), build, eval);
}

fn linear_nd<const D: usize>(case: &Case) {
    let (points, values, q) = (
        case.rows::<D>("points"),
        case.array("values"),
        case.rows::<D>("q"),
    );
    let eval = |built: &LinearNDInterpolator<D>, out: &mut [f64]| {
        batch(built, &q, out, f64::NAN);
    };
    if owned() {
        let make = || (points.clone(), values.clone());
        let build =
            |(points, values)| LinearNDInterpolator::from_owned(points, values).expect("valid");
        return run_given(case, OWNED, q.len(), make, build, eval);
    }
    let build = || LinearNDInterpolator::new(&points, &values).expect("valid");
    run(case, library(), q.len(), build, eval);
}

fn clough_tocher(case: &Case) {
    let (points, values, q) = (
        case.rows::<2>("points"),
        case.array("values"),
        case.rows::<2>("q"),
    );
    let eval = |built: &CloughTocher2DInterpolator, out: &mut [f64]| {
        batch(built, &q, out, f64::NAN);
    };
    if owned() {
        let make = || (points.clone(), values.clone());
        let build = |(points, values)| {
            CloughTocher2DInterpolator::from_owned(points, values).expect("valid")
        };
        return run_given(case, OWNED, q.len(), make, build, eval);
    }
    let build = || CloughTocher2DInterpolator::new(&points, &values).expect("valid");
    run(case, library(), q.len(), build, eval);
}

fn main() {
    let case = Case::from_args();
    if owned()
        && !matches!(
            case.family(),
            "grid" | "nearest_nd" | "linear_nd" | "clough_tocher" | "rbf"
        )
    {
        unsupported("built with `from_owned` for the types that keep their data");
    }
    match case.family() {
        "oned" => oned(&case),
        "ppoly" => {
            let (breakpoints, coefficients) =
                (case.array("breakpoints"), case.array("coefficients"));
            let degree = case.count("degree");
            scalar(&case, || {
                PPoly::new(&breakpoints, &coefficients, degree).expect("valid")
            });
        }
        "bspline" => {
            let (knots, control) = (case.array("knots"), case.array("control"));
            let degree = case.count("degree");
            scalar(&case, || {
                BSpline::new(&knots, &control, degree).expect("valid")
            });
        }
        "curve" => curve(&case),
        "grid" => match case.count("dim") {
            2 => grid::<2>(&case),
            3 => grid::<3>(&case),
            4 => grid::<4>(&case),
            other => unsupported(&format!("no grid of {other} dimensions")),
        },
        "nearest_nd" => nearest_nd(&case),
        "linear_nd" => match case.count_or_none("dim").unwrap_or(2) {
            2 => linear_nd::<2>(&case),
            3 => linear_nd::<3>(&case),
            other => unsupported(&format!("no points of {other} dimensions")),
        },
        "clough_tocher" => clough_tocher(&case),
        "rbf" => rbf(&case),
        other => unsupported(&format!("no family {other}")),
    }
}
