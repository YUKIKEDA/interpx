//! The `scirs2-interpolate` side of the comparison: the methods it shares
//! with SciPy's interface, on one variable, on grids, and the global RBF.
//!
//! Its evaluations return a new array; the values are copied into the place
//! the harness gives, inside the timed call (that is the interface there is).
//! Its grid interpolator takes the axes and the values by value, so
//! construction copies them, as interpx's does.

use perf_harness::{run, unsupported, Case};
use scirs2_core::ndarray::{Array1, ArrayD, ArrayView1, ArrayView2, IxDyn};
use scirs2_interpolate::interp1d::{ExtrapolateMode, Interp1d, InterpolationMethod};
use scirs2_interpolate::interpnd::{
    ExtrapolateMode as GridExtrapolate, InterpolationMethod as GridMethod, RegularGridInterpolator,
};
use scirs2_interpolate::{AkimaSpline, CubicSpline, PchipInterpolator, RBFInterpolator, RBFKernel};

const LIBRARY: &str = "scirs2";

fn copy(values: Array1<f64>, out: &mut [f64]) {
    for (slot, value) in out.iter_mut().zip(values.iter()) {
        *slot = *value;
    }
}

fn oned(case: &Case) {
    let (x, y, q) = (case.array("x"), case.array("y"), case.array("q"));
    let (x, y, q) = (ArrayView1::from(&x), ArrayView1::from(&y), ArrayView1::from(&q));
    let basic = |method| {
        run(
            case,
            LIBRARY,
            q.len(),
            || Interp1d::new(&x, &y, method, ExtrapolateMode::Error).expect("valid"),
            |built, out| copy(built.evaluate_array(&q).expect("in range"), out),
        );
    };
    match case.text("method") {
        "linear" => basic(InterpolationMethod::Linear),
        "nearest" => basic(InterpolationMethod::Nearest),
        "cubic_spline" => run(
            case,
            LIBRARY,
            q.len(),
            || CubicSpline::new_not_a_knot(&x, &y).expect("valid"),
            |built, out| copy(built.evaluate_array(&q).expect("in range"), out),
        ),
        "pchip" => run(
            case,
            LIBRARY,
            q.len(),
            || PchipInterpolator::new(&x, &y, false).expect("valid"),
            |built, out| copy(built.evaluate_array(&q).expect("in range"), out),
        ),
        "akima" => run(
            case,
            LIBRARY,
            q.len(),
            || AkimaSpline::new(&x, &y).expect("valid"),
            |built, out| copy(built.evaluate_array(&q).expect("in range"), out),
        ),
        other => unsupported(&format!("no method {other}")),
    }
}

fn grid(case: &Case) {
    let dim = case.count("dim");
    let axes: Vec<Array1<f64>> = (0..dim).map(|d| Array1::from(case.array(&format!("axis{d}")))).collect();
    if axes.iter().any(|axis| axis[1] < axis[0]) {
        unsupported("axes must increase");
    }
    let method = match case.text("method") {
        "linear" => GridMethod::Linear,
        "nearest" => GridMethod::Nearest,
        other => unsupported(&format!("no method {other}")),
    };
    let shape: Vec<usize> = axes.iter().map(Array1::len).collect();
    let values = ArrayD::from_shape_vec(IxDyn(&shape), case.array("values")).expect("the shape of the grid");
    let flat = case.array("q");
    let q = ArrayView2::from_shape((flat.len() / dim, dim), &flat).expect("rows of points");
    run(
        case,
        LIBRARY,
        flat.len() / dim,
        || RegularGridInterpolator::new(axes.clone(), values.clone(), method, GridExtrapolate::Error).expect("valid"),
        |built, out| copy(built.__call__(&q).expect("in range"), out),
    );
}

fn rbf(case: &Case) {
    if case.count_or_none("neighbors").is_some() {
        unsupported("no neighbors mode");
    }
    let kernel = match case.text("kernel") {
        "gaussian" => RBFKernel::Gaussian,
        "multiquadric" => RBFKernel::Multiquadric,
        "inverse_multiquadric" => RBFKernel::InverseMultiquadric,
        "thin_plate_spline" => RBFKernel::ThinPlateSpline,
        "linear" => RBFKernel::Linear,
        "cubic" => RBFKernel::Cubic,
        "quintic" => RBFKernel::Quintic,
        other => unsupported(&format!("no kernel {other}")),
    };
    let epsilon = case.number("epsilon");
    let (points, values, q) = (case.array("points"), case.array("values"), case.array("q"));
    let points = ArrayView2::from_shape((points.len() / 2, 2), &points).expect("rows of points");
    let q = ArrayView2::from_shape((q.len() / 2, 2), &q).expect("rows of points");
    let values = ArrayView1::from(&values);
    run(
        case,
        LIBRARY,
        q.nrows(),
        || RBFInterpolator::new(&points, &values, kernel, epsilon).expect("valid"),
        |built, out| copy(built.interpolate(&q).expect("valid"), out),
    );
}

fn main() {
    let case = Case::from_args();
    match case.family() {
        "oned" => oned(&case),
        "grid" => grid(&case),
        "rbf" => rbf(&case),
        other => unsupported(&format!("no family {other}")),
    }
}
