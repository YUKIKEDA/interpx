//! The `ndarray-interp` side of the comparison: linear and cubic-spline
//! interpolation of one variable, and bilinear interpolation on a grid of two.
//!
//! The interpolators are built on views of the data, the cheapest way the
//! crate offers.

use ndarray::{ArrayView1, ArrayView2, ArrayViewMut1};
use ndarray_interp::interp1d::cubic_spline::{BoundaryCondition, CubicSpline};
use ndarray_interp::interp1d::{Interp1DBuilder, Linear};
use ndarray_interp::interp2d::{Bilinear, Interp2DBuilder};
use perf_harness::{run, unsupported, Case};

const LIBRARY: &str = "ndarray_interp";

fn oned(case: &Case) {
    let (x, y, q) = (case.array("x"), case.array("y"), case.array("q"));
    let (x, y, q) = (ArrayView1::from(&x), ArrayView1::from(&y), ArrayView1::from(&q));
    match case.text("method") {
        "linear" => run(
            case,
            LIBRARY,
            q.len(),
            || Interp1DBuilder::new(y).x(x).strategy(Linear::new()).build().expect("valid"),
            |built, out| {
                built.interp_array_into(&q, ArrayViewMut1::from(out)).expect("in range");
            },
        ),
        "cubic_spline" => run(
            case,
            LIBRARY,
            q.len(),
            || {
                let strategy = CubicSpline::new().boundary(BoundaryCondition::NotAKnot);
                Interp1DBuilder::new(y).x(x).strategy(strategy).build().expect("valid")
            },
            |built, out| {
                built.interp_array_into(&q, ArrayViewMut1::from(out)).expect("in range");
            },
        ),
        other => unsupported(&format!("no method {other}")),
    }
}

fn grid(case: &Case) {
    if case.count("dim") != 2 || case.text("method") != "linear" {
        unsupported("bilinear on two axes only");
    }
    let (ax, ay, values) = (case.array("axis0"), case.array("axis1"), case.array("values"));
    if ax[1] < ax[0] {
        unsupported("axes must increase");
    }
    let rows = case.rows::<2>("q");
    let qx: Vec<f64> = rows.iter().map(|row| row[0]).collect();
    let qy: Vec<f64> = rows.iter().map(|row| row[1]).collect();
    let (qx, qy) = (ArrayView1::from(&qx), ArrayView1::from(&qy));
    let data = ArrayView2::from_shape((ax.len(), ay.len()), &values).expect("the shape of the grid");
    let (ax, ay) = (ArrayView1::from(&ax), ArrayView1::from(&ay));
    run(
        case,
        LIBRARY,
        rows.len(),
        || Interp2DBuilder::new(data).x(ax).y(ay).strategy(Bilinear::new()).build().expect("valid"),
        |built, out| {
            built.interp_array_into(&qx, &qy, ArrayViewMut1::from(out)).expect("in range");
        },
    );
}

fn main() {
    let case = Case::from_args();
    match case.family() {
        "oned" => oned(&case),
        "grid" => grid(&case),
        _ => unsupported("one variable, or a grid of two"),
    }
}
