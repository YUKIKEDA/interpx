//! The `interpn` side of the comparison: grids, with the methods it has.
//!
//! `interpn` borrows the axes and the values, so its construction is a check
//! of the sizes. Its queries are one array for each coordinate; they are laid
//! out that way before anything is timed. Its cubic method is a local cubic
//! through the four nearest grid points along each axis, not the spline SciPy
//! and interpx build, so the values of that method differ.

use interpn::{MulticubicRectilinear, MultilinearRectilinear, NearestRectilinear};
use perf_harness::{run, unsupported, Case};

const LIBRARY: &str = "interpn";

fn grid<const D: usize>(case: &Case) {
    let axes: [Vec<f64>; D] = core::array::from_fn(|d| case.array(&format!("axis{d}")));
    if axes.iter().any(|axis| axis[1] < axis[0]) {
        unsupported("axes must increase");
    }
    let values = case.array("values");
    let rows = case.rows::<D>("q");
    let columns: [Vec<f64>; D] = core::array::from_fn(|d| rows.iter().map(|row| row[d]).collect());
    let grids: [&[f64]; D] = core::array::from_fn(|d| axes[d].as_slice());
    let q: [&[f64]; D] = core::array::from_fn(|d| columns[d].as_slice());
    let count = rows.len();
    match case.text("method") {
        "linear" => run(
            case,
            LIBRARY,
            count,
            || MultilinearRectilinear::new(&grids, &values).expect("valid"),
            |built, out| built.interp(&q, out).expect("same length"),
        ),
        "nearest" => run(
            case,
            LIBRARY,
            count,
            || NearestRectilinear::new(&grids, &values).expect("valid"),
            |built, out| built.interp(&q, out).expect("same length"),
        ),
        "cubic" => run(
            case,
            LIBRARY,
            count,
            || MulticubicRectilinear::new(&grids, &values, false).expect("valid"),
            |built, out| built.interp(&q, out).expect("same length"),
        ),
        other => unsupported(&format!("no method {other}")),
    }
}

fn main() {
    let case = Case::from_args();
    if case.family() != "grid" {
        unsupported("grids only");
    }
    match case.count("dim") {
        2 => grid::<2>(&case),
        3 => grid::<3>(&case),
        4 => grid::<4>(&case),
        other => unsupported(&format!("no grid of {other} dimensions")),
    }
}
