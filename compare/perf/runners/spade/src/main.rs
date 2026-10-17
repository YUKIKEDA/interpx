//! The `spade` side of the comparison: linear interpolation on the Delaunay
//! triangulation of points of the plane (`LinearNDInterpolator` in two
//! dimensions).
//!
//! The triangulation is bulk-loaded from the points, each with its value. The
//! evaluation is `spade`'s barycentric interpolation, which has no value
//! outside the convex hull, as SciPy and interpx have none. Where four or more
//! points lie on one circle, `spade` splits the cell by its own rule, and the
//! values inside it may differ from SciPy's.

use perf_harness::{refused, run, unsupported, Case};
use spade::{DelaunayTriangulation, FloatTriangulation, HasPosition, Point2, Triangulation};

const LIBRARY: &str = "spade";

/// A point of the data with its value.
struct Site {
    position: Point2<f64>,
    value: f64,
}

impl HasPosition for Site {
    type Scalar = f64;

    fn position(&self) -> Point2<f64> {
        self.position
    }
}

fn linear_nd(case: &Case) {
    if case.count_or_none("dim").unwrap_or(2) != 2 {
        unsupported("points of the plane only");
    }
    let (points, values, q) = (
        case.rows::<2>("points"),
        case.array("values"),
        case.rows::<2>("q"),
    );
    run(
        case,
        LIBRARY,
        q.len(),
        || {
            let sites: Vec<Site> = points
                .iter()
                .zip(&values)
                .map(|(p, v)| Site {
                    position: Point2::new(p[0], p[1]),
                    value: *v,
                })
                .collect();
            DelaunayTriangulation::<Site>::bulk_load(sites).unwrap_or_else(|error| refused(error))
        },
        |built, out| {
            let barycentric = built.barycentric();
            for (query, slot) in q.iter().zip(out.iter_mut()) {
                let at = Point2::new(query[0], query[1]);
                *slot = barycentric
                    .interpolate(|vertex| vertex.data().value, at)
                    .unwrap_or(f64::NAN);
            }
        },
    );
}

fn main() {
    let case = Case::from_args();
    match case.family() {
        "linear_nd" => linear_nd(&case),
        other => unsupported(&format!("no family {other}")),
    }
}
