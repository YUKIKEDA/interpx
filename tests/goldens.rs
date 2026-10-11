//! Checks the golden files themselves: every file is well formed, the data sets
//! have the documented size, and the Rust copies of the analytic test functions
//! agree with the Python side.
//!
//! No interpolator exists yet. The rows that add one read the same files through
//! `common::load` and compare with the interpolator under test.

mod common;

use common::functions::{franke, genz, genz_a, rosenbrock};
use common::{assert_close, close, golden_dir, load, reals, GoldenFile, Real, DIRECT_F64};
use serde::Deserialize;
use serde_json::Value;
use std::collections::BTreeSet;

/// The files `just gen-goldens` writes. A generator that stops writing its file,
/// or a file nobody generates, fails `every_expected_file_exists_and_no_other`.
const FILES: [&str; 9] = [
    "calculus",
    "curve",
    "datasets",
    "functions",
    "grid",
    "oned",
    "rbf",
    "scattered",
    "tridiagonal",
];

/// docs/design.ja.md section 8.4 and `compare/common/golden.py`.
const MAX_FILE_BYTES: u64 = 1_000_000;
const MAX_TOTAL_BYTES: u64 = 4_000_000;

fn all() -> Result<Vec<(&'static str, GoldenFile<Value>)>, String> {
    FILES
        .iter()
        .map(|name| Ok((*name, load::<Value>(name)?)))
        .collect()
}

/// Every string found inside a numeric array, at any depth.
fn strings_in_arrays<'a>(v: &'a Value, inside_array: bool, out: &mut Vec<&'a str>) {
    match v {
        Value::Array(items) => items.iter().for_each(|w| strings_in_arrays(w, true, out)),
        Value::String(s) if inside_array => out.push(s),
        Value::Object(map) => map.values().for_each(|w| strings_in_arrays(w, false, out)),
        _ => {}
    }
}

#[test]
fn every_expected_file_exists_and_no_other() {
    let found: BTreeSet<String> = std::fs::read_dir(golden_dir())
        .expect("goldens directory")
        .map(|e| e.expect("directory entry").path())
        .filter(|p| p.extension().is_some_and(|x| x == "json"))
        .map(|p| {
            p.file_stem()
                .expect("file stem")
                .to_string_lossy()
                .into_owned()
        })
        .collect();
    let want: BTreeSet<String> = FILES.iter().map(|s| (*s).to_owned()).collect();
    assert_eq!(found, want);
}

#[test]
fn files_fit_the_size_budget() {
    let mut total = 0;
    for name in FILES {
        let len = std::fs::metadata(golden_dir().join(format!("{name}.json")))
            .expect("metadata")
            .len();
        assert!(len <= MAX_FILE_BYTES, "{name}.json is {len} bytes");
        total += len;
    }
    assert!(total <= MAX_TOTAL_BYTES, "goldens total {total} bytes");
}

#[test]
fn every_file_records_the_baseline_and_its_generator() {
    for (name, file) in all().expect("golden files") {
        assert!(
            file.meta.scipy.starts_with("1.18."),
            "{name}: SciPy baseline is 1.18, file says {}",
            file.meta.scipy
        );
        assert!(!file.meta.numpy.is_empty(), "{name}");
        let script = golden_dir()
            .parent()
            .expect("compare directory")
            .join(&file.meta.generator);
        assert!(
            script.is_file(),
            "{name}: no generator {}",
            script.display()
        );
    }
}

#[test]
fn every_case_is_well_formed() {
    for (name, file) in all().expect("golden files") {
        assert!(!file.cases.is_empty(), "{name}");
        let mut seen = BTreeSet::new();
        for case in &file.cases {
            let case_name = case["name"].as_str().expect("case name");
            let what = format!("{name}: {case_name}");
            assert!(seen.insert(case_name), "{what}: duplicate name");
            let reference = case["reference"].as_str().expect("reference");
            assert!(case["options"].is_object(), "{what}: options");

            // A data case holds inputs only. Any other case holds exactly one of
            // the result and the error class.
            let has_expected = case.get("expected").is_some();
            let has_error = case.get("expected_error").is_some();
            if reference == "data" {
                assert!(!has_expected && !has_error, "{what}");
            } else if reference.starts_with("scipy.interpolate.PPoly.from_spline") {
                assert!(case["coefficients"].is_array(), "{what}");
            } else if case_name.starts_with("calculus/") {
                // The results are per order: each entry has the values at the
                // queries.
                let derivatives = case["derivatives"].as_array().expect("derivatives");
                assert!(!derivatives.is_empty() && !has_expected, "{what}");
                let queries = case["queries"].as_array().expect("queries").len();
                for entry in derivatives {
                    assert!(entry["nu"].is_u64(), "{what}: nu");
                    let values = entry["values"].as_array().expect("values");
                    assert_eq!(values.len(), queries, "{what}: queries and values");
                }
                assert!(case["integrals"].is_array(), "{what}");
            } else {
                assert!(has_expected != has_error, "{what}: expected xor error");
            }

            if let (Some(q), Some(e)) = (case.get("queries"), case.get("expected")) {
                let (q, e) = (
                    q.as_array().expect("queries"),
                    e.as_array().expect("expected"),
                );
                assert_eq!(q.len(), e.len(), "{what}: queries and expected");
                assert!(!q.is_empty(), "{what}: no queries");
            }

            // Inside the numeric fields the only strings are the three spellings
            // of non-finite numbers. `options` is skipped: it holds names.
            let mut strings = Vec::new();
            for (key, value) in case.as_object().expect("case object") {
                if key != "options" {
                    strings_in_arrays(value, false, &mut strings);
                }
            }
            for s in strings {
                assert!(
                    matches!(s, "NaN" | "Infinity" | "-Infinity"),
                    "{what}: string {s:?} in a numeric field"
                );
            }
        }
    }
}

#[derive(Debug, Deserialize)]
struct Data {
    name: String,
    options: Value,
    x: Option<Vec<Real>>,
    y: Option<Vec<Real>>,
    points: Option<Vec<Vec<Real>>>,
    values: Option<Vec<Real>>,
    queries: Option<Vec<Vec<Real>>>,
    queries_values: Option<Vec<Real>>,
}

fn find_data(name: &str) -> Result<Data, String> {
    load::<Data>("datasets")?
        .cases
        .into_iter()
        .find(|c| c.name == name)
        .ok_or_else(|| format!("no data set {name}"))
}

#[test]
fn classic_data_sets_have_the_documented_size() {
    for (name, n) in [
        ("fritsch_carlson", 11),
        ("titanium", 49),
        ("square_with_near_point", 12),
    ] {
        let d = find_data(name).expect("data set");
        let (x, y) = (reals(&d.x.expect("x")), reals(&d.y.expect("y")));
        assert_eq!((x.len(), y.len()), (n, n), "{name}");
        assert!(
            x.windows(2).all(|w| w[0] < w[1]),
            "{name}: x not increasing"
        );
        assert!(x.iter().chain(&y).all(|v| v.is_finite()), "{name}");
        assert_eq!(d.options["kind"], "classic_1d", "{name}");
    }
    for (name, n) in [("franke_100", 100), ("franke_33", 33), ("lawson_25", 25)] {
        let d = find_data(name).expect("data set");
        let points = d.points.expect("points");
        assert_eq!(points.len(), n, "{name}");
        assert!(points.iter().all(|p| p.len() == 2), "{name}");
        let unique: BTreeSet<(u64, u64)> = points
            .iter()
            .map(|p| (p[0].0.to_bits(), p[1].0.to_bits()))
            .collect();
        assert_eq!(unique.len(), n, "{name}: duplicate nodes");
    }
}

#[test]
fn franke_100_reaches_outside_the_unit_square() {
    // The design notes this: queries on the unit square can fall outside the hull.
    let points = find_data("franke_100")
        .expect("data set")
        .points
        .expect("points");
    assert!(points.iter().any(|p| p[1].0 < 0.0));
}

#[test]
fn adversarial_sets_are_what_their_names_say() {
    let circle = find_data("cocircular_2d")
        .expect("data set")
        .points
        .expect("points");
    assert_eq!(circle.len(), 12);
    assert!(circle
        .iter()
        .all(|p| (p[0].0 * p[0].0 + p[1].0 * p[1].0 - 25.0).abs() < 1e-12));

    let sphere = find_data("cospherical_3d")
        .expect("data set")
        .points
        .expect("points");
    assert_eq!(sphere.len(), 30);
    assert!(sphere
        .iter()
        .all(|p| (p.iter().map(|c| c.0 * c.0).sum::<f64>() - 9.0).abs() < 1e-12));

    let line = find_data("collinear_2d")
        .expect("data set")
        .points
        .expect("points");
    assert!(line
        .iter()
        .all(|p| (p[1].0 - (2.0 * p[0].0 + 1.0)).abs() < 1e-12));

    let plane = find_data("coplanar_3d")
        .expect("data set")
        .points
        .expect("points");
    assert!(plane
        .iter()
        .all(|p| (p[2].0 - (p[0].0 + p[1].0)).abs() < 1e-12));

    let dup = find_data("duplicate_2d")
        .expect("data set")
        .points
        .expect("points");
    let bits: Vec<(u64, u64)> = dup
        .iter()
        .map(|p| (p[0].0.to_bits(), p[1].0.to_bits()))
        .collect();
    let unique: BTreeSet<_> = bits.iter().collect();
    assert_eq!(unique.len() + 1, bits.len(), "exactly one repeated point");

    let near = find_data("near_duplicate_2d")
        .expect("data set")
        .points
        .expect("points");
    let closest = near
        .iter()
        .enumerate()
        .flat_map(|(i, p)| near[i + 1..].iter().map(move |q| (p, q)))
        .map(|(p, q)| (p[0].0 - q[0].0).hypot(p[1].0 - q[1].0))
        .fold(f64::INFINITY, f64::min);
    assert!(closest > 0.0 && closest < 1e-11, "closest pair {closest:e}");
}

#[test]
fn adversarial_values_are_the_stated_affine_function() {
    // `LinearND` must reproduce an affine function exactly; the tests of that
    // rely on these values being affine in the points.
    let coef = [2.0, -3.0, 0.5];
    for name in [
        "cocircular_2d",
        "cospherical_3d",
        "collinear_2d",
        "coplanar_3d",
        "sliver_2d",
        "duplicate_2d",
        "near_duplicate_2d",
        "hull_exterior_2d",
    ] {
        let d = find_data(name).expect("data set");
        assert_eq!(d.options["kind"], "adversarial", "{name}");
        let affine = |p: &Vec<Real>| 1.0 + p.iter().zip(coef).map(|(c, k)| c.0 * k).sum::<f64>();
        let (points, values) = (d.points.expect("points"), d.values.expect("values"));
        assert_eq!(points.len(), values.len(), "{name}");
        for (p, v) in points.iter().zip(&values) {
            assert_close(affine(p), v.0, DIRECT_F64, name);
        }
        let (q, qv) = (
            d.queries.expect("queries"),
            d.queries_values.expect("values"),
        );
        assert_eq!(q.len(), qv.len(), "{name}");
        for (p, v) in q.iter().zip(&qv) {
            assert_close(affine(p), v.0, DIRECT_F64, name);
        }
    }
}

#[derive(Debug, Deserialize)]
struct FunctionCase {
    name: String,
    reference: String,
    options: Value,
    queries: Vec<Vec<Real>>,
    expected: Vec<Real>,
}

#[test]
fn rust_test_functions_match_the_python_side() {
    let file = load::<FunctionCase>("functions").expect("golden file functions");
    let mut kinds = BTreeSet::new();
    for case in &file.cases {
        kinds.insert(case.reference.as_str());
        for (k, (q, want)) in case.queries.iter().zip(&case.expected).enumerate() {
            let x = reals(q);
            let got = match case.reference.as_str() {
                "datasets.franke" => {
                    let which = case.options["k"].as_u64().expect("k");
                    franke(which, x[0], x[1]).expect("Franke function")
                }
                "datasets.genz" => {
                    let kind = case.options["kind"].as_str().expect("kind");
                    genz(kind, &x).expect("Genz family")
                }
                "datasets.rosenbrock" => rosenbrock(&x),
                other => panic!("{}: unknown reference {other}", case.name),
            };
            assert_close(got, want.0, DIRECT_F64, &format!("{} query {k}", case.name));
        }
        if case.reference == "datasets.genz" {
            let d = usize::try_from(case.options["d"].as_u64().expect("d")).expect("usize");
            let a: Vec<Real> = serde_json::from_value(case.options["a"].clone()).expect("a");
            for (got, want) in genz_a(d).iter().zip(reals(&a)) {
                assert_close(*got, want, DIRECT_F64, &case.name);
            }
        }
    }
    assert_eq!(
        kinds.len(),
        3,
        "Franke, Genz, and Rosenbrock are all present"
    );
    assert!(franke(7, 0.5, 0.5).is_none());
    assert!(genz("no_such_family", &[0.5, 0.5, 0.5]).is_none());
}

#[test]
fn franke_f6_is_the_sphere_of_radius_eight_ninths() {
    // The form with `/ 9` is the one that is a sphere centred at (1/2, 1/2, -1/2)
    // (docs/design.ja.md section 8.4.1). Without it the surface is not a sphere.
    for (x, y) in [(0.5, 0.5), (0.2, 0.7), (0.0, 0.0), (1.0, 0.3)] {
        let z = franke(6, x, y).expect("F6");
        let r2: f64 = (x - 0.5f64).powi(2) + (y - 0.5f64).powi(2) + (z + 0.5).powi(2);
        assert_close(r2, 64.0 / 81.0, DIRECT_F64, "F6 radius squared");
    }
    // F5 has the factor 1/3: its maximum is 1/3 at the centre.
    assert_close(
        franke(5, 0.5, 0.5).expect("F5"),
        1.0 / 3.0,
        DIRECT_F64,
        "F5 peak",
    );
    assert_close(
        franke(4, 0.5, 0.5).expect("F4"),
        1.0 / 3.0,
        DIRECT_F64,
        "F4 peak",
    );
}

/// Reads a JSON array of numbers (or the three non-finite spellings).
fn nums(v: &Value) -> Vec<f64> {
    reals(&serde_json::from_value::<Vec<Real>>(v.clone()).unwrap_or_default())
}

/// Reads a JSON matrix: one row per point.
fn rows(v: &Value) -> Vec<Vec<f64>> {
    serde_json::from_value::<Vec<Vec<Real>>>(v.clone())
        .unwrap_or_default()
        .iter()
        .map(|r| reals(r))
        .collect()
}

fn same_point(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(p, q)| p.to_bits() == q.to_bits())
}

/// The tolerance of the interpolation property at a data point. It is looser
/// than the comparison with SciPy: the reference itself solves a linear system.
const AT_DATA: common::Tol = common::Tol {
    rtol: 1e-9,
    atol: 1e-9,
};

#[test]
fn one_dimensional_goldens_reproduce_the_data_at_the_knots() {
    let file = load::<Value>("oned").expect("golden file oned");
    let mut checked = 0;
    for case in &file.cases {
        let name = case["name"].as_str().expect("name");
        let (Some(x), Some(y), Some(q), Some(e)) = (
            case.get("x"),
            case.get("y"),
            case.get("queries"),
            case.get("expected"),
        ) else {
            continue;
        };
        let (x, y, q, e) = (nums(x), nums(y), nums(q), nums(e));
        for (xi, yi) in x.iter().zip(&y) {
            let k = q
                .iter()
                .position(|v| v.to_bits() == xi.to_bits())
                .unwrap_or_else(|| panic!("{name}: knot {xi} is not among the queries"));
            assert_close(e[k], *yi, AT_DATA, &format!("{name} at knot {xi}"));
            checked += 1;
        }
    }
    assert!(checked > 1000, "only {checked} knots checked");
}

#[test]
fn curve_goldens_reproduce_the_points_at_the_parameters() {
    let file = load::<Value>("curve").expect("golden file curve");
    for case in &file.cases {
        let name = case["name"].as_str().expect("name");
        let (t, points) = (nums(&case["t"]), rows(&case["points"]));
        let (q, e) = (nums(&case["queries"]), rows(&case["expected"]));
        for (ti, p) in t.iter().zip(&points) {
            let k = q
                .iter()
                .position(|v| v.to_bits() == ti.to_bits())
                .unwrap_or_else(|| panic!("{name}: parameter {ti} is not among the queries"));
            for (got, want) in e[k].iter().zip(p) {
                assert_close(*got, *want, AT_DATA, &format!("{name} at t = {ti}"));
            }
        }
    }
}

#[test]
fn grid_goldens_reproduce_the_values_at_the_grid_nodes() {
    // This is the check that found SciPy's default solver for `cubic` and
    // `quintic` to be about 1e-5 away from the spline: `expected` must be the
    // interpolant, so it must return the data at a node for every method.
    let file = load::<Value>("grid").expect("golden file grid");
    let mut checked = 0;
    for case in &file.cases {
        let name = case["name"].as_str().expect("name");
        let axes = rows(&case["axes"]);
        let values = nums(&case["values"]);
        let (q, e) = (rows(&case["queries"]), nums(&case["expected"]));
        assert_eq!(
            values.len(),
            axes.iter().map(Vec::len).product::<usize>(),
            "{name}: values and shape"
        );
        for (point, got) in q.iter().zip(&e) {
            // The flat index in C order, when every coordinate is a node.
            let index = point.iter().zip(&axes).try_fold(0usize, |acc, (c, axis)| {
                let i = axis.iter().position(|a| a.to_bits() == c.to_bits())?;
                Some(acc * axis.len() + i)
            });
            if let Some(index) = index {
                assert_close(*got, values[index], AT_DATA, &format!("{name} at a node"));
                checked += 1;
            }
        }
    }
    assert!(checked > 100, "only {checked} nodes checked");
}

#[test]
fn grid_goldens_record_how_far_scipy_default_solver_is() {
    // docs/design.ja.md section 6.4: interpx solves the tensor-product spline
    // exactly, as SciPy does with `solver=spsolve` or the `*_legacy` methods.
    // SciPy's default (an iterative solver) differs from it: about 1e-5 on the
    // uniform 2D grids and up to about 1e-2 on the non-uniform and 3D grids
    // here. If the difference vanished, the note in the design would be stale.
    let file = load::<Value>("grid").expect("golden file grid");
    let mut worst: f64 = 0.0;
    for case in &file.cases {
        let name = case["name"].as_str().expect("name");
        let method = case["options"]["method"].as_str().expect("method");
        let has = case.get("expected_scipy_default").is_some();
        assert_eq!(has, matches!(method, "cubic" | "quintic"), "{name}");
        if has {
            let (a, b) = (
                nums(&case["expected"]),
                nums(&case["expected_scipy_default"]),
            );
            assert_eq!(a.len(), b.len(), "{name}");
            for (p, q) in a.iter().zip(&b) {
                assert_eq!(p.is_nan(), q.is_nan(), "{name}: NaN pattern");
                if !p.is_nan() {
                    worst = worst.max((p - q).abs() / p.abs().max(1.0));
                }
            }
        }
    }
    assert!(
        worst > 1e-3 && worst < 1e-1,
        "largest relative difference {worst:e}"
    );
}

#[test]
fn scattered_goldens_reproduce_the_values_at_the_data_points() {
    for (file_name, tol) in [("scattered", AT_DATA), ("rbf", AT_DATA)] {
        let file = load::<Value>(file_name).expect("golden file");
        let mut checked = 0;
        for case in &file.cases {
            let name = case["name"].as_str().expect("name");
            // A smoothing spline does not pass through the data.
            if case["options"]["smoothing"]
                .as_f64()
                .is_some_and(|s| s > 0.0)
            {
                continue;
            }
            let (points, values) = (rows(&case["points"]), nums(&case["values"]));
            let (q, e) = (rows(&case["queries"]), nums(&case["expected"]));
            for (query, got) in q.iter().zip(&e) {
                if let Some(i) = points.iter().position(|p| same_point(p, query)) {
                    assert_close(*got, values[i], tol, &format!("{file_name}: {name}"));
                    checked += 1;
                }
            }
        }
        assert!(
            checked > 50,
            "{file_name}: only {checked} data points checked"
        );
    }
}

#[test]
fn tolerance_rejects_a_tampered_value() {
    let tol = DIRECT_F64;
    assert!(close(1.0, 1.0, tol));
    assert!(close(1.0, 1.0 + 1e-13, tol));
    assert!(!close(1.0, 1.0 + 1e-9, tol));
}

#[test]
fn tolerance_handles_non_finite_values() {
    let tol = DIRECT_F64;
    assert!(close(f64::NAN, f64::NAN, tol));
    assert!(!close(f64::NAN, 0.0, tol));
    assert!(!close(0.0, f64::NAN, tol));
    assert!(close(f64::INFINITY, f64::INFINITY, tol));
    assert!(!close(f64::INFINITY, f64::NEG_INFINITY, tol));
    assert!(!close(f64::INFINITY, 1.0e308, tol));
}
