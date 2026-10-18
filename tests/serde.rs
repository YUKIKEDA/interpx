//! The `serde` feature: every interpolator through JSON and through a binary
//! format and back, evaluating bit for bit as before; and stored forms that
//! are broken, which are errors and not interpolators (docs/design.ja.md
//! section 6.5).

#![cfg(feature = "serde")]

use interpx::{
    Akima1D, AkimaMethod, BSpline, BSpline1D, BSplineCurve, BoundaryCondition, CubicHermiteSpline,
    CubicSpline, EndCondition, Extrapolate, GridMethod, Interpolator, Linear1D, Nearest1D,
    NearestKind, PPoly, Parametrization, PchipInterpolator, RegularGridInterpolator,
};
use serde::de::DeserializeOwned;
use serde::Serialize;
use serde_json::{json, Value as Json};
use std::fmt::Debug;

/// A small generator, so that the inputs are the same on every run.
struct Lcg(u64);

impl Lcg {
    /// The next number in [0, 1).
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1_u64 << 53) as f64
    }
}

/// Coordinates that increase by steps of no fixed length.
fn uneven(rng: &mut Lcg, n: usize) -> Vec<f64> {
    (0..n)
        .scan(0.0, |at, _| {
            Some(std::mem::replace(at, *at + 0.5 + rng.next()))
        })
        .collect()
}

/// A change to the JSON of a stored form that leaves no interpolator.
#[derive(Debug, Clone, Copy)]
enum Fault {
    /// The first two elements of the array change places.
    Swapped,
    /// The second element of the array becomes a copy of the first.
    Repeated,
    /// The array is turned around.
    Reversed,
}

/// The floats in a stored form, the first and the last of every array of
/// them and every one outside an array, zero left out (its bytes are found
/// anywhere).
fn floats(form: &Json, found: &mut Vec<f64>) {
    let mut take = |value: &Json| {
        if let Some(number) = value.as_f64().filter(|n| value.is_f64() && *n != 0.0) {
            found.push(number);
        }
    };
    match form {
        Json::Array(items) if items.iter().all(Json::is_number) => {
            items.first().into_iter().chain(items.last()).for_each(take);
        }
        Json::Array(items) => items.iter().for_each(|item| floats(item, found)),
        Json::Object(fields) => fields.values().for_each(|field| floats(field, found)),
        other => take(other),
    }
}

/// Every place of a stored form: the JSON pointers to its arrays that are
/// not empty, and to its versions.
fn places(form: &Json, at: &str, arrays: &mut Vec<String>, versions: &mut Vec<String>) {
    match form {
        Json::Array(items) => {
            if !items.is_empty() {
                arrays.push(at.to_owned());
            }
            for (i, item) in items.iter().enumerate() {
                places(item, &format!("{at}/{i}"), arrays, versions);
            }
        }
        Json::Object(fields) => {
            for (name, field) in fields {
                let place = format!("{at}/{name}");
                if name == "version" {
                    versions.push(place.clone());
                }
                places(field, &place, arrays, versions);
            }
        }
        _ => {}
    }
}

/// Where `needle` first is in `bytes`.
fn position(bytes: &[u8], needle: &[u8]) -> Option<usize> {
    bytes
        .windows(needle.len())
        .position(|window| window == needle)
}

/// Holds an interpolator against what comes back of it from JSON and from
/// `postcard`: every evaluation bit for bit the same (compared as the text
/// that tells `-0.0` from `0.0`), and the form of what came back the form
/// that was stored. Then breaks the stored form in every way that applies to
/// all types: each array one element short, each version another number,
/// each float not a number. `faults` are the ways of this type: coordinates
/// out of order, points twice.
///
/// `single` says that the scalars are `f32` (their bytes in the binary form
/// are four).
// Marked as test code so that `expect` is allowed here as in a `#[test]`
// function (clippy.toml allows it in tests only).
#[cfg(test)]
fn check<I>(name: &str, interp: &I, points: &[I::Point], faults: &[(&str, Fault)], single: bool)
where
    I: Interpolator + Serialize + DeserializeOwned,
    I::Output: Debug + Clone,
{
    let text = serde_json::to_string(interp).expect("stored as JSON");
    let bytes = postcard::to_stdvec(interp).expect("stored as postcard");
    let from_text: I = serde_json::from_str(&text).expect("read from JSON");
    let from_bytes: I = postcard::from_bytes(&bytes).expect("read from postcard");
    let fill = interp
        .eval_batch(points)
        .into_iter()
        .flatten()
        .next()
        .expect("a point with a value");
    for (format, back) in [("JSON", &from_text), ("postcard", &from_bytes)] {
        for point in points {
            let (before, after) = (interp.eval(point), back.eval(point));
            assert_eq!(
                format!("{before:?}"),
                format!("{after:?}"),
                "{name}, {format}: eval"
            );
        }
        let (before, after) = (interp.eval_batch(points), back.eval_batch(points));
        assert_eq!(
            format!("{before:?}"),
            format!("{after:?}"),
            "{name}, {format}: eval_batch"
        );
        let mut before = vec![fill.clone(); points.len()];
        let mut after = before.clone();
        assert_eq!(
            interp.eval_batch_into(points, &mut before, fill.clone()),
            Ok(())
        );
        assert_eq!(
            back.eval_batch_into(points, &mut after, fill.clone()),
            Ok(())
        );
        assert_eq!(
            format!("{before:?}"),
            format!("{after:?}"),
            "{name}, {format}: into"
        );
        assert_eq!(
            serde_json::to_string(back).expect("stored again"),
            text,
            "{name}, {format}"
        );
    }

    let form: Json = serde_json::from_str(&text).expect("JSON");
    let broken = |form: Json| serde_json::from_value::<I>(form).is_err();
    assert!(!broken(form.clone()), "{name}: the form as it is");
    let (mut arrays, mut versions) = (Vec::new(), Vec::new());
    places(&form, "", &mut arrays, &mut versions);
    assert!(!versions.is_empty(), "{name}: a version");
    for place in &versions {
        let mut other = form.clone();
        *other.pointer_mut(place).expect("the version") = json!(2);
        assert!(broken(other), "{name}: {place} is 2");
    }
    for place in &arrays {
        let mut short = form.clone();
        let array = short.pointer_mut(place).and_then(Json::as_array_mut);
        array.expect("the array").pop();
        assert!(broken(short), "{name}: {place} one short");
    }
    for (place, fault) in faults {
        let mut other = form.clone();
        let array = other.pointer_mut(place).and_then(Json::as_array_mut);
        let array = array.expect("the array of the fault");
        match fault {
            Fault::Swapped => array.swap(0, 1),
            Fault::Repeated => array[1] = array[0].clone(),
            Fault::Reversed => array.reverse(),
        }
        assert!(broken(other), "{name}: {place} {fault:?}");
    }

    // JSON has no NaN, so a float that is not finite is put into the binary
    // form, in place of the bytes of a float the form is known to hold.
    let mut held = Vec::new();
    floats(&form, &mut held);
    assert!(!held.is_empty(), "{name}: floats");
    for float in held {
        let (was, nan) = if single {
            (
                (float as f32).to_le_bytes().to_vec(),
                f32::NAN.to_le_bytes().to_vec(),
            )
        } else {
            (
                float.to_le_bytes().to_vec(),
                f64::NAN.to_le_bytes().to_vec(),
            )
        };
        let at = position(&bytes, &was).expect("the bytes of the float");
        let mut with_nan = bytes.clone();
        with_nan[at..at + nan.len()].copy_from_slice(&nan);
        assert!(
            postcard::from_bytes::<I>(&with_nan).is_err(),
            "{name}: {float} made NaN"
        );
    }
}

#[test]
fn every_type_of_one_variable_comes_back_and_a_broken_form_does_not() {
    let mut rng = Lcg(21);
    let x = uneven(&mut rng, 40);
    let y: Vec<f64> = x
        .iter()
        .map(|x| (0.3 * x).sin() + 0.01 * rng.next())
        .collect();
    let slopes: Vec<f64> = x.iter().map(|x| 0.3 * (0.3 * x).cos()).collect();
    let (low, high) = (x[0], x[x.len() - 1]);
    // Inside, on every data point, outside both ends, and no number.
    let mut points: Vec<f64> = (0..2000)
        .map(|_| low - 5.0 + rng.next() * (high - low + 10.0))
        .collect();
    points.extend(&x);
    let mut sorted = points.clone();
    sorted.sort_by(f64::total_cmp);
    points.extend(sorted);
    points.push(f64::NAN);
    let modes = [
        Extrapolate::Forbid,
        Extrapolate::Clamp,
        Extrapolate::Extend,
        Extrapolate::Periodic,
    ];
    let swapped = |place: &'static str| [(place, Fault::Swapped)];

    for mode in modes {
        check(
            &format!("Linear1D, {mode:?}"),
            &Linear1D::new(&x, &y).expect("valid").with_extrapolate(mode),
            &points,
            &swapped("/poly/breakpoints"),
            false,
        );
        for kind in [
            NearestKind::Nearest,
            NearestKind::NearestUp,
            NearestKind::Previous,
            NearestKind::Next,
        ] {
            let nearest = Nearest1D::new(&x, &y, kind).expect("valid");
            let back: Nearest1D<f64> = serde_json::from_value(json!(nearest)).expect("read");
            assert_eq!(
                (back.kind(), back.extrapolate()),
                (kind, Extrapolate::Forbid)
            );
            check(
                &format!("Nearest1D, {kind:?}, {mode:?}"),
                &nearest.with_extrapolate(mode),
                &points,
                &swapped("/x"),
                false,
            );
        }
        let cubic = CubicSpline::new(&x, &y, BoundaryCondition::NotAKnot).expect("valid");
        check(
            &format!("CubicSpline, {mode:?}"),
            &cubic.clone().with_extrapolate(mode),
            &points,
            &swapped("/poly/breakpoints"),
            false,
        );
        check(
            &format!("PchipInterpolator, {mode:?}"),
            &PchipInterpolator::new(&x, &y)
                .expect("valid")
                .with_extrapolate(mode),
            &points,
            &swapped("/poly/breakpoints"),
            false,
        );
        for method in [AkimaMethod::Akima, AkimaMethod::Makima] {
            let akima = Akima1D::new(&x, &y, method).expect("valid");
            let back: Akima1D<f64> = serde_json::from_value(json!(akima)).expect("read");
            assert_eq!(back.method(), method);
            check(
                &format!("Akima1D, {method:?}, {mode:?}"),
                &akima.with_extrapolate(mode),
                &points,
                &swapped("/poly/breakpoints"),
                false,
            );
        }
        check(
            &format!("CubicHermiteSpline, {mode:?}"),
            &CubicHermiteSpline::new(&x, &y, &slopes)
                .expect("valid")
                .with_extrapolate(mode),
            &points,
            &swapped("/poly/breakpoints"),
            false,
        );
        for degree in [1, 3, 5] {
            let spline = BSpline1D::new(&x, &y, degree, BoundaryCondition::NotAKnot);
            let spline = spline.expect("valid").with_extrapolate(mode);
            check(
                &format!("BSpline1D, degree {degree}, {mode:?}"),
                &spline,
                &points,
                &[("/spline/knots", Fault::Reversed)],
                false,
            );
            check(
                &format!("BSpline, degree {degree}, {mode:?}"),
                spline.as_bspline(),
                &points,
                &[("/knots", Fault::Reversed)],
                false,
            );
        }
        check(
            &format!("PPoly, {mode:?}"),
            &cubic.as_ppoly().clone().with_extrapolate(mode),
            &points,
            &swapped("/breakpoints"),
            false,
        );
    }

    // A periodic spline, and one with a slope given at an end.
    let mut closed = y.clone();
    closed[y.len() - 1] = y[0];
    check(
        "CubicSpline, periodic",
        &CubicSpline::new(&x, &closed, BoundaryCondition::Periodic).expect("valid"),
        &points,
        &swapped("/poly/breakpoints"),
        false,
    );
    let ends = BoundaryCondition::Ends(
        EndCondition::FirstDerivative(0.5),
        EndCondition::SecondDerivative(0.0),
    );
    check(
        "BSpline1D, ends",
        &BSpline1D::new(&x, &y, 3, ends).expect("valid"),
        &points,
        &[("/spline/knots", Fault::Reversed)],
        false,
    );

    // Values of two columns, a curve, and `f32`.
    let pairs: Vec<[f64; 2]> = y.iter().map(|v| [*v, 1.0 - v]).collect();
    check(
        "Linear1D, two columns",
        &Linear1D::new(&x, &pairs).expect("valid"),
        &points,
        &swapped("/poly/breakpoints"),
        false,
    );
    check(
        "CubicSpline, two columns",
        &CubicSpline::new(&x, &pairs, BoundaryCondition::NotAKnot).expect("valid"),
        &points,
        &swapped("/poly/breakpoints"),
        false,
    );
    for parametrization in [
        Parametrization::Uniform,
        Parametrization::ChordLength,
        Parametrization::Centripetal,
    ] {
        let curve = BSplineCurve::new(&pairs, 3, BoundaryCondition::NotAKnot, parametrization);
        let curve = curve.expect("valid").with_extrapolate(Extrapolate::Extend);
        let (start, end) = curve.domain();
        let mut along: Vec<f64> = (0..500)
            .map(|_| start - 0.1 + rng.next() * (end - start + 0.2))
            .collect();
        along.extend(curve.parameters());
        let back: BSplineCurve<f64, 2> = serde_json::from_value(json!(curve)).expect("read");
        assert_eq!(back.parameters(), curve.parameters());
        check(
            &format!("BSplineCurve, {parametrization:?}"),
            &curve,
            &along,
            &swapped("/parameters"),
            false,
        );
    }
    let x32: Vec<f32> = x.iter().map(|x| *x as f32).collect();
    let y32: Vec<f32> = y.iter().map(|y| *y as f32).collect();
    let points32: Vec<f32> = points.iter().map(|x| *x as f32).collect();
    check(
        "CubicSpline, f32",
        &CubicSpline::new(&x32, &y32, BoundaryCondition::NotAKnot).expect("valid"),
        &points32,
        &swapped("/poly/breakpoints"),
        true,
    );
    check(
        "Nearest1D, f32",
        &Nearest1D::new(&x32, &y32, NearestKind::Nearest).expect("valid"),
        &points32,
        &swapped("/x"),
        true,
    );
    check(
        "BSpline1D, f32",
        &BSpline1D::new(&x32, &y32, 3, BoundaryCondition::NotAKnot).expect("valid"),
        &points32,
        &[("/spline/knots", Fault::Reversed)],
        true,
    );
}

#[test]
fn a_type_does_not_read_the_polynomial_of_another_degree() {
    let x = [0.0, 1.0, 2.5, 4.0, 5.0, 7.0];
    let y = [1.0, 3.0, 2.0, 5.0, 4.0, 4.5];
    let line = Linear1D::new(&x, &y).expect("valid");
    let cubic = CubicSpline::new(&x, &y, BoundaryCondition::NotAKnot).expect("valid");
    let of_a_line = json!({ "version": 1, "poly": line.as_ppoly() });
    let of_a_cubic = json!({ "version": 1, "poly": cubic.as_ppoly() });

    // The form of each type, read back: the forms below are of that shape.
    assert!(serde_json::from_value::<CubicSpline<f64>>(of_a_cubic.clone()).is_ok());
    assert!(serde_json::from_value::<PchipInterpolator<f64>>(of_a_cubic.clone()).is_ok());
    assert!(serde_json::from_value::<CubicHermiteSpline<f64>>(of_a_cubic.clone()).is_ok());
    let mut akima = of_a_cubic.clone();
    akima["method"] = json!(AkimaMethod::Makima);
    assert!(serde_json::from_value::<Akima1D<f64>>(akima).is_ok());
    let mut with_last = of_a_line.clone();
    with_last["last"] = json!(4.5);
    assert!(serde_json::from_value::<Linear1D<f64>>(with_last).is_ok());

    assert!(serde_json::from_value::<CubicSpline<f64>>(of_a_line.clone()).is_err());
    assert!(serde_json::from_value::<PchipInterpolator<f64>>(of_a_line.clone()).is_err());
    assert!(serde_json::from_value::<CubicHermiteSpline<f64>>(of_a_line.clone()).is_err());
    let mut akima = of_a_line;
    akima["method"] = json!(AkimaMethod::Akima);
    assert!(serde_json::from_value::<Akima1D<f64>>(akima).is_err());
    let mut with_last = of_a_cubic;
    with_last["last"] = json!(4.5);
    assert!(serde_json::from_value::<Linear1D<f64>>(with_last).is_err());

    // `BSpline` has a degree of 0; `BSpline1D` does not.
    let steps = BSpline::new(&[0.0, 1.0, 2.0], &[5.0, 6.0], 0).expect("valid");
    let of_steps = json!({ "version": 1, "spline": steps });
    assert!(serde_json::from_value::<BSpline1D<f64>>(of_steps).is_err());
    let broken = BSpline1D::new(&x, &y, 1, BoundaryCondition::NotAKnot).expect("valid");
    let of_a_polygon = json!({ "version": 1, "spline": broken.as_bspline() });
    assert!(serde_json::from_value::<BSpline1D<f64>>(of_a_polygon).is_ok());
}

#[test]
fn the_stored_form_of_a_polynomial_is_this_text() {
    // The form of version 1, as a reader of another language would see it.
    let poly = PPoly::new(&[0.0, 1.0, 3.0], &[2.0, 1.0, -1.0, 3.0], 1).expect("valid");
    let text = serde_json::to_string(&poly.with_extrapolate(Extrapolate::Clamp)).expect("stored");
    assert_eq!(
        text,
        r#"{"version":1,"breakpoints":[0.0,1.0,3.0],"coefficients":[2.0,1.0,-1.0,3.0],"degree":1,"extrapolate":"Clamp"}"#
    );
    // Values of several columns are arrays within the array.
    let pairs = PPoly::new(&[0.0, 1.0], &[[1.0, -1.0], [0.0, 5.0]], 1).expect("valid");
    assert_eq!(
        serde_json::to_string(&pairs).expect("stored"),
        r#"{"version":1,"breakpoints":[0.0,1.0],"coefficients":[[1.0,-1.0],[0.0,5.0]],"degree":1,"extrapolate":"Extend"}"#
    );
}

#[test]
fn grids_come_back_and_a_broken_form_does_not() {
    let mut rng = Lcg(22);
    let axes = [
        uneven(&mut rng, 9),
        uneven(&mut rng, 7),
        uneven(&mut rng, 8),
    ];
    let values: Vec<f64> = (0..9 * 7 * 8_u32)
        .map(|n| f64::from(n * 37 % 101) / 17.0 + 0.001 * rng.next())
        .collect();
    let mut points: Vec<[f64; 3]> = (0..3000)
        .map(|_| core::array::from_fn(|d| -0.3 + rng.next() * (axes[d][axes[d].len() - 1] + 0.6)))
        .collect();
    // The points of the grid themselves, and a path.
    for i in 0..9 * 7 * 8 {
        points.push([axes[0][i / 56], axes[1][i / 8 % 7], axes[2][i % 8]]);
    }
    points.extend((0..1000).map(|i| {
        let share = f64::from(i) / 1000.0;
        core::array::from_fn(|d| share * axes[d][axes[d].len() - 1])
    }));
    points.push([f64::NAN, 1.0, 1.0]);
    let slices: [&[f64]; 3] = core::array::from_fn(|d| axes[d].as_slice());
    let faults = [
        ("/axes/0", Fault::Swapped),
        ("/axes/2", Fault::Repeated),
        // Turned around without its values: the form holds axes ascending.
        ("/axes/1", Fault::Reversed),
    ];
    for method in [
        GridMethod::Nearest,
        GridMethod::Linear,
        GridMethod::Cubic,
        GridMethod::Quintic,
        GridMethod::Pchip,
    ] {
        for mode in [Extrapolate::Forbid, Extrapolate::Clamp, Extrapolate::Extend] {
            let grid = RegularGridInterpolator::new(slices, &values, method).expect("valid");
            let grid = grid.with_extrapolate(mode);
            let back: RegularGridInterpolator<f64, 3> =
                serde_json::from_value(json!(grid)).expect("read");
            assert_eq!((back.method(), back.extrapolate()), (method, mode));
            assert_eq!(back.domain(), grid.domain());
            check(
                &format!("grid, {method:?}, {mode:?}"),
                &grid,
                &points,
                &faults,
                false,
            );
        }
        // What a spline method solved for belongs to no other method.
        let grid = RegularGridInterpolator::new(slices, &values, method).expect("valid");
        for other in [GridMethod::Linear, GridMethod::Cubic, GridMethod::Pchip] {
            let mut form = json!(grid);
            form["method"] = json!(other);
            let held = |method| !matches!(method, GridMethod::Nearest | GridMethod::Linear);
            let fits = other == method
                || held(method) == held(other) && method != GridMethod::Nearest
                || !held(method) && !held(other);
            let read = serde_json::from_value::<RegularGridInterpolator<f64, 3>>(form);
            assert_eq!(read.is_ok(), fits, "{method:?} read as {other:?}");
        }
    }

    // Without axes a grid is one value, whatever the method.
    for method in [GridMethod::Nearest, GridMethod::Linear, GridMethod::Pchip] {
        let point = RegularGridInterpolator::<f64, 0>::new([], &[5.0], method).expect("valid");
        let back: RegularGridInterpolator<f64, 0> =
            serde_json::from_value(json!(point)).expect("read");
        assert_eq!((back.method(), back.eval(&[])), (method, Some(5.0)));
        let mut form = json!(point);
        form["solved"] = json!([1.0]);
        assert!(serde_json::from_value::<RegularGridInterpolator<f64, 0>>(form).is_err());
    }

    // A descending axis is held ascending, and so stored; an axis of one
    // point; values of two columns; `f32`.
    let down: Vec<f64> = axes[1].iter().rev().copied().collect();
    let turned = [&axes[0][..], &down[..], &axes[2][..]];
    let grid = RegularGridInterpolator::new(turned, &values, GridMethod::Cubic).expect("valid");
    check(
        "grid, an axis given descending",
        &grid,
        &points,
        &faults,
        false,
    );
    let flat = [&axes[0][..], &axes[1][..1], &axes[2][..]];
    let on_a_plane: Vec<[f64; 3]> = points.iter().map(|p| [p[0], axes[1][0], p[2]]).collect();
    let grid = RegularGridInterpolator::new(flat, &values[..72], GridMethod::Linear);
    check(
        "grid, an axis of one point",
        &grid.expect("valid"),
        &on_a_plane,
        &[("/axes/0", Fault::Swapped)],
        false,
    );
    let pairs: Vec<[f64; 2]> = values.iter().map(|v| [*v, 1.0 - v]).collect();
    let grid = RegularGridInterpolator::new(slices, &pairs, GridMethod::Pchip).expect("valid");
    check("grid, two columns", &grid, &points, &faults, false);
    let axes32: [Vec<f32>; 3] =
        core::array::from_fn(|d| axes[d].iter().map(|c| *c as f32).collect());
    let values32: Vec<f32> = values.iter().map(|v| *v as f32).collect();
    let points32: Vec<[f32; 3]> = points.iter().map(|p| p.map(|c| c as f32)).collect();
    let grid = RegularGridInterpolator::from_owned(axes32, values32, GridMethod::Quintic);
    check("grid, f32", &grid.expect("valid"), &points32, &faults, true);
}

#[cfg(feature = "scattered")]
#[test]
fn scattered_types_come_back_and_a_broken_form_does_not() {
    use interpx::{
        LinearNDInterpolator, NearestNDInterpolator, RBFInterpolator, RbfDegree, RbfKernel,
        RbfOptions,
    };

    let mut rng = Lcg(23);
    let mut data: Vec<[f64; 2]> = (0..80).map(|_| [rng.next(), 5.0 * rng.next()]).collect();
    let values: Vec<f64> = data
        .iter()
        .map(|p| (3.0 * p[0]).sin() * p[1] + 0.3)
        .collect();
    let mut points: Vec<[f64; 2]> = (0..400).map(|_| [rng.next(), 5.0 * rng.next()]).collect();
    points.extend(&data);
    points.extend((0..100).map(|i| [f64::from(i) / 100.0, 2.0]));
    points.push([f64::INFINITY, 0.5]);

    let kernels = [
        RbfKernel::Linear,
        RbfKernel::ThinPlateSpline,
        RbfKernel::Cubic,
        RbfKernel::Quintic,
        RbfKernel::Multiquadric { epsilon: 1.5 },
        RbfKernel::InverseMultiquadric { epsilon: 2.5 },
        RbfKernel::InverseQuadratic { epsilon: 3.5 },
        RbfKernel::Gaussian { epsilon: 4.5 },
    ];
    let repeated = [("/points", Fault::Repeated)];
    for kernel in kernels {
        for options in [
            RbfOptions::new(),
            RbfOptions::new().with_smoothing(0.25),
            RbfOptions::new().with_degree(RbfDegree::Degree(2)),
            RbfOptions::new().with_neighbors(12),
            RbfOptions::new().with_neighbors(15).with_smoothing(0.125),
        ] {
            let rbf = RBFInterpolator::with_options(&data, &values, kernel, &options);
            let rbf = rbf.expect("valid");
            let back: RBFInterpolator<2> = serde_json::from_value(json!(rbf)).expect("read");
            assert_eq!((back.kernel(), back.len()), (kernel, 80));
            check(
                &format!("RBF, {kernel:?}, {options:?}"),
                &rbf,
                &points,
                &repeated,
                false,
            );
        }
    }
    let no_polynomial = RbfOptions::new().with_degree(RbfDegree::None);
    let rbf = RBFInterpolator::with_options(
        &data,
        &values,
        RbfKernel::Gaussian { epsilon: 3.0 },
        &no_polynomial,
    );
    check(
        "RBF, no polynomial",
        &rbf.expect("valid"),
        &points,
        &repeated,
        false,
    );
    let pairs: Vec<[f64; 2]> = values.iter().map(|v| [*v, 1.0 - v]).collect();
    let rbf = RBFInterpolator::new(&data, &pairs, RbfKernel::ThinPlateSpline).expect("valid");
    check("RBF, two columns", &rbf, &points, &repeated, false);

    // As many neighbors as points is one system; a form that says so of
    // neighbors is of no interpolant that was built.
    let mut form = json!(RBFInterpolator::<2>::with_options(
        &data,
        &values,
        RbfKernel::Cubic,
        &RbfOptions::new().with_neighbors(12)
    )
    .expect("valid"));
    form["neighbors"] = json!(80);
    assert!(serde_json::from_value::<RBFInterpolator<2>>(form.clone()).is_err());
    form["neighbors"] = json!(0);
    assert!(serde_json::from_value::<RBFInterpolator<2>>(form).is_err());

    // Without a polynomial there are as many coefficients as values: only
    // the number of neighbors tells the values of the one from the
    // coefficients of the other.
    let gaussian = RbfKernel::Gaussian { epsilon: 3.0 };
    let near = no_polynomial.with_neighbors(12);
    let near = RBFInterpolator::<2>::with_options(&data, &values, gaussian, &near);
    let mut form = json!(near.expect("valid"));
    form["neighbors"] = json!(80);
    assert!(serde_json::from_value::<RBFInterpolator<2>>(form).is_err());

    // The estimate of the condition is stored with one system (it is not
    // solved again to be read), and belongs to none with neighbors.
    let whole = RBFInterpolator::<2>::new(&data, &values, RbfKernel::Cubic).expect("valid");
    let estimate = whole.reciprocal_condition().expect("one system");
    assert!(estimate > 0.0 && estimate < 1.0, "{estimate}");
    let back: RBFInterpolator<2> = serde_json::from_value(json!(whole)).expect("read");
    assert_eq!(back.reciprocal_condition(), Some(estimate));
    let bytes = postcard::to_stdvec(&whole).expect("stored");
    let back: RBFInterpolator<2> = postcard::from_bytes(&bytes).expect("read");
    assert_eq!(back.reciprocal_condition(), Some(estimate));
    for wrong in [json!(1.5), json!(-0.25), json!(null)] {
        let mut form = json!(whole);
        form["reciprocal_condition"] = wrong.clone();
        assert!(
            serde_json::from_value::<RBFInterpolator<2>>(form).is_err(),
            "{wrong}"
        );
    }
    for edge in [0.0, 1.0] {
        let mut form = json!(whole);
        form["reciprocal_condition"] = json!(edge);
        let back: RBFInterpolator<2> = serde_json::from_value(form).expect("read");
        assert_eq!(back.reciprocal_condition(), Some(edge));
    }
    let options = RbfOptions::new().with_neighbors(12);
    let near = RBFInterpolator::<2>::with_options(&data, &values, RbfKernel::Cubic, &options);
    let near = near.expect("valid");
    assert_eq!(near.reciprocal_condition(), None);
    let back: RBFInterpolator<2> = serde_json::from_value(json!(near)).expect("read");
    assert_eq!(back.reciprocal_condition(), None);
    let mut form = json!(near);
    form["reciprocal_condition"] = json!(0.5);
    assert!(serde_json::from_value::<RBFInterpolator<2>>(form).is_err());

    // Points given twice: the tree has each once, and the first value. A
    // point of its own is given last, so that the last value is one the tree
    // points to (values after that one belong to no point of the form, and a
    // form without them is the form of fewer points given).
    data.extend_from_within(..20);
    data.push([0.5, 9.0]);
    let mut twice = values.clone();
    twice.extend((0..21).map(|i| -f64::from(i) - 1.0));
    let faults = [("/points", Fault::Repeated), ("/first", Fault::Swapped)];
    let nearest = NearestNDInterpolator::new(&data, &twice).expect("valid");
    check("NearestND", &nearest, &points, &faults, false);
    check(
        "NearestND, a bound",
        &nearest.clone().with_distance_bound(0.05),
        &points,
        &faults,
        false,
    );
    let rescaled = NearestNDInterpolator::rescaled(&data, &twice).expect("valid");
    let back: NearestNDInterpolator<2> = serde_json::from_value(json!(rescaled)).expect("read");
    assert_eq!(
        (back.len(), back.values()),
        (rescaled.len(), rescaled.values())
    );
    check(
        "NearestND, rescaled",
        &rescaled.with_distance_bound(0.25),
        &points,
        &faults,
        false,
    );
    let nearest = NearestNDInterpolator::new(&data, &pairs_of(&twice)).expect("valid");
    check("NearestND, two columns", &nearest, &points, &faults, false);

    // The places of the values of the points: within the values.
    let mut form = json!(NearestNDInterpolator::new(&data, &twice).expect("valid"));
    let last = form["first"].as_array().expect("places").len() - 1;
    form["first"][last] = json!(twice.len());
    assert!(serde_json::from_value::<NearestNDInterpolator<2>>(form).is_err());
    // A scale that is zero divides every coordinate by it.
    let mut form = json!(NearestNDInterpolator::rescaled(&data, &twice).expect("valid"));
    form["scaling"]["scale"][0] = json!(0.0);
    assert!(serde_json::from_value::<NearestNDInterpolator<2>>(form).is_err());

    // A bound of positive infinity is no bound: the interpolant holds none,
    // and comes back with none.
    let unbounded = NearestNDInterpolator::new(&data, &twice).expect("valid");
    let infinite = unbounded.clone().with_distance_bound(f64::INFINITY);
    assert_eq!(infinite.distance_bound(), None);
    let back: NearestNDInterpolator<2> = serde_json::from_value(json!(infinite)).expect("read");
    assert_eq!(back.distance_bound(), None);
    assert_eq!(back.eval_batch(&points), infinite.eval_batch(&points));
    let bytes = postcard::to_stdvec(&infinite).expect("stored");
    let back: NearestNDInterpolator<2> = postcard::from_bytes(&bytes).expect("read");
    assert_eq!(back.distance_bound(), None);
    // With a bound that is NaN or negative infinity no query has a value. A
    // text format would write `null` for either and read back no bound, with
    // which every query has one: neither is stored, in any format.
    for none in [f64::NAN, f64::NEG_INFINITY] {
        let never = unbounded.clone().with_distance_bound(none);
        assert_eq!(never.eval(&data[0]), None);
        assert!(serde_json::to_string(&never).is_err(), "{none}");
        assert!(postcard::to_stdvec(&never).is_err(), "{none}");
    }
    // A stored form that says so is not read either (a binary format can).
    let finite = unbounded.clone().with_distance_bound(0.375);
    let bytes = postcard::to_stdvec(&finite).expect("stored");
    let at = position(&bytes, &0.375_f64.to_le_bytes()).expect("the bound");
    for none in [f64::NAN, f64::NEG_INFINITY, f64::INFINITY] {
        let mut broken = bytes.clone();
        broken[at..at + 8].copy_from_slice(&none.to_le_bytes());
        assert!(
            postcard::from_bytes::<NearestNDInterpolator<2>>(&broken).is_err(),
            "{none}"
        );
    }

    // Without dimensions all points are one point.
    let nowhere = NearestNDInterpolator::new(&[[], []], &[7.0, 8.0]).expect("valid");
    let back: NearestNDInterpolator<0> = serde_json::from_value(json!(nowhere)).expect("read");
    assert_eq!(back.eval(&[]), Some(7.0));

    // LinearNDInterpolator keeps its triangulation: points given twice stay,
    // and no simplex names the second of them.
    let linear = LinearNDInterpolator::new(&data, &twice).expect("valid");
    check("LinearND", &linear, &points, &[], false);
    let rescaled = LinearNDInterpolator::rescaled(&data, &twice).expect("valid");
    check("LinearND, rescaled", &rescaled, &points, &[], false);
    let linear = LinearNDInterpolator::new(&data, &pairs_of(&twice)).expect("valid");
    check("LinearND, two columns", &linear, &points, &[], false);
    let back: LinearNDInterpolator<2, [f64; 2]> =
        serde_json::from_value(json!(linear)).expect("read");
    assert_eq!(
        (back.simplex_count(), back.points()),
        (linear.simplex_count(), linear.points())
    );
    // A simplex that names no point, a neighbor that names no simplex, a
    // scale of zero: each refused.
    let form = json!(LinearNDInterpolator::rescaled(&data, &twice).expect("valid"));
    let mut broken = form.clone();
    broken["simplices"][0] = json!(data.len());
    assert!(serde_json::from_value::<LinearNDInterpolator<2>>(broken).is_err());
    let mut broken = form.clone();
    let count = form["simplices"].as_array().expect("simplices").len() / 3;
    broken["neighbors"][0] = json!(count);
    assert!(serde_json::from_value::<LinearNDInterpolator<2>>(broken).is_err());
    // A neighbor taken away on one side only: the other simplex still names
    // this one (the review of #183).
    let first = form["neighbors"][0].as_u64().expect("a neighbor");
    let mut broken = form.clone();
    broken["neighbors"][0] = json!(if first == u64::from(u32::MAX) {
        0
    } else {
        u32::MAX
    });
    assert!(serde_json::from_value::<LinearNDInterpolator<2>>(broken).is_err());
    // A simplex that names a point twice.
    let mut broken = form.clone();
    broken["simplices"][1] = broken["simplices"][0].clone();
    assert!(serde_json::from_value::<LinearNDInterpolator<2>>(broken).is_err());
    let mut broken = form.clone();
    broken["scaling"]["scale"][1] = json!(0.0);
    assert!(serde_json::from_value::<LinearNDInterpolator<2>>(broken).is_err());
    let mut broken = form;
    broken["simplices"] = json!([]);
    broken["neighbors"] = json!([]);
    assert!(serde_json::from_value::<LinearNDInterpolator<2>>(broken).is_err());
    // Without dimensions: no simplex, and the first value.
    let nowhere = LinearNDInterpolator::new(&[[], []], &[7.0, 8.0]).expect("valid");
    let back: LinearNDInterpolator<0> = serde_json::from_value(json!(nowhere)).expect("read");
    assert_eq!(back.eval(&[]), Some(7.0));
}

/// Values of two columns from values of one.
#[cfg(feature = "scattered")]
fn pairs_of(values: &[f64]) -> Vec<[f64; 2]> {
    values.iter().map(|v| [*v, 1.0 - v]).collect()
}

#[test]
fn the_parameters_come_back() {
    /// Through JSON and through `postcard`.
    #[cfg(test)]
    fn back<P: Serialize + DeserializeOwned + PartialEq + Debug>(parameter: P) {
        let text = serde_json::to_string(&parameter).expect("stored");
        assert_eq!(serde_json::from_str::<P>(&text).expect("read"), parameter);
        let bytes = postcard::to_stdvec(&parameter).expect("stored");
        assert_eq!(postcard::from_bytes::<P>(&bytes).expect("read"), parameter);
    }

    for mode in [
        Extrapolate::Forbid,
        Extrapolate::Clamp,
        Extrapolate::Extend,
        Extrapolate::Periodic,
    ] {
        back(mode);
    }
    for method in [
        GridMethod::Nearest,
        GridMethod::Linear,
        GridMethod::Cubic,
        GridMethod::Quintic,
        GridMethod::Pchip,
    ] {
        back(method);
    }
    for kind in [
        NearestKind::Nearest,
        NearestKind::NearestUp,
        NearestKind::Previous,
        NearestKind::Next,
    ] {
        back(kind);
    }
    back(AkimaMethod::Akima);
    back(AkimaMethod::Makima);
    back(Parametrization::Uniform);
    back(Parametrization::ChordLength);
    back(Parametrization::Centripetal);
    back(EndCondition::<f64>::NotAKnot);
    back(EndCondition::FirstDerivative([1.5, -2.0]));
    back(EndCondition::SecondDerivative(0.25_f32));
    back(BoundaryCondition::<f64>::NotAKnot);
    back(BoundaryCondition::<[f64; 3]>::Periodic);
    back(BoundaryCondition::Ends(
        EndCondition::FirstDerivative(1.0),
        EndCondition::NotAKnot,
    ));
    #[cfg(feature = "scattered")]
    {
        use interpx::{RbfDegree, RbfKernel, RbfOptions};
        back(RbfKernel::Quintic);
        back(RbfKernel::Gaussian { epsilon: 0.5 });
        back(RbfDegree::Default);
        back(RbfDegree::None);
        back(RbfDegree::Degree(3));
        back(RbfOptions::new());
        back(
            RbfOptions::new()
                .with_smoothing(0.5)
                .with_degree(RbfDegree::Degree(1))
                .with_neighbors(9),
        );
    }
}
