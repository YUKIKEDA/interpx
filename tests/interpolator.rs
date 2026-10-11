//! The provided methods of `Interpolator`, exercised through small interpolators
//! defined here, and the thread-safety contract of the public types.

use interpx::{
    BoundaryCondition, EndCondition, Extrapolate, InterpolateError, Interpolator, Scalar, Value,
};

/// `2x` on `[0, 1]`, for either scalar type.
struct Double;

impl Interpolator for Double {
    type Point = f64;
    type Output = f64;

    fn eval(&self, x: &f64) -> Option<f64> {
        (0.0..=1.0).contains(x).then(|| 2.0 * x)
    }
}

/// A two-output interpolator on the unit square: `(x + y, x - y)`.
struct Pair;

impl Interpolator for Pair {
    type Point = [f32; 2];
    type Output = [f32; 2];

    fn eval(&self, p: &[f32; 2]) -> Option<[f32; 2]> {
        p.iter()
            .all(|c| (0.0..=1.0).contains(c))
            .then(|| [p[0] + p[1], p[0] - p[1]])
    }
}

fn bits(v: &[f64]) -> Vec<u64> {
    v.iter().map(|x| x.to_bits()).collect()
}

#[test]
fn eval_is_none_outside_the_domain_and_for_nan() {
    assert_eq!(
        Double.eval(&0.25).map(f64::to_bits),
        Some(0.5_f64.to_bits())
    );
    assert!(Double.eval(&1.5).is_none());
    assert!(Double.eval(&-0.5).is_none());
    assert!(Double.eval(&f64::NAN).is_none());
    assert!(Double.eval(&f64::INFINITY).is_none());
}

#[test]
fn eval_or_returns_the_fill_only_where_eval_is_none() {
    assert_eq!(Double.eval_or(&0.5, -1.0).to_bits(), 1.0_f64.to_bits());
    assert_eq!(Double.eval_or(&2.0, -1.0).to_bits(), (-1.0_f64).to_bits());
    assert_eq!(
        Double.eval_or(&f64::NAN, -1.0).to_bits(),
        (-1.0_f64).to_bits()
    );
    // A NaN fill marks the points without a value.
    assert!(Double.eval_or(&2.0, f64::NAN).is_nan());
}

#[test]
fn eval_batch_keeps_the_order_and_the_gaps() {
    let got = Double.eval_batch(&[0.0, 2.0, 0.5, f64::NAN, 1.0]);
    let want = [Some(0.0), None, Some(1.0), None, Some(2.0)];
    assert_eq!(got.len(), want.len());
    for (g, w) in got.iter().zip(want) {
        assert_eq!(g.map(f64::to_bits), w.map(f64::to_bits));
    }
    assert!(Double.eval_batch(&[]).is_empty());
}

#[test]
fn eval_batch_into_equals_the_loop_over_eval_or() {
    let points = [0.0, 2.0, 0.5, f64::NAN, 1.0, -3.0];
    let mut out = [f64::NAN; 6];
    assert_eq!(Double.eval_batch_into(&points, &mut out, 9.0), Ok(()));
    let looped: Vec<f64> = points.iter().map(|p| Double.eval_or(p, 9.0)).collect();
    assert_eq!(bits(&out), bits(&looped));
}

#[test]
fn eval_batch_into_rejects_a_length_mismatch_and_leaves_the_buffer() {
    // The buffer is longer, then shorter, than the points.
    for (n_points, n_out) in [(2, 3), (3, 2), (0, 1), (1, 0)] {
        let points = vec![0.5; n_points];
        let mut out = vec![7.0; n_out];
        assert_eq!(
            Double.eval_batch_into(&points, &mut out, 0.0),
            Err(InterpolateError::LengthMismatch {
                expected: n_points,
                actual: n_out,
            })
        );
        assert_eq!(bits(&out), bits(&vec![7.0; n_out]), "buffer was written");
    }
}

#[test]
fn eval_batch_into_accepts_an_empty_batch() {
    let mut out: [f64; 0] = [];
    assert_eq!(Double.eval_batch_into(&[], &mut out, 0.0), Ok(()));
}

#[test]
fn an_array_output_goes_through_the_same_methods() {
    let inside = [0.75_f32, 0.25];
    let outside = [0.5_f32, 2.0];
    let nan = [f32::NAN, 0.5];
    let fill = [-1.0_f32, -1.0];

    let v = Pair.eval(&inside).expect("inside the square");
    assert_eq!(v.map(f32::to_bits), [1.0_f32, 0.5].map(f32::to_bits));
    assert!(Pair.eval(&outside).is_none());
    assert!(Pair.eval(&nan).is_none());
    assert_eq!(
        Pair.eval_or(&nan, fill).map(f32::to_bits),
        fill.map(f32::to_bits)
    );

    let mut out = [[0.0_f32; 2]; 3];
    assert_eq!(
        Pair.eval_batch_into(&[inside, outside, nan], &mut out, fill),
        Ok(())
    );
    assert_eq!(out[0].map(f32::to_bits), [1.0_f32, 0.5].map(f32::to_bits));
    assert_eq!(out[1].map(f32::to_bits), fill.map(f32::to_bits));
    assert_eq!(out[2].map(f32::to_bits), fill.map(f32::to_bits));
    assert_eq!(Pair.eval_batch(&[inside, outside]).len(), 2);
}

/// A linear combination of two values, written against `Value` alone: what an
/// interpolator does with its data.
fn lerp<T: Scalar, V: Value<T>>(a: V, b: V, t: T) -> V {
    V::zero().axpy(T::one() - t, a).axpy(t, b)
}

#[test]
fn value_is_enough_to_write_an_interpolation() {
    assert_eq!(lerp(1.0_f64, 3.0, 0.25).to_bits(), 1.5_f64.to_bits());
    assert_eq!(lerp(1.0_f32, 3.0, 0.25).to_bits(), 1.5_f32.to_bits());
    let v = lerp([0.0_f64, 10.0], [4.0, 30.0], 0.5);
    assert_eq!(v.map(f64::to_bits), [2.0_f64, 20.0].map(f64::to_bits));
    assert!(v.all_finite());
}

fn assert_send_sync<T: Send + Sync>() {}

#[test]
fn public_types_are_send_and_sync() {
    // docs/design.ja.md section 3.4: every interpolator is Send + Sync. The types
    // interpolators are built from must be, too.
    assert_send_sync::<InterpolateError>();
    assert_send_sync::<Extrapolate>();
    assert_send_sync::<BoundaryCondition<f32>>();
    assert_send_sync::<BoundaryCondition<f64>>();
    assert_send_sync::<EndCondition<f64>>();

    fn scalars_are<T: Scalar>() {
        assert_send_sync::<T>();
    }
    scalars_are::<f32>();
    scalars_are::<f64>();

    fn values_are<T: Scalar, V: Value<T>>() {
        assert_send_sync::<V>();
    }
    values_are::<f64, f64>();
    values_are::<f32, [f32; 3]>();
}

#[test]
fn the_trait_can_be_used_as_a_bound_and_as_an_object() {
    fn total<I: Interpolator<Point = f64, Output = f64>>(f: &I, xs: &[f64]) -> f64 {
        xs.iter().map(|x| f.eval_or(x, 0.0)).sum()
    }
    assert_eq!(
        total(&Double, &[0.25, 0.5, 5.0]).to_bits(),
        1.5_f64.to_bits()
    );

    let boxed: Box<dyn Interpolator<Point = f64, Output = f64>> = Box::new(Double);
    assert_eq!(boxed.eval(&0.5).map(f64::to_bits), Some(1.0_f64.to_bits()));
}
