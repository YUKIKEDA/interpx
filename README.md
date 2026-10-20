English | [日本語](README.ja.md)

# interpx

Pure Rust interpolation: one variable, curves, regular grids, and scattered data. Inputs are checked at construction, and every interpolator owns its data and can be shared between threads.

```toml
[dependencies]
interpx = "0.1"
```

Features:

- `scattered` (on by default): `NearestNDInterpolator`, `LinearNDInterpolator`, `CloughTocher2DInterpolator`, and `RBFInterpolator`, and the dependencies only they need.
- `rayon`: batches evaluated on several threads (`ParInterpolator`), with the same results as the serial batches, bit for bit.
- `serde`: `Serialize` and `Deserialize` for every interpolator and the parameters it is built with.

## Example

```rust
use interpx::{BoundaryCondition, CubicSpline, Extrapolate, InterpolateError, Interpolator};

fn main() -> Result<(), InterpolateError> {
    let x = [0.0, 1.0, 2.0, 3.5, 5.0];
    let y = [1.0, 2.7, 5.8, 4.0, 2.2];

    // The not-a-knot cubic spline.
    let spline = CubicSpline::new(&x, &y, BoundaryCondition::NotAKnot)?;
    let at = spline.eval(&2.5);
    assert!(at.is_some());

    // No value outside the data, and a fill value at the call.
    let spline = spline.with_extrapolate(Extrapolate::Forbid);
    assert_eq!(spline.eval(&6.0), None);
    let mut out = [0.0; 3];
    spline.eval_batch_into(&[0.5, 4.0, 9.0], &mut out, f64::NAN)?;
    assert!(out[2].is_nan());
    Ok(())
}
```

The same with prints: `cargo run --example one_variable`. The other examples: `grid`, `scattered`, `parallel` (`--features rayon`), and `save_and_load` (`--features serde`).

## What a caller must know

- **The boundary condition is an argument.** `CubicSpline::new(x, y, bc)` and `BSpline1D::new(x, y, k, bc)` take it; `BoundaryCondition::default()` is `NotAKnot`.
- **The fill value is given at the call.** `eval` returns `Option`: `None` outside the domain when the interpolator does not extrapolate, and for a NaN query. `eval_or(point, fill)` and `eval_batch_into(points, out, fill)` take the value to fill.
- **What happens outside the domain is a mode, `Extrapolate`**, set with `with_extrapolate`. The default mode is set by the method of each interpolator:
  - `Extend` (the end pieces go on): `CubicSpline`, `PchipInterpolator`, `CubicHermiteSpline`, `BSpline1D`, `BSplineCurve`, `PPoly`, and `BSpline`. A spline built with a periodic boundary condition is `Periodic`.
  - `Forbid` (no value): `Akima1D`, `Nearest1D`, and `RegularGridInterpolator`.
  - `Clamp` (the end value): `Linear1D`.
- **Input is checked at construction.** In one variable the coordinates must be strictly increasing as given, and data must be finite everywhere: NaN or infinity in the data is an error (`InterpolateError::NonFinite`). A NaN of the data never ends up in the values.
- **Methods are chosen at construction**, not at each evaluation. `RBFInterpolator` takes its kernel as an argument.
- **Values are real.** A value can be a fixed-size array, `[T; M]`, for `M` columns interpolated at once.
- **Every interpolator owns its data** and is `Send + Sync + 'static`: it can be kept, moved to another thread, and shared. Constructors copy what they are given. The types that keep some of their data as it is given (`PPoly`, `BSpline`, `RegularGridInterpolator`, `NearestNDInterpolator`, `LinearNDInterpolator`, `CloughTocher2DInterpolator`, `RBFInterpolator`) also have `from_owned`, which takes that data as a `Vec` and keeps it without a copy.
- **Coordinates are `f32` or `f64`** in one variable, curves, and grids. Scattered data is `f64` only.
- **Grid values are in C order**: the last axis varies fastest.
- **NaN and infinity**: in the data they are an error at construction. A NaN query gives no value (`None`, or the fill value). An infinite query is outside the domain: no value under `Forbid` and `Periodic`, the end value under `Clamp`, and under `Extend` the value of the end piece there, which can be infinite or NaN.
- **`RBFInterpolator` over all points holds an `N` by `N` matrix while it is built**: 940 MB for 10 838 points. Over the nearest points (`RbfOptions::with_neighbors`) it solves a small system at each query instead. `reciprocal_condition` says how well the system over all points was solved.
- **Evaluation does not allocate**, except: `eval_batch` and `eval_batch_par`, which return a `Vec` (`eval_batch_par` one for each chunk too); `eval_batch_into_par`, where `rayon` allocates to share a batch out among threads; and an `RBFInterpolator` over neighbors, which allocates the working space a thread keeps at its first evaluation on that thread, and widens it when an interpolator with more neighbors is evaluated there.
- **0.x releases**: a minor version (0.1 to 0.2) may break the public API. The API is what the default features export.
- **Not yet in this crate**: derivatives of grid and scattered interpolants, roots of a piecewise polynomial, piecewise polynomials in several variables, barycentric and Krogh interpolation, smoothing splines, dimensions or output widths chosen at run time, `no_std`, and a fallback to a pseudo-inverse for an `RBFInterpolator` system that is singular.

## Usage

### `Interpolator`

Every interpolator implements `Interpolator`: `eval` at one point, `eval_or` with a fill value, `eval_batch` into a new `Vec`, `eval_batch_into` into a buffer of the caller's. A type of your own gets the batch methods by implementing `eval`.

```rust
use interpx::{InterpolateError, Interpolator};

/// 10 x on [0, 1], nothing outside.
struct Ramp;

impl Interpolator for Ramp {
    type Point = f64;
    type Output = f64;

    fn eval(&self, x: &f64) -> Option<f64> {
        (0.0..=1.0).contains(x).then(|| 10.0 * x)
    }
}

fn main() -> Result<(), InterpolateError> {
    assert_eq!(Ramp.eval(&0.5), Some(5.0));
    assert_eq!(Ramp.eval_or(&2.0, -1.0), -1.0);
    assert_eq!(Ramp.eval_batch(&[0.5, 2.0]), [Some(5.0), None]);
    let mut out = [0.0; 2];
    Ramp.eval_batch_into(&[0.5, 2.0], &mut out, -1.0)?;
    assert_eq!(out, [5.0, -1.0]);
    Ok(())
}
```

### `Hint`

`eval_hinted` remembers the interval of the last query in a `Hint`, so that queries in order find theirs in a step. The values are those of `eval`.

```rust
use interpx::{Hint, InterpolateError, Interpolator, Linear1D};

fn main() -> Result<(), InterpolateError> {
    let line = Linear1D::new(&[0.0, 1.0, 2.0, 3.0], &[0.0, 10.0, 20.0, 30.0])?;
    let mut hint = Hint::new();
    for x in [0.5, 1.5, 2.5] {
        assert_eq!(line.eval_hinted(x, &mut hint), line.eval(&x));
    }
    Ok(())
}
```

### `Extrapolate`

`Forbid` (no value), `Clamp` (the value at the nearest end), `Extend` (the end piece goes on), `Periodic` (the query is wrapped into the domain).

```rust
use interpx::{Extrapolate, InterpolateError, Interpolator, Linear1D};

fn main() -> Result<(), InterpolateError> {
    let line = Linear1D::new(&[0.0, 1.0, 2.0], &[0.0, 1.0, 0.0])?;
    assert_eq!(line.extrapolate(), Extrapolate::Clamp);
    assert_eq!(line.eval(&3.0), Some(0.0));
    assert_eq!(line.clone().with_extrapolate(Extrapolate::Forbid).eval(&3.0), None);
    assert_eq!(line.clone().with_extrapolate(Extrapolate::Extend).eval(&3.0), Some(-1.0));
    assert_eq!(line.with_extrapolate(Extrapolate::Periodic).eval(&2.5), Some(0.5));
    Ok(())
}
```

### Errors: `InterpolateError`, `AxisFault`, `DegeneracyKind`, `ParamReason`, `SolverError`

Constructors return `InterpolateError`, which says which input and where: `InsufficientPoints`, `LengthMismatch`, `NonFinite`, `NotSorted`, `DuplicatePoints`, `NotPeriodic`, `InvalidParameter` (with a `ParamReason`), `SolverFailed` (with a `SolverError`), `AllocationFailed`, `DegenerateGeometry` (with a `DegeneracyKind`, for scattered points that cannot be triangulated), and `GridAxis` (with an `AxisFault`, for one axis of a grid).

```rust
use interpx::{
    AxisFault, BSpline1D, BoundaryCondition, DegeneracyKind, GridMethod, InterpolateError,
    Linear1D, LinearNDInterpolator, ParamReason, RBFInterpolator, RbfKernel,
    RegularGridInterpolator, SolverError,
};

fn main() {
    let unsorted = Linear1D::new(&[0.0, 2.0, 1.0], &[0.0, 1.0, 2.0]);
    assert_eq!(unsorted.err(), Some(InterpolateError::NotSorted { index: 2 }));

    let degree = BSpline1D::new(&[0.0, 1.0, 2.0], &[0.0, 1.0, 0.0], 16, BoundaryCondition::NotAKnot);
    assert_eq!(
        degree.err(),
        Some(InterpolateError::InvalidParameter { name: "degree", reason: ParamReason::TooLarge { max: 15 } }),
    );

    let axis = RegularGridInterpolator::new([&[0.0, 1.0][..], &[0.0, 0.0]], &[0.0; 4], GridMethod::Linear);
    assert_eq!(
        axis.err(),
        Some(InterpolateError::GridAxis { axis: 1, fault: AxisFault::DuplicatePoints { first: 0, second: 1 } }),
    );

    // Points on a line, and a kernel whose polynomial has a term in the
    // second coordinate: that coefficient is not determined.
    let line = [[0.0, 3.0], [1.0, 3.0], [2.0, 3.0], [3.0, 3.0]];
    let singular = RBFInterpolator::new(&line, &[0.0, 1.0, 4.0, 9.0], RbfKernel::ThinPlateSpline);
    assert_eq!(singular.err(), Some(InterpolateError::SolverFailed(SolverError::Singular)));

    // The same points, on one line, have no triangles.
    let flat = LinearNDInterpolator::new(&line, &[0.0, 1.0, 4.0, 9.0]);
    assert_eq!(
        flat.err(),
        Some(InterpolateError::DegenerateGeometry(DegeneracyKind::LowerDimensional { affine_dim: 1 })),
    );
}
```

### `Scalar` and `Value`

`Scalar` is `f32` or `f64`. `Value<T>` is `T`, or `[T; M]` for `M` values at each point, interpolated column by column as `M` separate interpolants would be.

```rust
use interpx::{InterpolateError, Interpolator, Linear1D, Scalar, Value};

fn halfway<T: Scalar, V: Value<T>>(x: &[T], y: &[V]) -> Result<Option<V>, InterpolateError> {
    let line = Linear1D::new(x, y)?;
    Ok(line.eval(&((x[0] + x[x.len() - 1]) / (T::one() + T::one()))))
}

fn main() -> Result<(), InterpolateError> {
    assert_eq!(halfway(&[0.0_f32, 2.0], &[1.0_f32, 3.0])?, Some(2.0));
    assert_eq!(halfway(&[0.0, 2.0], &[[1.0, 10.0], [3.0, 30.0]])?, Some([2.0, 20.0]));
    Ok(())
}
```

### `PPoly`

A piecewise polynomial: breakpoints, and the coefficients of each interval, the highest power first. It has `derivative`, `antiderivative`, `eval_derivative`, and `integrate`, as `BSpline` has.

```rust
use interpx::{InterpolateError, Interpolator, PPoly};

fn main() -> Result<(), InterpolateError> {
    // x on [0, 1], then 1 + 2(x - 1) on [1, 2].
    let p = PPoly::new(&[0.0, 1.0, 2.0], &[1.0, 0.0, 2.0, 1.0], 1)?;
    assert_eq!(p.eval(&1.5), Some(2.0));
    assert_eq!(p.eval_derivative(1.5, 1), Some(2.0));
    assert_eq!(p.integrate(0.0, 2.0), Some(2.5));
    assert_eq!(p.derivative(1)?.eval(&0.5), Some(1.0));
    assert_eq!(p.antiderivative(1)?.eval(&1.0), Some(0.5));
    Ok(())
}
```

### `BSpline` and `BSPLINE_MAX_DEGREE`

A B-spline from knots and control points, of degree up to `BSPLINE_MAX_DEGREE` (15). It converts to and from `PPoly`.

```rust
use interpx::{BSpline, InterpolateError, Interpolator, PPoly, BSPLINE_MAX_DEGREE};

fn main() -> Result<(), InterpolateError> {
    // A quadratic on [0, 1]: 2 x (1 - x).
    let b = BSpline::new(&[0.0, 0.0, 0.0, 1.0, 1.0, 1.0], &[0.0, 1.0, 0.0], 2)?;
    assert_eq!(b.eval(&0.5), Some(0.5));
    let area = b.integrate(0.0, 1.0).unwrap_or(f64::NAN);
    assert!((area - 1.0 / 3.0).abs() < 1e-15);

    let p = PPoly::from_bspline(&b)?;
    assert_eq!(p.eval(&0.5), Some(0.5));
    assert!(BSpline::from_ppoly(&p)?.degree() <= BSPLINE_MAX_DEGREE);
    Ok(())
}
```

### `Linear1D`

Linear interpolation. At every data point the value is the data value exactly.

```rust
use interpx::{InterpolateError, Interpolator, Linear1D};

fn main() -> Result<(), InterpolateError> {
    let line = Linear1D::new(&[0.0, 1.0, 3.0], &[0.0, 2.0, 6.0])?;
    assert_eq!(line.eval(&2.0), Some(4.0));
    assert_eq!(line.domain(), (0.0, 3.0));
    Ok(())
}
```

### `Nearest1D` and `NearestKind`

The data value of a neighbor: the nearest (ties to the left or to the right), the previous, or the next.

```rust
use interpx::{InterpolateError, Interpolator, Nearest1D, NearestKind};

fn main() -> Result<(), InterpolateError> {
    let x = [0.0, 1.0, 2.0];
    let y = [10.0, 20.0, 30.0];
    assert_eq!(Nearest1D::new(&x, &y, NearestKind::Nearest)?.eval(&0.5), Some(10.0));
    assert_eq!(Nearest1D::new(&x, &y, NearestKind::NearestUp)?.eval(&0.5), Some(20.0));
    assert_eq!(Nearest1D::new(&x, &y, NearestKind::Previous)?.eval(&0.9), Some(10.0));
    assert_eq!(Nearest1D::new(&x, &y, NearestKind::Next)?.eval(&0.1), Some(20.0));
    Ok(())
}
```

### `CubicSpline`, `BoundaryCondition`, and `EndCondition`

The C² cubic spline. `BoundaryCondition` is `NotAKnot`, `Periodic`, or `Ends` with an `EndCondition` at each end (`NotAKnot`, `FirstDerivative`, `SecondDerivative`); `natural()` is a zero second derivative at both ends, and `clamped()` a zero first derivative.

```rust
use interpx::{BoundaryCondition, CubicSpline, EndCondition, InterpolateError, Interpolator};

fn main() -> Result<(), InterpolateError> {
    let x = [0.0, 1.0, 2.0, 3.0];
    let y = [0.0, 1.0, 0.0, 1.0];
    let natural = CubicSpline::new(&x, &y, BoundaryCondition::natural())?;
    let clamped = CubicSpline::new(&x, &y, BoundaryCondition::clamped())?;
    let slopes = BoundaryCondition::Ends(EndCondition::FirstDerivative(1.0), EndCondition::SecondDerivative(0.0));
    let mixed = CubicSpline::new(&x, &y, slopes)?;
    assert_eq!(natural.eval(&1.0), Some(1.0));
    assert_eq!(clamped.eval(&2.0), Some(0.0));
    assert!(mixed.eval(&1.5).is_some());

    // Periodic: the first and the last values agree.
    let periodic = CubicSpline::new(&x, &[0.0, 1.0, 0.5, 0.0], BoundaryCondition::Periodic)?;
    assert!(periodic.eval(&4.5).is_some());
    Ok(())
}
```

### `PchipInterpolator`

The shape-preserving C¹ cubic (PCHIP): monotone where the data is, and never past the range of the data.

```rust
use interpx::{InterpolateError, Interpolator, PchipInterpolator};

fn main() -> Result<(), InterpolateError> {
    let pchip = PchipInterpolator::new(&[0.0, 1.0, 2.0, 3.0], &[0.0, 0.0, 1.0, 1.0])?;
    let v = pchip.eval(&1.5).unwrap_or(f64::NAN);
    assert!((0.0..=1.0).contains(&v));
    Ok(())
}
```

### `Akima1D` and `AkimaMethod`

Akima's C¹ cubic, with Akima's weights or the modified ones (makima): each slope from the five nearest points.

```rust
use interpx::{Akima1D, AkimaMethod, InterpolateError, Interpolator};

fn main() -> Result<(), InterpolateError> {
    let x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    let y = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    let makima = Akima1D::new(&x, &y, AkimaMethod::Makima)?;
    // A flat stretch stays flat.
    assert_eq!(makima.eval(&0.5), Some(0.0));
    assert!(Akima1D::new(&x, &y, AkimaMethod::Akima)?.eval(&2.5).is_some());
    Ok(())
}
```

### `CubicHermiteSpline`

A cubic through given values with given slopes.

```rust
use interpx::{CubicHermiteSpline, InterpolateError, Interpolator};

fn main() -> Result<(), InterpolateError> {
    // x³ from its values and slopes at -1, 0, 1.
    let cube = CubicHermiteSpline::new(&[-1.0, 0.0, 1.0], &[-1.0, 0.0, 1.0], &[3.0, 0.0, 3.0])?;
    assert_eq!(cube.eval(&0.5), Some(0.125));
    Ok(())
}
```

### `BSpline1D`

An interpolating B-spline of degree 1 to 15. `as_bspline` lends the `BSpline` it is.

```rust
use interpx::{BSpline1D, BoundaryCondition, InterpolateError, Interpolator};

fn main() -> Result<(), InterpolateError> {
    let x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    let y = [0.0, 1.0, 4.0, 9.0, 16.0, 25.0];
    let quintic = BSpline1D::new(&x, &y, 5, BoundaryCondition::NotAKnot)?;
    assert_eq!(quintic.degree(), 5);
    let v = quintic.eval(&2.5).unwrap_or(f64::NAN);
    assert!((v - 6.25).abs() < 1e-12);
    assert_eq!(quintic.as_bspline().degree(), 5);
    Ok(())
}
```

### `BSplineCurve` and `Parametrization`

A parametric curve through points of `D` dimensions. `new` chooses the parameters (`Uniform`, `ChordLength`, or `Centripetal`); `with_parameters` takes them from the caller.

```rust
use interpx::{BSplineCurve, BoundaryCondition, InterpolateError, Interpolator, Parametrization};

fn main() -> Result<(), InterpolateError> {
    let points = [[0.0, 0.0], [1.0, 2.0], [3.0, 3.0], [4.0, 1.0]];
    let curve = BSplineCurve::new(&points, 3, BoundaryCondition::NotAKnot, Parametrization::Centripetal)?;
    let (start, end) = curve.domain();
    assert_eq!(curve.eval(&start), Some(points[0]));
    assert!(curve.eval(&end).is_some());

    let t = [0.0, 1.0, 2.0, 3.0];
    let given = BSplineCurve::with_parameters(&t, &points, 3, BoundaryCondition::NotAKnot)?;
    assert_eq!(given.parameters(), &t);
    Ok(())
}
```

### `RegularGridInterpolator` and `GridMethod`

Values on a rectilinear grid in `D` dimensions: `Nearest`, `Linear`, `Cubic`, `Quintic`, or `Pchip`. Axes may be unevenly spaced, and may be given in descending order.

```rust
use interpx::{Extrapolate, GridMethod, InterpolateError, Interpolator, RegularGridInterpolator};

fn main() -> Result<(), InterpolateError> {
    let x = [0.0, 1.0, 2.0, 3.0];
    let y = [0.0, 2.0, 4.0, 6.0];
    // f(x, y) = x + y, in C order: values[i * 4 + j] = x[i] + y[j].
    let values: Vec<f64> = x.iter().flat_map(|a| y.iter().map(move |b| a + b)).collect();
    let grid = RegularGridInterpolator::new([&x, &y], &values, GridMethod::Linear)?;
    assert_eq!(grid.eval(&[1.5, 3.0]), Some(4.5));
    assert_eq!(grid.shape(), [4, 4]);
    assert_eq!(grid.eval(&[4.0, 0.0]), None);
    assert_eq!(grid.with_extrapolate(Extrapolate::Extend).eval(&[4.0, 0.0]), Some(4.0));

    let cubic = RegularGridInterpolator::new([&x, &y], &values, GridMethod::Cubic)?;
    let v = cubic.eval(&[1.5, 3.0]).unwrap_or(f64::NAN);
    assert!((v - 4.5).abs() < 1e-12);
    Ok(())
}
```

### `NearestNDInterpolator`

The value of the nearest data point, in `D` dimensions. `rescaled` measures distances after each coordinate is brought to the same range; `with_distance_bound` gives no value beyond a distance.

```rust
use interpx::{InterpolateError, Interpolator, NearestNDInterpolator};

fn main() -> Result<(), InterpolateError> {
    let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
    let nearest = NearestNDInterpolator::new(&points, &[1.0, 2.0, 3.0])?;
    assert_eq!(nearest.eval(&[0.9, 0.2]), Some(2.0));

    let bounded = nearest.with_distance_bound(0.5);
    assert_eq!(bounded.eval(&[5.0, 5.0]), None);

    let rescaled = NearestNDInterpolator::rescaled(&points, &[1.0, 2.0, 3.0])?;
    assert_eq!(rescaled.len(), 3);
    Ok(())
}
```

### `LinearNDInterpolator`

Linear interpolation on the Delaunay triangulation of the points, in `D` dimensions. No value outside the convex hull of the points. `rescaled` triangulates after each coordinate is brought to the same range. `with_simplex_tolerance` widens what counts as inside a simplex.

```rust
use interpx::{InterpolateError, Interpolator, LinearNDInterpolator};

fn main() -> Result<(), InterpolateError> {
    // f(x, y) = 1 + 2x + 3y at the corners of the unit square.
    let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]];
    let linear = LinearNDInterpolator::new(&points, &[1.0, 3.0, 4.0, 6.0])?;
    let v = linear.eval(&[0.25, 0.5]).unwrap_or(f64::NAN);
    assert!((v - 3.0).abs() < 1e-12);
    assert_eq!(linear.eval(&[1.5, 0.5]), None);

    let rescaled = LinearNDInterpolator::rescaled(&points, &[1.0, 3.0, 4.0, 6.0])?;
    assert_eq!(rescaled.simplex_count(), 2);
    Ok(())
}
```

### `CloughTocher2DInterpolator`

A piecewise cubic that is continuous with its first derivatives, on the Delaunay triangulation of points of the plane. The gradients at the points are estimated from the data so that the surface bends least along the edges; `gradients()` returns them. No value outside the convex hull. `rescaled`, `with_simplex_tolerance`, and the search for the triangle of a query are those of `LinearNDInterpolator`.

```rust
use interpx::{CloughTocher2DInterpolator, InterpolateError, Interpolator};

fn main() -> Result<(), InterpolateError> {
    // f(x, y) = x^2 + y at a few points of the unit square.
    let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.5, 0.5], [0.3, 0.8]];
    let values: Vec<f64> = points.iter().map(|p| p[0] * p[0] + p[1]).collect();
    let cubic = CloughTocher2DInterpolator::new(&points, &values)?;
    // At a data point, the value given.
    let v = cubic.eval(&[0.5, 0.5]).unwrap_or(f64::NAN);
    assert!((v - 0.75).abs() < 1e-12);
    assert_eq!(cubic.eval(&[1.5, 0.5]), None);
    assert_eq!(cubic.gradients().len(), points.len());
    Ok(())
}
```

### `RBFInterpolator`, `RbfKernel`, `RbfDegree`, and `RbfOptions`

Radial basis function interpolation in `D` dimensions, with eight kernels; the four with a shape parameter carry their `epsilon`. `RbfOptions` sets the smoothing, the degree of the polynomial terms (`RbfDegree`), and the number of neighbors.

```rust
use interpx::{InterpolateError, Interpolator, RBFInterpolator, RbfDegree, RbfKernel, RbfOptions};

fn main() -> Result<(), InterpolateError> {
    let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.5, 0.5]];
    let values = [0.0, 1.0, 1.0, 2.0, 1.2];

    let tps = RBFInterpolator::new(&points, &values, RbfKernel::ThinPlateSpline)?;
    let v = tps.eval(&[1.0, 0.0]).unwrap_or(f64::NAN);
    assert!((v - 1.0).abs() < 1e-12);
    assert!(tps.reciprocal_condition().is_some());

    let options = RbfOptions::new()
        .with_smoothing(0.1)
        .with_degree(RbfDegree::Degree(1))
        .with_neighbors(4);
    let gaussian = RbfKernel::Gaussian { epsilon: 2.0 };
    let local = RBFInterpolator::with_options(&points, &values, gaussian, &options)?;
    assert_eq!(local.neighbors(), Some(4));
    assert!(local.eval(&[0.3, 0.6]).is_some());
    Ok(())
}
```

### `ParInterpolator` (feature `rayon`)

`eval_batch_par` and `eval_batch_into_par`, for every interpolator that is `Sync`. A batch times its first points and shares out the rest only when it is worth it; the results are those of the serial batch, bit for bit.

```rust
use interpx::{InterpolateError, Interpolator, Linear1D, ParInterpolator};

fn main() -> Result<(), InterpolateError> {
    let line = Linear1D::new(&[0.0, 1.0], &[0.0, 10.0])?;
    let queries: Vec<f64> = (0..10_000).map(|i| f64::from(i) / 10_000.0).collect();
    assert_eq!(line.eval_batch_par(&queries), line.eval_batch(&queries));
    let mut out = vec![0.0; queries.len()];
    line.eval_batch_into_par(&queries, &mut out, f64::NAN)?;
    Ok(())
}
```

### Saving and reading back (feature `serde`)

Every interpolator is written by any format of `serde` and read back without being built again: what is stored is what it was built into (the coefficients of a spline, not the data), with a version of the form. Reading runs the checks of the constructor, so a broken form is an error, never an interpolator that panics. What comes back evaluates bit for bit as what was stored when the format gives every number back as it was (binary formats; `serde_json` with its `float_roundtrip` feature).

```rust
use interpx::{BoundaryCondition, CubicSpline, Interpolator};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let spline = CubicSpline::new(&[0.0, 1.0, 2.5, 4.0], &[1.0, 3.0, 2.0, 5.0], BoundaryCondition::NotAKnot)?;
    let text = serde_json::to_string(&spline)?;
    let back: CubicSpline<f64> = serde_json::from_str(&text)?;
    assert_eq!(back.eval(&1.7), spline.eval(&1.7));
    Ok(())
}
```

## Verification

The results are checked against SciPy 1.18: the tests compare every interpolator with values SciPy 1.18 produced (NumPy 2.5's `numpy.interp` for `Linear1D`), committed in the package (`compare/goldens/`), within fixed tolerances. Where SciPy 1.18 computes something else, the results differ as follows.

- **Grid `Cubic` and `Quintic`**: interpx solves the spline exactly, and matches SciPy's `method="cubic_legacy"` and `"quintic_legacy"` to rounding. SciPy's `cubic` and `quintic` use an iterative solver by default, whose values were 1e-5 to 1e-2 of the largest absolute value away from the spline on smooth data (2-D and 3-D grids), and as far as the values themselves on random values.
- **`CloughTocher2DInterpolator`**: interpx solves for the gradients until the residual is below 1e-14 of the right-hand side. SciPy stops at `tol=1e-6` or `maxiter=400`; its values were up to 3e-8 of the largest absolute value away from the solution on the 11 sets compared.
- **`LinearNDInterpolator` and `CloughTocher2DInterpolator`**: where four or more points lie on one circle (five on one sphere, ...), the region they bound is divided by a rule of the points' indices, and SciPy's triangulation (Qhull) divides it otherwise; points on a grid are such points. Inside such a region the linear interpolant differs; the Clough-Tocher interpolant differs everywhere, since its gradients come from one system over all edges. Of points with the same coordinates, interpx takes the value of the one given first, and Qhull keeps one, not always the first. `LinearNDInterpolator` takes one-dimensional points, which SciPy refuses.
- **`RBFInterpolator`**: the system is solved by an LU decomposition, and SciPy uses a symmetric indefinite one; on the test systems the values differ by less than 1e-12 of the largest absolute value. SciPy has a default kernel (the thin plate spline) and takes `epsilon` for every kernel; in interpx the kernel is an argument, and the four kernels with a shape parameter carry `epsilon`. SciPy warns for a polynomial degree below the minimum of the kernel and takes duplicate points with smoothing; interpx refuses both. A query that is not finite has no value (SciPy: NaN), and with neighbors, a query whose neighbors give a singular system has no value (SciPy raises). Of neighbors at the same distance, interpx takes the one given first, which can be another point than SciPy takes.
- **One variable**: SciPy's `interp1d` sorts the coordinates, and `numpy.interp` reads them as given unchecked; interpx refuses coordinates that are not strictly increasing, and needs two points at least (SciPy and NumPy take one). `Nearest1D` gives no value at a NaN query (SciPy: the last data value), and its `Clamp` and `Extend` give the end value on both sides for every kind (SciPy's `"previous"` and `"next"` give NaN on one side). `Akima1D` interpolates each column of array values on its own; SciPy decides where the weights vanish over all columns together.
- **`NearestNDInterpolator`**: a query that is not finite has no value (SciPy raises). Of points with the same coordinates, interpx takes the one given first, and SciPy's choice depends on its tree. `rescaled` sums the mean of each coordinate in another order, so that a query at exactly the distance bound can be decided otherwise.
- **`LinearNDInterpolator`**: a query that rounding puts just beyond a face of the hull is looked for in every simplex, and barycentric coordinates are solved by an LU decomposition, where SciPy multiplies by the inverse: a linear function in a tetrahedron of condition number about 1e4 was off by 1.8e-15 here and 1.8e-12 there.
- **`BSpline1D`**: degree 0 is refused (the step function is `Nearest1D` with `NearestKind::Previous`), and degree 3 takes not-a-knot at one end and a derivative at the other, which SciPy refuses. Degree 1 with `Periodic` extrapolates periodically, where SciPy does not. For an even degree, a knot halfway between two data points that rounds onto a data point is `DuplicatePoints` (SciPy: a singular matrix or NaN), and a halfway point whose sum would overflow is computed without overflow (SciPy: infinite knots).
- **NaN or infinity in the data**: refused at construction. SciPy's `RegularGridInterpolator` takes NaN and returns NaN near it, and `numpy.interp` takes infinite values.
- **`RegularGridInterpolator` outside the grid**: no value under `Forbid`, where SciPy raises by default (`bounds_error=True`).
- **`PPoly` and `BSpline`**: `integrate` with a limit outside the domain under `Forbid` gives `None` (SciPy: NaN for `PPoly`, the part inside the domain for `BSpline`), and a NaN limit gives `None` (SciPy: `PPoly` raises, `BSpline` gives NaN). A `BSpline` derivative of an order above the degree is the zero spline, and a `BSpline` with inner knots repeated `degree + 1` times is differentiated piece by piece; SciPy's `BSpline` raises for both. `BSpline::antiderivative` is zero at the left end of the domain for any knots; SciPy's is so when the first knot is repeated `degree + 1` times. A `PPoly` built from a `BSpline` has no interval of zero width at a repeated knot, where SciPy's keeps one.
- **Grid `Pchip`** needs two points an axis; SciPy asks for four.

## Performance

The comparison with SciPy 1.18 and the other pure-Rust interpolation crates, time and memory, with the machine, the versions, the data, and the command: [docs/performance.md](https://github.com/YUKIKEDA/interpx/blob/v0.1.0/docs/performance.md).

## License

Licensed under either of [Apache License, Version 2.0](LICENSE-APACHE) or [MIT license](LICENSE-MIT) at your option.
