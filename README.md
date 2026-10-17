English | [日本語](README.ja.md)

# interpx

Pure Rust interpolation covering the main classes of SciPy 1.18's `scipy.interpolate`: one variable, curves, regular grids, and scattered data. The values are SciPy's to rounding wherever SciPy computes the same interpolant, and the tests hold them against values SciPy 1.18 produced.

```toml
[dependencies]
interpx = "0.1"
```

Features:

- `scattered` (on by default): `NearestNDInterpolator` and `RBFInterpolator`, and the dependencies only they need.
- `rayon`: batches evaluated on several threads (`ParInterpolator`), with the same results as the serial batches, bit for bit.
- `serde`: `Serialize` and `Deserialize` for every interpolator and the parameters it is built with.

## Coming from SciPy

What differs from SciPy, before the first example:

- **The boundary condition is an argument.** `CubicSpline::new(x, y, bc)` and `BSpline1D::new(x, y, k, bc)` take it, and SciPy's default is `BoundaryCondition::NotAKnot` (also `BoundaryCondition::default()`).
- **There is no `fill_value` at construction.** `eval` returns `Option`: `None` outside the domain when the interpolator does not extrapolate, and for a NaN query. The fill value is given at the call: `eval_or(point, fill)`, or `eval_batch_into(points, out, fill)`.
- **What happens outside the domain is a mode, `Extrapolate`**, set with `with_extrapolate`. The defaults follow SciPy:
  - `Extend` (the end pieces go on): `CubicSpline`, `PchipInterpolator`, `CubicHermiteSpline`, `BSpline1D`, `BSplineCurve`, `PPoly`, and `BSpline`. A spline built with a periodic boundary condition is `Periodic`.
  - `Forbid` (no value): `Akima1D`, `Nearest1D`, and `RegularGridInterpolator`.
  - `Clamp` (the end value): `Linear1D`, as `numpy.interp`.
- **Input is checked, and rejected where SciPy would compute something else.** In one variable the coordinates must be increasing and distinct (SciPy's `interp1d` sorts them), and data must be finite everywhere: NaN or infinity in the data is an error (`InterpolateError::NonFinite`). A NaN of the data never ends up in the values.
- **The grid's `cubic` and `quintic` are solved exactly.** SciPy's default solves them with an iterative solver, whose values are 1e-5 to 1e-3 of the largest value away from the spline. interpx gives the spline itself: SciPy's `method="cubic_legacy"` to rounding.
- **`RBFInterpolator` has no default kernel**, solves its system by an LU decomposition (SciPy: a symmetric one; the values differ by less than 1e-12), and refuses duplicate points.
- **Values are real.** No complex dtypes. A value can be a fixed-size array, `[T; M]`, for `M` columns interpolated at once.
- **Methods are chosen at construction**, not at each evaluation.

## Example

```rust
use interpx::{BoundaryCondition, CubicSpline, Extrapolate, InterpolateError, Interpolator};

fn main() -> Result<(), InterpolateError> {
    let x = [0.0, 1.0, 2.0, 3.5, 5.0];
    let y = [1.0, 2.7, 5.8, 4.0, 2.2];

    // SciPy's `CubicSpline(x, y)`.
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

- **Every interpolator owns its data** and is `Send + Sync + 'static`: it can be kept, moved to another thread, and shared. Constructors copy what they are given. The types that keep some of their data as it is given (`PPoly`, `BSpline`, `RegularGridInterpolator`, `NearestNDInterpolator`, `RBFInterpolator`) also have `from_owned`, which takes that data as a `Vec` and keeps it without a copy.
- **Coordinates are `f32` or `f64`** in one variable, curves, and grids. Scattered data is `f64` only.
- **Grid values are in C order**: the last axis fastest, as NumPy lays out `values[i, j, k]`.
- **NaN and infinity**: in the data they are an error at construction. A NaN query gives no value (`None`, or the fill value). An infinite query is outside the domain: no value under `Forbid` and `Periodic`, the end value under `Clamp`, and under `Extend` the value of the end piece there, which can be infinite or NaN.
- **`RBFInterpolator` over all points holds an `N` by `N` matrix while it is built**: 940 MB for 10 838 points. Over the nearest points (`RbfOptions::with_neighbors`) it solves a small system at each query instead. `reciprocal_condition` says how well the system over all points was solved.
- **Evaluation does not allocate**, except `eval_batch`, which returns a `Vec`, and the first evaluation on a thread of an `RBFInterpolator` over neighbors, which allocates the working space the thread keeps.
- **0.x releases**: a minor version (0.1 to 0.2) may break the public API. The API is what the default features export.
- **Not in this crate**: `interp2d` (removed from SciPy) and `interp1d` as a type (`Linear1D`, `Nearest1D`, and `BSpline1D` cover its kinds). Not yet: derivatives of grid and scattered interpolants, `roots`, `NdPPoly`, `BarycentricInterpolator` and `KroghInterpolator`, the FITPACK family (`UnivariateSpline`, `splrep`, `RectBivariateSpline`), and dimensions or output widths chosen at run time.

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

### Errors: `InterpolateError`, `AxisFault`, `ParamReason`, `SolverError`

Constructors return `InterpolateError`, which says which input and where: `InsufficientPoints`, `LengthMismatch`, `NonFinite`, `NotSorted`, `DuplicatePoints`, `NotPeriodic`, `InvalidParameter` (with a `ParamReason`), `SolverFailed` (with a `SolverError`), `AllocationFailed`, and `GridAxis` (with an `AxisFault`, for one axis of a grid).

```rust
use interpx::{
    AxisFault, BSpline1D, BoundaryCondition, GridMethod, InterpolateError, Linear1D,
    ParamReason, RBFInterpolator, RbfKernel, RegularGridInterpolator, SolverError,
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

A piecewise polynomial (SciPy's `PPoly`): breakpoints, and the coefficients of each interval, the highest power first. It has `derivative`, `antiderivative`, `eval_derivative`, and `integrate`, as `BSpline` has.

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

A B-spline from knots and control points (SciPy's `BSpline`), of degree up to `BSPLINE_MAX_DEGREE` (15). It converts to and from `PPoly`.

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

Linear interpolation (`numpy.interp`, `interp1d(kind="linear")`). At every data point the value is the data value exactly.

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

The data value of a neighbor (`interp1d` with `kind` `"nearest"`, `"nearest-up"`, `"previous"`, or `"next"`).

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

The C² cubic spline (SciPy's `CubicSpline`). `BoundaryCondition` is `NotAKnot`, `Periodic`, or `Ends` with an `EndCondition` at each end (`NotAKnot`, `FirstDerivative`, `SecondDerivative`); `natural()` and `clamped()` are SciPy's `"natural"` and `"clamped"`.

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

The shape-preserving C¹ cubic (SciPy's `PchipInterpolator`): monotone where the data is, and never past the range of the data.

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

Akima's C¹ cubic (SciPy's `Akima1DInterpolator`, `method="akima"` or `"makima"`): each slope from the five nearest points.

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

A cubic through given values with given slopes (SciPy's `CubicHermiteSpline`).

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

An interpolating B-spline of degree 1 to 15 (SciPy's `make_interp_spline`). `as_bspline` lends the `BSpline` it is.

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

A parametric curve through points of `D` dimensions. `new` chooses the parameters (`Uniform`, `ChordLength`, or `Centripetal`); `with_parameters` takes them, as `make_interp_spline(t, points)` does.

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

Values on a rectilinear grid in `D` dimensions (SciPy's `RegularGridInterpolator`): `Nearest`, `Linear`, `Cubic`, `Quintic`, or `Pchip`. Axes may be unevenly spaced, and may be given in descending order.

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

The value of the nearest data point, in `D` dimensions (SciPy's `NearestNDInterpolator`). `rescaled` measures distances after each coordinate is brought to the same range (`rescale=True`); `with_distance_bound` gives no value beyond a distance.

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

### `RBFInterpolator`, `RbfKernel`, `RbfDegree`, and `RbfOptions`

Radial basis function interpolation in `D` dimensions (SciPy's `RBFInterpolator`). The kernels are those of SciPy; the four with a shape parameter carry their `epsilon`. `RbfOptions` sets the smoothing, the degree of the polynomial terms (`RbfDegree`), and the number of neighbors.

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

## Performance

The comparison with SciPy 1.18 and the other pure-Rust interpolation crates, time and memory, with the machine, the versions, the data, and the command: [docs/performance.md](https://github.com/YUKIKEDA/interpx/blob/main/docs/performance.md).

## License

Licensed under either of [Apache License, Version 2.0](LICENSE-APACHE) or [MIT license](LICENSE-MIT) at your option.
