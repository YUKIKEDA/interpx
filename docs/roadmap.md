# Roadmap

Order and status of work. A row holds ID, kind, title, Issue, and status only. Acceptance text lives on the Issue. The spec is `docs/design.ja.md` (milestones: §9.2).

`v0.1.0` is published once, after M1 to M5 are all `Done`. Milestones are internal. Status values: `Todo`, `Set after Grill on #n` (do not start), `In progress`, `Done`.

## Current work

M1-7.

## Dependencies

- M4 starts after `convx` is published with a license (M5-9 tracks the `Cargo.toml` dependency switch).
- M1 to M3 do not depend on `convx`.

## M0: Bootstrap

| ID   | Kind | Title                                              | Issue | Status      |
| ---- | ---- | -------------------------------------------------- | ----- | ----------- |
| M0-1 | Task | Repository rules, docs, ADRs, CI, and build config | #1    | Done        |

## M1: Core, 1D, and curve

| ID    | Kind | Title                                                                                                                                            | Issue | Status |
| ----- | ---- | ------------------------------------------------------------------------------------------------------------------------------------------------ | ----- | ------ |
| M1-1  | Task | Golden pipeline and data sets: pinned SciPy 1.18, layer-A goldens for every family, `tests/common`                                               | #2    | Done   |
| M1-2  | Feat | `Scalar`, `Value`, `Interpolator`, `InterpolateError`, `Extrapolate`, `BoundaryCondition`, `tests/alloc.rs` harness, restore the CI coverage job | #11   | Done   |
| M1-3  | Feat | `PPoly` (interval-major layout, Horner evaluation, `Hint`) and the shared input checks in their fixed order                                      | #14   | Done   |
| M1-4  | Feat | `BSpline` (de Boor) and `PPoly` ⇄ `BSpline` conversion                                                                                           | #18   | Done   |
| M1-5  | Feat | `derivative`, `antiderivative`, `integrate`, `eval_derivative` on `PPoly` and `BSpline`                                                          | #23   | Done   |
| M1-6  | Feat | Tridiagonal solver (Thomas, periodic via Sherman-Morrison)                                                                                       | #28   | Done   |
| M1-7  | Feat | `Linear1D` and `Nearest1D` (nearest, previous, next)                                                                                             |       | Todo   |
| M1-8  | Feat | `CubicSpline` (NotAKnot, Periodic, first or second derivative at each end)                                                                       |       | Todo   |
| M1-9  | Feat | `PchipInterpolator`                                                                                                                              |       | Todo   |
| M1-10 | Feat | `Akima1D` (including makima)                                                                                                                     |       | Todo   |
| M1-11 | Feat | `CubicHermiteSpline`                                                                                                                             |       | Todo   |
| M1-12 | Feat | `BSpline1D` (any degree `k`), and the banded LU with partial pivoting it solves with                                                             |       | Todo   |
| M1-13 | Feat | `BSplineCurve<T, D>`                                                                                                                             |       | Todo   |

## M2: Regular grid

| ID   | Kind | Title                                                                                                  | Issue | Status |
| ---- | ---- | ------------------------------------------------------------------------------------------------------ | ----- | ------ |
| M2-1 | Feat | `RegularGridInterpolator`: shape, C-order values, `nearest`, `linear`, `Extrapolate`                   |       | Todo   |
| M2-2 | Feat | `RegularGridInterpolator`: `slinear`, `cubic`, `quintic` (tensor-product coefficients at construction) |       | Todo   |
| M2-3 | Feat | `RegularGridInterpolator`: `pchip` (lazy per-axis convolution)                                         |       | Todo   |

## M3: Scattered, without `convx`

| ID   | Kind | Title                                                                              | Issue | Status |
| ---- | ---- | ---------------------------------------------------------------------------------- | ----- | ------ |
| M3-1 | Feat | `NearestNDInterpolator<D>` (`kiddo`)                                               |       | Todo   |
| M3-2 | Feat | `RBFInterpolator<D>`: 8 kernels, polynomial terms, shift and scale                 |       | Todo   |
| M3-3 | Feat | `RBFInterpolator<D>`: global solve (`faer` partial-pivot LU, singular is an error) |       | Todo   |
| M3-4 | Feat | `RBFInterpolator<D>`: local `neighbors` mode                                       |       | Todo   |

## M4: Scattered, with `convx`

| ID   | Kind | Title                                                                                | Issue | Status |
| ---- | ---- | ------------------------------------------------------------------------------------ | ----- | ------ |
| M4-1 | Feat | `LinearNDInterpolator<D>`: Delaunay, Jump-and-Walk, tolerances, brute-force fallback |       | Todo   |
| M4-2 | Feat | `CloughTocher2DInterpolator`: global gradient estimate, Bernstein-Bézier patch       |       | Todo   |

## M5: Integration and quality

| ID    | Kind | Title                                                                       | Issue | Status |
| ----- | ---- | --------------------------------------------------------------------------- | ----- | ------ |
| M5-1  | Feat | `rayon` features: `eval_batch_par`, `eval_batch_into_par`                   |       | Todo   |
| M5-2  | Feat | `ndarray` interoperability                                                  |       | Todo   |
| M5-3  | Feat | `serde` for all interpolators (topology stored, format version, validation) |       | Todo   |
| M5-4  | Test | Property tests (`proptest`) and the SciPy 1.18 golden suite for every type  |       | Todo   |
| M5-5  | Task | `compare/perf`: SciPy and other pure-Rust crates, time and memory           |       | Todo   |
| M5-6  | Task | `README.md`, `README.ja.md`, and `examples/`                                |       | Todo   |
| M5-7  | Task | English `docs/design.md`                                                    |       | Todo   |
| M5-8  | Task | Run `publish-review`, fix rustdoc and README failures                       |       | Todo   |
| M5-9  | Task | Switch `convx` to the published version; set license and publish fields     |       | Todo   |
| M5-10 | Task | Release PR for v0.1.0                                                       |       | Todo   |

## After v0.1.0

Named so that no deferral is anonymous. Each becomes a roadmap row with an Issue when its turn comes.

| ID  | Kind | Title                                                                                       | Issue | Status |
| --- | ---- | ------------------------------------------------------------------------------------------- | ----- | ------ |
| A-1 | Feat | `roots` / `solve` for `PPoly` and `BSpline`                                                 |       | Todo   |
| A-2 | Feat | `NdPPoly`                                                                                   |       | Todo   |
| A-3 | Feat | Derivatives for grid and scattered interpolators                                            |       | Todo   |
| A-4 | Feat | Runtime dimension and runtime output width for grid and scattered                           |       | Todo   |
| A-5 | Feat | RBF pseudo-inverse fallback                                                                 |       | Todo   |
| A-6 | Feat | `BarycentricInterpolator` and `KroghInterpolator`                                           |       | Todo   |
| A-7 | Feat | FITPACK family: `UnivariateSpline`, `splrep`, `splprep`, `RectBivariateSpline`, and related |       | Todo   |
| A-8 | Task | `no_std` (`alloc`) feasibility                                                              |       | Todo   |
| A-9 | Feat | `BSpline::from_ppoly`: lower the knot multiplicity where the pieces join smoothly           |       | Todo   |

## Intentionally out of scope

- `interp2d`: removed from SciPy.
- `interp1d` as a type: replaced by `Nearest1D` (nearest, previous, next) and `BSpline1D` (`slinear`, `quadratic`).
