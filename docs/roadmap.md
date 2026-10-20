# Roadmap

Order and status of work. A row holds ID, kind, title, Issue, and status only. Acceptance text lives on the Issue. The spec is `docs/design.ja.md` (milestones: §9.2).

`v0.1.0` is published once, after M1 to M5 are all `Done`. Milestones are internal. Status values: `Todo`, `Set after Grill on #n` (do not start), `In progress`, `Done`.

## Current work

The rows of M5 that do not wait for `convx` are done, and M4 is done: M4-1 (#180) and M4-2 (#185). P-21 (#186) took `convx` 0.1.1. The rows the comparison left open are done too, before the release (the owner, 2026-10-19): P-20 (#184), P-10 (#138), P-11 (#144), and P-19 (#181), which ran the whole comparison again for `docs/performance.md`. M5-9 (#201) checked the publish fields, and M5-12 (#203) ran `publish-review` again and fixed what it found. M5-10 (#207) prepared v0.1.0: the version, `CHANGELOG.md`, the `Release` workflow for the versions after it, and the links at the tag `v0.1.0`. The owner published 0.1.0 on crates.io on 2026-10-19 (UTC), by hand, and #213 dated it in the changelog. Every later version goes through the `Release` workflow. V-1 (#215) decided where the documentation names SciPy (ADR 0011). V-2 (#217) brought the README, the rustdoc, the description, and the keywords to it; they reach crates.io with the next published version. Next a row of After v0.1.0, which the owner chooses.

## Dependencies

- M4 depends on `convx` from crates.io (published with `MIT OR Apache-2.0`); M4-1 added the dependency at 0.1.0, P-21 takes 0.1.1.
- M1 to M3 do not depend on `convx`.

## M0: Bootstrap

| ID   | Kind | Title                                              | Issue | Status |
| ---- | ---- | -------------------------------------------------- | ----- | ------ |
| M0-1 | Task | Repository rules, docs, ADRs, CI, and build config | #1    | Done   |

## M1: Core, 1D, and curve

| ID    | Kind | Title                                                                                                                                            | Issue    | Status |
| ----- | ---- | ------------------------------------------------------------------------------------------------------------------------------------------------ | -------- | ------ |
| M1-1  | Task | Golden pipeline and data sets: pinned SciPy 1.18, layer-A goldens for every family, `tests/common`                                               | #2       | Done   |
| M1-2  | Feat | `Scalar`, `Value`, `Interpolator`, `InterpolateError`, `Extrapolate`, `BoundaryCondition`, `tests/alloc.rs` harness, restore the CI coverage job | #11      | Done   |
| M1-3  | Feat | `PPoly` (interval-major layout, Horner evaluation, `Hint`) and the shared input checks in their fixed order                                      | #14      | Done   |
| M1-4  | Feat | `BSpline` (de Boor) and `PPoly` ⇄ `BSpline` conversion                                                                                           | #18      | Done   |
| M1-5  | Feat | `derivative`, `antiderivative`, `integrate`, `eval_derivative` on `PPoly` and `BSpline`                                                          | #23      | Done   |
| M1-6  | Feat | Tridiagonal solver (Thomas, periodic via Sherman-Morrison)                                                                                       | #28      | Done   |
| M1-7  | Feat | `Linear1D` and `Nearest1D` (nearest, nearest-up, previous, next)                                                                                 | #35      | Done   |
| M1-8  | Feat | `CubicSpline` (NotAKnot, Periodic, first or second derivative at each end)                                                                       | #40      | Done   |
| M1-9  | Feat | `PchipInterpolator`                                                                                                                              | #45      | Done   |
| M1-10 | Feat | `Akima1D` (including makima)                                                                                                                     | #50      | Done   |
| M1-11 | Feat | `CubicHermiteSpline`                                                                                                                             | #55      | Done   |
| M1-12 | Feat | `BSpline1D` (degrees 1 to 15), and the banded LU with partial pivoting it solves with                                                            | #64, #60 | Done   |
| M1-13 | Feat | `BSplineCurve<T, D>`                                                                                                                             | #67      | Done   |
| M1-14 | Task | Golden data budget: 8 MB for the directory, the goldens stay in the package (ADR 0010)                                                           | #70      | Done   |

## M2: Regular grid

| ID   | Kind | Title                                                                                            | Issue | Status |
| ---- | ---- | ------------------------------------------------------------------------------------------------ | ----- | ------ |
| M2-1 | Feat | `RegularGridInterpolator`: shape, C-order values, `nearest`, `linear`, `Extrapolate`             | #74   | Done   |
| M2-2 | Feat | `RegularGridInterpolator`: `cubic`, `quintic` (tensor-product coefficients at construction)      | #79   | Done   |
| M2-3 | Feat | `RegularGridInterpolator`: `pchip` (last-axis slopes kept, the other axes folded at evaluation)  | #87   | Done   |
| M2-4 | Task | Grid `eval_batch_into` on a path: close the gap to a loop over `eval_hinted` (does not block M3) | #90   | Done   |

## M3: Scattered, without `convx`

| ID   | Kind | Title                                                                              | Issue | Status |
| ---- | ---- | ---------------------------------------------------------------------------------- | ----- | ------ |
| M3-1 | Feat | `NearestNDInterpolator<D>` (`kiddo`)                                               | #95   | Done   |
| M3-2 | Feat | `RBFInterpolator<D>`: 8 kernels, polynomial terms, shift and scale                 | #99   | Done   |
| M3-3 | Feat | `RBFInterpolator<D>`: global solve (`faer` partial-pivot LU, singular is an error) | #103  | Done   |
| M3-4 | Feat | `RBFInterpolator<D>`: local `neighbors` mode                                       | #107  | Done   |
| M3-5 | Task | RBF neighbors: one factorization for queries of a batch with the same neighbors    | #108  | Done   |

## M4: Scattered, with `convx`

| ID   | Kind | Title                                                                                         | Issue | Status |
| ---- | ---- | --------------------------------------------------------------------------------------------- | ----- | ------ |
| M4-1 | Feat | `LinearNDInterpolator<D>`: Delaunay, a walk from cells, tolerances, brute-force fallback, `serde` | #180  | Done   |
| M4-2 | Feat | `CloughTocher2DInterpolator`: global gradient estimate, Bernstein-Bézier patch, `serde`       | #185  | Done   |

## P: What the comparison found (before M5-1)

The comparison (`docs/performance.md`) is run again for the cases of a row when the row is done. A row that M4 adds to the comparison and that loses is added here.

| ID   | Kind   | Title                                                                                                                                                              | Issue | Status |
| ---- | ------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ----- | ------ |
| P-1  | Design | Construction: what it copies and what it checks (grid 180x and more, 1-D 1.3x to 5x)                                                                               | #122  | Done   |
| P-2  | Task   | 1-D evaluation of queries in no order on 10 000 points (2.4x against `ndarray-interp`)                                                                             | #123  | Done   |
| P-3  | Task   | Grid `linear` and `nearest` evaluation (1.1x to 1.6x against `interpn`)                                                                                            | #124  | Done   |
| P-4  | Design | Grid: an axis with one point (SciPy builds it; eight cases of SciPy's suite do not run)                                                                            | #125  | Done   |
| P-5  | Task   | RBF `quintic` at 10 points (1.12x against SciPy, one case)                                                                                                         | #126  | Done   |
| P-6  | Feat   | `from_owned` for `NearestNDInterpolator` and `RBFInterpolator` (the types that keep data)                                                                          | #128  | Done   |
| P-7  | Task   | `BSpline1D` of degree 1 on 10 000 points: evaluation in order (1.14x against SciPy)                                                                                | #130  | Done   |
| P-8  | Task   | Construction on 10 000 points: `Linear1D` 1.9x, `Nearest1D` 1.1x, `BSpline1D` degree 1 2.9x                                                                        | #131  | Done   |
| P-9  | Task   | Search: the time for points in no order moves twofold with where it is built in                                                                                    | #132  | Done   |
| P-10 | Task   | Search: the cause of P-9 in the disassembly, the five cases it made slower, a two-place check                                                                      | #138  | Done   |
| P-11 | Task   | Grid: the search on two and three axes against `interpn` (up to 1.5x, part of it unexplained)                                                                      | #144  | Done   |
| P-12 | Bug    | RBF, Gaussian kernel: on an ill-conditioned system the error is three times SciPy's, outside SciPy's own range (it was inside: one order of the points among many) | #159  | Done   |
| P-13 | Task   | Grid splines: a direct against an iterative solve, inside SciPy and against interpx; whether the solver is the user's to choose                                    | #162  | Done   |
| P-14 | Bug    | RBF on real station data: the cubic and quintic kernels give interpolants off by tens to thousands of degrees, far from SciPy's (stations that nearly coincide; the interpolant now says how badly its system is conditioned) | #164  | Done   |
| P-15 | Bug    | compare: eight runs reported as failures of other libraries with no reason; all eight were our runners' own panic | #165  | Done   |
| P-16 | Bug    | compare: SciPy's memory in `docs/performance.md` is that of a 4 MB launcher, not of the interpreter | #166  | Done   |
| P-17 | Task   | compare: a construction of minutes is built twice, and one run sits at the time limit | #168  | Done   |
| P-18 | Task   | Some constructions of interpx take one of two times depending on the process (`CubicSpline` of 1000 points: 15 or 38 µs), so whether the comparison counts them as slower depends on the run | #172  | Done   |
| P-19 | Task   | After the warm-up of P-18, some constructions at 10 000 points take 1.3 to 1.6 times as long in some rounds (`Akima1D`, `Pchip`: medians 83 to 138 µs, fastest 81 to 91 µs) | #181  | Done   |
| P-20 | Task   | `LinearNDInterpolator`: what each evaluation always pays, profiled (at the pixels of an image it is level with SciPy, 1.01 times; the walk is 1.8 steps) | #184  | Done   |
| P-21 | Task   | convx 0.1.1 (its 2-D Delaunay build level with `spade`): take it, and measure again what `LinearND` and `CloughTocher2D` pay on top of it; before M5-10 | #186  | Done   |

## M5: Integration and quality

| ID    | Kind | Title                                                                                                                                                           | Issue | Status |
| ----- | ---- | --------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----- | ------ |
| M5-1  | Feat | `rayon` features: `eval_batch_par`, `eval_batch_into_par`                                                                                                       | #150  | Done   |
| M5-3  | Feat | `serde` for all interpolators (topology stored, format version, validation)                                                                                     | #154  | Done   |
| M5-4  | Test | Property tests (`proptest`) and the SciPy 1.18 golden suite for every type                                                                                      | #157  | Done   |
| M5-5  | Task | `compare/perf`: SciPy and other pure-Rust crates, time and memory                                                                                               | #120  | Done   |
| M5-6  | Task | `README.md`, `README.ja.md`, and `examples/`                                                                                                                    | #174  | Done   |
| M5-7  | Task | English `docs/design.md`                                                                                                                                        | #176  | Done   |
| M5-8  | Task | Run `publish-review`, fix rustdoc and README failures                                                                                                           | #178  | Done   |
| M5-9  | Task | Check the publish fields of `Cargo.toml` (the `convx` dependency is added from crates.io by M4-1; the license is set)                                          | #201  | Done   |
| M5-12 | Task | Run `publish-review` again on the release candidate, fix what fails (before M5-10)                                                                         | #203  | Done   |
| M5-10 | Task | Release PR for v0.1.0                                                                                                                                           | #207  | Done   |
| M5-11 | Task | `compare/perf`: the real-world data sets of design section 8.4.7 (left out of M5-5) | #160  | Done   |

## V: The documentation after 0.1.0

Where the documentation names SciPy (ADR 0011). The pages reach crates.io with the next published version.

| ID  | Kind   | Title                                                                                                 | Issue | Status      |
| --- | ------ | ----------------------------------------------------------------------------------------------------- | ----- | ----------- |
| V-1 | Design | Where the documentation names SciPy: ADR 0011, the rules, the first line of the design                 | #215  | Done        |
| V-2 | Docs   | The README, the rustdoc, the description, and the keywords, as ADR 0011 says                          | #217  | Done        |

## After v0.1.0

Named so that no deferral is anonymous. Each becomes a roadmap row with an Issue when its turn comes.

| ID   | Kind | Title                                                                                                              | Issue | Status |
| ---- | ---- | ------------------------------------------------------------------------------------------------------------------ | ----- | ------ |
| A-1  | Feat | `roots` / `solve` for `PPoly` and `BSpline`                                                                        |       | Todo   |
| A-2  | Feat | `NdPPoly`                                                                                                          |       | Todo   |
| A-3  | Feat | Derivatives for grid and scattered interpolators                                                                   |       | Todo   |
| A-4  | Feat | Runtime dimension and runtime output width for grid and scattered                                                  |       | Todo   |
| A-5  | Feat | RBF pseudo-inverse fallback                                                                                        |       | Todo   |
| A-6  | Feat | `BarycentricInterpolator` and `KroghInterpolator`                                                                  |       | Todo   |
| A-7  | Feat | FITPACK family: `UnivariateSpline`, `splrep`, `splprep`, `RectBivariateSpline`, and related                        |       | Todo   |
| A-8  | Task | `no_std` (`alloc`) feasibility                                                                                     |       | Todo   |
| A-9  | Feat | `BSpline::from_ppoly`: lower the knot multiplicity where the pieces join smoothly                                  |       | Todo   |
| A-10 | Feat | `BSpline1D`: several derivatives at an end (`bc_type` as lists of pairs), and knots given by the caller            |       | Todo   |
| A-11 | Feat | `RBFInterpolator`: a smoothing value for each point                                                                |       | Todo   |
| A-12 | Feat | `ndarray` interoperability: constructors and batches that take `ArrayView` (was M5-2; needs data that is borrowed) |       | Todo   |
| A-13 | Feat | Grid: a cubic that solves nothing (centered-difference Hermite, as `interpn` has), to choose beside the spline     | #163  | Todo   |

## Intentionally out of scope

- `interp2d`: removed from SciPy.
- `interp1d` as a type: replaced by `Nearest1D` (nearest, nearest-up, previous, next) and `BSpline1D` (`slinear`, `quadratic`).
- `griddata` as a function: it builds `NearestNDInterpolator`, `LinearNDInterpolator`, or `CloughTocher2DInterpolator` by its `method`, evaluates, and drops it. The types do the same, and can be evaluated again (the owner's decision, #179).
