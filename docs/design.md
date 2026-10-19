English | [日本語](design.ja.md)

# `interpx` Architecture & Design Specification

> **An interpolation library in pure Rust that covers the main interpolation classes of `scipy.interpolate`.**

This file is the English translation of [design.ja.md](design.ja.md), and both carry the whole specification. When the two differ, follow `design.ja.md`.

The SciPy version of reference is fixed at **1.18** (the numerical agreement tests and the coverage table are all defined against this version).

---

## 1. Vision & Design Principles

### 1.1 Background and purpose
The largest problem of Rust's scientific computing ecosystem is fragmentation. Interpolation in one dimension, on grids, and of scattered points is split over separate crates, with APIs and type systems that do not agree. Multidimensional interpolation of irregular scattered data in particular (`LinearNDInterpolator`, `griddata`) is practically a blank.

`interpx` aims to be **the all-in-one de facto standard covering the main interpolation classes of `scipy.interpolate` (the coverage table of §5)**. What is outside the coverage is listed in §5.3 and handled by the roadmap (§9).

### 1.2 Core design principles
1. **Pure Rust (no C or Fortran dependency)**: it builds anywhere with `cargo build` alone (cross compilation included).
2. **Static type safety**: dimensions optimized at compile time with `const D: usize`, and stack allocation, come first.
3. **Rust idioms and high throughput together**: safe out-of-domain results with `Option`, reuse of memory (zero allocation), and parallel batch evaluation with `rayon`.
4. **Modular**: the heavy dependencies (`convx`, `kiddo`, `faer`) are kept behind the `scattered` feature, so that 1D / curve / grid build light.
5. **`std` required**: v0.1.0 assumes `std`. `no_std` is out of scope (roadmap).

---

## 2. Technology Stack

| Role                  | Crate            | Use and reason for the choice                                                                                  |
| :-------------------- | :--------------- | :------------------------------------------------------------------------------------------------------------- |
| **Geometry engine**   | **`convx`**      | Pure Rust $N$-dimensional Delaunay simplicial decomposition and neighbor topology. `scattered` only.           |
| **Linear algebra**    | **`faer`**       | Only the dense LU factorization of RBF. `default-features = false, features = ["std"]`, the version the same as `convx`'s (now 0.24). |
| **Spatial neighbors** | **`kiddo`**      | k-d tree (const generic dimension). `NearestND`, the local neighbors of RBF. `scattered` only.           |
| **Scalar abstraction**| **`num-traits`** | The `Float` generics of 1D / curve / grid (f32 / f64).                                                         |
| **Parallelism (opt)** | **`rayon`**      | Parallel batch evaluation.                                                                                     |
| **Serialization (opt)** | **`serde`**    | Saving and restoring interpolators (§6.5).                                                                     |
| **Errors**            | **`thiserror`**  | The error types.                                                                                               |

* Tridiagonal systems (CubicSpline and others) are solved by **an O(n) Thomas algorithm of the crate's own** (periodic boundaries as a cyclic tridiagonal system by Sherman–Morrison). Banded systems (`BSpline1D`, the grid's `cubic` / `quintic`) are solved by a banded LU with partial pivoting of the crate's own (§6.1). Neither uses `faer`.
* `convx` is a dependency on 0.1.1 from crates.io (`MIT OR Apache-2.0`) (M4-1 added it at 0.1.0; P-21 took 0.1.1, whose 2-D Delaunay build is faster).
* The MSRV is **1.89**, that of `convx` and `kiddo`.

### 2.1 Development dependencies (dev-dependencies)

For tests and benchmarks only; they do not appear among the dependencies of the published package. Nothing other than what is named here is added (an addition is a design change first. `.cursor/rules/layout.mdc`).

| Crate            | Use                                                                                        | Row that first needs it              |
| :--------------- | :----------------------------------------------------------------------------------------- | :----------------------------------- |
| **`serde_json`** | Reading the golden JSON (`compare/goldens/*.json`). Used by the loader of `tests/common`.  | M1-1                                 |
| **`serde`**      | `derive`, to read the JSON into typed structs in the loader above (only when `serde` is needed as a dev dependency). | M1-1  |
| **`proptest`**   | Property tests (§8.1).                                                                     | The row that writes the first property test (M1-3 or later). |
| **`criterion`**  | Timing in `benches/` (§8.3). Not run in CI.                                                | The row that adds the first benchmark |
| **`stats_alloc`** | The counting allocator of `tests/alloc.rs` (§8.1; also used for the memory measurements of §8.3). | M1-2                          |

* The counting allocator of the zero-allocation check (`tests/alloc.rs`, §8.1) **uses the external crate `stats_alloc`**. Implementing `GlobalAlloc` needs an `unsafe impl`, and `unsafe_code = "forbid"` applies to `interpx`'s own code (tests included), so the crate cannot implement it itself. The `unsafe` lives in the dependency (`gprx` uses the same crate).
* The version and features of each dependency are decided in the row that adds it.

---

## 3. Type System & APIs

### 3.1 Scalar type, value type, and number of dimensions
* **Scalar `Scalar`**: a sealed trait, implemented for `f32` and `f64` only. It includes `num_traits::Float + Debug + Send + Sync`. The bound of an interpolator's type parameter `T` is this one trait.
* **1D / curve / grid**: generic over `T: Scalar` (`f32` / `f64`).
* **scattered**: **`f64` only** (`convx` is f64 only). The caller converts f32 input. There is no implicit cast.
* **Value type `Value`**: the type of the value interpolated (one row of SciPy's `values`). A sealed trait, implemented for a scalar `T` and for a fixed-size array `[T; M]`. The operations it needs are zero, a scalar multiple with an addition, and a test of finiteness, nothing else.

  ```rust
  pub trait Value<T: Scalar>: Copy + Debug + Send + Sync {
      fn zero() -> Self;
      /// self + k * x
      fn axpy(self, k: T, x: Self) -> Self;
      fn all_finite(&self) -> bool;
  }
  ```

  * **Column-wise operations inside the crate**: the public operations are only the three above (what a linear combination needs). Some interpolators need nonlinear operations column by column (the slopes of PCHIP and Akima are decided by comparisons of signs and harmonic means per column; the periodicity check of `CubicSpline` looks at the size of differences per column). These are given to the crate-private trait used for sealing, and are not public. It has "combine several values column by column with a function" (for two values and for any number; an Akima slope is decided from five: four secants and a threshold) and "a column-wise test of closeness". Users do not add value types, so the public `Value` need not grow.
  * Several outputs (SciPy's `values` of shape `(N, m)`) are a value type `[T; M]`. For example: `LinearNDInterpolator<const D: usize, V = f64>`, `RegularGridInterpolator<T, const D: usize, V = T>`, `BSplineCurve<T, D>` (whose value type is `[T; D]`).
  * `[T; M]` is contiguous, so it agrees with the coefficient layout of `PPoly`, `[interval][power][output column]` (§6.1).
  * A number of output columns `m` decided at run time is out of scope for v0.1.0 (roadmap A-4; the same as a dynamic dimension).
* The dimension is **`const D: usize` only** (grid and scattered alike). A wrapper with a dynamic dimension is out of scope for v0.1.0 (roadmap). `kiddo` is meant for low dimensions (fewer than about 10), and the doc states the recommended dimensions.

### 3.2 The common trait `Interpolator`
```rust
pub trait Interpolator {
    /// The type of a query point. T in 1D, [T; D] in N-D.
    type Point;
    /// The type of the output. T, or [T; M].
    type Output;

    /// Evaluates at one point. None outside the domain (`Extrapolate::Forbid`, outside the convex hull) and for a NaN query.
    fn eval(&self, point: &Self::Point) -> Option<Self::Output>;

    /// Returns `fill` outside the domain and for a NaN query.
    fn eval_or(&self, point: &Self::Point, fill: Self::Output) -> Self::Output {
        self.eval(point).unwrap_or(fill)
    }

    /// Returns the results in a new Vec (allocates).
    fn eval_batch(&self, points: &[Self::Point]) -> Vec<Option<Self::Output>> {
        points.iter().map(|p| self.eval(p)).collect()
    }

    /// Writes the results to `out` (does not allocate). Outside the domain and for a NaN query, `fill`.
    /// `points.len() != out.len()` is Err(LengthMismatch), and `out` is not touched.
    fn eval_batch_into(
        &self,
        points: &[Self::Point],
        out: &mut [Self::Output],
        fill: Self::Output,
    ) -> Result<(), InterpolateError>
    where
        Self::Output: Clone;
}
```
* The only required method is `eval`; the others have default implementations. Each interpolator overrides `eval_batch_into` where it can be made faster (reusing the `Hint` of the interval search within a batch, and so on).
* `Point` is `Sized` (v0.1.0 has `const D` only). The trait is not sealed (users can implement it for interpolators of their own).
* There is no `fill_value` constructor argument as in SciPy. It is always the `fill` of the call (to fill with NaN, pass `T::nan()`).
* `eval`, `eval_or`, and `eval_batch_into` do not allocate (the zero-allocation check of §8.1). `eval_batch` returns a `Vec` and so allocates.
  * There is one exception. The neighbors mode of `RBFInterpolator` (§6.3) solves a small linear system at every evaluation. Its working space is kept per thread and allocated at the thread's first evaluation (widened when the number of neighbors grows). Evaluations after that do not allocate.

**Parallel forms** (the `rayon` feature. M5-1, #148): `eval_batch_par` / `eval_batch_into_par` are provided as methods of a trait separate from `Interpolator`.

```rust
/// Exists only with the `rayon` feature.
pub trait ParInterpolator: Interpolator + Sync
where
    Self::Point: Sync,
    Self::Output: Send + Sync + Clone,
{
    fn eval_batch_par(&self, points: &[Self::Point]) -> Vec<Option<Self::Output>>;

    fn eval_batch_into_par(
        &self,
        points: &[Self::Point],
        out: &mut [Self::Output],
        fill: Self::Output,
    ) -> Result<(), InterpolateError>;
}

impl<I> ParInterpolator for I where I: Interpolator + Sync + ?Sized, /* the conditions above */ {}
```

* **Where it lives**: adding feature-gated methods to `Interpolator` is not taken (the shape of the trait would change with a feature, and the `Sync` conditions would be written per method). The extension trait is implemented once for every `Interpolator + Sync`, so interpolators that users implement get the parallel forms too. The methods cannot be overridden (there is one implementation).
* **The unit of parallelism is the serial batch**: the points are split into chunks, and each chunk calls the serial `eval_batch_into` (`eval_batch` for `eval_batch_par`). What each type does for its batches (hints, reusing the factorization of the same set of neighbors) works as it is. The result of a serial batch is bit-identical to repeated `eval`, so the parallel result is bit-identical to the serial one whatever the chunks and the number of threads. Calling `eval` per point is not taken (the batch work would be lost).
* **Whether to split, and the length of a chunk, are decided by timing the first few points.** A prototype was compared with chunks of fixed lengths (16 threads; speed over the serial `eval_batch_into`; under 1 is slower):

  | Evaluation (serial, per point) | Points | Chunk 256 | Chunk 1024 | Chunk 4096 | Timed |
  | :--- | ---: | ---: | ---: | ---: | ---: |
  | `Linear1D`, in no order (5 ns) | $10^3$ | 0.40 | 0.99 | 1.00 | 0.97 |
  | | $10^4$ | 2.51 | 2.12 | 1.23 | 1.04 |
  | | $10^5$ | 5.61 | 5.91 | 5.08 | 4.40 |
  | | $10^6$ | 6.55 | 6.67 | 6.55 | 6.52 |
  | Grid, 3-D (25 ns) | $10^3$ | 0.83 | 0.91 | 0.99 | 0.99 |
  | | $10^4$ | 4.78 | 2.62 | 1.93 | 2.92 |
  | | $10^5$ | 6.56 | 6.14 | 5.29 | 6.08 |
  | `NearestND` (150 ns) | $10^3$ | 1.24 | 0.92 | 0.92 | 1.55 |
  | | $10^4$ | 6.18 | 3.90 | 2.37 | 7.33 |
  | | $10^5$ | 9.02 | 8.10 | 8.09 | 9.04 |
  | RBF over 20 neighbors (5.5 µs) | $10^3$ | 3.42 | 0.99 | 0.98 | 4.39 |
  | | $10^4$ | 7.44 | 5.70 | 2.26 | 7.92 |

  * Every fixed length loses on one side: 256 points falls to 0.4 times on $10^3$ points of a cheap evaluation (the cost of waking threads exceeds the 5 µs of work in all), and 1024 points or more gains nothing on $10^3$ points of an expensive evaluation (one chunk). Whether splitting is worth it is decided by "the time of the work", not by "the number of points".
  * **Procedure**: evaluate the first few points by a serial batch and time them (the results go to the output as they are, so nothing is wasted). If the time measured is too short against the step of the clock, time again with more points. If the estimate of the rest (time per point × points left) does not reach "the time worth waking threads", evaluate the rest serially as well. If it does, split the rest into chunks of the length that makes "the time of the work of one chunk" and hand them to `rayon`.
  * **The values decided** (the implementation of M5-1, #150. `benches/par.rs`, 16 threads): 16 points are timed first. If that takes less than 2 µs (20 times the clock step of 100 ns), time again with 16 times as many points (256, 4096, …). But if the points left are no more than the number to be timed next, the rest is evaluated serially without timing (the timing is done once for 100 points and twice for 1000). The time worth waking threads is 100 µs, the time of the work of one chunk is 10 µs, and a chunk is a multiple of 16 points of at least 16.
    * With a worth of 30 µs, 1000 points of a grid (26 µs in all) fell to 0.71 times. The first points of a batch are slower than the later ones, with nothing in the cache. The rest of the work is overestimated, and split. With 100 µs it is not split.
    * With 64 points timed first, 100 points of the RBF over 20 neighbors (0.5 ms) stayed at 1.2 times (the first 64 points are evaluated on one thread). With 16 points it is 2.0 times.
    * The result with these values (speed over the serial `eval_batch_into`):

      | Points | `Linear1D`, in no order | Grid, 3-D | RBF over 20 neighbors |
      | ---: | ---: | ---: | ---: |
      | $10^2$ | 0.88 | 0.98 | 2.03 |
      | $10^3$ | 0.96 | 1.08 | 6.81 |
      | $10^4$ | 1.02 | 3.08 | 9.08 |
      | $10^5$ | 3.99 | 6.10 | — |
      | $10^6$ | 6.88 | 7.67 | — |

      Below 1 is the reading of the clock in a batch that is not split (0.06 µs on 100 points of `Linear1D`). 10 000 points of a cheap evaluation (50 µs in all) do not reach the worth, and are not split.
    * A batch of 16 points or fewer calls the serial batch as it is, without timing (every point would be evaluated to be timed, and nothing would be left to decide).
    * The result over all 277 cases of the comparison (§8.3) (8 threads; speed over serial interpx; grouped by the serial time):

      | Serial time | Cases | Median | Min | Max | Largest time added |
      | ---: | ---: | ---: | ---: | ---: | ---: |
      | under 1 µs | 51 | 0.89 | 0.60 | 1.07 | 0.14 µs |
      | 1–10 µs | 57 | 0.95 | 0.77 | 1.05 | 0.74 µs |
      | 10–100 µs | 63 | 0.99 | 0.85 | 1.10 | 6.75 µs |
      | 0.1–1 ms | 67 | 2.88 | 1.35 | 5.36 | — |
      | 1 ms or more | 39 | 5.20 | 3.90 | 6.28 | — |

      Batches that are not split (under 100 µs) are 0.60 to 1.1 times the serial one (above 1 is the spread of the measurement; there is no reason for a batch that is not split to be faster than the serial one). The time added is the reading of the clock to time, and splitting the points timed and the rest into two serial batches. For batches known to be small, use the serial forms.
  * What is timed is time, but the results do not depend on time (as above, they are bit-identical whatever the chunks).
  * A chunk length per type (returned by a trait method) is not taken. It would need a table per type to keep up, would not apply to users' interpolators, and would not catch differences between machines.
* **Thread pool**: `rayon`'s global pool. A caller who wants another pool calls from inside `ThreadPool::install`. A form that takes a pool as an argument is not added (`rayon`'s types would appear in the API). Inside a parallel region nothing is parallelized further.
* **Allocation**: `eval_batch_into_par` allocates nothing of its own (what `rayon` allocates inside is not counted). `eval_batch_par` allocates, besides the result `Vec`, one `Vec` per chunk, which the serial `eval_batch` returns (to keep the serial batch's work. A chunk is one "time of the work of one chunk", so the number of allocations is proportional to the amount of work). The neighbors mode of `RBFInterpolator` allocates the per-thread working space at the thread's first evaluation (the exception above).
* **Use the serial forms for small batches**: the parallel forms are slower than the serial forms of the same work by the reading of the clock (about 0.1 µs a time). On 100 points of a cheap evaluation (0.4 µs) that is a fifth.
* **Tests (the implementation PR)**: for every type, the parallel result is bit-identical to the serial one (with a pool of one thread and of several, with numbers of points that are not split and that are). The error of a length mismatch is the same as the serial one's.

### 3.3 Extrapolation and boundary modes
```rust
#[non_exhaustive]
pub enum Extrapolate { Forbid, Clamp, Extend, Periodic }   // 1D / curve / grid

pub enum BoundaryCondition<V> {                            // CubicSpline / BSpline1D. V is the value type
    NotAKnot,                                 // the default (as SciPy)
    Periodic,
    Ends(EndCondition<V>, EndCondition<V>),   // left end, right end
}
pub enum EndCondition<V> { NotAKnot, FirstDerivative(V), SecondDerivative(V) }
```
* **`Extrapolate`**: `Forbid` gives `None` outside the domain. `Clamp` gives the value at the end (the behavior of `numpy.interp`). `Extend` continues the end polynomial as it is (SciPy's `extrapolate=True`). `Periodic` wraps into the domain.
* **The defaults follow SciPy 1.18 type by type.**

  | Type                                                                    | Default     |
  | :---------------------------------------------------------------------- | :---------- |
  | `Linear1D`                                                              | `Clamp`     |
  | `CubicSpline`, `PchipInterpolator`, `CubicHermiteSpline`, `BSpline1D`, `BSplineCurve` | `Extend`    |
  | `CubicSpline` and `BSpline1D` with the boundary condition `Periodic`    | `Periodic`  |
  | `Akima1D`, `Nearest1D`, `RegularGridInterpolator`                       | `Forbid`    |

* **`BoundaryCondition`**: a condition per end, so that the left and the right can differ (SciPy's `bc_type=((1, a), (2, b))`). `Natural` is `SecondDerivative(0)` at both ends and `Clamped` is `FirstDerivative(0)` at both ends, made with `BoundaryCondition::natural()` / `clamped()`. Periodic is given for both ends together; a state with one end periodic cannot be expressed by the type. The type parameter is **the value type** (settled in M1-8). The derivative at an end has the shape of a value, so a spline whose values are `[T; M]` takes `FirstDerivative([..; M])`, column by column (SciPy also takes an array per column).
* The scattered types have no extrapolation mode: outside the domain (outside the convex hull) is `None` (replaced by `fill`) only.

### 3.4 Ownership and thread safety
* Every interpolator is **owned** and has no lifetime parameter. No `T: 'static` bound is imposed.
* Constructors provide `::new(&[...])` (built from borrowed data). The types that keep the data they are given as it is also provide `::from_owned(...)` (a zero-copy move): `PPoly`, `BSpline`, `RegularGridInterpolator`, `NearestNDInterpolator`, `RBFInterpolator` (these two added in P-6, #128), `LinearNDInterpolator` (M4-1; moves the points and the values). What is moved is only the arrays the type keeps as they are: `NearestNDInterpolator::from_owned(points: &[..], values: Vec<V>)` moves the values and borrows the points (the points go into the tree, and the array given is not kept). `RBFInterpolator::from_owned(points: Vec<..>, values: Vec<V>, kernel, &options)` moves the points, and in the neighbors mode the values too (over all points the values go into the system solved and are not kept). For these two types the time of construction does not change (building the tree and solving take it: 9.2 ms and 8.5 ms for a `NearestNDInterpolator` of $10^5$ points). What they gain is not holding two large arrays. The one-dimensional types and the curve make coefficients from the data and keep nothing but the coordinates, so they have no `from_owned` (all it would save is one copy of the coordinates).
  * **What construction copies and checks** (P-1, #122. The owner's decision): in the comparison with other libraries (§8.3), construction lost most, because the others (`interpn`, `ndarray-interp`, SciPy) only borrow the data and do not check it. It was decided as follows.
    * **No entry that skips the checks** (#146. A decision confirmed after the full measurement): grid construction, even without the copy by `from_owned`, is slower than SciPy in 36 of 87 cases and than `interpn` in 43 of 43 (on a grid of 200 MB: `new` 56 ms, `from_owned` 18 ms, SciPy 0.009 ms). What is left is the check of the values. An entry without the checks (something like `from_owned_unchecked`) would bring it to the same order, but it would be an entry where data with NaN makes evaluation return NaN silently. Construction matters only to uses that rebuild a large grid at short intervals, so the promise is kept. In the comparison's tables, the difference stays as "the cost of the checks".
    * **No borrowing form.** An interpolator owns its data and has no lifetime (it can be carried as `Send + Sync + 'static`). A caller who wants to avoid the copy uses `from_owned`. Types with a lifetime or an entry with `Arc` are not added (the types would split into two families, and there would be three entries).
    * **The check of finite values stays.** Returning `NonFinite { index }` for NaN and infinity at construction is a promise of this crate (SciPy does not check, and NaN comes out in evaluation).
    * **Tighten without changing the API.** The checks look at 64 values at a time (a loop that does not exit midway lets the compiler look at several at once). The check of finiteness goes through one function for every type (`validate::first_failing`): the one-dimensional types and the grid (P-1), the points and values of scattered data, the knots and control points of `BSpline`, the points of a curve, the slopes and coefficients of the grid (aligned in the review of P-6). What remains is only the check of the order of `BSpline`'s knots, which compares two neighbors and is a different form. The slopes of `Linear1D`, the slopes of `Akima1D`, and the reordering of a descending last axis were changed in the same spirit.
  * **The differences that remain, and why** (`docs/performance.md`):
    * **The copy**: `new` copies its input. On a grid of 200 MB this takes 30 ms, and the gap to the others, who only borrow (close to 0), does not close. `from_owned` does not copy.
    * **The first write to new memory**: from 10 000 points (some hundred KB) on, the time of writing to memory just allocated (the assignment of pages) takes most of construction. On this machine (Windows), copying 10 000 points of $x$ and of $y$ alone takes 27 to 47 µs (0.2 µs for 1000 points). Types that keep coefficients (4 per interval) take several times that.
    * **What construction of 10 000 points is made of** (P-8, #131. `Linear1D` and `Nearest1D`, this machine, medians):

      | | `Linear1D` | `Nearest1D` |
      | :--- | ---: | ---: |
      | A program that repeats only construction (5 of 7 runs) | 20–22 µs | 11–12 µs |
      | The same program (the other 2 runs) | 68–71 µs | 45–47 µs |
      | The comparison's runner (all 6 times) | 109–113 µs | 76–78 µs |
      | SciPy (the comparison's runner, 4 times) | 47–51 µs | — |

      The time of the same code differs by 5 times from process to process. The comparison runner's row measured the slow time after the runner has read its data and evaluated once (P-18, §8.3's "A warm-up before construction is timed"). In the comparison with the warm-up, `Linear1D` is 19.6 µs and `Nearest1D` 11.9 µs, the same as the faster of the program that only repeats the construction, and faster than SciPy (44.0 µs and 53.0 µs). The work the type does is small: the three passes of the checks that read the coordinates and values take 5 µs, and writing 30 000 numbers to memory already there takes 3.5 µs. The rest is decided by whether the memory the allocator hands out is "pages not yet backed" (when allocations of 80 to 160 KB are made and freed repeatedly, there is a state in which the memory handed out is a new page each time and one in which it is not. What decides the state has not been looked into). In the state of new pages, the time is proportional to the amount written: `Linear1D` writes 3 numbers per point (the coordinate and 2 coefficients per interval), `Nearest1D` and SciPy 2.
      * **Not changed.** If `Linear1D` kept the values as they are and found the slope at every evaluation, it would write 2, but evaluation would gain a division, and the public `as_ppoly()` would go. The choice made for uses where evaluations far outnumber constructions is kept.
    * **Making coefficients**: the one-dimensional types make the coefficients of each interval's polynomial at construction. Construction is slower than the other that keeps only slopes (`scirs2-interpolate`) (1000 points of `Akima1D`: 9.7 µs against 4.9 µs; `PchipInterpolator`: 8.9 µs against 4.6 µs), and evaluation is faster (15 µs against 1.5 ms on the same cases). A choice made for uses where evaluations far outnumber constructions.
    * **The checks**: even against the other that copies (`scirs2-interpolate`), grid construction is 1.2 to 2.5 times slower (14 µs against 5.6 µs on a grid of 40 × 20 × 10 × 5). The other does not check the values.
    * **Descending axes**: the values are reordered into ascending order and kept so. The last axis by reversing each row, the axes before it by swapping blocks; each axis goes over the whole once.
* **Every interpolator is `Send + Sync`** is a public contract (the bounds of `Scalar` and `Value` guarantee it). There is no shared mutable state such as `Cell` / `RefCell` / `Mutex`; evaluation through `&self` is pure.
* The state of the hinted search (the optimization that remembers the previous interval) is in a `Hint` value the caller holds, or local to a batch.

### 3.5 Input validation
* NaN / Inf in the input data (coordinates and values) is `Err(NonFinite)` at construction.
* A NaN query point does not panic: `eval` returns `None` and `eval_or` returns `fill`.
* **One-dimensional coordinates must be strictly increasing** (the same as SciPy's `CubicSpline`). Only the grid accepts descending axes.
* **The order of the checks is fixed.** For input with several faults, the first error met in the following order is returned.
  1. `LengthMismatch` (the lengths of coordinates and values)
  2. `InsufficientPoints`
  3. `NonFinite` (coordinates first, then values. The first index)
  4. `DuplicatePoints` (neighboring coordinates are equal) or `NotSorted` (a decrease), whichever comes first
* **Types with point-wise input besides the values** (the slopes of `CubicHermiteSpline`) check it after the values at each of the steps above. Step 1 is the length of the values, then of the slopes. Step 3 is coordinates, values, slopes, in that order. The order of the steps does not change (with NaN in the slopes and unsorted coordinates, it returns `NonFinite`).

---

## 4. Error Handling

```rust
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum InterpolateError {
    #[error("insufficient points: required {required}, got {actual}")]
    InsufficientPoints { required: usize, actual: usize },

    #[error("length mismatch: expected {expected}, got {actual}")]
    LengthMismatch { expected: usize, actual: usize },

    #[error("non-finite value at index {index}")]
    NonFinite { index: usize },

    #[error("coordinates are not strictly increasing at index {index}")]
    NotSorted { index: usize },

    #[error("duplicate points at indices {first} and {second}")]
    DuplicatePoints { first: usize, second: usize },

    #[error("periodic boundary needs equal values at both ends")]
    NotPeriodic,
}
```
* It is `#[non_exhaustive]`, so a variant is added in the row that first needs it (not defined before it is used). The above are what M1 (1D, curve) uses. `NotPeriodic` is defined in M1-2 and first returned in M1-8 (`CubicSpline`).
* M1-8 adds `NonFinite` to `ParamReason` (an argument that is NaN or infinite: the derivatives of a boundary condition).
* M1-12 adds `TooSmall { min: usize }` (below a lower bound: degree 0 of `BSpline1D`) and `Unsupported` (the value itself is valid but cannot be used with the other arguments: derivatives at an end given to a `BSpline1D` that is not cubic).
* Variants added in later rows:

  | Variant                                                    | Row that adds it   | Use                                                    |
  | :--------------------------------------------------------- | :----------------- | :----------------------------------------------------- |
  | `InvalidParameter { name: &'static str, reason: ParamReason }` | M1-4 (`BSpline`) | An argument out of range. M1-4 puts in `ParamReason::TooLarge { max: usize }` (the upper bound of the degree). A `ParamReason` such as a required argument missing or a degree under the minimum is added by the row that returns it |
  | `GridAxis { axis: usize, fault: AxisFault }` | M2-1 (`RegularGridInterpolator`) | A fault of the grid's axis `axis`. `AxisFault` is `InsufficientPoints { required, actual }`, `NonFinite { index }`, `NotMonotonic { index }`, `DuplicatePoints { first, second }` (positions within that axis). The 1D variants cannot carry which axis |
  | `SolverFailed(SolverError)`                                | M3 (RBF)           | A singular linear system (`SolverError::Singular`)     |
  | `AllocationFailed { bytes: usize }`                        | #110               | An allocation far larger than the input failed for lack of memory ("The size of allocations and their failure" below) |
  | `DegenerateGeometry(DegeneracyKind)`                       | M4 (`LinearND`)    | `LowerDimensional { affine_dim }`, `NumericalLimit`    |
  | `SolverError::NotConverged`                                | M4-2 (`CloughTocher2D`) | the system of the gradient estimate was not solved within the limit on steps (§6.3) |

* The `index` of `NonFinite { index }` is the position of the coordinate when a coordinate is not finite, and the position of the value only when every coordinate is finite (the coordinates are checked first. §3.5). The variant does not say which array.
* No string payloads (only typed values that can be matched). `ParamReason`, `AxisFault`, `SolverError`, and `DegeneracyKind` are `#[non_exhaustive]` too.
* There is no variant for a mismatch of dimensions (with `const D` only, the type cannot express it).
* The `index` of `NotSorted { index }` is the first element smaller than the coordinate before it. The `expected` of `LengthMismatch` is the length of the coordinates (or queries), and `actual` the length of the values (or the output buffer).
* **`convx::Error` is not public** (so that `convx`'s version and SemVer are not tied to this crate's). In M4-1 it is converted inside as follows. `convx::Error` is `#[non_exhaustive]`, so there is a place for kinds not known:
  * `InsufficientPoints` → `InsufficientPoints` (`actual` is the number of distinct points after merging)
  * `DegenerateDimension{actual_dim,..}` → `DegenerateGeometry(LowerDimensional{affine_dim})`
  * `ExactEvaluationExhausted` → `DegenerateGeometry(NumericalLimit)`
  * Converted so that none is lost, though the checks of the constructor come first and they do not happen: `NonFiniteCoordinate{index}` → `NonFinite{index}`, `LengthMismatch` → `LengthMismatch`, `TooManyPoints` → `InvalidParameter` (`TooLarge { max: u32::MAX }`)
  * What the Delaunay build does not return (`NonPositiveDimension` (`convx` is not called with zero dimensions), `NonFiniteFacetPlane`, `NonFiniteCircumcenter`), and kinds a later version adds → `DegenerateGeometry(NumericalLimit)`
* **The size of allocations and their failure**: the library does not fail in a way the caller cannot handle (a panic, the end of the process). For allocations it holds to this in two steps.
  1. **The limit of a `Vec`**: when the caller's arguments decide the size of an allocation, and it exceeds what a `Vec` can hold (`isize::MAX` bytes), it is checked before the allocation and `InvalidParameter { reason: TooLarge { max } }` is returned (allocating it as it is would panic). `max` means "up to this value it fits in a `Vec`".
  2. **Lack of memory**: within the limit there are still allocations that do not fit in the machine's memory. **An allocation that can exceed a constant multiple of the size of the input** is made in the form that returns its failure as an error (`Vec::try_reserve_exact`), and a failure returns `AllocationFailed { bytes }` (`bytes` the size it tried to allocate). Of the places where a small input can cause a huge allocation, four apply now.

     | Place | Size | Example |
     | :--- | :--- | :--- |
     | The global matrix of `RBFInterpolator` (§6.3) | the square of the number of points | 80 GB for $10^5$ points (an input of 1.6 MB) |
     | `PPoly::antiderivative` (§6.1) | intervals × (degree + $\nu$) | grows without bound with $\nu$ |
     | The working space of the neighbors mode of `RBFInterpolator` (§6.3) | the square of the number of neighbors | 80 GB for $k = 10^5$ |
     | The copy of the simplices of `LinearNDInterpolator`, and the transform of each (§6.2) | the number of simplices × ($2D^2 + D$); in high dimensions the simplices can be far more than the points | 79 simplices for 25 points in 22 dimensions (measured); the bound on the number of simplices is of the order of the number of points to the power $\lceil D/2 \rceil$ |

     A later row that adds an allocation meeting this condition makes it in the same form and adds it to this table.
  * **What is not covered**: allocations only proportional to the input (copies of the arrays given, coefficients of the same size, the `Vec` that `eval_batch` returns) are made as before. If the input is in memory, a constant multiple of it nearly is too. Covering every allocation would mean replacing some tens of places on the construction paths, and would protect nothing more.
  * **What cannot be covered** (written in the rustdoc too):
    * Allocations inside the dependencies (`faer`, `kiddo`, `convx`). If one fails, the process ends. The three above are all the crate's own allocations (the matrix of `RBFInterpolator` is allocated by the crate, and `faer` factors it in place).
    * Environments that promise more memory than there is (Linux by default). The allocation succeeds, and the process may be stopped when it writes. This works in environments that report a failed allocation (Windows, Linux with overcommit off). macOS too may let a large allocation succeed and then get stuck in the writes (not confirmed on a machine; what CI confirms on macOS is only that an allocation too large for the address space fails).
  * **History**: in M1-5 (#24) it had been decided not to catch lack of memory. For two reasons: (1) returning it as `TooLarge { max }` would stop `max` meaning "up to this value it passes", and (2) the inside of the dependencies cannot be caught, so making only part of it an error would not be a rule. In M3-3 (#104) it turned out that 100 000 points end the process with an allocation of 80 GB, and leaving no failure the caller cannot handle was preferred (#110). (1) is solved by a variant of its own. (2) remains, and is written above under "What cannot be covered".
* Cospherical point sets are not an error (`convx` decomposes them into simplices). The degeneracy errors are only "all points on a lower-dimensional subspace" and "fewer than D+1 distinct points after merging".

---

## 5. Module Structure & SciPy Coverage (Coverage Matrix)

### 5.1 Module structure
```text
interpx/
├── oned/             # interpolation in one dimension
│   ├── linear        # numpy.interp
│   ├── cubic_spline  # CubicSpline
│   ├── pchip         # PchipInterpolator
│   ├── akima         # Akima1DInterpolator (makima included)
│   ├── cubic_hermite # CubicHermiteSpline
│   ├── hermite       # the common procedure that makes a cubic PPoly from values and slopes (crate-private. §6.1)
│   ├── bspline1d     # make_interp_spline (degree 1 to 15)
│   └── nearest       # Nearest1D (nearest / nearest-up / previous / next)
├── curve/            # N-dimensional parametric curves
│   └── bspline_curve # BSplineCurve, Parametrization
├── grid/             # interpolation on N-dimensional rectilinear grids
│   └── regular_grid  # RegularGridInterpolator
└── scattered/        # interpolation of N-dimensional irregular scattered data
    ├── nearest_nd    # NearestNDInterpolator (kiddo)
    ├── mesh          # the triangulation and point location shared by LinearND and CloughTocher2D
    ├── linear_nd     # LinearNDInterpolator (convx + walk from cells)
    ├── clough_tocher # CloughTocher2DInterpolator (Mesh + gradient estimate)
    └── rbf           # RBFInterpolator (faer + kiddo)
```
(The module is named `oned`, avoiding Rust's constraint on identifiers.)

* This is where the source lives. Modules such as `oned` are **private**, and every public item is re-exported from `lib.rs`. Public paths are at the crate root, like `interpx::Linear1D` (settled in M1-7).
* The foundations the interpolators share (`PPoly`, `BSpline`, the interval search, the input checks, the tridiagonal solver) are outside this tree, directly under `src/`.

### 5.2 Coverage table

| Category          | SciPy 1.18                                      | `interpx` struct                              | Core engine inside                     |
| :---------------- | :---------------------------------------------- | :-------------------------------------------- | :------------------------------------- |
| **1D**            | `numpy.interp`                                  | `Linear1D`                                    | `PPoly`                                |
| **1D**            | `CubicSpline`                                   | `CubicSpline`                                 | Thomas algorithm + `PPoly`             |
| **1D**            | `PchipInterpolator`                             | `PchipInterpolator`                           | `PPoly` (derivatives limited for monotonicity) |
| **1D**            | `Akima1DInterpolator` (with `method='makima'`)  | `Akima1D`                                     | `PPoly`                                |
| **1D**            | `CubicHermiteSpline`                            | `CubicHermiteSpline`                          | `PPoly`                                |
| **1D**            | `make_interp_spline` (degree 1 to 15)           | `BSpline1D`                                   | `BSpline` (de Boor)                    |
| **1D**            | nearest / previous / next of `interp1d`         | `Nearest1D`                                   | binary search                          |
| **1D**            | derivatives and integrals of `PPoly` / `BSpline` | `derivative` / `antiderivative` / `integrate` | conversion `PPoly` ⇄ `BSpline`         |
| **N-D Curve**     | `make_interp_spline` (N-D y)                    | `BSplineCurve<T, D>`                          | `BSpline` (vector control points)      |
| **N-D Grid**      | `RegularGridInterpolator` `nearest`             | `RegularGridInterpolator`                     | rounding of the index                  |
| **N-D Grid**      | `linear` / `slinear` (the same function)        | the same (`GridMethod::Linear`)               | tensor-product linear                  |
| **N-D Grid**      | `cubic` / `quintic`                             | the same                                      | tensor-product B-spline coefficients at construction (§6.4) |
| **N-D Grid**      | `pchip`                                         | the same                                      | PCHIP along each axis, from the last axis (§6.4) |
| **N-D Scattered** | `NearestNDInterpolator`                         | `NearestNDInterpolator<D>`                    | `kiddo`                                |
| **N-D Scattered** | `LinearNDInterpolator`                          | `LinearNDInterpolator<D>`                     | `convx` + walk from cells              |
| **N-D Scattered** | `CloughTocher2DInterpolator`                    | `CloughTocher2DInterpolator`                  | `convx` + gradient estimate + Bernstein-Bézier |
| **N-D Scattered** | `RBFInterpolator`                               | `RBFInterpolator<D>`                          | `faer` LU + `kiddo`                    |

* `interp1d` is a legacy API in SciPy. `kind='quadratic'` / `'slinear'` are replaced by `BSpline1D`. `interp2d` has been removed from SciPy and is out of scope.
* `griddata` is not provided as a function (the owner's decision, #179). `griddata` builds `NearestNDInterpolator`, `LinearNDInterpolator`, or `CloughTocher2DInterpolator` according to its `method`, evaluates, and drops it. The same is written with the three types, and a type can be evaluated again as often as needed. A function would add no capability and make two entries to the same thing.

### 5.3 Out of scope (not implemented in v0.1.0. §9 roadmap)
* The FITPACK family: `UnivariateSpline`, `splrep`, `splprep`, `RectBivariateSpline` and the other FITPACK 2D classes
* `BarycentricInterpolator`, `KroghInterpolator`
* `NdPPoly`
* `roots` / `solve` (`PPoly` / `BSpline`)
* Derivatives on grids and of scattered data
* A wrapper with a dynamic dimension, `no_std`, a pseudo-inverse fallback for RBF

---

## 6. Internal Algorithms

### 6.1 Spline foundations (`PPoly` and `BSpline`)
* **`PPoly`**: keeps the coefficients in **a flat array contiguous per interval** (`[interval][power][output column]`, the output column innermost). This differs from SciPy's `(k+1)×(N-1)` layout, for locality of evaluation and vectorization. `eval` finds the interval by binary search (or a Hint) and evaluates by Horner's method.
  * **Type**: `PPoly<T: Scalar, V: Value<T> = T>`. The degree $k$ is a run-time value. It keeps breakpoints $x_0 < x_1 < \dots < x_n$ and $k+1$ coefficients per interval.
  * **The order of the coefficients**: the value on interval $i$ is $\sum_{j=0}^{k} c_{i,j}\,(x - x_i)^{k-j}$, **the highest power first** within an interval. Horner's method reads the memory front to back. SciPy's `c[power][interval]` also has the highest power first, so these are the same numbers reordered per interval.
  * **Construction**: `PPoly::new(&[T], &[V], degree)` (a copy) and `from_owned(Vec<T>, Vec<V>, degree)` (a move). The degree is an explicit argument and is not inferred from the length of the coefficients. The checks are in the order of §3.5: a number of coefficients other than `intervals × (degree + 1)` is `LengthMismatch`, and fewer than 2 breakpoints `InsufficientPoints`.
  * **Which interval a point is in**: intervals are closed on the left and open on the right, $[x_i, x_{i+1})$, and only the last one contains its right end (the same as SciPy).
  * **Outside the domain**: the four modes of `Extrapolate`. The default is `Extend` (the same as SciPy's `PPoly`), changed with `with_extrapolate`. `Periodic` wraps to $x_0 + ((x - x_0) \bmod (x_n - x_0))$. The right end $x_n$ wraps to $x_0$ too (one period contains only its left end; the same as SciPy). For a periodic function this makes no difference, but where the values or derivatives differ at the two ends, it gives the value of the left end's interval. A NaN query is `None`. An infinite query is `None` under `Forbid` and `Periodic`, the end value under `Clamp`, and under `Extend` the value of the end polynomial evaluated as it is (in general infinite or NaN).
  * **`Hint`**: `eval_hinted(x, &mut Hint)` looks at the hint's interval and its two neighbors first, and goes back to binary search when they miss. The result is bit-identical to `eval`; only the time changes.
  * **When a hint pays**: the hinted search starts from the interval of the previous point, so each point waits for the previous point's result. When close points follow each other, the wait is short and it is faster by the binary search saved (about 2 times on 10 000 ascending points). On points in no order, the wait is a whole binary search, and it is about 2 times slower than a plain search where points do not wait for each other (measured by `benches/ppoly.rs`). `eval_hinted` is for queries in order; use `eval` for points in no order.
  * **Batches**: `eval_batch_into` splits the points into runs of 64 and decides per run whether to use the hint. It does when the run is monotone (ascending or descending) and spans no more intervals than it has points. A run of fewer than 16 points is not tested and is searched plainly (the fixed cost of the test is larger). On every path the interval of each point is that of `eval`, and the results are bit-identical.
* **`BSpline`**: keeps a knot vector and control points, and evaluates as a linear combination of the nonzero basis functions and control points.
  * **Type**: `BSpline<T: Scalar, V: Value<T> = T>`. The degree $k$ is a run-time value. It keeps control points $c_0, \dots, c_{n-1}$ and knots $t_0 \le t_1 \le \dots \le t_{n+k}$ (repeats allowed).
  * **Construction**: `BSpline::new(&[T], &[V], degree)` (a copy) and `from_owned(Vec<T>, Vec<V>, degree)` (a move). The checks are in this order, and the first error met is returned.
    1. `InvalidParameter { name: "degree", reason: TooLarge { max: 15 } }` (the degree exceeds the upper bound)
    2. `LengthMismatch` (the number of knots is not `control points + degree + 1`)
    3. `InsufficientPoints` (fewer than `degree + 1` control points)
    4. `NonFinite` (knots first, then control points)
    5. `NotSorted` (the knots decrease)
    6. `DuplicatePoints { first: k, second: n }` (a domain of zero width, that is $t_k = t_n$)
  * **Domain, and which interval a point is in**: the domain is $[t_k, t_n]$. The interval of a point $x$ is the knot interval of positive width with $t_\ell \le x < t_{\ell+1}$, and the right end of the domain $t_n$ belongs to the last interval of positive width (closed on the left and open on the right like `PPoly`, closed only at the end). An interval of zero width made by repeated knots is never the interval of any point.
  * **Evaluation**: on interval $\ell$, the values of the nonzero basis functions $B_{\ell-k}, \dots, B_\ell$ ($k+1$ of them) are found in an array on the stack by the de Boor–Cox recurrence, and $\sum_j B_{\ell-k+j}(x)\,c_{\ell-k+j}$ is summed with `Value::axpy`. The interval $\ell$ evaluated has positive width, so every divisor of the recurrence ($t_{\ell+n} - t_{\ell+n-j}$) spans interval $\ell$ and is positive even with repeated knots. The branch of SciPy's implementation (`_deBoor_D`) that "skips the terms of equal knots" is not needed.
  * **The form of the recurrence**: values are found in the triangular form that keeps the distances from $x$ to the knots on both sides and updates the values in place (Piegl–Tiller's A2.2). SciPy's form, which "copies the working array per degree", is used only for derivatives (the Taylor coefficients of `from_bspline`, and M1-5). Computing values in the latter was about 2 times slower for cubics (measured by `benches/bspline.rs`: 40.8 µs against 21.4 µs on 1000 points). Taking the reciprocals of the divisors first made a difference within the noise, and is not taken.
  * **The upper bound of the degree**: evaluation does not allocate (§3.2). The working arrays are on the stack with a fixed length, so the degree is at most 15. What SciPy's `make_interp_spline` uses in practice is degree 5 at most. The bound is public as the crate-level constant `BSPLINE_MAX_DEGREE` (as an associated constant of `BSpline` it could not be written without type parameters).
  * **Outside the domain**: the four modes of `Extrapolate`, the default `Extend` (the same as SciPy's `BSpline`). `Extend` continues the polynomial of the first or the last interval of positive width. `Periodic`, NaN, and infinite queries are as for `PPoly`.
  * **The shared interval search**: binary search, the hinted search, and the per-run decision of a batch to use the hint are in the crate module `search`, shared by `PPoly` and `BSpline`. The `Hint` type is common, and `BSpline` has `eval_hinted` too.
* **Conversion between `PPoly` and `BSpline`**: both are exact conversions (but for rounding errors), and carry the extrapolation mode over.
  * **`PPoly::from_bspline(&BSpline) -> Result<PPoly, _>`** (SciPy's `PPoly.from_spline`): the breakpoints are the distinct knots within the domain $[t_k, t_n]$. The coefficients of interval $[u_j, u_{j+1}]$ are the values at $u_j$ of the spline's $m$-th derivatives divided by $m!$ (the Taylor coefficients about $u_j$). SciPy keeps every knot as a breakpoint, so repeated knots become intervals of zero width; the breakpoints of `PPoly` are strictly increasing, so no interval of zero width is made.
  * **`BSpline::from_ppoly(&PPoly) -> Result<BSpline, _>`**: a knot of multiplicity $k+1$ is put at each breakpoint, and on each interval the coefficients in the power basis are turned into coefficients in the Bernstein basis (Bézier control points). It assumes no continuity between intervals, so a discontinuous piecewise polynomial is represented exactly too. There are `intervals × (k + 1)` control points. Examining the continuity to lower the multiplicity of knots is out of scope for v0.1.0 (roadmap A-9). A `PPoly` of degree over 15 returns `InvalidParameter`. **Precision outside the domain**: inside the breakpoints the two representations agree to rounding. When an end interval is short, continuing its Bézier representation far outside the domain loses digits that `PPoly` keeps ($r$ times the width of the end interval away, the rounding error grows about $r^k$ times: a cubic interval of width $10^{-3}$ evaluated 700 widths away gave a relative $4 \times 10^{-8}$). Evaluate far extrapolation with `PPoly`.
  * The round trips (`PPoly` → `BSpline` → `PPoly`, and the reverse) are tested to give values that agree within a tolerance.
* **Derivatives and integrals (M1-5)**: four functions of the same form on `PPoly` and `BSpline`. SciPy 1.18 behaves differently for the two types in three places; here they are made the same ("Differences from SciPy" below).
  * **Functions**:
    * `derivative(nu: usize) -> Result<Self, InterpolateError>`: the $\nu$-th derivative.
    * `antiderivative(nu: usize) -> Result<Self, InterpolateError>`: the $\nu$-th antiderivative.
    * `integrate(a: T, b: T) -> Option<V>`: the definite integral from $a$ to $b$.
    * `eval_derivative(x: T, nu: usize) -> Option<V>`: the value of the $\nu$-th derivative at $x$.
  * `derivative` and `antiderivative` make a new object (they allocate). `nu = 0` returns a copy. Coefficients that overflow return `NonFinite`. `integrate` and `eval_derivative` are evaluations and do not allocate (§3.2; covered by `tests/alloc.rs`). A `nu` whose `antiderivative` does not fit in a `Vec` is `InvalidParameter { name: "nu", reason: TooLarge { max } }`, and one that fits but not in memory `AllocationFailed { bytes }` ("The size of allocations and their failure" of §4).
  * **`derivative`**:
    * `PPoly`: the breakpoints stay, the degree is $k - \nu$. The coefficient of $(x - x_i)^m$ is multiplied by $m(m-1)\cdots(m-\nu+1)$, and the $\nu$ lowest are dropped.
    * `BSpline`: per order, one knot is dropped at each end, and the control points become $c'_i = k\,\dfrac{c_{i+1} - c_i}{t_{i+k+1} - t_{i+1}}$ ($n - 1$ of them). The domain $[t_k, t_n]$ does not change. The denominator is 0 when $t_{i+1} = \dots = t_{i+k+1}$; that basis function has an empty support, so its control point is 0. So a spline with knots of multiplicity $k+1$, as `from_ppoly` makes, can be differentiated too.
    * For $\nu > k$ both types return the zero function. `PPoly` is of degree 0 with the breakpoints as they are, `BSpline` of degree 0 with knots only at the two ends of the domain.
  * **`antiderivative`**:
    * `PPoly`: the degree is $k + \nu$. Per order, the coefficient of $(x - x_i)^m$ is divided by $m + 1$ and its power raised, and the constant terms are summed up from the left so that the antiderivative is continuous at the breakpoints. `PPoly` has no upper bound on the degree, but a `nu` whose resulting coefficients exceed what a `Vec` can hold (`isize::MAX` bytes) returns `InvalidParameter { name: "nu", reason: TooLarge { max } }` (`max` the largest order that fits within that bound. Allocation as in §4).
    * `BSpline`: per order, an end knot is added at each end, and the control points are $C_0 = 0$, $C_{j+1} = C_j + c_j\,\dfrac{t_{j+k+1} - t_j}{k+1}$ ($n + 1$ of them). A result of degree over 15 returns `InvalidParameter { name: "nu", reason: TooLarge { max } }` (`max` is $15 - k$).
    * **The constant of integration**: for both types, 0 at the left end of the domain. `BSpline` subtracts from the control points above their value at the left end of the domain (the basis functions sum to 1, so the same value is subtracted from every control point). With end knots of multiplicity $k+1$ the value subtracted is exactly 0, which agrees with SciPy.
  * **The extrapolation mode carried over**:

    | Original mode | `derivative` | `antiderivative` |
    | :------------ | :----------- | :--------------- |
    | `Extend`      | `Extend`     | `Extend`         |
    | `Forbid`      | `Forbid`     | `Forbid`         |
    | `Periodic`    | `Periodic`   | `Forbid`         |
    | `Clamp`       | `Forbid`     | `Forbid`         |

    The antiderivative of a periodic function is in general not periodic (SciPy too changes it to "no extrapolation"). `Clamp` continues outside the domain with a constant, so the derivative is 0 outside and the antiderivative linear outside, neither expressible in the four modes. Where they cannot be expressed, the mode is `Forbid`, which gives no value. When the values outside of `Clamp` are needed, the original object's `eval_derivative` and `integrate` return them.
  * **`eval_derivative`**: returns the $\nu$-th derivative of the function that `eval` defines.
    * The interval is chosen as by `eval`. At a breakpoint where the derivative jumps it gives the value of the interval on the right, and only at the right end of the domain the value of the interval on the left (the same as SciPy).
    * Outside the domain, `Extend` gives the derivative of the end polynomial, `Periodic` the value at the wrapped point, `Forbid` `None`. `Clamp` gives 0 for $\nu \ge 1$ (it is constant there). NaN is `None`.
    * $\nu = 0$ is bit-identical to `eval`. $\nu > k$ is 0 at points that have a value.
    * The evaluation is Horner's method multiplying the coefficients by ratios of factorials for `PPoly`, and the recurrence for derivatives (the "form of the recurrence" above) for `BSpline`.
    * **No batch form and no hinted form**. For many points, use `eval_batch_into` on the object `derivative(nu)` makes. For `BSpline` this path is faster even per point (the recurrence for values is faster). `derivative(nu)` does not carry `Clamp` over, so points outside `Clamp` are found with `eval_derivative`.
  * **`integrate`**: returns the definite integral of the function that `eval` defines.
    * $a > b$ gives the value with its sign reversed, $a = b$ 0. If $a$ or $b$ is not finite, `None`.
    * `Forbid`: `None` if $a$ or $b$ is outside the domain.
    * `Extend`: integrates the end polynomials continued.
    * `Clamp`: the part outside is the end value × the length.
    * `Periodic`: with period $P$, $b - a$ is split into a whole multiple of $P$ and a remainder; the multiple is that multiple of the integral over one period, and the remainder is found on the wrapped interval (split in two where it crosses the right end of the domain) (the same procedure as SciPy).
    * `PPoly`: sums, over the intervals spanned, the difference of the polynomial's antiderivative between their ends. The cost is proportional to the number of intervals spanned.
    * `BSpline`: finds $F(b) - F(a)$ directly, without making the control points of the antiderivative $F$ in an array. On the interval $\ell_a$ of $a$ and $\ell_b$ of $b$, the values of the basis functions of degree $k+1$ ($k+2$ each) are found on the stack, and $C_j$ is summed up from $j = \ell_a - k$ to $\ell_b + 1$, each multiplied by its weight (the basis value at $b$ − the basis value at $a$) and added. The basis functions sum to 1, so starting from $C_{\ell_a - k} = 0$ does not change the difference. The cost is proportional to the number of knot intervals spanned. The basis of degree $k+1$ can be computed on the original knot vector as it is (the knots added at both ends do not appear in the recurrence of the intervals of the domain). The working array has 17 entries, so a spline of degree 15 can be integrated too.
  * **Differences from SciPy** (to make `PPoly` and `BSpline` behave alike):

    | Case | SciPy 1.18 | `interpx` |
    | :--- | :--------- | :-------- |
    | `derivative` of a `BSpline` with $\nu > k$ | `ValueError` | the zero function (as `PPoly`) |
    | `derivative` of a `BSpline` with an inner knot of multiplicity $k+1$ | `ValueError` | differentiable |
    | Antiderivative of a `BSpline` whose end knots are not of multiplicity $k+1$ | not 0 at the left end of the domain | 0 at the left end of the domain (as `PPoly`) |
    | `integrate` leaving the domain without extrapolation | NaN for `PPoly`; `BSpline` integrates clipped to the domain | `None` for both |
    | A NaN end of integration | an exception for `PPoly`, NaN for `BSpline` | `None` for both |

    Each was confirmed by running SciPy 1.18.1.
* **Tridiagonal solver (M1-6)**: for the linear system that the construction of `CubicSpline` solves; kept inside the crate and not public (`scipy.interpolate` has no corresponding public function either).
  * **How the matrix is held**: three arrays of length $n$, `lower`, `diag`, `upper`: `lower[i]` is $A_{i,i-1}$, `diag[i]` is $A_{i,i}$, `upper[i]` is $A_{i,i+1}$. `lower[0]` and `upper[n-1]` are outside the tridiagonal band, and the non-periodic form does not read them.
  * **Right-hand side**: an array of `V: Value<T>`. With values `[T; M]`, the $M$ columns are solved in one elimination. Scalars are `T: Scalar`.
  * **No allocation**: it rewrites the arrays given. `diag` becomes the reciprocals of the pivots after elimination, and the right-hand side becomes the solution. `lower` and `upper` are only read. The allocations of construction are gathered in one place in the caller.
  * **Method**: the Thomas algorithm without pivoting (forward elimination and back substitution, $O(n)$). One pass down does the elimination and the forward substitution, one pass up the back substitution. The cyclic form has two right-hand sides (the data's, and Sherman–Morrison's $u$), and carries both on the same two passes.
  * **Keeping the chain of dependence short**: this method is a recurrence, and each row waits for the result of the row before. What decides its speed is less the number of operations than the length of the sequence of operations from one row's result to the next row's (the chain of dependence). The pivot is written $p_i = d_i - (l_i u_{i-1})/p_{i-1}$, with the product $l_i u_{i-1}$ found off the chain, so that the chain is one division and one subtraction. The substitution multiplies by the reciprocal of the pivot first and then subtracts the neighboring row, so that its chain is one multiplication and one addition. Separate passes for elimination and substitution would add the times of their chains, so they are one pass. The form that "reduces divisions to one per row" (keeping only the reciprocal of the pivot in the recurrence) has fewer divisions but a longer chain, and was not faster (measured by `benches/cubic_spline.rs`: 21 µs for the construction of 1000 points, unchanged; the form with the short chain 16 µs).
  * **Why no pivoting**: SciPy's `CubicSpline` solves with a banded solver with partial pivoting (`solve_banded`), but the matrices solved here are solved stably without it. The matrices of the natural, clamped (first derivative), and second-derivative boundary conditions are strictly diagonally dominant. Not-a-knot has first and last rows that are not diagonally dominant, but after the first row $[\,h_1,\; h_0 + h_1\,]$ eliminates the second $[\,h_1,\; 2(h_0 + h_1),\; h_0\,]$, the pivot is $h_0 + h_1 > 0$ ($h_i = x_{i+1} - x_i$). The pivot of the last row $[\,h_{n-3} + h_{n-2},\; h_{n-3}\,]$ is $h_{n-3}\,\bigl(1 - (h_{n-3} + h_{n-2})/p\bigr)$, and the pivot $p$ of the row before is larger than $2h_{n-3} + h_{n-2}$, so this is positive too. Agreement with SciPy is confirmed on the goldens of `CubicSpline` in M1-8.
  * **Breakdown**: when a pivot becomes 0 or not finite, a crate-private error naming that row is returned (no panic). Which public variant it maps to is decided by M1-8, which first returns it. The caller checks whether the solution is finite (the construction of `PPoly` checks the coefficients).
  * **Cyclic tridiagonal (periodic boundary)**: a system with one entry each at the upper right $A_{0,n-1}$ and the lower left $A_{n-1,0}$ besides the tridiagonal band, solved by the Sherman–Morrison formula. The upper right is put in `lower[0]` and the lower left in `upper[n-1]` (the two entries the non-periodic form does not read). With $\gamma = -A_{0,0}$, it is split as $A = A' + u v^{\mathsf T}$, $u = (\gamma, 0, \dots, 0, A_{n-1,0})^{\mathsf T}$, $v = (1, 0, \dots, 0, A_{0,n-1}/\gamma)^{\mathsf T}$. $A'$ is tridiagonal: $A$ with the two ends of the diagonal replaced by $A_{0,0} - \gamma$ and $A_{n-1,n-1} - A_{n-1,0}A_{0,n-1}/\gamma$. It solves $A' y = b$ and $A' z = u$, and sets $x = y - z\,\dfrac{v^{\mathsf T} y}{1 + v^{\mathsf T} z}$. One elimination and two substitutions; the caller passes the working array of length $n$ for $z$.
  * **Sizes**: the non-periodic form takes $n \ge 1$, the cyclic form $n \ge 3$ (for $n \le 2$ the cyclic form is not determined). Periodic splines with few points are handled separately by `CubicSpline`, as SciPy does (M1-8).
* **Banded solver (M1-12)**: the collocation matrix of `BSpline1D` and the per-axis systems of the grid's `cubic` / `quintic` (§6.4) have wider bands than tridiagonal (for degree $k$, with a not-a-knot knot vector, both the lower and the upper bandwidth are at most $k$). They are solved by **a banded LU with partial pivoting** (SciPy uses LAPACK's `gbsv`). M1-12 implements it inside the crate, and M2-2 uses it. The tridiagonal case is not served by the banded LU because the construction of `CubicSpline` is the most used path, and the dedicated Thomas algorithm is faster.
  * **Place and form**: the crate module `banded`. The matrix is `T` and the right-hand side `V: Value<T>`, so values `[T; M]` are solved column by column independently (the choice of pivots depends on the matrix alone). Factorization and solution are separate: one factorization solves any number of right-hand sides of different types (a periodic spline solves a right-hand side of `T` and one of `V` with the same factorization).
  * **Storage**: an $n \times n$ matrix with $k_l$ lower and $k_u$ upper diagonals is kept in the banded form of LAPACK's `gbsv` ($2k_l + k_u + 1$ entries per column; the top $k_l$ rows are room for the fill that row interchanges make). With $k_l = k_u = n - 1$ the same function solves a dense matrix (used for the small systems of periodic splines).
  * **Elimination**: per column, the entry of largest absolute value on or below the diagonal is chosen as the pivot and the rows are interchanged (partial pivoting). A pivot that is 0 or not finite returns `Breakdown { row }` (the same form as the tridiagonal solver; no panic). The caller checks whether the solution is finite.
  * **Allocation**: the matrix and the record of interchanges are in space the caller passes. The solver allocates nothing.
  * **Cyclic banded matrices**: a matrix with corner blocks at the upper right and lower left outside the band (periodic splines) is split into the banded matrix $A$ without the corners and a low-rank correction $U V^{\mathsf T}$, and solved by the Woodbury formula (the same as SciPy's `_woodbury_algorithm`): solve $A Z = U$ and $A y = b$, and $x = y - Z\,(I + V^{\mathsf T} Z)^{-1} V^{\mathsf T} y$. The rank of the correction is $k_l + k_u$, and $I + V^{\mathsf T} Z$ is a dense matrix of that size (solved by the same solver). The corners must not overlap the band ($n > k_l + k_u$); smaller systems are solved as dense matrices.
* **Interpolators that are piecewise polynomials, and `PPoly`**: `Linear1D`, `CubicSpline`, `PchipInterpolator`, `Akima1D`, and `CubicHermiteSpline` keep a `PPoly` inside and lend it with `as_ppoly(&self) -> &PPoly<T, V>` (no allocation). Derivatives and integrals (`derivative`, `antiderivative`, `integrate`, `eval_derivative`) are those of this `PPoly`. The interpolator types do not carry the same four functions again. The extrapolation mode of the `PPoly` is the interpolator's. `BSpline1D` (M1-12) keeps a `BSpline` by the same rule and lends it with `as_bspline(&self) -> &BSpline<T, V>`.
* **`Linear1D<T: Scalar, V: Value<T> = T>` (M1-7, `numpy.interp`)**: a polyline.
  * **Construction**: `Linear1D::new(&[T], &[V])`. The checks are in the order of §3.5, with at least 2 points (`InsufficientPoints { required: 2 }`).
  * **What it keeps**: a `PPoly` of degree 1 (the coefficients of interval $i$ are the slope $s_i = (y_{i+1} - y_i)/(x_{i+1} - x_i)$ and the left value $y_i$), and the last data value $y_{n-1}$.
  * **Evaluation**: $y_i + s_i\,(x - x_i)$ on interval $i$ (the same formula as `numpy.interp`). At the left end $x = x_i$ it returns $y_i$ exactly. The only right end of an interval $x = x_{i+1}$ that is met is the right end of the domain, and there $y_{n-1}$ is returned as it is, without the formula (the formula can be 1 ulp off). So the data value is returned exactly at every node. The `eval` of `as_ppoly()` is the formula, so it can differ by 1 ulp at the right end of the domain.
  * **Outside the domain**: the four modes. The default is `Clamp` (the end data value, as `numpy.interp`; on the right, $y_{n-1}$ exactly). `Extend` continues the end segment. `Periodic` wraps into the domain, as for the other types.
  * **The arguments of `numpy.interp`**: there are no fill values outside, `left` / `right`. Use `Forbid` and `eval_or(x, fill)` (for different values on the two sides, the caller takes the `None` and compares with the domain). There is no `period` either. `numpy.interp`'s `period` wraps and reorders the nodes too by the period, which means something other than `Periodic`.
  * **`eval_hinted` and batches**: as for `PPoly` (through `search`).
  * **Evaluation at low degrees** (P-7, #130): the basis function values are found in an array of 4 up to degree 3, and in an array of 16 (up to degree 15) above that. It is the same recurrence, so the values are the same bits. When three arrays of 16 were filled with zeros every time at every degree, that took 20 to 30% of the evaluation at low degrees (evaluating 40 000 ordered points on 10 000 points: degree 1 316 → 262 µs, degree 2 422 → 329 µs, degree 3 593 → 400 µs. SciPy: 304, 768, 1386 µs).
  * **Evaluation of points in no order is slower than the piecewise polynomial types** (657 µs on 10 000 points at degree 1; `Linear1D` 191 µs). Measured with stages cut out, the recurrence of the basis functions (at least one division per point) takes 5.5 ns for points in no order and 2 ns for ordered points. This comes from the different design of `PPoly`, which keeps coefficients made beforehand, and is not changed. SciPy takes 49 ms on the same case.
  * **The construction of `BSpline1D` at degree 1**: the control points are the values themselves, and the knots are the coordinates with the ends repeated, so there is nothing to solve. `BSpline` checking these two checked parts once more was stopped (`from_checked`. 77 → 69 µs on 10 000 points, 0.32 → 0.21 µs on 100). On 10 000 points it stays 2.9 times SciPy (24 µs). The reason is the one in §3's "The differences that remain, and why": writing two newly allocated arrays (knots and control points, 160 KB) takes most of the time. SciPy does not copy the values.
* **`Nearest1D<T: Scalar, V: Value<T> = T>` (M1-7, `nearest` / `nearest-up` / `previous` / `next` of `interp1d`)**: returns the value of the nearest (or the previous or next) data point.
  * **Construction**: `Nearest1D::new(&[T], &[V], NearestKind)`. `NearestKind::{Nearest, NearestUp, Previous, Next}` (`#[non_exhaustive]`). The checks are in the order of §3.5, with at least 2 points.
  * **What it keeps**: the coordinates, values, and kind only. No `PPoly`, and no `as_ppoly()` (it is piecewise constant, but which interval a point belongs to differs from `PPoly`).
  * **Which point's value it returns**: for $x$ in the domain on the interval $[x_i, x_{i+1}]$, the value of this point.

    | Kind | Point returned |
    | :--- | :------------- |
    | `Nearest` | with the midpoint $m_i$, $x_i$ if $x \le m_i$, otherwise $x_{i+1}$ (exactly at the midpoint, the left) |
    | `NearestUp` | $x_i$ if $x < m_i$, otherwise $x_{i+1}$ (exactly at the midpoint, the right) |
    | `Previous` | $x_i$ if $x < x_{i+1}$, $x_{i+1}$ if $x = x_{i+1}$ (the largest node not above $x$) |
    | `Next` | $x_i$ if $x = x_i$, otherwise $x_{i+1}$ (the smallest node not below $x$) |

    The midpoint is found by SciPy's formula $m_i = x_i/2 + x_{i+1}/2$ (the sum is not taken first, so it does not overflow). It chooses the same point as SciPy, rounding included.
  * **Outside the domain**: the four modes. The default is `Forbid` (SciPy's `interp1d` raises by default). `Clamp` and `Extend` both return the end data value (continuing a piecewise constant gives the same value). `Periodic` wraps. A NaN query is `None`.
  * **`eval_hinted` and batches**: the interval is found by `search` and the point chosen by the table above. It is the same path as `PPoly`'s, so the results of hints and batches are bit-identical to `eval`.
* **Differences from NumPy and SciPy (M1-7)**:

  | Case | NumPy 2.5 / SciPy 1.18 | `interpx` |
  | :--- | :--------------------- | :-------- |
  | One data point | `numpy.interp` gives a constant, and `interp1d` builds too | `InsufficientPoints` (the interval search assumes at least 2 points) |
  | Repeated $x$, $x$ out of order | `numpy.interp` does not check. `interp1d` sorts, and accepts repeats | `DuplicatePoints` / `NotSorted` (§3.5) |
  | Infinite values | `numpy.interp` accepts them | `NonFinite` (§3.5) |
  | A NaN query to `Nearest1D` | the last data value | `None` |
  | Left of `previous`, right of `next`, with `fill_value="extrapolate"` | NaN | `Clamp` / `Extend` give the end data value |

  Each was confirmed by running NumPy 2.5.4 and SciPy 1.18.1.
* **`CubicSpline<T: Scalar, V: Value<T> = T>` (M1-8, SciPy's `CubicSpline`)**: the $C^2$ cubic spline.
  * **Construction**: `CubicSpline::new(&[T], &[V], BoundaryCondition<V>)`. The boundary condition is a required argument; the default not-a-knot is written `BoundaryCondition::default()`. At least 2 points.
  * **Allocations at construction**: 3 (the copy of the coordinates, the coefficients, the working space). The slopes are solved into the coefficient array, and the coefficients are written over it from the last interval forward (interval $i$ reads the slopes $s_i$, $s_{i+1}$ and writes from position $4i$ on, so slopes still needed are not overwritten). The working space holds the matrix's three diagonals, the reciprocals of the interval widths (used by both the secants and the coefficients, so found once), and the cyclic form's working array, in one array. The input is checked once, and the crate-private construction path that makes a `PPoly` from checked parts (`PPoly::from_checked`) is used. `Linear1D` uses the same path.
  * **What it keeps**: a cubic `PPoly`, lent with `as_ppoly()`. Evaluation, `eval_hinted`, and batches are those of the `PPoly` as they are (unlike `Linear1D`, there is no special case at the right end. SciPy too returns the polynomial's value at the right end).
  * **How it is found**: the slopes $s_i$ at the data points are found from a linear system, and the coefficients of interval $i$ made from the slopes. Let $h_i = x_{i+1} - x_i$ and $\delta_i = (y_{i+1} - y_i)/h_i$.
    * Inner rows ($1 \le i \le n-2$): $h_i\,s_{i-1} + 2(h_{i-1} + h_i)\,s_i + h_{i-1}\,s_{i+1} = 3\,(h_i\,\delta_{i-1} + h_{i-1}\,\delta_i)$.
    * End rows (the left end; the right end is symmetric):

      | End condition | Row | Right-hand side |
      | :------------ | :-- | :-------------- |
      | `NotAKnot` | $h_1\,s_0 + (h_0 + h_1)\,s_1$ | $\bigl((h_0 + 2d)\,h_1\,\delta_0 + h_0^2\,\delta_1\bigr)/d$, $d = h_0 + h_1$ |
      | `FirstDerivative(v)` | $s_0$ | $v$ |
      | `SecondDerivative(v)` | $2h_0\,s_0 + h_0\,s_1$ | $-\tfrac12 v\,h_0^2 + 3\,(y_1 - y_0)$ |

    * Coefficients: with $t = (s_i + s_{i+1} - 2\delta_i)/h_i$, interval $i$ is $\dfrac{t}{h_i}(x - x_i)^3 + \Bigl(\dfrac{\delta_i - s_i}{h_i} - t\Bigr)(x - x_i)^2 + s_i\,(x - x_i) + y_i$.
    * The system is solved by the tridiagonal solver (the Thomas algorithm, no pivoting). The matrix and the right-hand side are assembled as SciPy (`_cubic.py`) assembles them.
  * **Periodic boundary**: $s_{n-1} = s_0$, so the unknowns are the $n-1$ values $s_0, \dots, s_{n-2}$. Row $i$ is the formula of the inner rows, with the indices taken modulo $n-1$. This is a cyclic tridiagonal system with one entry each at the upper right and the lower left, solved by the cyclic solver (Sherman–Morrison).
  * **Few points** (treated as SciPy treats them):

    | Points | Boundary condition | Result |
    | :----- | :----------------- | :----- |
    | 2 | `NotAKnot` (at either end) | that end's condition becomes `FirstDerivative(δ_0)`. At both ends, a line |
    | 2 | `Periodic` | both ends' conditions become `FirstDerivative(δ_0)`. $y_0 = y_1$, so it is a constant |
    | 3 | `NotAKnot` at both ends | the parabola through the 3 points. The rows are $[1, 1, 0]$, $[h_1, 2(h_0 + h_1), h_0]$, $[0, 1, 1]$, the right-hand side $2\delta_0$, $3(h_0\delta_1 + h_1\delta_0)$, $2\delta_1$ (tridiagonal, so the same solver solves it) |
    | 3 | `Periodic` | 2 unknowns, outside the range of the cyclic solver. $s_0 = s_1 = s_2 = \dfrac{\delta_0/h_0 + \delta_1/h_1}{1/h_0 + 1/h_1}$ |
    | 4 or more | `Periodic` | the cyclic solver (3 rows or more) |

    Other cases (3 points with `NotAKnot` at one end only, and so on) use the general formulas.
  * **The periodicity check**: with `Periodic`, $y_0$ and $y_{n-1}$ that do not agree return `NotPeriodic`. Agreement is judged with SciPy's tolerance, column by column: $|y_0 - y_{n-1}| \le \tau\,(1 + |y_{n-1}|)$, $\tau = 4.5\,\varepsilon$ (about $10^{-15}$ for `f64`; SciPy uses `rtol = atol = 1e-15`). Exact agreement is not required, so that rounding differences such as those of `sin(0)` and `sin(2π)` are not refused. Once the check is passed, $y_{n-1}$ is used as given (the same as SciPy; it goes into $\delta_{n-2}$ of the last interval).
  * **Outside the domain**: the default is `Extend`, and `Periodic` when the boundary condition is `Periodic` (§3.3, the same as SciPy). `with_extrapolate` sets any of the four modes.
  * **Errors**: the order of the checks is the common checks of §3.5 → the values of the boundary condition's derivatives (NaN or infinity is `InvalidParameter { name: "bc", reason: ParamReason::NonFinite }`) → periodicity (`NotPeriodic`) → the solve. A breakdown of the tridiagonal solver and overflow of the slopes or coefficients return `NonFinite { index }` (`index` the data point concerned). For valid input the matrix is nonsingular, so a breakdown happens only by overflow. No variant of its own is added (treated like the overflow of `Linear1D`'s slopes).
  * **Differences from SciPy**: only those of the common checks (repeated or unsorted $x$, values that are not finite). The minimum number of points (2), the handling of few points, and the tolerance of periodicity are SciPy's.
* **Making a cubic `PPoly` from values and slopes (the cubic Hermite form; separated in M1-9)**: once the value $y_i$ and slope $s_i$ at each data point are decided, `CubicSpline`, `PchipInterpolator`, `Akima1D`, and `CubicHermiteSpline` become a cubic `PPoly` by the same procedure. They differ only in how the slopes are decided. The procedure is one crate-private function, shared by the four types.
  * The coefficients of interval $i$ are, with $\delta_i = (y_{i+1} - y_i)/h_i$ and $t = (s_i + s_{i+1} - 2\delta_i)/h_i$, $\dfrac{t}{h_i}(x - x_i)^3 + \Bigl(\dfrac{\delta_i - s_i}{h_i} - t\Bigr)(x - x_i)^2 + s_i\,(x - x_i) + y_i$ (the same formula as SciPy's `CubicHermiteSpline`).
  * It is called with the slopes at the front of the coefficient array (4 per interval), and writes the coefficients over it from the last interval forward (the form put into `CubicSpline` in M1-8). It adds no allocation. The caller has done the input checks, and the `PPoly` is made with `PPoly::from_checked`.
  * Coefficients that are not finite give `NonFinite { index }` (`index` the first interval that overflowed).
* **`PchipInterpolator<T: Scalar, V: Value<T> = T>` (M1-9, SciPy's `PchipInterpolator`)**: a shape-preserving $C^1$ cubic interpolation. On intervals where the data is monotone the interpolant is monotone too, and it does not overshoot the range of the data.
  * **Construction**: `PchipInterpolator::new(&[T], &[V])`. The checks are in the order of §3.5, with at least 2 points.
  * **What it keeps**: a cubic `PPoly`, lent with `as_ppoly()`. Evaluation, `eval_hinted`, and batches are the `PPoly`'s.
  * **Slopes** (the same as SciPy. $h_i = x_{i+1} - x_i$, $m_i = (y_{i+1} - y_i)/h_i$):
    * **Inner points** $x_i$: if $m_{i-1}$ and $m_i$ differ in sign or either is 0, $s_i = 0$. Otherwise the weighted harmonic mean $\dfrac{1}{s_i} = \dfrac{1}{w_1 + w_2}\Bigl(\dfrac{w_1}{m_{i-1}} + \dfrac{w_2}{m_i}\Bigr)$, $w_1 = 2h_i + h_{i-1}$, $w_2 = h_i + 2h_{i-1}$ (Fritsch–Butland).
    * **Ends** (the left end; the right end is symmetric): the one-sided three-point formula $d = \dfrac{(2h_0 + h_1)\,m_0 - h_0\,m_1}{h_0 + h_1}$ is found and limited to preserve the shape. If $d$ differs in sign from $m_0$, $s_0 = 0$. Otherwise, if $m_0$ and $m_1$ differ in sign and $|d| > 3|m_0|$, $s_0 = 3m_0$. Otherwise $s_0 = d$ (Moler's `pchiptx`).
    * **2 points**: $s_0 = s_1 = m_0$ (a line).
    * Signs are compared with three values, positive, negative, and 0, as SciPy's `sign`.
  * **With values `[T; M]`**: the rules above are applied to each column independently (signs can differ by column), with the "column-wise operations inside the crate" of §3.1.
  * **Outside the domain**: the default is `Extend` (§3.3, the same as SciPy). `with_extrapolate` sets any of the four modes. Monotonicity and "no overshoot" are properties within the domain and do not extend to the end cubics continued.
  * **Errors**: the common checks of §3.5. Overflow of the secants, slopes, or coefficients is `NonFinite { index }`.
  * **Differences from SciPy**: only those of the common checks (repeated or unsorted $x$, values that are not finite).
* **`Akima1D<T: Scalar, V: Value<T> = T>` (M1-10, SciPy's `Akima1DInterpolator`)**: a $C^1$ cubic interpolation. The slope at each point is decided from the five nearby points only, so one outlying point does not reach far, and the waviness of splines is less likely.
  * **Construction**: `Akima1D::new(&[T], &[V], AkimaMethod)`. `AkimaMethod::{Akima, Makima}` (`#[non_exhaustive]`; SciPy's `method="akima"` / `"makima"`). The checks are in the order of §3.5, with at least 2 points.
  * **What it keeps**: a cubic `PPoly`, lent with `as_ppoly()`. The coefficients are made by the common procedure (the cubic Hermite form).
  * **Slopes** (the same as SciPy): the secants $m_i = (y_{i+1} - y_i)/h_i$ ($i = 0, \dots, n-2$) are extended by two lines beyond each end: $m_{-1} = 2m_0 - m_1$, $m_{-2} = 2m_{-1} - m_0$, $m_{n-1} = 2m_{n-2} - m_{n-3}$, $m_n = 2m_{n-1} - m_{n-2}$. The slope at point $x_i$ is decided from the two secants before it, $m_{i-2}, m_{i-1}$, and the two after it, $m_i, m_{i+1}$.
    * Weights: $f_1 = w(m_i, m_{i+1})$, $f_2 = w(m_{i-2}, m_{i-1})$. `Akima` has $w(a, b) = |b - a|$, `Makima` has $w(a, b) = |b - a| + \tfrac12 |a + b|$.
    * Points where $f_1 + f_2$ exceeds the cut-off threshold: $s_i = m_{i-1} + \dfrac{f_2}{f_1 + f_2}\,(m_i - m_{i-1})$.
    * Points where it does not (the two secants before are equal, and the two after are equal; the formula does not decide the slope): $s_i = \tfrac12\,(m_{i-2} + m_{i+1})$.
    * The threshold: $10^{-9}$ times the largest $f_1 + f_2$ over all points (SciPy's `break_mult`).
    * **2 points**: $s_0 = s_1 = m_0$ (a line).
  * **A property of `Makima`**: the slopes are 0 at both ends of a flat stretch (3 or more equal values in a row), and a flat stretch stays flat. `Akima` has a jump in its treatment of points where two secants are equal (the threshold only moves where it is, and cannot remove it), and can overshoot next to a flat stretch.
  * **With values `[T; M]`**: the slopes and the threshold are decided per column, and each column is interpolated independently. The result is bit-identical to building one column at a time.
  * **Outside the domain**: the default is `Forbid` (§3.3; SciPy's Akima does not extrapolate by default). `with_extrapolate` sets any of the four modes.
  * **Errors**: the common checks of §3.5. Overflow of the secants, slopes, or coefficients is `NonFinite { index }`.
  * **Differences from SciPy**: besides those of the common checks (repeated or unsorted $x$, values that are not finite), there is one. With several columns of values, SciPy takes the largest value used for the threshold **over all columns together**, as one (`xp.max(f12)`). So the result of one column depends on the size of the values of the other columns, and can differ from interpolating one column at a time. `interpx` takes the largest value per column. The result changes only at points where a column's $f_1 + f_2$ exceeds $10^{-9}$ times that column's largest value but not $10^{-9}$ times the largest over all columns. For scalar values it is the same as SciPy.
* **`CubicHermiteSpline<T: Scalar, V: Value<T> = T>` (M1-11, SciPy's `CubicHermiteSpline`)**: a $C^1$ cubic interpolation given both the value and the first derivative (slope) at each data point. It has no rule for the slopes and uses the slopes given as they are. For data whose derivatives are known (the solution of a differential equation and the like), or to decide the slopes by a rule of one's own.
  * **Construction**: `CubicHermiteSpline::new(x: &[T], y: &[V], dydx: &[V])`. The slopes have the value type `V`; with values `[T; M]` they are slopes per column. At least 2 points (the same as SciPy; 2 points give one cubic).
  * **What it keeps**: a cubic `PPoly`, lent with `as_ppoly()`. Evaluation, `eval_hinted`, and batches are the `PPoly`'s. The coefficients are made by the common procedure (the cubic Hermite form). The slopes are copied to the front of the coefficient array and passed, so construction allocates 3: the copy of the coordinates, the coefficients, and the reciprocals of the interval widths.
  * **Reproducing the values and slopes given**: the constant and linear coefficients of interval $i$ are $y_i$ and $s_i$ themselves (bit for bit). So at the left end of each interval the value equals $y_i$ and the first derivative $s_i$. At the right end (the left end of the next interval too), the value and derivative of the left interval's cubic agree with them to rounding.
  * **With values `[T; M]`**: each column is interpolated independently. The result is bit-identical to building one column at a time.
  * **Outside the domain**: the default is `Extend` (§3.3, the same as SciPy). `with_extrapolate` sets any of the four modes.
  * **Periodicity is not checked**: `Periodic` is an extrapolation mode, and this type has no boundary condition. Putting `Periodic` on data whose values or slopes do not agree at the two ends makes the values or slopes jump at the seam of the period (the same as SciPy's `extrapolate='periodic'`, `PchipInterpolator`, and `Akima1D`). Only `CubicSpline`, which chooses the period as a boundary condition, returns `NotPeriodic`.
  * **Errors**: the checks of §3.5, in the order that includes the slopes. Slopes of another length are `LengthMismatch`, and NaN or infinity in the slopes `NonFinite { index }` (`index` the position in `dydx`). `NonFinite` does not carry which array, so whether it is coordinates, values, or slopes follows from the order of the checks (values if every coordinate is finite, slopes if every value is finite too). No variant to tell the arrays apart is added. Overflow of the secants or coefficients is `NonFinite { index }`.
  * **Differences from SciPy**: only those of the common checks (repeated or unsorted $x$). SciPy too refuses values and slopes that are not finite.
* **`BSpline1D<T: Scalar, V: Value<T> = T>` (M1-12, SciPy's `make_interp_spline`)**: the B-spline of degree $k$ through the data points. Its derivatives up to order $k - 1$ are continuous.
  * **Construction**: `BSpline1D::new(x: &[T], y: &[V], degree: usize, bc: BoundaryCondition<V>)`.
  * **What it keeps**: a `BSpline` (knots and control points), lent with `as_bspline()`. Evaluation, `eval_hinted`, batches, derivatives, and integrals are the `BSpline`'s.
  * **Degree**: 1 to 15 (`BSPLINE_MAX_DEGREE`). 16 or more is `InvalidParameter { name: "degree", reason: TooLarge { max: 15 } }`, and 0 is `InvalidParameter { name: "degree", reason: TooSmall { min: 1 } }`.
    * **Why there is no degree 0**: SciPy's degree 0 (knots $x_0, \dots, x_{n-1}, x_{n-1}$, control points $y_i$) puts the value of the last data point on an interval of zero width. `BSpline` does not evaluate intervals of zero width ("Domain, and which interval a point is in" above), so it would return $y_{n-2}$ and not $y_{n-1}$ at the right end, and not pass through the data. The staircase that returns the value of the point on the left is `Previous` of `Nearest1D` (with $y_{n-1}$ at the right end).
  * **Combinations of boundary condition and degree** (SciPy's `bc_type`; in SciPy too, `"natural"` and `"clamped"` can be used for cubics only):

    | `bc` | Degrees | SciPy |
    | :--- | :--- | :--- |
    | `NotAKnot` (`Ends(NotAKnot, NotAKnot)` is the same) | 1–15 | `None` / `"not-a-knot"` |
    | `Periodic` | 1–15 | `"periodic"` |
    | `Ends` with `FirstDerivative` / `SecondDerivative` at either end | 3 only | `"natural"`, `"clamped"`, `([(1, a)], [(2, b)])` |

    Derivatives at an end for a degree other than 3 are `InvalidParameter { name: "bc", reason: Unsupported }`. Giving several derivatives per end (SciPy's `([(1, a), (2, b)], ...)`; a quintic needs 2 per end) and giving the knot vector (SciPy's `t`) are out of scope for v0.1.0 (roadmap A-10).
  * **Knot vector** (the same as SciPy. The domain is $[x_0, x_{n-1}]$, and the knots at both ends have multiplicity $k + 1$):
    * **Degree 1**: the inner knots are $x_1, \dots, x_{n-2}$. The control points are the $y_i$ themselves, and no linear system is solved (a polyline).
    * **`NotAKnot`, odd degree**: the inner knots are the data points without $(k+1)/2$ at each end ($x_{(k+1)/2}, \dots, x_{n-1-(k+1)/2}$).
    * **`NotAKnot`, even degree**: the inner knots are the midpoints of neighboring data points $(x_i + x_{i+1})/2$ ($n - 1$ of them), without $k/2$ at each end. The midpoint is found as $x_i/2 + x_{i+1}/2$ (the sum does not overflow; where it would not, the result is the same).
    * **An end given derivatives, at degree 3**: the data point next to that end stays a knot too (with derivatives at both ends, the inner knots are $x_1, \dots, x_{n-2}$). At a `NotAKnot` end, the data point next to it ($x_1$ or $x_{n-2}$) is left out.
    * **`Periodic`, degree 2 or more**: the knots within the domain are $x_0, \dots, x_{n-1}$ for odd degrees. For even degrees the two ends stay, and the inner $x_i$ ($1 \le i \le n-2$) are replaced by $(x_{i-1} + x_i)/2$. Beyond them, $k$ on each side, the spacing of the knots is repeated periodically (SciPy's `_periodic_knots`). The end knots are not repeated.
    * **Even degrees, with two neighboring points too close**: even degrees put a knot between two points. When two points are so close that they are neighbors in floating point, the midpoint is rounded to one of them and the knot falls on a data point (the system can be singular; the interval of the coinciding knots has zero width and its basis functions cannot be evaluated). If there is an interval where the knot does not fall strictly between its two points, `DuplicatePoints { first: i, second: i + 1 }` is returned (two points this scheme cannot tell apart). With `NotAKnot` every interval is checked, including those whose midpoint is not used as a knot (to keep one rule). With `Periodic` the intervals that get a knot (all but the last) are checked. Odd degrees and degree 1 have the data points themselves as knots, so this does not happen.
    * **Knots that overflow**: periodic knots extend outside the domain, so knots can overflow even with finite data. Knots that are not finite, including where the difference of two neighboring points overflows, are `NonFinite { index }` (`index` the nearer data point).
  * **The linear system** (degree 2 or more): the unknowns are the control points. The rows are the interpolation conditions at each data point, $\sum_j B_j(x_i)\,c_j = y_i$, and the conditions of the derivatives at the ends, $\sum_j B_j^{(\nu)}(x_0)\,c_j = v$. The rows of the left end's derivatives come first, the right end's last. In this order the matrix is banded with lower and upper bandwidths within $k$, and the banded LU above solves it.
    * **The width of the band**: without rows of derivatives, a band of $k - 1$ on each side is enough. A point lies in an interval of $k + 1$ basis functions, but only the first and the last rows reach the outermost diagonals, and the entries there are basis functions that are 0 at the end of the interval (exactly 0 by the recurrence). Leaving those two out, the band is declared as $k - 1$. With rows of derivatives (cubic) it stays $k$. The cost of elimination is proportional to the square of the bandwidth, so declaring it narrower is faster (`benches/bspline1d.rs`).
  * **The linear system of `Periodic`**: on a periodic knot vector the control points are periodic too ($c_{j + n - 1} = c_j$). Reducing the unknowns to $n - 1$ and lining up the interpolation conditions at $x_0, \dots, x_{n-2}$ gives a cyclic banded matrix (lower and upper bandwidths $\lfloor k/2 \rfloor$). For odd degrees $x_i$ is a knot, so the last of the $k + 1$ basis functions of the interval is exactly 0, and leaving it out makes the band symmetric (cyclic tridiagonal for cubics). If $n - 1 > 2\lfloor k/2 \rfloor$, it is solved by the Woodbury formula above; otherwise the band wraps around onto itself, and it is solved as a dense matrix with the overlapping entries added.
    * For $n \le k$, SciPy solves another dense matrix ($n + k - 1$ unknowns) that lines up conditions of continuity of the derivatives. The solution is the same (the PR of M1-12 checks every combination of degrees 1 to 15 and 2 points or more against SciPy). To keep one way of solving, only the form above is used here.
  * **The minimum number of points**: the interpolation conditions and the derivative conditions together must be the $k + 1$ that decide one polynomial of degree $k$. That is, $n \ge k + 1 - (\text{the number of ends given derivatives})$, and at least 2 points. `NotAKnot` needs $k + 1$ points (exactly that gives one polynomial of degree $k$), a cubic with derivatives at both ends 2 points, and `Periodic` 2 points whatever the degree (2 points give a constant). Too few is `InsufficientPoints { required, actual }`.
    * `CubicSpline` makes 3 points with not-a-knot a parabola and 2 points a line (the same as SciPy's `CubicSpline`). `BSpline1D` does not accept them, as SciPy's `make_interp_spline` does not.
  * **A cubic `BSpline1D` and `CubicSpline`**: with the same data and the same boundary condition (satisfying the numbers of points above), they give the same curve (to rounding). They differ in the representation (B-spline and piecewise polynomial) and in how they solve (banded LU and the Thomas algorithm), and `CubicSpline` builds faster. The combination with `NotAKnot` at one end only is not in SciPy's `make_interp_spline`, but is accepted so that the same `BoundaryCondition` as `CubicSpline`'s is accepted with the same meaning.
  * **The periodicity check**: with `Periodic`, the agreement of $y_0$ and $y_{n-1}$ is checked with `CubicSpline`'s tolerance, and `NotPeriodic` if they do not agree.
  * **With values `[T; M]`**: the matrix does not depend on the column, so the result is bit-identical to solving each column independently.
  * **Outside the domain**: the default is `Extend`, and `Periodic` when the boundary condition is `Periodic` (§3.3). `with_extrapolate` sets any of the four modes.
  * **Errors**: the order of the checks is step 1 of §3.5 (lengths) → the degree (`InvalidParameter`) → the combination of boundary condition and degree (`InvalidParameter`) → step 2 (the number of points; the minimum as above) → steps 3 and 4 (not finite, the order of coordinates) → the values of the boundary condition's derivatives (`InvalidParameter { name: "bc", reason: NonFinite }`) → periodicity (`NotPeriodic`) → the knots (`DuplicatePoints` for even degrees, `NonFinite` for overflow) → the solve. A breakdown of the solver and overflow of the control points are `NonFinite { index }` (`index` the row concerned; the same treatment as `CubicSpline`). Even with finite data, coordinates so far apart that their differences overflow make the basis functions impossible to evaluate, and this error comes as a breakdown of the solver.
    * **Is the matrix nonsingular**: the matrices of `NotAKnot` and of derivatives at the ends are nonsingular because the knots and data points satisfy the Schoenberg–Whitney conditions. The banded matrix $A$ solved by the Woodbury formula is a square submatrix of the collocation matrix and satisfies the same conditions. For the whole periodic matrix (that $I + V^{\mathsf T} Z$ is nonsingular), SciPy too assumes it, but no proof is held here. The condition numbers of 20 000 random data sets of degrees 2 to 8, 2 to 13 points, and ratios of neighboring spacings up to $10^3$ were examined, and no case approaching singularity was found (the worst $2.5 \times 10^9$, for degree 8; it grows smoothly with the ratio of the spacings). If a singular case is found, this is changed to return `SolverFailed(SolverError::Singular)` (§4; a variant planned for M3).
  * **Differences from SciPy**: besides those of the common checks (repeated or unsorted $x$, values that are not finite), the following.

    | Item | SciPy 1.18 | `interpx` |
    | :--- | :--- | :--- |
    | Degree 0 | returns a staircase | `InvalidParameter` (the reason above) |
    | `Periodic` at degree 1 | extrapolation is `True` (not periodic) | extrapolation is `Periodic` (as for degree 2 or more) |
    | One point at degree 1 | accepted (the domain has zero width) | `InsufficientPoints` (2 points or more) |
    | Not-a-knot at one end only, at degree 3 | `make_interp_spline` does not accept it | accepted (the same curve as `CubicSpline`) |
    | Even degree, with the midpoint of two neighboring points rounded to one of them | an exception for a singular matrix, or NaN | `DuplicatePoints` |
    | Even degree, with the sum for a midpoint overflowing | the knot becomes infinite | found in a form that does not overflow |
    | Several derivatives per end, a given knot vector | accepted | not provided (A-10) |
* **`BSplineCurve<T: Scalar, const D: usize>` (M1-13)**: the parametric curve $u \mapsto [T; D]$ through points of $D$ dimensions. It is a `BSpline1D` with values `[T; D]`, plus a way of deciding the parameters.
  * **The relation to `BSpline1D`**: SciPy has no class for curves, and `make_interp_spline(t, points)` is given a two-dimensional `y`. That is the same as `BSpline1D<T, [T; D]>`, with the caller providing the parameter $t_i$ of each point. A curve usually has only the points and no $t_i$, so `BSplineCurve` decides $t_i$ from the points. It keeps a `BSpline1D<T, [T; D]>` inside, and the degree, boundary condition, knot vector, linear system, and errors are those of `BSpline1D`.
  * **Construction**:
    * `BSplineCurve::new(points: &[[T; D]], degree: usize, bc: BoundaryCondition<[T; D]>, parametrization: Parametrization)`: decides the parameters from the points.
    * `BSplineCurve::with_parameters(t: &[T], points: &[[T; D]], degree: usize, bc: BoundaryCondition<[T; D]>)`: is given the parameters (SciPy's `make_interp_spline(t, points)`). The result is bit-identical to a `BSpline1D` with the same arguments.
  * **`Parametrization`** (`#[non_exhaustive]`; it has no default and is an explicit argument): with the distance between neighboring points $d_i = \lVert p_{i+1} - p_i \rVert_2$,

    | Value | $t_i$ (before normalization) | Property |
    | :--- | :--- | :--- |
    | `Uniform` | $i$ | ignores the spacing of the points. Wavy on unevenly spaced points |
    | `ChordLength` | $\sum_{j < i} d_j$ | approximates the length along the curve. The same as the default of FITPACK's `splprep` |
    | `Centripetal` | $\sum_{j < i} \sqrt{d_j}$ | overshoots less at sharp turns |

    Each is normalized to $[0, 1]$ by dividing by the last value ($t_0 = 0$, $t_{n-1} = 1$). The sums are taken front to back.
  * **Computing the distances**: $d_i$ divides the difference by its component of largest absolute value before taking the sum of squares (the squares in between neither overflow nor underflow). So $d_i = 0$ only when the two points are exactly the same.
  * **What it keeps**: the `BSpline1D<T, [T; D]>`, and the parameter of each data point, returned by `parameters() -> &[T]` (when they are decided automatically, there is no other way to know which $u$ is which point). `as_bspline()` lends the `BSpline<T, [T; D]>`. Tangents (first derivatives) and second derivatives are taken with `as_bspline().eval_derivative(u, nu)`; the curve type does not carry the same function again.
  * **Closed curves**: `BoundaryCondition::Periodic`. A first and last point that do not agree (with `BSpline1D`'s tolerance) are `NotPeriodic`.
  * **Outside the domain**: as for `BSpline1D` (the default `Extend`, `Periodic` when periodic). The domain of the parameter is $[0, 1]$ for `new` and $[t_0, t_{n-1}]$ for `with_parameters`.
  * **Errors**: the order of the checks of `new` is the degree → the combination of boundary condition and degree → the number of points (the minimum as for `BSpline1D`) → points not finite (`NonFinite { index }`) → the parameters → then `BSpline1D`'s.
    * With `ChordLength` and `Centripetal`, two neighboring points that are the same are `DuplicatePoints { first: i, second: i + 1 }`. With `Uniform` the parameters do not coincide, so repeated points pass.
    * Even when two points are not the same, if they are too close compared with the whole length, the normalized parameters can coincide by rounding. This too is `DuplicatePoints` (by the check of `BSpline1D`'s coordinates).
    * Even with finite points, points so far apart that differences or sums overflow make the parameters not finite, `NonFinite { index }`.
    * The errors of `with_parameters` are those of `BSpline1D::new`.
  * **What it does not have**: arc length, reparametrization by arc length, the search for the nearest point. SciPy has nothing corresponding, and they are outside this design. Fitting with smoothing (`splprep` with $s > 0$) belongs to the FITPACK family (roadmap A-7).
  * **Differences from SciPy**: `with_parameters` is as `BSpline1D` (the table of §6.1). `make_interp_spline` has no function corresponding to `new`. The formulas of the parameters are as above; for the goldens, the generating script computes $t_i$ by the same formulas and passes them to `make_interp_spline`.

### 6.2 Geometric search for scattered data (with `convx`)
* **Delaunay decomposition**: `[[f64; D]]` is flattened to `&[f64]` (row-major) without a copy and passed to `convx`.
* **Duplicate points**: `convx`'s `representative()` merges identical coordinates into the smallest index. `NearestND` / `LinearND` / `CloughTocher2D` **accept them without an error**, and use the value of the point with the smallest index (stated in the doc). `RBFInterpolator` returns `Err(DuplicatePoints)`.
* **How a point is located** (M4-1, #180; how the walk starts is the owner's decision):
  1. **Where the walk starts**: the box of the points is cut into $p$ cells along each axis, $p$ the largest integer with $p^D$ no more than twice the number of points, $2n$: there are no more than twice as many cells as points ($p = 1$ in high dimensions). It was $\lceil n^{1/D} \rceil$ at first, which made $p \ge 2$ however few the points, and $2^D$ cells: a million cells for 23 points in 20 dimensions, 0.37 s to build (the review of #183 measured 3.5 s at 22 dimensions; fixed, 0.010 s at 20 and 0.017 s at 22). With $n$ as the bound, the cells at 100 points in 3-D went from 125 to 64 and evaluation was 13% slower, so the bound is $2n$ (the suggestion of the review of #183; measured in three trees in turn, three rounds: against $n$, evaluation 5 to 15% faster, construction 0 to 5% slower; 125 cells at 100 points in 3-D; the owner's decision). At construction, a walk goes to the middle of each cell from the start of the cell before it, and the last simplex of that walk (the one that holds the middle, or one on the hull toward it) is the start of the cell. The cells are visited in the order of a snake: the last axis fastest, each axis turning back at its end instead of jumping to its start; jumping back walked across the triangulation once a row. In P-21 this took construction from 57.3 to 54.5 ms at 100 000 points and from 445 to 410 µs at 1000; evaluations are the same bit for bit. A query walks from the start of the cell it is in (the nearest cell, just outside the box). The start depends on the query's coordinates alone, so `eval`, a batch, and a parallel batch give the same results bit for bit.
     * The design first had a jump: the nearest point found with `kiddo`, and a simplex at it. In the comparison, evaluation at the pixels of an image took 1.6 to 3.6 times SciPy's (SciPy walks from the simplex it found last, and hardly walks from one pixel to the next). The time was the search of the tree for every query (about 100 ns at 1000 points). With cells, evaluation is 5 to 9 times faster (at 100 000 points, 245 → 28 ns a query in 2-D, 480 → 79 ns in 3-D), and construction did not change (70.7 ms and 69.7 ms at 100 000 points in 2-D). The losses in evaluation in the comparison are gone.
     * As many cells as simplices make evaluation a little faster (23 ns at 100 000 points in 2-D), but made construction in 3-D a fifth longer, so there are as many as points.
     * Walking from the last simplex, as SciPy does, is not taken: a point on a face is found in either simplex, so how a batch is cut would change the simplex found, and the bits of the result.
  2. **The walk**: across the face opposite the first vertex whose barycentric coordinate is below `-eps` (SciPy's `_find_simplex_directed`). At most as many steps as there are simplices. When a coordinate is NaN or above `1 + eps` (a degenerate simplex in the way), or the steps run out, every simplex is looked at.
  3. **Outside the hull**: a walk that would cross a face of the hull is outside it, unless the coordinate is `-eps_broad` or above (within rounding): then every simplex is looked at. In the property test, a data point itself was measured at $-2.8 \times 10^{-14}$ from one simplex, and taken as outside (the rounding of the transform). On the same 34 points SciPy gave every data point a value (another triangulation and other transforms do not round the same way; the same after any other point).
  4. **Barycentric coordinates and tolerances** (SciPy's rules for what is inside): for each simplex, the LU decomposition with partial pivoting $PA = LU$ of the matrix $A$ of the other vertices less the last one is made at construction, and $M = L^{-1}P$ and $U$ (its diagonal replaced by the reciprocals) are kept. A query solves $U c = M (x - v_D)$ by back substitution. SciPy keeps $A^{-1}$ and multiplies by it. $A^{-1}$ holds the large numbers of $U^{-1}$, so on a thin simplex its rounding grows with the condition number (at 20 000 points inside a tetrahedron with a condition number of about $10^4$, the error of an affine function is $1.8 \times 10^{-12}$ with $A^{-1}$ and $1.8 \times 10^{-15}$ with the solve; the value at a data point at its vertex was $1.5 \times 10^{-12}$ off, which the property test found on CI; the owner's decision). It is 1.5 times the multiplications of $A^{-1}$, with no division and no reordering; evaluation takes 1.3 to 1.75 times as long as with $A^{-1}$ ("Measured speed" below). Inside is `-eps <= c_i <= 1 + eps`, `eps = 100 * f64::EPSILON`. A simplex whose reciprocal condition number in the 1-norm is below `1000 * f64::EPSILON` is degenerate (SciPy uses LAPACK's estimate; here it is computed from the inverse found from the LU at construction). In the look at every simplex, next to a degenerate simplex, a leeway of `eps_broad = sqrt(f64::EPSILON)` is allowed toward it (SciPy's `_find_simplex_bruteforce`). A point outside the box of the points by more than `eps` is outside without a walk.
  5. **The value**: the vertices' values summed with the barycentric coordinates, in SciPy's order (from zero, the vertices in turn).
* **Boundary sentinel**: `NO_NEIGHBOR = u32::MAX`. `convx`'s rustdoc promises this value, so it is defined in one place inside `interpx`. `convx` is not asked for a public constant (the owner's decision, #180).
* The order of `convx`'s `simplices` is not guaranteed across versions, so nothing depends on re-triangulation (§6.5).

* **`LinearNDInterpolator<const D: usize, V: Value<f64> = f64>` (M4-1, #180, SciPy's `LinearNDInterpolator`)**: linear on each simplex of the Delaunay triangulation of the points. `Point = [f64; D]`, `Output = V`.
  * **Construction**: `new(points, values)`, `from_owned(points: Vec<..>, values: Vec<V>)` (moves the points and the values), `rescaled` (SciPy's `rescale=True`; triangulates after the change of `NearestND`), `rescaled_from_owned` (changes the moved points in place). It keeps the points and the values. `points()` returns the points as triangulated (after the change, as SciPy's `points`).
  * **`with_simplex_tolerance(t)`** (SciPy 1.18's `simplex_tolerance`; the owner's decision in M4-1): the tolerance of what is inside a simplex becomes `100ε × t`. SciPy takes it at each call; `Interpolator::eval` has no place for an argument, so the interpolant holds it, as `NearestND` holds its distance bound. `eps_broad` does not change (as in SciPy). A `t` that is not finite and above zero is `InvalidParameter` (SciPy: `ValueError`). In the goldens (36 queries moved out of the edges of the hull by 1e-14 to 1e-11, with the default and 100 times it), which queries have a value is SciPy's (1 with the default, 15 with 100 times). It is stored.
  * **The order of the checks**: length → number of points ($D + 1$ or more) → the limit of the number (`u32`) → coordinates → values → (`rescaled` only) mean and range → the triangulation (fewer than $D + 1$ distinct points, a lower dimension, the working space of the exact arithmetic).
  * **Zero and one dimension**: with zero dimensions all points are one point, and the first value is returned (`convx` is not called; as `NearestND`). With one, it is linear interpolation between neighboring points. SciPy refuses one-dimensional points (`ValueError: Need at least 2-D data`).
  * **Agreement with SciPy** (21 golden cases, `compare/goldens/scattered.json`; 2, 3, and 4 dimensions): the largest relative difference is $2.6 \times 10^{-15}$ (the tolerance $10^{-13}$). In the case of a grid with points taken out, four or more points lie on one circle, and 20 of 64 queries differ (regions divided otherwise). There too an affine function comes back exactly from both (a test).
  * **Differences from SciPy**: how points on one circle (sphere) are split (above). Of points with the same coordinates, SciPy (Qhull) keeps one that is not always the first (6 of 25 sets in the golden case of duplicates); interpx takes the value of the first. One dimension is taken. A point beyond a face of the hull by rounding is looked for in every simplex (3 above). The barycentric coordinates are solved with an LU decomposition, not $A^{-1}$ (4 above; three orders more accurate inside a thin simplex).
  * **Measured speed** (criterion, this machine, uniform points of the unit square and cube, 1000 queries): evaluation 16 to 27 ns a query in 2-D, 31 to 102 ns in 3-D (100 to 100 000 points; one run after P-20. One criterion run moves a lot from run to run on this machine, so it is not compared with the figures before; what P-20 gained is told by the measurement in turn below). Construction at 100 000 points: 56 ms in 2-D, 512 ms in 3-D (after P-21's `convx` 0.1.1 and snake order of the cells; before, 79 ms and 566 ms). With $A^{-1}$, evaluation was 16 to 28 ns in 2-D and 25 to 79 ns in 3-D, and construction 71 ms and 484 ms.
    * **The walk is short**: on the comparison's cases at the pixels of an image, 1.76 to 1.84 steps a query (1 is found in the starting simplex). Walking from the last simplex, as SciPy does, would save at most 0.8 of a step.
    * **A form tried and not taken**: deciding the walk by $A^{-1}$ and solving with the LU only in the simplex found. Measured in two trees in turn, three rounds: 0 to 15% faster (about the same at 1000 points in 2-D, a tenth at 100 points in 3-D), with 30 numbers a simplex in 3-D instead of 21. Not taken (the owner's decision). That most of an evaluation's time is outside the walk was read, and P-20's profile proved it wrong (below).
    * **Where an evaluation's time goes** (P-20, #184; Visual Studio's CPU profiler, evaluations repeated on the data of the comparison's cases): at 1000 points and 250 000 pixels in 2-D (22 ns a query) the walk took 58%, `locate` 10%, and the rest (the check of the coordinates, the box, the cell, the sum of the values) about 30%. At 100 000 points in 3-D with 10 000 queries in no order (231 ns) the walk took 89%. The walks are 1.63 and 2.73 steps on average (counted with a temporary counter). A step in 3-D takes about 75 ns, long for the barycentric coordinates (about 15 multiplications), and most of it is read to be the wait to read each step's new simplex's LU decomposition ($M$ and $U$) from memory (the cache misses were not counted).
      * **The change kept**: the `floor` that finds the cell's index was left as a call of a library function (2% of the profile). It is replaced by a conversion to an integer, which saturates (below zero to zero, too large to the largest), so with the clamp to the last cell the start is the one before. Evaluations are the same bit for bit (checked on 300 000 queries in 2-D, 100 000 in 3-D, and the data points). Against the version before in turn, 3 rounds, one core, best against best: `LinearND` 22.7 → 21.1 ns in 2-D, 227 → 216 ns in 3-D, `CloughTocher2D` 48.1 → 44.5 ns. These 5 to 7% are this machine's (Windows, the default target). On the machine of #193's review (Linux on x86-64, the default target) it was 15 to 20% faster in 2-D and 10 to 20% in 3-D. With the default target (no SSE4.1) `floor` is a call of a function, and where it becomes the one instruction `roundsd` (`target-cpu=native`, for one) the gap should be smaller (not measured). The gain is larger than `floor`'s share of the profile (2% in 2-D, 0.3% in 3-D). There may be a part that a share does not show, such as what the call does to the order of the instructions around it, but the cause has not been confirmed.
      * **Forms tried and not taken**: inlining `walk` and `locate` made it slower, 21.5 → 24.2 ns in 2-D and 221 → 253 ns in 3-D. Running the check of the walk by an index instead of `chain` made 2-D faster, 21.5 → 19.0 ns, and 3-D slower, 221 → 236 ns (the same way in all 7 rounds); a form that makes 3-D slower is not kept. Sorting the queries in the order of space hardly changed the time in 3-D, 228 → 221 ns (10 000 queries against about 670 000 simplices do not read a simplex again however they are sorted, so this does not tell the share of the cache).
  * **70% of construction is `convx`'s triangulation, and `convx` 0.1.1 is about level with `spade`'s construction** (P-21; 2-D, uniform points, one process, one core, medians of 5 runs; `convx` with the copy):

    | points | interpx, whole | of it `convx` | `spade` (`bulk_load`) |
    | ---: | ---: | ---: | ---: |
    | 10 | 2.35 µs | 1.69 µs | 1.69 µs |
    | 100 | 25.9 µs | 20.3 µs | 18.5 µs |
    | 1000 | 427 µs | 306 µs | 278 µs |
    | 10 000 | 5.12 ms | 3.72 ms | 3.51 ms |
    | 100 000 | 55.7 ms | 39.0 ms | 43.1 ms |

    * With `convx` 0.1.0 it was 4.4 µs at 10 points, 475 µs at 1000, and 59.1 ms at 100 000, and `convx` alone was longer than `spade`'s whole construction. `convx`'s R2-1 (YUKIKEDA/convx#667) made it faster.
    * **The rest (interpx's own part) is about 17 ms at 100 000 points**: a temporary build with a clock around each stage gave the LU decomposition and condition number of every simplex about 7 ms (about 200 000 simplices) and the starts of the cells about 11 ms (about 200 000 cells). The snake order of the cells took about 3 ms off the whole construction (the stages were not timed again). This part is what construction loses to `spade` in the comparison. Both are tens of nanoseconds apiece: the LU decomposition is there for the accuracy of evaluation and the cells for its speed (both the owner's decisions). The policy is to fit a use in which evaluations far outnumber constructions, as in §3's "Making coefficients" (confirmed with the owner in P-21). Against solving by Cramer's rule from the three vertices at every evaluation, as `spade` does, the break-even was estimated at queries a few times the number of points (what solving in advance costs was measured, about 7 ms at 100 000 points; the added cost of solving at each evaluation was not). Another way of computing the condition number would change which simplices are taken as degenerate, so it is not changed.


* **`NearestNDInterpolator<const D: usize, V: Value<f64> = f64>` (M3-1, SciPy's `NearestNDInterpolator`)**: returns the value of the data point nearest to the query point. `Point = [f64; D]`, `Output = V`. The coordinates are `f64` (the rule of §3.1 for scattered data).
  * **Construction**: `new(points: &[[f64; D]], values: &[V])` and `rescaled(points, values)`. The points and values are copied and kept.
    * `rescaled` corresponds to SciPy's `rescale=True`. Per axis, the mean is subtracted from the coordinates and the result divided by the range (largest − smallest) before distances are measured. An axis of range 0 is divided by 1 (the same as SciPy). It is for data whose axes have different units (latitude and elevation, and so on). The same change is applied to query points. The distance of `with_distance_bound` is also the one after the change (the same as SciPy).
  * **Tree**: `ImmutableKdTree<f64, D>` of `kiddo` 6.3 (a tree that is not changed once built; 32 points per leaf). The tree that can be added to (`MutableKdTree`) is not taken. Measured with a prototype at design time (random points, 10 000 queries, best values):

    | | $D = 2$, $10^3$ points | $D = 2$, $10^5$ points | $D = 3$, $10^5$ points | $D = 5$, $10^5$ points |
    | :--- | :--- | :--- | :--- | :--- |
    | Construction (per point): immutable / mutable tree | 20 ns / 37 ns | 47 ns / 75 ns | 51 ns / 89 ns | 57 ns / 104 ns |
    | Nearest neighbor (once): immutable / mutable tree | 96 ns / 98 ns | 141 ns / 169 ns | 220 ns / 277 ns | 774 ns / 1126 ns |

    The mutable tree cannot be built when points with the same coordinate on one axis exceed the size of a leaf (with 100 identical points, the 33rd addition was an error). Grid-like data and duplicate points do that, so it cannot be used.
  * **Duplicate points**: points whose coordinates are all equal are merged into one at construction, using the value of the smallest index (the rule of §6.2). Only the points after merging go into the tree. SciPy puts duplicate points into the tree as they are, and which it returns depends on how the tree was built.
    * Coordinates are compared as numbers (`0.0` and `-0.0` are the same point).
    * Duplicates are found with a hash set of the coordinates seen (the standard library's hash). Compared with sorting the points (measured in M3-1; random points in two dimensions, per point): 33 ns and 32 ns for 100 points, 25 ns and 43 ns for 1000, 45 ns and 97 ns for $10^5$, 52 ns and 231 ns for $10^6$. Sorting wins only on points already in order (12 ns and 27 ns on a grid ordered row by row). A lighter hash (a multiplication and a rotation) was 2 times faster on random points but 5 times slower on a grid, whose coordinates differ in few bits, and is not taken.
  * **Equidistant points** (not duplicates; different points): the value of one of them is returned. Which one is fixed for a given set of points, but is not defined as a rule (SciPy is the same: it depends on how the tree was built).
    * "The smallest index" was measured and not taken. To see every equidistant point with `kiddo`, a second query for the points within the nearest distance is needed after finding it. In the prototype one evaluation took 2 times as long ($D = 2$, $10^5$ points: 141 ns became 279 ns; $D = 5$: 774 ns became 1525 ns), with one allocation every time to return the results. That breaks the requirement that evaluation does not allocate (§3.2).
    * Using the tree as it is, 141 of the 361 cell centers of a 20×20 grid returned a point that was not the smallest index.
  * **Distance bound**: with `with_distance_bound(bound)`, a value is returned only when the nearest distance is **less than** `bound`, and `None` otherwise (SciPy's `distance_upper_bound`; SciPy too uses less than, and a point exactly at `bound` is excluded). The default is no bound. `distance_bound()` reads it.
  * **Distance**: Euclidean only (SciPy's `p=2`). The search is exact (SciPy's `eps=0`).
  * **Evaluation**: one query of the tree for the nearest point, and one look-up of the value. No allocation (the prototype confirmed 0 allocations over 1000 queries). A batch repeats `eval`. There is no hint as on the grid (a tree search has no way to use the result of the previous point).
  * **Queries that are not finite**: a coordinate that is NaN or infinite gives `None` (SciPy raises).
  * **Overflow**: a query so far away that the square of the distance does not fit in `f64` gives `None` (when the nearest distance is not finite. SciPy gets an infinite distance and returns NaN).
  * **Dimension 0** (`D = 0`): there are no coordinates, so every point is the same single point. No tree is built, and the first value is returned (`kiddo` panics when a tree of one point in dimension 0 is queried; found while implementing M3-1).
  * **What it keeps**: the tree (of the merged points), the mapping from the merged points to the original indices, and the original values. With `rescaled`, the mean and range per axis. It returns `len()` (the number of points given) and `values()`.
  * **Errors**: in the order of §3.5.
    1. `LengthMismatch` (the numbers of points and values)
    2. `InsufficientPoints { required: 1, actual: 0 }` (it builds from 1 point; the same as SciPy)
    3. More points than `u32::MAX` is `InvalidParameter { name: "points", reason: TooLarge { max } }` (the tree keeps point numbers as `u32`. The limit depends on the input alone, so it is returned as an error. Counting is enough to know it, so it comes before the next check, which looks at every point)
    4. `NonFinite { index }` (coordinates first, then values; `index` the point's number)
    5. `rescaled` only: if the sum of an axis's coordinates (for the mean) or its range overflows, `NonFinite { index }` (`index` the point where the sum overflowed, or the point with the largest coordinate on that axis). Going on would make the changed coordinates the same value or infinite for every point, and no query would have a value. SciPy raises in building the tree when the mean overflows
  * **Differences from SciPy**:

    | Item | SciPy 1.18 | `interpx` |
    | :--- | :--- | :--- |
    | Duplicate points | depends on how the tree was built | the value of the smallest index |
    | Equidistant points | depends on how the tree was built | the same (not defined). Not necessarily the same point as SciPy |
    | Distance | `p` (Minkowski) can be chosen | Euclidean only |
    | Approximate search (`eps`) | yes | no |
    | Queries that are not finite | an exception | `None` |
    | Distance bound | passed with every query | set after construction with `with_distance_bound` |
    | Parallelism (`workers`) | yes | M5-1 (`rayon`) |
    | Coordinate type | `float64` | `f64` |
    | How `rescale` sums for the mean | pairwise (NumPy's `mean`) | in order. The last few bits of the mean can differ. The order of distances does not change, so the point returned is the same. With a bound, a query exactly at the bound's distance can be decided differently |
  * The goldens compare values only for queries whose nearest point is unique. For queries with equidistant points, the test checks that the value returned is "the value of one of the equidistant points".

### 6.3 `CloughTocher2DInterpolator` and `RBFInterpolator`
**CloughTocher2D** (M4-2, #185, SciPy's `CloughTocher2DInterpolator`)

`CloughTocher2DInterpolator<V: Value<f64> = f64>`. C1 interpolation (the value and the first derivatives are continuous) on the Delaunay triangulation of points of the plane: each triangle is cut at its centroid into three, and on each part the interpolant is a cubic Bernstein–Bézier patch. `Point = [f64; 2]`, `Output = V`.

* **Construction**: `new`, `from_owned`, `rescaled`, `rescaled_from_owned`, `with_simplex_tolerance`, of the same shape as `LinearNDInterpolator`. The triangulation, the way a point is located, the change of coordinates of `rescaled`, the tolerance, and the checks of the stored form are held by `Mesh` (`src/scattered/mesh.rs`), which `LinearND` shares (owner's decision). Getters: `len`, `values`, `points`, `gradients` (SciPy's `grad`), `simplex_count`, `simplex_tolerance`.
* **Gradients at the points** (the global estimate of Nielson 1983 and Renka and Cline 1984; the equations of SciPy's `estimate_gradients_2d_global`): the gradients $G_i$ that minimize the sum over the edges of the integral of the square of the second derivative along the edge. Setting the derivative to zero gives, at each vertex $i$, with the edges $E = p_j - p_i$ and $L = \lvert E \rvert$,
  $$
  \sum_j \frac{4 E E^{\mathsf T}}{L^3} G_i + \sum_j \frac{2 E E^{\mathsf T}}{L^3} G_j = -\sum_j \frac{6 (f_i - f_j) E}{L^3}.
  $$
  The whole is a symmetric positive definite system.
  * **How it is solved**: conjugate gradients with the 2 × 2 block of each vertex as the preconditioner, all columns of the values at once (the inner products and the coefficients are per column). A column that is solved stops there. Every column comes out bit for bit as it does when solved alone (a test).
  * **SciPy does not solve it through**: SciPy stops its Gauss–Seidel iteration over the vertices at `tol` (default $10^{-6}$) or `maxiter` (default 400), and only warns when it has not converged. Its gradients, and its values, are that far from the solution. interpx solves until the residual is below $10^{-14}$ of the right-hand side (owner's decision). There are no `tol` and `maxiter`. A system not solved within the limit ($10 n + 100$ steps) returns `SolverFailed(SolverError::NotConverged)`.
  * **Convergence**: about a quarter per step, and below $10^{-16}$ within 26 steps on every set tried (random points up to $10^5$, coordinates whose ranges are a thousand times apart, clustered points, two points $10^{-12}$ apart, triangles $10^{-9}$ thin). No input that gives `NotConverged` has been found.
  * **The threshold**: at $10^{-14}$ the values are within $4.5 \times 10^{-13}$ of the largest value from SciPy iterated to convergence. At $10^{-12}$ they were $2.9 \times 10^{-11}$ apart (`units/plain`).
  * A point that is no vertex of a triangle (the second and later of points with the same coordinates) has the gradient zero.
  * **Independent of size** (the review of #187): the test of convergence compares sums of squares of the residual and the right-hand side, so with values near 1e±160 or coordinates near 1e±100 the squares overflowed (or became zero), and the gradients came back zero with `Ok`, or the system was not solved. Each column of the values and the span of the points are brought near one by powers of two before the system is assembled, and the solution is scaled back by the same powers. Multiplying by a power of two rounds nothing, square roots included, so for input of the usual sizes the gradients are the same bit for bit as before (checked on 5 sets of points, up to 10 000, with and without `rescaled`). Values times 1e±300, and the coordinates of both axes times the same 1e±300, give the same interpolant (a test). Two limits remain (from the re-review of #187, checked here too). Values above about 6e307 in magnitude overflow the sums of the patch, and evaluations give NaN (40 points: finite up to 6e307, and at 1e308 168 of 182 queries NaN; SciPy 1.18 is finite up to 1e300 and its gradients are NaN at 1e306 to 1e307). Axes whose spans are about 1e200 apart (x times 1e100 and y times 1e-100, for one) overflow the system, `NotConverged`; `rescaled` solves it (written in the rustdoc). When the sum of squares of the residual is no longer finite (two points 1e-100 apart in a span of 1, for one), it returns `NotConverged`. An input that does not converge runs to the limit of $10 n + 100$ steps, so the error takes time of the order of $n^2$ (written in the rustdoc).
* **The coefficients of the patch**: SciPy's `_clough_tocher_2d_single`, copied down to the order of operations. The direction of the derivative across an edge (SciPy's `g[k]`) comes from the centroid of the triangle across the edge, in the barycentric coordinates of this one ($-1/2$ on an edge of the hull). It depends on the positions of the points only, so it is computed for each triangle at construction and kept (owner's decision; SciPy computes it at every evaluation).
* **Evaluation**: `Mesh` finds the triangle and the barycentric coordinates (§6.2; the search and tolerance of `LinearND`), and the patch is evaluated. A query outside the convex hull, or with a coordinate that is not finite, is `None`. It does not allocate (`tests/alloc.rs`).
* **Order of the checks**: that of `LinearNDInterpolator` in two dimensions (lengths → number of points (3 or more) → the limit on the number of points → coordinates → values → (`rescaled` only) mean and range → triangulation), and then `NotConverged`.
* **Agreement with SciPy** (golden, `compare/goldens/scattered.json`):
  * 10 cases with SciPy iterated to convergence (`tol=1e-14`, `maxiter=10^6`, warnings as errors): within $4.5 \times 10^{-13}$ of the largest absolute value (allowed $10^{-12}$). Four sets of points on Franke's first function, coordinates in different units (with and without `rescaled`), two columns of values, points with the same coordinates, queries just outside an edge of the hull (tolerance 100), 2000 points.
  * 11 cases at SciPy's defaults: within $2.8 \times 10^{-8}$ (allowed $10^{-7}$), the part SciPy leaves unsolved.
  * The 2 cases with four or more points on one circle (`franke_33`, and a lattice with points taken out) are not compared. The triangulation differs from Qhull's, and since the gradients are the solution of one system through all edges, one edge that differs moves the values everywhere (`franke_33`: one edge of 88 differs, $6 \times 10^{-3}$ of the largest value). For `LinearND` only the values in that region differed.
* **Properties** (tests): a linear function is reproduced exactly, with its gradient. A quadratic is reproduced only when the gradients at the points are replaced by the exact ones (a property of the patch; from the estimated gradients the error stays at $3 \times 10^{-3}$, in SciPy too). At a data point, the value given. The difference of the slopes across an edge goes down tenfold or more when the step does (C1; for linear interpolation it stays).
* **Order of convergence**: §8.4.2.
* **Measured speed** (criterion, this machine, uniform points of the unit square, 1000 queries): construction 61 µs for 100 points, 0.80 ms for 1000, 125 ms for 100 000 (after P-21's `convx` 0.1.1 and snake order of the cells; `LinearND` in 2-D: 56 ms for 100 000). At the end of M4-2 it was 71 µs, 0.99 ms, and 146 ms. Evaluation 40 ns a query (100 points), 52 ns (1000), 104 ns (100 000). These are after the fix of "Independent of size" above; before it, construction was 72 µs, 1.03 ms, and 152 ms. Timed first in a process without a warm-up, the construction of 100 points was 20 to 30% slower after the fix (about 60 µs → 70 to 80 µs; 3 sets of 3 rounds in turn). A version with only the values scaled, and one with only the coordinates scaled, were as slow; one with only the check for non-finite sums was not. With a warm-up (criterion) there is no difference, and the cause is not known. The comparison's runners warm construction up for 1 second before they time it.
  * **Where construction went** (a temporary build with a clock around each stage, 3 to 4 constructions in one process; before the three changes below): at 100 000 points, the triangulation and the search's preparation (`Mesh`) 94 to 98 ms, the gradient estimate 144 to 178 ms, `g[k]` 7 to 10 ms. 70 to 80% of the estimate was the conjugate gradient steps (19); collecting the edges and assembling the system were the rest (at 1000 points the steps were 40 to 80%). The figures moved by up to a factor of two from one construction to the next. In a profile (Visual Studio's CPU profiler), the gradient estimate took about 49%, `convx`'s triangulation about 28%, the rest of `Mesh` about 11%, and the sort of the edges about 7%.
  * **Three changes that made the estimate faster** (a workload of 15 constructions of 100 000 points, against the version before in turn, 3 rounds, 6 for (3), whose difference was small; one core, best against best): (1) the edges are taken once each from the neighbors of the triangles instead of collected and sorted (184.0 → 175.1 ms). (2) The product of the matrix and a vector runs over the edges instead of the rows of the vertices: each edge's weight is read once and added to both ends, and nothing the steps do not use is kept (170.3 → 150.7 ms). (3) The inverse of each block of the preconditioner is computed once, before the steps, and the updates of the solution, the residual, and the preconditioned residual are made in one loop with the two dot products (152.8 → 144.6 ms; faster in 5 rounds of 6). (3) keeps the order of the sums, so the results are the same bit for bit; (1) and (2) change it. The distance from the goldens is unchanged (4.4e-13 from SciPy iterated to convergence). Most of what is left is `convx`'s triangulation (YUKIKEDA/convx#667).
* **Differences from SciPy**: the gradients are solved through (above). The division of points on one circle (above). Of points with the same coordinates, the value of the first is used (as in `LinearND`).
* **Stored form**: what `LinearND` stores (points, values, vertices of the simplices, their neighbors, change of coordinates, tolerance) and the gradients. Reading does not solve the system again. `g[k]` is rebuilt (§6.5).

**RBF**
* The linear system:
  $$
  \begin{bmatrix} A & P \\ P^T & 0 \end{bmatrix}
  \begin{bmatrix} w \\ c \end{bmatrix}
  =
  \begin{bmatrix} y \\ 0 \end{bmatrix},
  \quad A_{ij} = \phi(\|x_i - x_j\|) + \lambda \delta_{ij}
  $$
  ($\lambda$ is `smoothing`, added to the diagonal. The same sign as SciPy 1.18; confirmed in M3-2)
* **Eight kernels** (the enum `RbfKernel`): `linear`, `thin_plate_spline`, `cubic`, `quintic`, `gaussian`, `multiquadric`, `inverse_multiquadric`, `inverse_quadratic`.
* How `epsilon` and `degree` are held, the polynomial terms, and the change of coordinates are as in "The parts of RBF (M3-2)" below.
* The public type, checks, solve, evaluation, and memory are as in "Public, and the global solve (M3-3)" below.
* The form that solves over neighbors only is as in "The neighbors mode (M3-4)" below.

**The parts of RBF (M3-2)**

M3-2 puts the parts up to assembling the linear system into the crate. It exposes no public type. This is so as not to publish a type that cannot solve first; `RBFInterpolator` and `RbfKernel` are published by M3-3, which brings the solve. The parts are in `src/scattered/rbf.rs`, and are checked against the output of SciPy's internal functions.

* **Kernels**: functions of $r = \varepsilon \lVert x - y \rVert$ (the Euclidean distance multiplied by `epsilon`). The formulas and signs are SciPy's.

  | Kernel | $\phi(r)$ | `epsilon` | Minimum degree | Default degree |
  | :--- | :--- | :--- | :--- | :--- |
  | `Linear` | $-r$ | no | 0 | 0 |
  | `ThinPlateSpline` | $r^2 \log r$ (0 at $r = 0$) | no | 1 | 1 |
  | `Cubic` | $r^3$ | no | 1 | 1 |
  | `Quintic` | $-r^5$ | no | 2 | 2 |
  | `Multiquadric` | $-\sqrt{r^2 + 1}$ | yes | 0 | 0 |
  | `InverseMultiquadric` | $1/\sqrt{r^2 + 1}$ | yes | none | 0 |
  | `InverseQuadratic` | $1/(r^2 + 1)$ | yes | none | 0 |
  | `Gaussian` | $e^{-r^2}$ | yes | none | 0 |

  * Some have a negative sign to make them conditionally positive definite (the same as SciPy; Fasshauer, *Meshfree Approximation Methods with MATLAB*, chapter 8).
  * Distances are found in SciPy's order: the coordinates are multiplied by `epsilon`, then the difference taken, then the square root of the sum of squares.
* **`epsilon`**: the kernels that need it (the four marked "yes") carry it in the variant (in the form `RbfKernel::Gaussian { epsilon }`). The four that do not need it do not carry it. The state "needed but not given" cannot be written, so there is no error for it. Only finite positive values are accepted (the check and the variant of its reason are brought by M3-3, which exposes the public type).
  * The four without it are kernels that do not depend on scale (changing `epsilon` does not change the interpolant without smoothing; confirmed with SciPy 1.18: 30 points, `epsilon` 1 and 2.5, differences of at most $10^{-13}$). SciPy lets `epsilon` be given to these four too (default 1). With smoothing the result changes, but that is the same as changing the strength of the smoothing. $\phi(\varepsilon r)$ is $\varepsilon^p \phi(r)$ (for `ThinPlateSpline`, plus what the polynomial terms absorb), so SciPy's $(\varepsilon, \lambda)$ corresponds to $\lambda / \varepsilon^p$ without `epsilon` ($p$ is 1 for `Linear`, 2 for `ThinPlateSpline`, 3 for `Cubic`, 5 for `Quintic`). The goldens of M3-3 check this.
* **The degree of the polynomial**: held as `Option<usize>`. `Some(k)` adds the monomials up to degree $k$. `None` adds no polynomial terms (SciPy's `degree=-1`).
  * The default is the "default degree" of the table (the same as SciPy; kernels with no minimum degree still default to degree 0, adding a constant term).
  * `Some(k)` with $k$ below the minimum degree is `InvalidParameter { name: "degree", reason: TooSmall { min } }` (SciPy only warns and goes on; the solution may not be unique).
  * `None` is accepted for every kernel (SciPy does not warn for `-1` either). Adding no terms for a kernel with a minimum degree can make the system singular; the solve of M3-3 then returns an error.
* **Monomials**: there are $\binom{D + k}{D}$ monomials of degree up to $k$ in $D$ variables. Their order is SciPy's: by increasing degree, and within a degree "the multisets of variables in lexicographic order" ($1,\ x,\ y,\ x^2,\ xy,\ y^2$ for degree 2 in 2 variables).
  * A number that does not fit in `usize` is `InvalidParameter { name: "degree", reason: TooLarge { max } }` (`max` the largest degree whose number fits in that dimension).
  * Fewer points than monomials is `InsufficientPoints { required, actual }` (`required` the number of monomials. The matrix of the polynomial cannot have full column rank. SciPy errors on the same condition).
* **The change of coordinates** (for the polynomial terms only): the coordinates put into the monomials are $\hat{x}_d = (x_d - s_d) / c_d$ per axis. $s_d$ is the midpoint of the largest and the smallest, $c_d$ half the range, and an axis of range 0 has $c_d = 1$ (the same as SciPy). The data's coordinates fall in $[-1, 1]$, so the values of monomials of high degree do not grow too large. The coordinates put into the kernel are not changed (only multiplied by `epsilon`).
  * This is a different thing from the change of `rescaled` of `NearestNDInterpolator` (mean and range). Both are kept as SciPy's formulas.
* **System**: with $N$ points and $R$ monomials, an $(N + R) \times (N + R)$ symmetric matrix is made in a flat row-major array. The upper left is $A_{ij} = \phi(\varepsilon \lVert y_i - y_j \rVert) + \lambda\,\delta_{ij}$, the upper right $P_{ij}$ (the value of monomial $j$ at point $i$), the lower left $P^T$, the lower right 0. The right-hand side is the $N$ values followed by $R$ zeros.
  * An allocation size $(N + R)^2$ beyond what a `Vec` can hold is `InvalidParameter { name: "points", reason: TooLarge { max } }`. Within it, an allocation that cannot be made is `AllocationFailed { bytes }` ("The size of allocations and their failure" of §4).
  * The kernel value is found once per pair of points and written in two places of the matrix (it is symmetric). It goes tile by tile of 64 points square, so that the two places written are in the cache. Compared with finding every entry (measured in M3-2; two dimensions, per entry; the first figure of each pair is the form adopted): `ThinPlateSpline`, with its logarithm, 3.4 ns and 6.3 ns on 100 points, 5.3 ns and 7.2 ns on 3000. `Gaussian`, with its exponential, 3.2 ns and 6.1 ns, 5.0 ns and 6.8 ns. The light `Linear`, `Cubic`, `Quintic` are 1.0–1.5 ns and 1.6–2.1 ns on 100 points, and 2.7–3.1 ns and 2.2–2.3 ns on 3000: slower on large matrices (writing the second place costs more than the computation). Assembly grows with the square of the number of points and the solve after it with the cube, so assembly weighs less the larger the matrix. The form fast on small matrices and the default kernel is taken (its share in the whole construction is measured in M3-3, which brings the solve). Copying the lower half afterwards in one pass was slower still on the light kernels (3.1–3.5 ns on 3000 points).
* **The evaluation row**: for a query point $x$, the $N$ values $\phi(\varepsilon \lVert x - y_i \rVert)$ followed by the $R$ values of the monomials at $\hat{x}$. Its inner product with the coefficients is the interpolated value (M3-3). It is written to space the caller passes, without allocation.
* **Checks**: the goldens are the output of SciPy's internal functions of `_rbfinterp_pythran`: `_build_system` (matrix, right-hand side, change of coordinates), `_kernel_matrix`, `_polynomial_matrix`, `_build_evaluation_coefficients` (the evaluation row), and `_monomial_powers` (the order of the monomials).

**Public, and the global solve of RBF (M3-3)**

* **Type**: `RBFInterpolator<const D: usize, V: Value<f64> = f64>`. `Point = [f64; D]`, `Output = V`. It solves one linear system over all the data points (SciPy's `neighbors=None`). The form over neighbors only is M3-4.
* **Construction**: `new(points: &[[f64; D]], values: &[V], kernel: RbfKernel)` and `with_options(points, values, kernel, &RbfOptions)`.
  * `RbfKernel` (`#[non_exhaustive]`): the eight of the table of M3-2. The four variants that need `epsilon` carry it. There is no default kernel (SciPy's default is `thin_plate_spline`; it is always chosen as an argument of `new`).
  * `RbfOptions`: `RbfOptions::new()` makes the defaults (smoothing 0, the kernel's default degree), changed with `with_smoothing(f64)` and `with_degree(RbfDegree)`. The fields are not public (so that adding an option breaks nothing).
  * `RbfDegree`: `Default` (the kernel's default degree), `None` (no polynomial terms), `Degree(usize)`. It corresponds to the `Option<usize>` of M3-2.
* **Smoothing**: a scalar $\lambda \ge 0$. 0 interpolates (passes through the data points), positive smooths (does not pass). SciPy also accepts an array per point; that is not in v0.1.0 (roadmap A-11). Negative values pass in SciPy but are refused here (what the system represents is not determined).
* **The order of the checks** (this type's after the order of §3.5):
  1. `LengthMismatch` (the numbers of points and values)
  2. `InsufficientPoints { required: 1, actual: 0 }`
  3. `NonFinite { index }` (coordinates, then values)
  4. `InvalidParameter { name: "epsilon", reason }`: `NonFinite` if not finite, `NotPositive` if 0 or less
  5. `InvalidParameter { name: "smoothing", reason }`: `NonFinite` if not finite, `Negative` if negative
  6. `InvalidParameter { name: "degree", reason }` (below the minimum degree, monomials too many to count. M3-2)
  7. `InsufficientPoints { required, actual }` (fewer points than monomials. M3-2)
  8. `DuplicatePoints { first, second }` (points with the same coordinates. `first` is the point given first. Coordinates are compared as numbers, `0.0` and `-0.0` being the same)
  9. `InvalidParameter { name: "points", reason: TooLarge { max } }` (the matrix does not fit in a `Vec`. M3-2). If it fits but not in memory, `AllocationFailed { bytes }` (§4)
  10. `SolverFailed(SolverError::Singular)` ("Singular systems" below)

  The neighbors mode (`with_neighbors`) adds the checks of `neighbors` to this order and replaces 9 and 10. Where they are added is written, by these numbers, under "Checks" of "The neighbors mode" below.
* **The values of monomials** (P-5, #126): a term $x^p$ of the polynomial is found on the spot by multiplications that repeat squaring (`f64::powi` is a function call). The powers of monomials are mostly 0, 1, 2, whose values are exactly $1$, $x$, $x \cdot x$. From 3 on, they can differ by one rounding from the result of `powi` (both are products of the same number of factors). What `powi` does can differ by environment, so the values of powers of 3 or more could differ by one rounding by environment before. Now the order of the product is decided by the code, and the value is the same in every environment. When `powi` was called, more than half of the evaluation of 10 points of a two-dimensional `quintic` (6 monomials) was this call (75 → 24 µs for an evaluation of 1000 points; SciPy 71 µs). Both the global and the neighbors forms use these values at every evaluation.
* **Duplicate points**: an error with or without smoothing. With smoothing the system can be solved even with duplicates (SciPy accepts them), but one rule is kept. They are found with the same hash set as `NearestNDInterpolator`, sharing the part.
* **Solve**: the matrix is factored by the partially pivoted LU factorization of `faer` 0.24, and the forward and back substitutions are the crate's own. The factorization works in place on the assembled matrix (`faer`'s low-level `lu_in_place`; the matrix is symmetric, so the row-major array can be passed as it is as column-major). `faer`'s high-level `partial_piv_lu()` copies the matrix and then splits $L$ and $U$ into separate matrices, using the memory of three matrices. With values `[f64; M]` the right-hand side has $M$ columns. `Value` has no operation that takes out a column, so `faer`'s $L$, $U$, and permutation are taken out, and the substitutions are done by linear combinations of values (`axpy`) (the same form as the banded LU of M1-12). Substitution is $O(N^2)$, small beside the $O(N^3)$ of the factorization.
  * **Why LU** (measured with a prototype at design time): SciPy solves by a symmetric indefinite factorization (LAPACK `dsysv`). `faer` has the same kind of factorization too (`lblt`). Both solved the 72 global cases of `rbf.json`, and were compared with SciPy's values (the difference as a ratio to the largest absolute value).

    | | LU (partial pivoting) | Symmetric indefinite (`lblt`) |
    | :--- | :--- | :--- |
    | Largest difference from SciPy | $9.0 \times 10^{-13}$ | $5.5 \times 10^{-13}$ |
    | Time, 100 points | 0.138 ms | 0.135 ms |
    | Time, 1000 points | 20.1 ms | 18.0 ms |
    | Time, 3000 points | 448 ms | 358 ms |

    The precision is about the same, and the symmetric factorization is 2 to 20% faster. LU is still taken because it refuses singular systems reliably. On a system with the same point twice (singular), the LU's pivot became 0, and the symmetric factorization returned a finite solution (a wrong one). The symmetric factorization has no way to check its pivots from outside. The fully pivoted LU took 111 ms on 1000 points, 5 times slower, with the same precision, and is not taken.
    * **"About the same precision" is about small well-conditioned systems** (P-14, #164). On real station data (10 838 points; §8.3), comparing the residual of the coefficients solved put back into the same matrix (600 rows, compensated sums), LU's is larger than the symmetric factorization's:

      | Kernel (the system's reciprocal condition number) | interpx (`faer`'s LU) | SciPy | LAPACK `gesv` (LU) | LAPACK `sysv` (symmetric indefinite) |
      | :--- | ---: | ---: | ---: | ---: |
      | `linear` ($10^{-10}$) | $8.2 \times 10^{-8}$ | $9.1 \times 10^{-9}$ | $3.4 \times 10^{-8}$ | $8.7 \times 10^{-9}$ |
      | `cubic` ($10^{-20}$) | 21.0 | 1.97 | 9.11 | 2.00 |

      SciPy's values agree with `sysv`'s. LU, even LAPACK's, is about 4 times that, and interpx about 2.4 times more (the reason for this 2.4 is not known). **But this difference is not a difference in the correctness of the solution**: the residual $8.2 \times 10^{-8}$ of `linear` is $6 \times 10^{-18}$ of the size of the matrix and solution (about $1.4 \times 10^{10}$), below the rounding of double precision. Both are solved stably, and the precision of the solution guaranteed for a system of condition number $10^{10}$ (about 6 digits) is the same for both. The reason for taking LU (it refuses singular systems reliably) does not change.
    * **Correcting the solution by its residual (iterative refinement) was built, measured, and not put in** (P-14. The owner's decision). The residual is found with compensated sums, solved again with the same factorization, and added to the solution. On the `linear` system above, the residual went from $8.2 \times 10^{-8}$ to $9.4 \times 10^{-10}$ (a tenth of SciPy's). But the interpolated values moved by $10^{-7}$ °C at most (the values are −52 to 44 °C), and the error at the held-out points was the same to the 10th digit. In exchange, construction gets slower by making the matrix's entries once more: always correcting, 2.0 times on 100 points (slower than SciPy), 1.5 times on 1000, 1.07 times on 10 000. The form that decides by the condition number whether to correct (when the reciprocal condition number is below $\sqrt{\varepsilon}$ and above $\varepsilon$) was also 1.55 times on 1000 points. And on a system not solved in double precision (`quintic`, reciprocal condition number $5 \times 10^{-30}$), correcting made the residual 170 times worse.
* **Estimating the condition number** (P-14, #164. The owner's decision): after the system over all points is solved, an estimate of its reciprocal condition number (1-norm) is kept and returned by `reciprocal_condition()` (`None` for neighbors). **Under the $\varepsilon$ of double precision ($2.2 \times 10^{-16}$), the system is not solved in double precision**: the interpolator is built, but its values mean little. Neither interpx nor SciPy refuses such a system. This number is the way to know.
  * **Not an error**: it would refuse input for which SciPy returns values, and which range would be refused has not been measured (the 1000 points of the Gaussian on synthetic data are $5 \times 10^{-20}$ too). The number is returned, and the judgment left to the user. The rustdoc states the numbers of the station data below and what to do (merge close points, use `linear`, use neighbors).
  * **A number that does not depend on the unit of the coordinates** (the review of #170): the condition number of the system as it is solved changes with the unit of the coordinates. The monomials are made after the points are moved to $[-1, 1]$, but the kernels are of the distances as they are, so multiplying the coordinates by $s$ multiplies the `cubic` kernel by $s^3$ and leaves the monomials. On 100 points of Franke's first function, taking $s$ from 1 to $10^6$ took the reciprocal condition number of `cubic` from $3.2 \times 10^{-7}$ to $1.9 \times 10^{-38}$, and that of `quintic` to $7.9 \times 10^{-63}$. How far the interpolant misses its data points stays $10^{-13}$ and $4 \times 10^{-12}$ for every $s$. With nothing wrong, writing in metres alone gives a number that reads as "not solved". So the condition number of $DAD$, balanced on both sides by a diagonal matrix $D$, is estimated. $D$ puts $a = 1/\sqrt{\max |K_{ij}|}$ on the kernel rows and columns and $1/a$ on the monomial rows and columns. The kernel block of $DAD$ has entries of size 1 at most, the blocks joining kernels and monomials do not change, and the block of monomials stays 0. For kernels that are powers of the distance, $DAD$ is the same matrix whatever $s$ (`linear` $2.1 \times 10^{-4}$, `cubic` $1.4 \times 10^{-7}$, `quintic` $8.9 \times 10^{-10}$, for every $s$). `thin_plate_spline` is $r^2 \log r$ and changes shape with the unit, so over a factor of $10^6$ it moves from $1.8 \times 10^{-5}$ to $3.5 \times 10^{-7}$. For kernels with the shape parameter `epsilon` (the four such as `multiquadric`), the number moves with the unit. $\varepsilon r$ changes with the unit and it becomes another problem, so this is right (a Gaussian with `epsilon` 3, the same 100 points: $6.2 \times 10^{-13}$ at $s = 1$, under $10^{-19}$ at $s = 10^{-4}$ ($5.0 \times 10^{-20}$ in the review's measurement of #170, $7.6 \times 10^{-20}$ here; far below $\varepsilon$, an order of magnitude only). `multiquadric`: $7.3 \times 10^{-10}$ at $s = 1$, $2.8 \times 10^{-19}$ at $s = 10^{-4}$). The rustdoc says so too. The system solved does not change (the solution is not balanced). Solving with $DAD$ is applying $D^{-1} A^{-1} D^{-1}$ with the factored $A$, so no extra factorization is needed.
  * **How it is estimated**: the norm of the matrix is found from the column sums before the factorization. The norm of the inverse is estimated by a few solves with the factored system, taking the larger of two methods. (1) Hager's search (the same as LAPACK's `dgecon`). (2) A vector whose entries alternate in sign and grow is put through the inverse once, what comes out is scaled to size 1, and put through twice more. The estimate is a lower bound of the norm of the inverse, so the value returned is not below the true value.
  * **Why (2) is needed**: the search starts from a vector of equal entries, so it cannot find an inverse that is large only across two nearly equal rows (two close points). On a `linear` system with two points $10^{-7}$ apart, the search alone read $1.6 \times 10^{-2}$ (correctly $5.6 \times 10^{-9}$). That close points worsen the condition is what this row found, so an estimate loose there is no use. LAPACK's `dgecon` itself is 4 times loose on a system of the same form ($5.5 \times 10^{-4}$ as $2.1 \times 10^{-3}$).
  * **The precision measured** (202 symmetric matrices: Hilbert, random, with graded eigenvalues, and systems of the four kernels with close and far points added; up to 200 points. Estimate / the true norm of the inverse; 1 is exact, smaller is looser):

    | Configuration | Worst |
    | :--- | ---: |
    | The search alone | 0.001 |
    | The search + LAPACK's final check (the alternating vector through once, measured by 2/3 of a mean) | 0.004 |
    | (2) alone | 0.26 |
    | The search + (2) (the one taken) | 0.40 |
    | The search + (2) through once only | 0.15 |
    | The search + (2) starting from a vector of equal entries | 0.05 |

    The form taken is within 2.5 times the true value. LAPACK's final check changed no result when added, and is not put in. Putting the vector through 3 times did not change the worst either. A starting vector without alternating signs (only growing) has the same worst of 0.40, with 0.61 becoming 0.49 for one kind of matrix only: nearly the same, and the form matching LAPACK was kept. On the four systems of the station data (balanced), interpx's estimate and LAPACK's `dgecon` are as follows. On the systems before balancing, the two agreed ($1.04 \times 10^{-10}$, $1.1 \times 10^{-16}$, $1.9 \times 10^{-20}$, $5.4 \times 10^{-30}$).

    | Kernel | interpx | LAPACK `dgecon` |
    | :--- | ---: | ---: |
    | `linear` | $6.89 \times 10^{-11}$ | $6.89 \times 10^{-11}$ |
    | `thin_plate_spline` | $5.03 \times 10^{-17}$ | $5.03 \times 10^{-17}$ |
    | `cubic` | $6.84 \times 10^{-21}$ | $4.50 \times 10^{-21}$ |
    | `quintic` | $3.56 \times 10^{-24}$ | $3.00 \times 10^{-24}$ |

    `cubic` and `quintic` are 5 orders or more below $\varepsilon$, where the factorization and substitutions themselves are rounding. The estimate there is an order of magnitude, and the difference between the two values (up to 1.5 times) means nothing.
  * **The cost**: construction gets slower. In full runs of the comparison, before the estimate (`056f117`) and now (`c5fe29c`, balancing included) were compared (P-17; medians of the rounds):

    | Case | Before | Now | Ratio | SciPy |
    | :--- | ---: | ---: | ---: | ---: |
    | `thin_plate_spline`, 10 points | 1.92 µs | 3.29 µs | 1.71 | 12.3 µs |
    | `linear`, 100 points | 49.5 µs | 75.2 µs | 1.52 | 83.0 µs |
    | `thin_plate_spline`, 100 points | 72.4 µs | 97.3 µs | 1.35 | 117 µs |
    | `thin_plate_spline`, 1000 points | 18.9 ms | 23.9 ms | 1.27 | 66.7 ms |
    | `thin_plate_spline`, 3000 points | 389 ms | 444 ms | 1.14 | 1.94 s |
    | `linear`, the 10 838 stations | 16.6 s | 17.8 s | 1.07 | 116 s |

    Over the eight kernels at 10 to 3000 points, the ratio is 1.7–2.3 at 10 points, 1.35–1.63 at 100, 1.25–1.33 at 1000. At 100 points it is level with SciPy (the median ratio interpx / SciPy 0.83–1.05; 1.05 is 5% slower in the median and does not reach the threshold of 1.10. §8.3). The 10 838 stations are measured once and move from 0.93 to 1.07 over the four kernels. (The PR of #170 said 7 to 78%, from runs of the comparison's runner once each.)

    (2) alone would be lighter (1.26 times for `linear` at 100 points, 1.14 times for `thin_plate_spline`), but the estimate becomes up to 4 times looser, and whether it stays in that range on systems of 10 000 points has not been measured. The search was kept.

    Of this, the part of the form that does not depend on the unit of the coordinates (balancing) was measured separately. Balancing adds a scan for the largest kernel value and the multiplications in computing the norm. Before balancing (the `rbf.rs` of d9d43a1) and now were measured from two source trees built apart, in turn for 3 rounds, comparing the best value of each (criterion, pinned to one core. P-17):

    | Case | Before | Now | Ratio | Range over the rounds (before / now) |
    | :--- | ---: | ---: | ---: | :--- |
    | `thin_plate_spline`, 10 points | 3.64 µs | 3.36 µs | 0.92 | 3.64–3.96 / 3.36–3.47 |
    | `linear`, 10 points | 2.78 µs | 2.93 µs | 1.05 | 2.78–3.22 / 2.93–3.06 |
    | `thin_plate_spline`, 100 points | 93.1 µs | 100.5 µs | 1.08 | 93.1–98.1 / 100.5–102.0 |
    | `linear`, 100 points | 64.5 µs | 72.3 µs | 1.12 | 64.5–67.7 / 72.3–73.4 |
    | `thin_plate_spline`, 1000 points | 21.1 ms | 21.6 ms | 1.02 | 21.1–22.1 / 21.6–22.5 |
    | `linear`, 1000 points | 19.0 ms | 20.0 ms | 1.05 | 19.0–19.5 / 20.0–21.2 |

    At 100 points now is slower in every round, and the cost is 8 to 12%. At 10 points and for `thin_plate_spline` at 1000 the ranges over the rounds overlap, and no difference is measured. `linear` at 1000 points has ranges apart, 2 to 5%. (#170 said 6 to 8% at 100 points from criterion run once each; it had not been measured by the procedure of `bench.mdc`.)

    Reusing the working vectors inside the estimate (an optional suggestion of the review of #170: not making a `Vec` at every solve) was measured by the same procedure, as a third source tree. Its ratio to now is 0.95–0.98 at 10 points, 1.00–1.02 at 100, 0.99–1.01 at 1000, and in every case its range over the rounds overlaps now's (`linear` at 10 points: 2.77–2.95 µs and 2.93–3.06 µs). No difference is measured, so it is not put in (`bench.mdc`: a gain within the spread is not put in). #170 said "7 to 10% faster at 10 points, 3 to 6% slower at 100" from criterion run once each; by this procedure that did not reproduce.
  * The estimate goes into the stored form (it cannot be found without solving again. §6.5). Reading checks that it is within 0 to 1.
  * `faer` is `default-features = false, features = ["std"]` (parallelism in M5-1).
* **Singular systems**: when the solution has a value that is not finite, `SolverFailed(SolverError::Singular)`. A pivot of 0 in the factorization makes the values divided by it infinite or NaN, so this one check catches both singular systems and systems whose numbers do not fit (a check of the pivots on their own was confirmed by a mutation check to change no result when removed, and is not kept). Anything else is solved and returned. Systems only near singular (points on a slanted line with linear polynomial terms added, and so on) can return values of meaningless size. SciPy is the same: it returned values of $10^{13}$ on the same input. Points on a line parallel to a coordinate axis make a pivot exactly 0, and both give an error. Refusing by a threshold on the size of pivots is not taken (it would refuse even systems for which SciPy returns the same values despite bad conditioning, such as a `Gaussian` with a small `epsilon`).
  * Solving again by a pseudo-inverse is roadmap A-5.
* **What it keeps**: the points ($N \times D$), the coefficients ($N + R$ values), the kernel, and the polynomial terms (exponents, shift, and scale). It does not keep the matrix or the factorization. During construction one $(N + R)^2$ matrix is allocated (about 800 MB for $N = 10^4$). Measured, the peak of allocation is 1.001 times the matrix (3000 points), and what remains after construction is 24 bytes per point (two dimensions, values `f64`). Besides, `faer` allocates a working space for matrix products once per thread and keeps it (about 4.2 MB; measured on a construction of 1000 points). No upper bound is set (only sizes that do not fit in a `Vec` are refused). For large $N$, use the neighbors mode of M3-4. Written in the rustdoc.
* **Evaluation**: $f(x) = \sum_i w_i\,\phi(\varepsilon \lVert x - y_i \rVert) + \sum_j c_j\,p_j(\hat{x})$ is found by adding the terms one at a time. No allocation (the evaluation row is not made as an array; the `evaluation_row` of M3-2 is changed to this way of summing). The cost is proportional to the number of data points (about 6 µs per evaluation for 1000 points and `ThinPlateSpline`, from the measurements of M3-2).
  * There is no domain. It returns a value at every point (away from the data, a value decided by the polynomial terms and the shape of the kernel). It returns `None` only when the coordinates are not finite and when the result is not finite.
  * A batch repeats `eval`. SciPy splits the evaluation to fit memory and computes it by matrices, but the values are the same.
* **Getters**: `len()`, `kernel()`, `degree()` (`Option<usize>`, the degree actually used), `smoothing()`.
* **Differences from SciPy**:

  | Item | SciPy 1.18 | `interpx` |
  | :--- | :--- | :--- |
  | Solve | symmetric indefinite factorization (`dsysv`) | LU with partial pivoting. The values differ by less than $10^{-12}$ (the table above) |
  | Default kernel | `thin_plate_spline` | none (always chosen) |
  | `epsilon` | can be given to any kernel; negative values pass | carried only by the four that need it; finite and positive. `epsilon` on the four that do not depend on scale can be reproduced by rescaling the smoothing (M3-2) |
  | `degree` below the minimum | a warning | an error |
  | Smoothing | also an array per point; negative values pass | a scalar, 0 or more (A-11) |
  | Duplicate points | solvable with smoothing; "Singular matrix" without | always an error (`DuplicatePoints`) |
  | Value type | complex too | real (`f64`, `[f64; M]`) |
  | Queries that are not finite | returns NaN | `None` |

**The neighbors mode of RBF (M3-4)**

* **How it is chosen**: `RbfOptions::with_neighbors(k)` (SciPy's `neighbors=k`). Each query point is interpolated over its $k$ nearest points only. The global solve needs time cubic and memory quadratic in the number of points, so this is the form for many points (above several thousand). It is read with `RbfOptions::neighbors()` and `RBFInterpolator::neighbors()` (`Option<usize>`).
  * With $k$ at least the number of points, it becomes the global solve (the result is the same, so it is solved once at construction; SciPy too clips $k$ to the number of points). `RBFInterpolator::neighbors()` then returns `None`.
* **Checks** (added to the order of M3-3):
  * After 2: more points than `u32::MAX` is `InvalidParameter { name: "points", reason: TooLarge { max } }` (the tree keeps point numbers as `u32`. The same check as `NearestNDInterpolator`'s; counting is enough to know it, so it comes before the checks that look at each point. In the neighbors mode only)
  * After 6: `InvalidParameter { name: "neighbors", reason: TooSmall { min: 1 } }` ($k = 0$)
  * 7 reads: when there are fewer points, or a $k$ smaller, than the monomials, `InsufficientPoints { required, actual }` (`actual` the smaller of the number of points and $k$. SciPy too makes the same check with the number of neighbors as the number of points)
  * 9 reads (in the same place, after the duplicate points of 8 and before building the tree): a $k$ whose neighbor system $(k + R)^2$ does not fit in a `Vec` is `InvalidParameter { name: "neighbors", reason: TooLarge { max } }` at construction. A case that fits but not in memory cannot be known at construction (the working space is allocated per evaluating thread), so evaluation returns `None` ("Working space" below).
  * 10 (singular) does not happen in the construction of the neighbors mode (no system is made at construction).
  * Duplicate points are `DuplicatePoints` at construction, as in the global form.
* **Construction**: the points go into a k-d tree. The values are kept as they are. No system is made. The time is $O(N \log N)$ (9.1 ms for $10^5$ points), and the memory only that of the points and values (44 bytes per point in two dimensions: the points in the order given, the points in the tree's order, the numbers, the values. Measured in M3-4).
* **The search for the $k$ nearest points**: a small k-d tree inside the crate (`src/scattered/kdtree.rs`). `kiddo` is not used: `kiddo`'s search for $k$ neighbors allocates the `Vec` of the results every time and has no form that does not allocate (in the prototype at design time, one allocation per query for every $k$).
  * The tree splits at the median, taking the point in the middle of the range as the node and changing axes in turn. Ranges of 16 points or fewer are searched in order. The points are kept reordered in the tree's order.
  * The search goes down the near side first, and looks at the far side only when "the distance from the query point to that box" is nearer than the current $k$-th. The candidates are in a heap of size $k$. The heap and the search's stack are in the working space, with no allocation.
  * Points at the same distance count the one given earlier as nearer (to make the result deterministic. SciPy depends on how the tree was built).
  * Compared with `kiddo` (measured with a prototype at design time; random points, one query. The points chosen were the same as `kiddo`'s in every condition):

    | | $k = 1$ | $k = 20$ | $k = 50$ | $k = 100$ |
    | :--- | :--- | :--- | :--- | :--- |
    | $D = 2$, $10^5$ points: `kiddo` / own | 145 ns / 267 ns | 0.87 µs / 1.27 µs | 1.98 µs / 2.69 µs | 4.24 µs / 5.76 µs |
    | $D = 3$, $10^5$ points | 235 ns / 502 ns | 1.63 µs / 2.21 µs | 3.54 µs / 4.86 µs | 7.06 µs / 8.72 µs |
    | $D = 5$, $10^5$ points | 745 ns / 1169 ns | 3.87 µs / 5.76 µs | 7.68 µs / 10.7 µs | 13.1 µs / 18.6 µs |

    The crate's own tree is 1.2 to 1.5 times slower for $k$ neighbors, and 1.6 to 2.1 times for the single nearest point. Construction is about the same (47–59 ns and 52–70 ns per point for $10^5$ points). In the neighbors mode, solving the system of $k$ points after the search costs more, so this difference is small in the whole evaluation (the ratio is measured in M3-4). For the single nearest point `kiddo` is clearly faster, so `NearestNDInterpolator` stays on `kiddo`.
  * The distance to a far box is found exactly by carrying the overhang of each axis (compared with estimating from one axis's distance alone, 13.1 µs became 10.7 µs at five dimensions with $k = 50$; no change at two and three dimensions).
* **Evaluation**: for each query point $x$,
  1. Find the $k$ nearest points and sort them into the order given (SciPy too sorts them by number).
  2. With those $k$ points, make the change of coordinates (midpoint and half range) and the matrix $A$ of $(k + R)$ square (the formulas of M3-2. The change is found per neighborhood, as SciPy does).
  3. Make the evaluation row $e$ of $x$ ($k$ kernels and $R$ monomials) and solve $A\,\lambda = e$ (the LU of M3-3).
  4. The value is $\sum_{i < k} \lambda_i\, d_i$ ($d_i$ the values of the neighbors).
  * SciPy solves $A\,c = (d, 0)$ and finds $e \cdot c$. $A$ is symmetric, so $e \cdot A^{-1} (d, 0) = (A^{-1} e) \cdot (d, 0)$, the same value. Solving for the row is taken because what is solved is one column of `f64`, one solve even when the values are arrays, and the working space does not depend on the value type.
  * A singular neighbor system (a solution that is not finite) gives `None`. SciPy raises. It cannot be known at construction (the system is made at every evaluation).
  * The values can be discontinuous where the set of neighbors changes (the same as SciPy). Written in the rustdoc.
  * **Cost** (measured in M3-4; two dimensions, $10^5$ points, one evaluation):

    | | Search and sorting | Assembling system and row | Factorization | Substitution | Total |
    | :--- | :--- | :--- | :--- | :--- | :--- |
    | `ThinPlateSpline`, $k = 10$ | 0.93 µs | 0.85 µs | 0.70 µs | 0.13 µs | 2.6 µs |
    | `ThinPlateSpline`, $k = 20$ | 1.51 µs | 2.35 µs | 1.79 µs | 0.24 µs | 5.9 µs |
    | `ThinPlateSpline`, $k = 50$ | 3.61 µs | 11.1 µs | 8.88 µs | 0.66 µs | 24.2 µs |
    | `Linear`, $k = 20$ | 1.47 µs | 1.01 µs | 1.40 µs | 0.22 µs | 4.1 µs |

    The search takes a third at $k = 10$ and a seventh at $k = 50$. The larger $k$, the more the $k^2$ kernels and the $k^3$ factorization dominate. Unlike the global evaluation (6.5 µs for 1000 points), it hardly depends on the number of points (5.7 µs and 6.3 µs for $k = 20$ on $10^4$ and $10^5$ points).
  * The size of a leaf (16 points), compared at 4, 8, 16, and 32, made differences within 10% (about the spread of the measurement), so it stays 16.
* **Working space** (the exception of §3.2): the matrix, the row, the permutation, `faer`'s working space, and the heap are kept once per thread (`thread_local`). They are allocated at the thread's first evaluation, and widened the first time an interpolator with a larger $k$ is evaluated. After that nothing is allocated. `tests/alloc.rs` checks for this type that "evaluations after one evaluation do not allocate".
  * The size is 5.8 KB at $k = 20$ and 20.2 KB at $k = 50$ (`ThinPlateSpline`, two dimensions; the 5.7 KB and 20 KB measured in M3-4 plus the $k$ neighbor numbers (4 bytes each) added in M3-5). The first evaluation allocates 10 times (`tests/alloc.rs` prints the count).
  * An evaluation that could not allocate its working space returns `None` ("The size of allocations and their failure" of §4; `eval` returns `Option`, so it goes in the same place as a singular neighborhood). The working space keeps what it could allocate as it is.
  * Forms not taken: on the stack (no allocation at all, but $k + R$ needs an upper bound, a bound SciPy does not have). The caller passes the working space (there would be two, an `eval_with` that does not allocate and the trait's `eval` that does).
* **Batches**: `eval_batch_into` reuses, for a point whose set of neighbors is the same as the previous point's, the system factored for the previous point as it is (M3-5, #108). Only the row and the substitutions are done per point. It uses the same numbers, so the result is bit-identical to repeated `eval`. The factorization is carried only within one batch, not to `eval` or to another interpolator (the working space is shared by the thread, so it goes back to "not held" at the start of a batch and at every `eval`).
  * It pays when the queries are denser than the data (pixels of an image, a path). Measured with `thin_plate_spline`, random points in two dimensions, a row of 1000 pixels across the unit square (per 1000 points, repeated `eval` → batch):

    | | $k = 10$ | $k = 20$ | $k = 50$ |
    | :--- | :--- | :--- | :--- |
    | $10^4$ points (10 pixels per spacing of points) | 2.30 → 1.46 ms (1.6 times) | 5.51 → 3.76 ms (1.5 times) | 24.1 → 17.9 ms (1.3 times) |
    | $10^5$ points (3 pixels per spacing of points) | 2.54 → 2.28 ms | 5.91 → 5.54 ms | 24.4 → 25.4 ms |
    | 1000 points in no order ($10^4$ points) | 2.57 → 2.52 ms | 5.90 → 5.81 ms | 24.3 → 24.5 ms |

    The spread from run to run was 2 to 8%, so the row on $10^5$ points and the points in no order are "no difference". A sequence where sets do not repeat saves nothing, but is not slower either (comparing a set looks at $k$ numbers only).
  * SciPy gathers the queries with the same set over the whole batch (they need not be consecutive). Here only the consecutive ones are used. Gathering needs an allocation proportional to the number of queries, which breaks the requirement that evaluation does not allocate (§3.2).
* **Differences from SciPy** (added to the table of M3-3):

  | Item | SciPy 1.18 | `interpx` |
  | :--- | :--- | :--- |
  | How the neighbor system is solved | solves for the coefficients, then the inner product with the row | solves for the row, then the inner product with the values. The same values (to rounding) |
  | Points at the same distance | depends on how the tree was built | the one given earlier. When points at the $k$-th distance tie, the points chosen can differ from SciPy's |
  | Singular neighbor systems | an exception | `None` |
  | Duplicate points | accepted (singular when duplicates fall in a neighborhood) | an error at construction |

### 6.4 `RegularGridInterpolator`
* **Type**: `RegularGridInterpolator<T: Scalar, const D: usize, V: Value<T> = T>`. It interpolates values on a $D$-dimensional rectilinear grid (the Cartesian product of the coordinates of each axis). `Point = [T; D]`, `Output = V`.
* **Construction**: `new(axes: [&[T]; D], values: &[V], method: GridMethod)`, and `from_owned(axes: [Vec<T>; D], values: Vec<V>, method: GridMethod)`, which takes ownership.
  * `values` is a flat array **in C order (row-major, the last axis fastest)** (the standard layout of `ndarray`, and numpy's). The value of point $(i_0, \dots, i_{D-1})$ is `values[((i_0 n_1 + i_1) n_2 + \dots) + i_{D-1}]`.
  * The shape follows from the lengths of the axes $n_d$ and is not an argument. A `values.len()` other than $\prod_d n_d$ is `LengthMismatch`. `shape() -> [usize; D]` returns it.
* **`GridMethod`** (`#[non_exhaustive]`): decided at construction and not changed afterwards (SciPy can change it at every evaluation, but the spline methods solve their coefficients at construction). M2-1 adds `Nearest` and `Linear`, M2-2 `Cubic` and `Quintic`, M2-3 `Pchip`. There is no variant for SciPy's `slinear` ("`slinear`" below).
* **Axes**: strictly ascending or strictly descending (the same as SciPy). The points per axis are at least 1 for `Nearest` and `Linear` ("Axes of one point" below), at least 4 for `Cubic`, at least 6 for `Quintic` (degree + 1; the same as SciPy), and at least 2 for `Pchip` (the same as the one-dimensional `PchipInterpolator`; SciPy's grid asks for at least 4).
  * A descending axis is reordered into ascending order at construction and kept so (the axis's coordinates reversed, and the values reversed along that axis). Evaluation is always on an ascending grid. `axes()` and `values()` return them after the reordering.
* **Interval and outside the domain (per axis)**: on each axis the interval $i_d$ of the coordinate $x_d$ ($g_{i_d} \le x_d < g_{i_d + 1}$, the top end in the last interval) is found by binary search. Outside the domain is handled by the four modes of `Extrapolate`, **the same mode on every axis, applied to each axis independently**. The rules are those of `PPoly` (§6.1) (the code of the interval search is shared). The default is `Forbid`.

  | Mode | Per axis | Result of `Linear` | Result of `Nearest` |
  | :--- | :--- | :--- | :--- |
  | `Forbid` | no value if any coordinate is outside | `None` | `None` |
  | `Extend` | the end interval, with the coordinate as it is | the formula of the end cell continued (SciPy's `fill_value=None`) | the value of the end point (the same as SciPy) |
  | `Clamp` | the coordinate clipped to $[g_0, g_{n-1}]$ | the value on the surface of the box | the same as `Extend` |
  | `Periodic` | the coordinate wrapped into $[g_0, g_{n-1})$ | the value at the wrapped point | the same |

  `Periodic` does not check the periodicity of the values (it is an extrapolation mode, not a boundary condition). A point with NaN is `None`. Infinite coordinates are as for `PPoly` (`None` under `Forbid` and `Periodic`, the end under `Clamp`, the formula evaluated as it is under `Extend`). There is no form that gives a different mode per axis (SciPy has none either).
* **`Linear`**: with the normalized distance along axis $d$, $y_d = (x_d - g_{i_d}) / (g_{i_d + 1} - g_{i_d})$, the values of the cell's $2^D$ corners are summed with the weights $\prod_d w_d$ ($w_d = 1 - y_d$ if the corner is on the lower side, $y_d$ if on the upper). The corners are summed in the order that counts axis 0 outermost (the same order as SciPy). No allocation.
* **`Nearest`**: per axis, the lower point if $y_d \le 0.5$ and otherwise the upper, returning the value of that grid point. Exactly halfway is the lower one (the same as SciPy). On an axis given descending, it is the lower one after the reordering into ascending order (SciPy too decides after reordering).
* **Hints**: `eval_hinted(&[T; D], &mut [Hint; D])`. One hint per axis, and the result bit-identical to `eval`. `eval_batch_into` looks at the points in runs of 64 and uses the hints only when every axis moves in one direction, by at most one cell per point (paths, rows of an image). Other sequences are searched point by point (on points in no order hints are slower; the same test as `PPoly`'s batches, measured by `benches/grid.rs`). Either way the result is bit-identical to `eval`.
  * **Batches along a path, against repeated `eval_hinted`** (M2-4, #90): the added cost of a batch is only the test per run. With the loop put in the caller's function and compared, the difference is 1 to 2 ns per point for `Linear` and `Nearest` (2, 3, 4 dimensions), and within the run-to-run spread for `Cubic` and `Pchip`. In the disassembly too, both call the same `eval_hinted` per point. The numbers of `benches/grid.rs` make the batch look slower (13 → 19 µs for `Linear` in 3 dimensions, 7 to 13% for `Cubic` and `Pchip`), but that is because a loop written directly in the benchmark is faster, and there is nothing to fix on the batch's side (why has not been followed).
* **Evenly spaced axes**: the interval search (`search` of §3.4, shared with the one-dimensional types) guesses, on an evenly spaced grid, the interval a point should be in by a division, and uses it after checking that the interval contains the point. If it does not, it goes back to binary search (P-2, #123). Only one interval contains the point, so a hit is the binary search's answer, and the rounding of the division does not affect the result. A shortcut that uses the guess without checking is not put in (rounding can give another interval).
  * **Whether the grid is evenly spaced is checked once at construction, and the type remembers it** (`PPoly`, `BSpline`, `Nearest1D`, the grid per axis). If every spacing between neighboring points agrees with the mean spacing within 1/16 of it, the grid is evenly spaced. This test decides only how the search starts, and does not affect the result. Comparing each point with "its position if evenly spaced" is slower by turning the point's index into a number (2 ns per point; the construction of a `Linear1D` of 100 points 292 → 491 ns), and was not taken.
  * **Guessing on every grid was not taken**: it did that at first, but evaluation on grids with unevenly spaced axes (10 to 40 points per axis) became up to 2 times slower (242 → 453 µs for 10 000 points on a 40 × 20 × 10 grid). On axes with few points, binary search ends in a few steps, and the division and check of a missed guess add on top as they are.
  * The hinted search (`eval_hinted`, and the runs where a batch uses the hints) does not guess. Starting from the previous point's interval is faster.
* **What the evaluation of `Linear` is made of, and the gap to `interpn`** (P-3, #124. A grid with unevenly spaced axes of 40, 20, 10, 5 points, 10 000 points in no order, per point. Measured by taking out one stage at a time):

  | | 2-D | 3-D | 4-D |
  | :--- | ---: | ---: | ---: |
  | Total (before P-3) | 15.1 ns | 21.3 ns | 40.4 ns |
  | of which the interval search | 9.3 ns | 13.1 ns | 16 ns |
  | of which the sum of the corners (weights and positions included) | 2 ns | 3.7 ns | 20 ns |
  | Total (after P-3) | 15.1 ns | 21.3 ns | 29 ns |

  * **The weights and positions of the corners are made by doubling per axis** (up to 4 dimensions): each axis added splits every corner there is into a lower and an upper one. Against multiplying all axes again per corner (64 multiplications in 4 dimensions), it takes 30, and no corner waits for the previous one. The weights multiply the same products in the same order, so the values do not change by a bit, and the order of summing stays SciPy's. It became 1.35 times faster in 4 dimensions, level with `interpn` (about 30 ns in the same measurement). 5 dimensions and more find them per corner as before.
  * **The times in 2 and 3 dimensions do not change with this** (4 and 8 corners, the sum being 10 to 20% of the whole to begin with; 15.1 → 15.1 ns, 21.3 → 21.3 ns in the table above). That the ratio to `interpn` in 2 dimensions moved from 1.12 to 1.16 in the comparison is run-to-run difference.
  * **In 2 and 3 dimensions the search is 60%**, and the gap remains there (1.1 to 1.5 times `interpn` in the comparison). The rest of the investigation is P-11 (#144). What is known is the check of the domain: calling the bare search (without the range test, the extrapolation mode, and the handling of no value) is about 1 ns faster per axis (15.1 → 12.6 ns in 2 dimensions). `interpn` does not make this check, and puts points outside the range in the end interval. The reason for the rest of the gap (1 to 2 ns per axis) is not known.
  * **Forms tried and not taken**: folding the corners two at a time (`interpn`'s form; fewer multiplications, but the same time, 1.00 times, and the order of summing changes from SciPy's). On short grids, "counting the points not above `x`" instead of binary search (64 points or fewer; 2 to 7% slower). Removing the branches for an axis of one point and for even spacing (no change).
* **A step of binary search is not a branch** (P-9, #132): with queries in no order, which half a point is in cannot be predicted. The search is a short loop of the crate's own, with that one step written with `core::hint::select_unpredictable`. While it was left to the standard `partition_point`, the speed of the search differed by up to 2.3 times depending on where it was inlined (the same interpolator, data, and queries evaluated from two places in one program: 328 µs and 707 µs; after the fix, both about 307 µs). Over the whole comparison too, the grid evaluations where interpx is fastest rose from 58 of 87 to 67.
  * **What the disassembly showed** (P-10, #138): #132 guessed that `partition_point` became a branch in one place and a conditional move in another. In the comparison's runner built from the code just before #137 with the rustc of the time, 1.99.0, all 34 loops of binary search are conditional moves (`ucomisd` and `cmovb`), and none is a branch. The guess was wrong, and the cause of #132's factor of two is not confirmed. In a program that evaluates one interpolator from two places (`benches/search_places.rs`; 100 points not evenly spaced, 40 000 queries in no order), the search before #137 gave the two places a gap of a tenth and more (`Linear1D` 276 / 306 µs), the search now within 3% (304 / 304 µs).
  * **A step of the search now still has the bounds check of the index** (P-10): the loop is `mov shr lea cmp jae ucomisd cmovae sub cmp ja`, and `cmp jae` is the check of `grid[middle]`. The standard `partition_point` leaves it out with `unsafe` inside, and this crate forbids `unsafe`. Measured again against the search before #137, in turn, 3 to 5 rounds, of the 6 cases that #137 made slower in the comparison, 4 are as fast as before or faster and 1 (`oned/nearest/n10/shuffled`) is within the spread. `oned/linear/n10000/shuffled-uneven` is slower: in a run of 5 rounds, 587.5 µs before #137 and 616.5 µs now, 5%. A version for the experiment with the check removed by `unsafe` (not kept), measured in turn with the version now in one run, took 556 µs against 618 µs; on 100 points, 2678 ns with the check and 2367 ns without. So the gap is the bounds check. Two forms that remove the check without `unsafe` were tried, and both were slower and are not taken: clamping the index with `min`, 772 µs against 620 µs for the version now (the `min` joins the chain of dependence of each step), and choosing between the halves of a slice cut by `split_at`, 1076 µs against 677 µs. The version now reads differently in each pair because each pair is a run of its own.
* **Axes of one point** (P-4, #125): `Nearest` and `Linear` accept an axis with only one point (the same as SciPy; SciPy's own benchmark uses a 10 × 5 × 2 × 1 grid). Along that axis the grid is constant.
  * If the query's coordinate equals that point, the value of the grid of the other axes is returned. If not, it is outside the axis: `Forbid` gives no value, `Clamp` and `Extend` the value of that point (one point has no slope to continue, so continuing is constant; the same as SciPy's `fill_value=None`), `Periodic` no value (one point is no period). NaN gives no value under every extrapolation.
  * It has no direction (ascending or descending). The low and high ends of `domain()` are the same value.
  * `Pchip` (2 points for slopes) and the splines (degree + 1 points) still do not accept it.
  * Implementation: an axis of one point does not go through the interval search (the search assumes at least 2 points). The upper side of the weights points to the same point with weight 0.
* **Errors**: the order of the checks follows the steps of §3.5.
  1. `LengthMismatch` (`values.len()` and $\prod_d n_d$. A product that does not fit in `usize` matches no array, so it is the same error)
  2. `GridAxis { axis, fault: InsufficientPoints { required, actual } }` (from axis 0 on. `required` the fewest points for the method)
  3. `GridAxis { axis, fault: NonFinite { index } }` (from axis 0 on), then `NonFinite { index }` for the values (`index` the position in the flat array)
  4. The direction of the axes (from axis 0 on): if the first 2 points are equal, `DuplicatePoints { first: 0, second: 1 }`. After that, equal to a neighbor is `DuplicatePoints`, and going against the direction the first 2 points set is `NotMonotonic { index }` (the position of the point that went against it).
  * Positions are in the order given (before reordering). The same for positions within an axis and in the flat array of values; the positions of overflow found during construction (the spline coefficients, the slopes of `Pchip`, the collocation matrix of an axis) are mapped back over descending axes before being returned.
* **Differences from SciPy**:

  | Item | SciPy 1.18 | `interpx` |
  | :--- | :--- | :--- |
  | Default outside the domain | an exception (`bounds_error=True`) | `None` (`Forbid`). `eval_or(point, fill)` corresponds to `fill_value` |
  | NaN in the values | accepted (evaluations touching NaN give NaN) | refused with `NonFinite` (§3.5. The spline methods cannot solve systems with NaN, so there is one rule whatever the method). Fill the gaps of a grid with missing values before passing it |
  | Method | can change at every evaluation | decided at construction |
  | Axes of one point | `nearest` and `linear` accept them | the same (P-4, #125. At first 2 points or more were required) |
  | Axes of one point, periodic extrapolation | no periodic extrapolation | no value off the point's coordinate |
  | Axes of one point, an infinite coordinate under `Extend` (`fill_value=None`) | the point's value | the same. On axes of 2 points or more, infinity goes into the formula continued (infinite or NaN for `Linear`, and no value). An axis of one point is constant when continued, so the result differs by axis |
  | `slinear` | a name separate from `linear` | corresponds to `Linear`. No variant of its own |
  | Coefficients of `cubic` and `quintic` | an iterative solver by default | a direct solve. Agrees with `*_legacy` |
  | Partial derivatives (`nu`) | available for `cubic`, `quintic`, and others | not provided (A-3) |
  | `Clamp` and `Periodic` | none | provided |
  | `linear` in 2 dimensions | a path of its own (another order of summing) | the same formula in every dimension. Differences within rounding |
| Points per axis of `pchip` | at least 4 | at least 2 (the same as the one-dimensional `PchipInterpolator`) |
| Evaluation of `pchip` | builds one-dimensional interpolants of whole columns at every evaluation | the slopes along the last axis are found at construction, the rest only from the values around the point. The same result |
* **`slinear`**: SciPy's `slinear` is the tensor-product spline of degree 1, the same function as `linear` (the continuation outside the domain too. Only the path of computation and the rounding differ: $7 \times 10^{-16}$ on an uneven two-dimensional grid, measured at design time of M2-2). The same function does not get two names; `GridMethod::Linear` corresponds to it. The golden cases of `slinear` are checked against `Linear`.
* **`Cubic` / `Quintic` (M2-2)**: the tensor-product B-spline of degree $k$ (3 or 5) along each axis. It corresponds to the `NdBSpline` that SciPy's `make_ndbspl` makes.
  * **Knot vectors**: per axis, the same as the not-a-knot of `BSpline1D` (§6.1; SciPy's `_not_a_knot`). Odd degrees only, so the inner knots are the data points without $(k+1)/2$ at each end.
  * **Coefficients**: found once at construction. The interpolation operators of the axes commute, so it is enough to apply, from axis 0 on, one-dimensional interpolation (solving with the collocation matrix) to every column along that axis. The collocation matrix is factored once per axis (the banded LU of M1-12), and every column is solved with the same factorization. The part of `BSpline1D` that "makes the knot vector and factors the collocation matrix" is a crate-private piece, shared.
  * **What it keeps**: the coefficients (the same shape as the grid), the knot vector of each axis, and the original values (returned by `values()`). The memory is proportional to twice the number of grid points. `Nearest` and `Linear` keep only the values.
  * **Evaluation**: per axis, the values of the $k + 1$ basis functions nonzero on the interval are found on the stack (the same recurrence as `BSpline`), and the $(k+1)^D$ coefficients summed with the products of the weights. No allocation. The interval and outside the domain are as in the table above; `Extend` continues the polynomial of the end cell (the same as SciPy's `fill_value=None`).
  * **`Periodic`**: only wraps the coordinate, and does not make the spline periodic (it is not a periodic boundary condition). At the seam the value is that of the wrapped point, and the derivatives generally jump. There is no grid spline with a periodic boundary condition (SciPy has none either).
  * **Derivatives**: for this method only, SciPy can take the order of a partial derivative (`nu`) at evaluation. Not in v0.1.0 (roadmap A-3).
  * **Errors**: the same order as M2-1. The minimum points per axis becomes degree + 1. When an axis's collocation matrix cannot be factored, `GridAxis { axis, fault: NonFinite { index } }` (`index` the row. The collocation matrix of an axis is nonsingular, so this happens only on an axis whose coordinates are so far apart that their differences overflow and entries are not finite. The cause is in the axis, so it says which axis). Overflow of the coefficients is `NonFinite { index }` (`index` the position in the flat array of coefficients).
* **The difference from SciPy 1.18's default (`cubic` / `quintic`)**: SciPy finds the coefficients of these two methods by an iterative solver (`gcrotmk`) by default. Its result is away from the tensor-product spline by about $10^{-5}$ for the values of a smooth function on an evenly spaced 2D grid, up to about $10^{-2}$ on unevenly spaced or 3D grids (measured when the goldens were made), and up to the size of the values themselves for random values (on an uneven 7×9 grid, 0.67 for `cubic` and 6.5 for `quintic`; measured at design time of M2-2). `interpx` finds the coefficients by a direct solve, so it does not agree with SciPy's default output.
  * **What it agrees with**: SciPy's `cubic_legacy` / `quintic_legacy` (an implementation that fits one-dimensional splines along each axis), and the `NdBSpline` made by applying `make_interp_spline` along each axis. Both are the same tensor-product spline, and agree with each other to rounding.
  * **Why `spsolve` is not the reference**: `solver=scipy.sparse.linalg.spsolve` gives the same spline too, but it solves one sparse matrix over all the axes together, so its condition number is the product of the axes' condition numbers. Compared with exact rational arithmetic (an uneven 7×9, random values), the error of the coefficients, relative to the largest coefficient, was $8 \times 10^{-16}$ (`cubic`) and $1 \times 10^{-14}$ (`quintic`) for the axis-by-axis solve, and $2 \times 10^{-15}$ and $5 \times 10^{-13}$ for `spsolve`. For `quintic` `spsolve` is about 40 times coarser.
  * The `expected` of the goldens is the value of `*_legacy`, and the default's value is given beside it as `expected_scipy_default` (§8.4.1). This difference is written in the rustdoc too.
* **The ways of finding the coefficients compared** (P-13, #162): five ways were measured on the same grids, values, and queries (`compare/perf/solvers.py`, commit `08dbe10`; a commit with the same content is `0b462c2` of #173. As in the comparison: one core, one thread, memory from outside, the median of 3 rounds). The values are a smooth function on grids with uneven spacing, and the elevation tile (3601 × 2401; §8.4.7). The queries are 10 000 points within the grid (the 25.91 million held-out points for the tile). The difference is the largest difference from the values of the axis-by-axis direct solve, divided by their largest value.
  * The five ways: SciPy's default (one sparse matrix over all the axes, solved iteratively by `gcrotmk`), `solver=spsolve` solving the same system directly (SuperLU), `*_legacy` (solves nothing at construction, and fits one-dimensional splines along the axes at every evaluation), the axis-by-axis direct solve (`make_interp_spline` applied to each axis in turn, making an `NdBSpline`; interpx's way of posing the problem done in SciPy, and the reference of the differences), and interpx.
  * Extracts (time of construction / memory added while building or evaluating / difference):

    | Grid | Iterative (default) | `spsolve` | Axis by axis (SciPy) | interpx |
    | :--- | ---: | ---: | ---: | ---: |
    | `cubic`, 2-D 10 × 10 | 844 µs / 1.2 MB / $6 \times 10^{-5}$ | 342 µs / 1.0 MB / $6 \times 10^{-16}$ | 136 µs / 0.8 MB | 1.59 µs / 1.4 MB / $5 \times 10^{-16}$ |
    | `cubic`, 2-D 100 × 100 | 4.42 ms / 5.8 MB / $1 \times 10^{-5}$ | 38 ms / 20 MB | 409 µs / 1.0 MB | 98 µs / 1.8 MB |
    | `cubic`, 2-D 1000 × 1000 | 748 ms / 354 MB / $2 \times 10^{-5}$ | 61.9 s / 3055 MB | 42.5 ms / 30 MB | 8.28 ms / 16 MB |
    | `cubic`, 3-D 40³ | 62.7 ms / 78 MB / $1 \times 10^{-5}$ | 70.5 s / 2322 MB | 3.55 ms / 3.5 MB | 304 µs / 2.4 MB |
    | `cubic`, the elevation tile | 13.5 s / 3362 MB / $7 \times 10^{-5}$ | out of memory (SuperLU) | 569 ms / 264 MB | 76.4 ms / 324 MB |
    | `quintic`, 2-D 1000 × 1000 | 1.02 s / 676 MB / $2 \times 10^{-5}$ | 509 s / 10377 MB | 43.5 ms / 31 MB | 10.6 ms / 16 MB |
    | `quintic`, 3-D 40³ | 174 ms / 226 MB / $6 \times 10^{-5}$ | 716 s / 7243 MB | 3.93 ms / 2.9 MB | 547 µs / 3.0 MB |
    | `quintic`, the elevation tile | 19.6 s / 6733 MB / $7 \times 10^{-4}$ | out of memory (SuperLU) | 594 ms / 265 MB | 91.9 ms / 307 MB |

    The differences of `spsolve`, `*_legacy`, and interpx are at most $1.4 \times 10^{-15}$ for `cubic` and $1.2 \times 10^{-14}$ for `quintic`: the same tensor-product spline. The whole table is `compare/perf/out/solvers/table.md` (written by `solvers.py`).
  * **Inside SciPy, the iterative solver is faster than the direct one**. In 2 dimensions from 40 points a side (at 10 points `spsolve` is 2.5 to 3 times faster), in 3 dimensions from 10 points a side. The gap widens with the grid: 83 times (`cubic`) and 500 times (`quintic`) at 2-D 1000 × 1000, 1100 and 4100 times at 3-D 40³. The memory is also 9 to 32 times less on the large grids (2-D 1000 × 1000, 3-D 40³). On the elevation tile, `spsolve` ended with SuperLU unable to allocate memory ("Out of memory during gssv"). **The price is precision**: the iterative solver's values are $1 \times 10^{-5}$ to $3 \times 10^{-4}$ of the largest value away from the tensor-product spline, and $7 \times 10^{-4}$ for `quintic` on the elevation tile.
  * **Solving axis by axis is faster than both, uses less memory, and is exact to rounding**. Even within SciPy, the axis-by-axis direct solve is 6 to 45 times faster than the default iterative solver (18 and 23 times at 2-D 1000 × 1000, 24 and 33 times on the elevation tile), and its memory is 1/12 to 1/78 on the large grids and the elevation tile. The reason is in the next item.
  * **One system, and systems per axis**: the coefficients $C$ of the tensor-product spline through every point of the grid (an array of the grid's shape) satisfy, for the values $Y$,
    $$
    (A_0 \otimes A_1 \otimes \cdots \otimes A_{D-1}) \, \mathrm{vec}(C) = \mathrm{vec}(Y)
    $$
    where $A_d$ is the collocation matrix of the one-dimensional spline of axis $d$ ($n_d \times n_d$, with a bandwidth about the degree $k$) and $\otimes$ is the Kronecker product. SciPy's default and `spsolve` make the matrix on the left ($N = \prod_d n_d$ rows, $(k+1)^D$ nonzeros per row) as one sparse matrix and solve it. But the inverse of a Kronecker product is $(A_0 \otimes \cdots \otimes A_{D-1})^{-1} = A_0^{-1} \otimes \cdots \otimes A_{D-1}^{-1}$, so the same $C$ is obtained by solving each column of $Y$ along axis $d$ with $A_d$, for $d = 0, 1, \ldots, D-1$ in turn. That is the axis-by-axis direct solve (interpx; SciPy's `make_interp_spline` applied to each axis; what `*_legacy` does at every evaluation). Only banded matrices of $n_d$ unknowns are solved, and a banded LU costs in proportion to $n_d$ per column, so the whole cost and memory are proportional to $N$. Solved as one sparse matrix, this structure goes unused. The direct solve gains nonzeros during the factorization, and its memory and time grow (3 GB at 2-D 1000 × 1000 in the table above; impossible to allocate on the elevation tile). The iterative solver multiplies by $N \times (k+1)^D$ nonzeros per iteration (16 times the number of points for a 2-D `cubic`, 64 times in 3-D). The number of iterations was not measured.
  * **The conclusion of SciPy's PR and the conclusion of this section do not contradict each other**: what the PR compared was solving one sparse matrix directly or iteratively, and in that comparison the measurement above too finds the iterative solver faster (83 and 500 times at 2-D 1000 × 1000). What this section added is the comparison of solving one system over everything against solving axis by axis. Split by axis, the direct solve (banded LU) becomes faster and more exact than either way of solving the one system. Most of the gap is the posing of the problem, not the kind of solver.
  * **What the gap in construction between interpx and SciPy is made of** (`cubic` at 2-D 1000 × 1000): the default iterative solver 748 ms → the axis-by-axis direct solve 42.5 ms (the posing, 18 times) → interpx 8.28 ms (the implementation, 5.1 times). Changing only the kind of solver (iterative or direct) makes it slower, keeping the one system (61.9 s). Most of the gap is not the solver but solving all the axes together or axis by axis. On small grids the implementation's part is large (at 10 × 10, 1.59 µs against SciPy's 136 µs axis by axis; the cost of Python's calls).
  * **Evaluation**: the four ways that keep coefficients take about the same time (SciPy's three evaluate the same `NdBSpline`). `*_legacy` solves at every evaluation, so 10 000 points take 0.45 to 2 seconds (120 to 450 times the forms with coefficients). On the elevation tile it ended with `MemoryError` (no message).
  * **Why SciPy made the iterative solver the default**: from the description of the SciPy PR that put `NdBSpline` into `RegularGridInterpolator` (scipy/scipy#19633, January 2024). It makes the design matrix as one sparse matrix over all the axes and solves it. The direct `spsolve` (SuperLU) is "extremely slow and memory hungry for high dimensions and large degrees", so the iterative solver was made the default as "the classic answer". On choosing `gcrotmk` it says "No specific preference from me". It may not converge, so the user can change the solver and its arguments (`solver`, `solver_args`). No mention was found, in the PR or its discussion, of considering solving axis by axis (using the structure of the Kronecker product). The measurement above confirms the slowness and memory of `spsolve` that this PR describes, and confirms that the iterative solver is faster than it. Beyond that, solving axis by axis needs neither.
    * Why SciPy made it one system is not written in the PR. The following is inference without a source: the machinery that makes and solves the matrix over all points can also be used as it is for data not on a grid and for weighted least-squares fitting, where the form of a Kronecker product does not hold. In SciPy, least-squares multidimensional splines have since been advancing on the same machinery (scipy/scipy#25472, #25501). Restricted to interpolation through every point of a grid, solving axis by axis is better, and that interpolation is all interpx handles.
  * **The solver is not the user's to choose** (the owner's decision. The proposal came from the owner: that the iterative solver or the direct solve over all axes might be selectable too). On every grid measured, the axis-by-axis direct solve (interpx) is 45 to 720 times faster than SciPy's default iterative solver, faster still than the direct solve over all axes, uses less memory, and is exact to rounding. There was no size at which a choice would gain the user anything. SciPy's default values (the iterative solver's, about $10^{-5}$ away from the tensor-product spline) are not reproduced ("The difference from SciPy 1.18's default" above in §6.4).
    * **Conditions not measured**: 4 dimensions or more, grids whose axes have very different numbers of points (such as 10000 × 10), an iterative solver with a looser tolerance (`solver_args`), and solving again with only the values changed on the same grid (the previous solution could be the starting point; interpx has no entry to swap the values). One iteration alone costs a multiplication by the matrix over all axes of points × $(k+1)^d$ (16 times the number of points for a 2-D `cubic`), and the cost of the axis-by-axis direct solve is points × axes × a small constant, so being overtaken in these is unlikely. But this is inference, not measurement.
* **`Pchip` (M2-3)**: the one-dimensional PCHIP (`PchipInterpolator` of §6.1) applied along each axis in turn. It corresponds to SciPy's `pchip`.
  * **The order of the axes**: from the last axis to axis 0. First, each column along the last axis is interpolated by one-dimensional PCHIP to find its value at the coordinate $x_{D-1}$. Those values, with one axis fewer, are interpolated the same way along axis $D-2$, and so on to axis 0. PCHIP's slopes are not linear in the values, so another order of the axes gives another interpolant (on an uneven 6×7 grid of random values, swapping the axes changes the values by $6.5 \times 10^{-3}$; measured at design time). The order is SciPy's.
  * **Decided by the neighborhood of the point only**: one-dimensional PCHIP on the interval $[x_i, x_{i+1}]$ is decided by the 4 values $y_{i-1}, \dots, y_{i+2}$ only (3 on an end interval, 2 on an axis of 2 points), because the slopes at its two ends look only as far as the neighboring points. So evaluating at one point needs at most 4 values per axis, $4^D$ in all. SciPy builds the interpolation of whole columns at every evaluation, but the same result comes from these values alone.
  * **What it keeps**: the values, and the slopes along the last axis (the same shape as the grid). The memory is proportional to twice the number of grid points (the same as `Cubic`). The first fold (the $4^{D-1}$ columns along the last axis) is over the original values, so its slopes do not depend on the point evaluated and can be found at construction. The folds of the axes before it work on folded values that differ by point, so they are found at evaluation.
    * The grounds for the decision (measured with a prototype at design time; the grids of `benches/grid.rs`, the time to evaluate 1000 points; best of 3 runs each, in turn; path / in no order):

      | | 2-D (200 points × 2) | 3-D (40 points × 3) | 4-D (10 points × 4) |
      | :--- | :--- | :--- | :--- |
      | Without slopes (values only) | 90.8 µs / 93.6 µs | 267 µs / 294 µs | 843 µs / 903 µs |
      | With the slopes along the last axis | 56.3 µs / 59.3 µs | 145 µs / 154 µs | 427 µs / 457 µs |
      | Ratio | 0.62 / 0.63 | 0.54 / 0.52 | 0.51 / 0.51 |

      The form without slopes finds two slopes again for each of the $4^{D-1}$ columns at every evaluation.
    * Keeping even the coefficients of the cubics along the last axis (4 per interval) is not taken. The memory would be 5 times the values. What it would save is only making the coefficients of the cubics along the last axis, and the whole computation along the last axis is about a third of an evaluation (about 53 ns of 156 ns in 3 dimensions).
  * **Evaluation**: the interval is found per axis ("Interval and outside the domain" above), and on the $4^{D-1}$ columns the interval's cubic is made from the slopes and values kept and evaluated. Those results are folded from axis $D-2$ to axis 0, 4 at a time per axis (3 or 2 at the ends). A fold finds the secants, the slopes at both ends, and the cubic from those 4 values. The working arrays are on the stack, with no allocation.
    * The formula of the slopes, the cubic of an interval, and its evaluation are those of `PchipInterpolator`, made a crate-private piece and shared. The result is bit-identical to `PchipInterpolator` applied along each axis in the same order (the tests say so too).
    * The formula of the slopes stays in SciPy's form (the weighted harmonic mean with 4 divisions, branching on the sign test). A form with 2 divisions and a form with a select instead of the branch were tried, but none was faster in every dimension (measured at design time; 1000 points in no order, best of 3):

      | Formula of the slopes | 2-D | 3-D | 4-D |
      | :--- | :--- | :--- | :--- |
      | Branch, 4 divisions (SciPy's form) | 61.5 µs | 171 µs | 480 µs |
      | Branch, 2 divisions | 56.2 µs | 157 µs | 457 µs |
      | Select, 4 divisions | 52.4 µs | 172 µs | 559 µs |
      | Select, 2 divisions | 51.9 µs | 163 µs | 521 µs |

      With a select, the divisions are done even at points whose slope is 0 (about half for random values), so it is slower in 4 dimensions with many columns. The difference of the 2-division forms (5 to 9%) is within the spread of repeating the same version (10% or more in the measurements of the same table). There is no reason to change a form that matches SciPy down to the order of operations.
    * On smooth values the slopes do not become 0, and the divisions are done at every point. That was measured too (measured in M2-3; `pchip_smooth` of `benches/grid.rs`, 1000 points in no order, best of 3). The 4-division form took 55.7 µs / 163 µs / 522 µs, the 3-division form ($(w_1 + w_2) / (w_1/m_{i-1} + w_2/m_i)$) 55.1 µs / 151 µs / 484 µs, the 2-division form 56.1 µs / 158 µs / 489 µs: within 7%. Construction (400×400) was within 7% too. The formula stays.
    * At construction the reciprocals of the interval widths are found once per axis (finding them again per column would add as many divisions as those of the slopes). On a 400×400 grid, 2.09 ms became 1.65 ms. The construction of `PchipInterpolator`, which uses the same piece, became 9 to 19% faster too (9.3 µs became 7.5 µs on 1000 points).
  * **It preserves shape**: within the domain, the result lies between the smallest and largest of the values at the $2^D$ corners of the cell containing the point (but for rounding). One-dimensional PCHIP stays between the values at the two ends of its interval, and every fold keeps that. It is flat where the values are flat. `Cubic` and `Quintic` do not have this property. It does not extend beyond the domain (`Extend`).
  * **Outside the domain**: the default is `Forbid`. `Extend` continues the cubic of the end interval per axis (the same as SciPy's `fill_value=None`). `Clamp` and `Periodic` map the coordinate, then use the same formulas.
  * **Overflow**: what is found at construction is only the slopes along the last axis, so only those can be checked. A secant or slope along the last axis that is not finite is `NonFinite { index }` (`index` the position in the flat array of values). The secants, slopes, and coefficients of the axes before it are found at evaluation, so an overflow there makes the result a value that is not finite (returned as IEEE 754 gives it; no check is added). It happens only when the differences of values are close to the largest value of `T`.
  * **Errors**: the same order as M2-1. The minimum points per axis is 2. After that, the `NonFinite` above.
  * **What it does not have**: partial derivatives (SciPy cannot take them for `pchip` either). Complex values (not in `Value`; deprecated in SciPy too).

### 6.5 Serialization (the `serde` feature)
* The vertex array, neighbor array, values, and parameters are all serialized (no re-triangulation), so that saved data is not tied to `convx`'s version.
* The stored form carries a format version. Reading checks lengths and index ranges, and returns `Err` if they are wrong. For `LinearNDInterpolator` it also checks that the vertices of a simplex are distinct and that the neighbors are symmetric (the neighbor $t$ of a simplex $s$ has $s$ among its neighbors; the review of #183). The Delaunay property of the stored topology is not checked again.
* **The stored forms of the types there are** (M5-3, #154. Decided by the owner in the Grill of 2026-10-13):
  * **What is stored is the computed state.** Reading does not solve the linear systems again. Each type has a private struct for storing (a version `version` and the content), and the interpolator reads and writes through it.
  * **What is cheap to rebuild is not stored.** At reading it is rebuilt from the stored content. What is not stored cannot disagree with the content.

    | Type | What is stored | What is rebuilt at reading |
    | :--- | :--- | :--- |
    | `PPoly` | breakpoints, coefficients, degree, extrapolation | whether evenly spaced |
    | `BSpline` | knots, control points, degree, extrapolation | the knot intervals of the domain, whether evenly spaced |
    | `CubicSpline`, `PchipInterpolator`, `CubicHermiteSpline` | the piecewise cubic `PPoly` | (as for `PPoly`) |
    | `Akima1D` | the piecewise cubic `PPoly`, the method | the same |
    | `Linear1D` | the piecewise linear `PPoly`, the last data value | the same |
    | `Nearest1D` | coordinates, values, kind, extrapolation | whether evenly spaced |
    | `BSpline1D` | the `BSpline` | (as for `BSpline`) |
    | `BSplineCurve` | the `BSpline1D`, the parameters of the data points | the same |
    | `RegularGridInterpolator` | the axes (ascending, as kept), values, method, extrapolation, what the method solved (coefficients for `cubic` and `quintic`, slopes for `pchip`) | strides, the spline knots, whether each axis is evenly spaced |
    | `NearestNDInterpolator` | the tree's points (one per coordinate; after the change of coordinates for `rescaled`), the position of each point's value, the values, the change of coordinates, the distance bound | the k-d tree |
    | `LinearNDInterpolator` | the points (after the change of coordinates for `rescaled`), the values, the vertices of the simplices, their neighbors, the change of coordinates | the LU decomposition of every simplex's matrix, the box of the points, the starts of the cells |
    | `CloughTocher2DInterpolator` | what `LinearNDInterpolator` stores, and the gradient at each point | what `LinearNDInterpolator` rebuilds, and `g[k]` of every triangle (§6.3) |
    | `RBFInterpolator` | points, kernel, smoothing, degree of the polynomial, number of neighbors; when solved over all points the coefficients and the estimate of the reciprocal condition number, with neighbors the values | the powers of the monomials, the change of coordinates of the polynomial, the k-d tree of the neighbors |

    * The k-d tree is not stored so that the stored form is not tied to the internal order of a version of `kiddo` (the same reason as for `convx`). A tree rebuilt from the same points is the same tree.
    * The last data value of `Linear1D` comes out of the segment's coefficients off by rounding, so it is stored.
  * **Reading makes every check of the constructor.** Lengths, finiteness, increasing order, agreement with the degree, duplicate points, the ranges of parameters. The same as the decision not to add entries without checks (§3). If wrong, it is an error of the deserializer (no panic). Checks added per type:
    * The types that wrap a `PPoly` do not read one of another degree than theirs (3; 1 for `Linear1D`). `BSpline1D` does not read degree 0.
    * The parameters of `BSplineCurve` are finite and increasing, with their ends at the ends of the spline's domain.
    * The grid's axes are ascending (an axis given descending is turned over, values and all, at construction and kept ascending). What the method solved is checked for length and finiteness.
    * For `NearestNDInterpolator`, the tree's points are distinct, the positions of values increase and fall within the values, the change of coordinates is finite with positive scales, and the distance bound is finite.
    * For `RBFInterpolator`, the number of neighbors is less than the number of points (equal or more makes construction one system over all points), and the number of coefficients is the number of points plus the number of monomials.
  * **What is not checked**: that the stored coefficients were solved from some data (without solving again it cannot be known). What is checked is that the interpolator can be evaluated.
  * **Version**: `version: u32` per type (now 1). An unknown version is an error. The form of another version is not guessed.
  * **How scalars and values are written**: values (scalars or fixed-size arrays) and points are written with a private trait (`stored::Stored`), not `serde`'s `Serialize`. `serde` cannot handle arrays of any length (only up to 32). Also, making `Serialize` a bound of `Scalar` and `Value` would change the public traits with the feature. `Stored` is a supertrait of the sealing trait, and the definitions of `Scalar` and `Value` are the same with or without `serde`. Fixed-size arrays are written as tuples without their length (`[[1.0,-1.0],[0.0,5.0]]` in JSON).
  * **Numbers that are not finite are not stored.** Every floating-point number in a stored form is finite (checked at reading). JSON has neither NaN nor infinity. Only the distance bound of `NearestNDInterpolator` can be a value that is not finite in the interpolator. Positive infinity is the same as "no bound", so `with_distance_bound` keeps it as no bound (`distance_bound()` returns `None`. It is not made a special case on the storing side so that two interpolators that return the same values are the same thing). NaN and negative infinity make an interpolator with no value for any query, and storing it is an error (JSON writes both as `null`, which reads back as "no bound", that is, another interpolator with a value for every query). The first version stored infinity of either sign as "no bound", and an interpolator with negative infinity began returning values once stored (the review of #155).
  * **JSON and floating point**: by default `serde_json` can read a floating-point number off by the last bit. To read back an interpolator that is the same to the bit, `serde_json`'s `float_roundtrip` feature is needed (the tests enable it; written in the crate documentation). Binary formats (such as `postcard`) come back as they were.
  * **Scope**: every interpolator there is, and their public parameters (`Extrapolate`, `GridMethod`, `NearestKind`, `AkimaMethod`, `Parametrization`, `EndCondition`, `BoundaryCondition`, `RbfKernel`, `RbfDegree`, `RbfOptions`). The error types and `Hint` are not included. The types of M4 (`LinearNDInterpolator`, `CloughTocher2DInterpolator`) get `serde` in their own rows.
  * **Tests** (`tests/serde.rs`): every type goes round trip through JSON and `postcard`, checking that `eval`, `eval_batch`, and `eval_batch_into` are the same to the bit, and that the stored form of what was read back is the same string as the original. Broken stored forms, by breakages common to every type (every array one shorter, every version 2, every floating-point number NaN) and per type (coordinates swapped, points repeated, an axis turned over), are each checked to be errors. NaN cannot be written in JSON, so the bytes of that number in `postcard`'s byte string are replaced.

---

## 7. Crate Features

```toml
[features]
default = ["scattered"]

# 1D / curve / grid are always built (no heavy dependencies)
scattered = ["dep:convx", "dep:kiddo", "dep:faer"]   # without it, the scattered types and heavy dependencies go away

# opt-in
rayon   = ["dep:rayon"]
serde   = ["dep:serde"]
```
* **No `ndarray` feature in v0.1.0** (#146. The owner's decision). The API takes slices, so an `ndarray` array can be passed as it is with `as_slice()` (the order of the grid's values is the standard layout of `ndarray`. Results can be written into the user's array with `eval_batch_into`). A dedicated entry would add only the convenience of taking an `ArrayView` directly, and does nothing for speed. Taking a view without copying would need interpolators that borrow their data, which changes the decision on ownership (§3); it is decided once a use that needs it appears (roadmap A-12).
* There is no `std` feature (always `std`). `convx` / `faer` / `kiddo` are not exposed as separate public features.
* What `scattered` contains is added by the rows that use the dependencies (the above is the form of v0.1.0). In M3-1 it was `scattered = ["dep:kiddo"]`; `faer` is added by M3-3 (`["dep:kiddo", "dep:faer"]`), and `convx` by M4.
* `kiddo` is `default-features = false` (dropping its default `fixed` (fixed-point coordinates) and `multi-threaded` (`rayon`). Parallelism comes in M5-1 from `interpx`'s `rayon` feature). In this form the transitive dependencies are 29 crates (at design time of M3-1).
* `faer` is `default-features = false, features = ["std"]`.

---

## 8. Quality Assurance (Testing, Precision, Benchmarks)

### 8.1 Tests
* **Golden data**: generated by the Python scripts of `compare/` (`compare/generate_*.py`, with SciPy 1.18 pinned by `uv` in `pyproject.toml` / `uv.lock`) and committed to `compare/goldens/*.json`. `cargo test` only reads them and does not call Python; CI does not run SciPy either. They are regenerated with `just gen-goldens` (added in M1-1). `compare/goldens/` is included in the published package (`include` of `Cargo.toml`).
* **Zero-allocation check**: `tests/alloc.rs` checks with a counting allocator that `eval` / `eval_or` / `eval_batch_into` allocate 0 times (part of `just test`). The counter counts over the whole process, so this file runs without the test harness (`harness = false`, a `main` function). The harness has threads of its own, and allocations made there during a measurement would be counted too. The parallel `eval_batch_par` is tested to be bit-identical to the serial form.
* **`proptest` property tests**: agreement at the data points, the monotonicity of PCHIP, exact reproduction of linear functions by `LinearND`, the continuity of CubicSpline's derivatives, and so on.

### 8.2 Numerical precision criteria
| Subject                          | `rtol` / `atol`   |
| :------------------------------- | :---------------- |
| Direct methods (f64)             | `1e-12` / `1e-12` |
| CloughTocher (goldens of SciPy iterated to convergence / at SciPy's defaults) | `1e-12` (goldens iterated to convergence), `1e-7` (goldens at SciPy's defaults); both relative to the largest absolute value (§6.3) |
| RBF (the global solve)           | the test's constant (below) |
| f32                              | `1e-5`            |

* The values of RBF are the solution of a linear system, solved another way than SciPy (§6.3). The difference is proportional to the system's condition number, so `1e-12` is not applied across the board. In the prototype at design time, the differences over the 72 global cases of `rbf.json` were at most $9 \times 10^{-13}$ relative to the largest absolute value (`Gaussian`, 64 points). The tests of M3-3 compare as a ratio to the size of the values, with the tolerance a named constant with its measurement. If a case does not meet it, exact rational arithmetic decides which is right (the same way as for the coefficients of M2-2).

### 8.3 Benchmarks
* The root crate is a single crate (no workspace). The materials for comparison are in `compare/` at the root (`pyproject.toml`, `uv.lock`, `generate_*.py`, `goldens/`, `perf/`). The Rust runners of `compare/perf/` are independent Cargo packages that depend on other crates (each with its own `[workspace]`), so that the libraries compared do not mix into `interpx`'s own dependencies.
* The comparison is run by hand and is not a gate of CI (the raw results in `compare/perf/out/`, ignored by git). The summary is `docs/performance.md`, written by `compare/perf/report.py` (not edited by hand). The `criterion` benchmarks are in `benches/`.
* **The place of the comparison** (M5-5, #120): a library slower than the existing ones has no reason to exist. The comparison comes before the features that sit on top of the types (`rayon`, `serde`), and where interpx loses to a library computing the same thing, it is fixed before those features (what cannot be fixed is written here with its numbers and reason). When types are added (M4), their row adds them to the same machinery.
* **How it runs**: one process per case and library, pinned to one core, on one thread. Every library reads the same data files. The libraries take turns case by case, and the whole is gone through 3 rounds. Each process goes through the same phases (read the data → build once → evaluate once and write the values → time construction and evaluation, each repeated).
* **How long it takes**: decided not by the work measured but by the number of processes (cases × libraries × rounds), the time given to each timed loop (50 ms each), and the start of Python and SciPy (0.7 s a process). At first each loop was given 250 ms, there was a wait of 60 ms after every phase of every case, and pairs of a case and a library that does not support it were started too: all cases took 2 hours (fixed in P-1).
* **Speed**: the time of one call (construction, and the evaluation of all the case's queries, separately). interpx's serial `eval_batch_into`, the others' call with an array. A call shorter than 200 µs is repeated within one sample and divided (the clock has steps of 100 ns). The median of a process's samples is taken, and the median of its 3 rounds reported.
* **A warm-up before construction is timed** (P-18, #172): before construction is timed by repeating it, constructions are made, untimed, for 1 second (`WARM` in the shared Rust part, and the SciPy runner). For a while after a runner has read its case and evaluated once, a construction is slower. Measured in windows of 25 ms, one construction of the grid's 3-D `linear` at `max40` (8000 values) went down 25 → 17 → 4 µs and then stayed at 3 µs. Repeating only the construction, without reading the queries, gives 3 µs from the start. In the slow time a construction has about 9 page faults; once it has settled, 0.07. The timed loop of 50 ms fell within that time, by a length that differed from process to process, so the same construction took one of two times depending on the process (`CubicSpline` at 1000 points: 15 µs and 38 µs). How long the slow time lasts depends on the case. Measured in the runner in windows of 50 ms for 3 seconds, it was about 150 ms for the grid's 2-D `linear` at `max200` and for `CubicSpline` at 1000 points, about 250 ms for `Nearest1D` at 10 000 points, and about 450 ms for `Linear1D` at 10 000 points (going down from 99 µs to 21 µs). The warm-up of 250 ms put in first did not get `Linear1D` at 10 000 points past it (in the run of the whole comparison `7c7bbbd`, medians of the rounds 44 to 63 µs, fastest 21 to 24 µs), and 1 second was chosen as more than twice the longest time. With a warm-up of 250 ms in the runner, over 6 processes each, `CubicSpline` at 1000 points went from 34.5 to 41.5 µs to 15.5 to 17.0 µs, the grid from 15.1 to 22.8 µs to 2.9 to 3.2 µs, and `Akima1D` at 1000 points from 11.5 to 20.2 µs to 8.3 to 9.1 µs. The cause of the slow time is not in interpx: with fewer allocations in interpx (a prototype in which a thread keeps the working space of the construction), `CubicSpline` after the warm-up was 15.5 to 17.1 µs, unchanged. Why the page faults happen (whether the Windows heap gives freed space back and takes it again) was not checked. Evaluation is not warmed up (it is timed after construction, past that time).
* **Memory**: the resident set size (RSS) of the process is compared as **a series over time**, read from outside every 1 ms. With the times of the ends of the phases laid over it, the point the data was read, the highest during construction, after construction, and the highest during evaluation are tabulated and drawn. It covers the large cases (`MEMORY` of `cases.py`), and only for those cases there is a wait of 60 ms after each phase, so that the phase shows in the samples (small cases are a few MB for every library, with no difference). A single number for the peak does not tell when and for how long the memory was held, and cannot compare Python and Rust by the same measurement. The counting allocator is used for the check "does it allocate" (§8.1), not for comparison between libraries.
* **Checking the values**: the values each library writes are compared with SciPy's. A library that differs beyond rounding computes another interpolant under the same name (another end condition, a local cubic instead of a spline, and so on). Its time is in the tables, with a note that it is not counted as a loss.
* **The Rust crates compared** (looked up on crates.io when M5-5 started: maintained, pure Rust, computing the same things as this crate):

  | Crate | Version | Scope | Why chosen |
  | :--- | :--- | :--- | :--- |
  | `scirs2-interpolate` | 0.6.5 | one variable, grids, RBF | aims to port SciPy; the widest Rust crate |
  | `ndarray-interp` | 0.6.0 | one variable (linear, cubic spline), bilinear in two dimensions | the most used for one variable |
  | `interpn` | 0.11.2 | rectilinear grids (linear, nearest, cubic) | dedicated to grids, and claims speed |
  | `spade` | 2.15.1 | the Delaunay triangulation of the plane, and interpolation on it | computes what `LinearNDInterpolator` does in 2-D (barycentric interpolation on the Delaunay triangulation; looked up in M4-1). `scirs2-interpolate` computes it too, with `triangulation_interp::make_linear_triangulation` (2-D, Bowyer–Watson, a tolerance of 1e-10), and is compared for `LinearND` in 2-D (its values are SciPy's to 1e-15 relative). At first only its `scattered_2d`, which is not public, was looked at, and it was left out (found and put right in the work of #183). The linear method of its `griddata` is an approximation by natural neighbors, another thing |

  * **The current results** (all 317 cases measured in one run at the last commit of M4-2's code, `a2498de`, which adds `CloughTocher2DInterpolator`, moves the location of points into `Mesh`, shared with `LinearND`, inlines its `find`, and makes the gradient estimate faster. Construction is warmed up for 1 second before it is timed (P-18). 3 rounds; a pair whose call took 2 seconds or more and was timed by one call, 1 round. `docs/performance.md` is the result of this run. The check of the driver (`just compare-check`) was passed before the run):

    | Area | Operation | Cases | interpx fastest | Slower than SciPy | Slower than a Rust crate |
    | :--- | :--- | ---: | ---: | ---: | ---: |
    | One variable | evaluation | 112 | 112 | 0 | 0 |
    | One variable | construction | 112 | 76 | 0 | 36 |
    | `PPoly`, `BSpline`, curve | evaluation and construction | 10 | all | 0 | 0 |
    | Grid | evaluation | 92 | 74 | 0 | 18 |
    | Grid | construction | 92 | 29 | 38 | 46 |
    | `NearestNDInterpolator` | evaluation and construction | 16 | all | 0 | 0 |
    | `LinearNDInterpolator` | evaluation | 19 | 19 | 0 | 0 |
    | `LinearNDInterpolator` | construction | 19 | 3 | 0 | 16 |
    | `CloughTocher2DInterpolator` | evaluation | 7 | 7 | 0 | 0 |
    | `CloughTocher2DInterpolator` | construction | 7 | 5 | 0 | 2 (a library with other values; not counted) |
    | `RBFInterpolator` | evaluation and construction | 61 | all | 0 | 0 |

    * **No evaluation is slower than SciPy** (by the threshold). The speed of evaluation against SciPy (SciPy's time / interpx's time, the median of the per-case ratios): one variable 5.6 times, `PPoly` 5.7 times, `BSpline` 3.3 times, curve 3.3 times, grid 7.9 times, `NearestND` 3.1 times, `LinearND` 1.8 times, `CloughTocher2D` 2.0 times, RBF 1.9 times. `LinearND`'s queries at the pixels of an image: 1.20 times at 1000 points and 250 000 pixels (interpx 5.85 ms, SciPy 7.04 ms; in the previous run `58c8390` 1.02 times, interpx 6.50 ms. Between the two runs `LinearND`'s evaluation moved by a median of 0.995 over its cases, 0.90 to 1.12, and this case is at the faster end). For queries in no order the gap is wide (a 27th of SciPy's time at 100 000 points in 3-D). One variable going from 6.0 to 5.6 times against the previous run is not a change of the cases' times (the median per-case ratio is 1.006 for interpx and 1.008 for SciPy): the ratios of the 112 cases are sparse near their median (5.5 to 6.5 times), so the median moved.
    * **The remaining losses in evaluation are 18 cases (18 comparisons)**, all on the grid against `interpn` (1.11 to 1.49 times; P-11, #144). In the cases on real data there is no loss in evaluation.
    * **Losses in construction remain** (36 cases in one variable, 63 on the grid, 16 of `LinearND`). No construction in one variable is slower than SciPy. The reasons are as in §3's "The differences that remain, and why". `LinearND` loses to `spade` in every 2-D case, by 1.58 to 3.32 times, and in the 3 cases of 10 points to `scirs2-interpolate` as well, by 2.83 to 2.85 times (5.70 to 5.75 µs against 2.00 to 2.03 µs): 70 to 90% of its construction is `convx`'s triangulation (§6.2). `convx`'s R2-1 (YUKIKEDA/convx#667) made its 2-D build about level with `spade`, to be published as 0.1.1; taking it and measuring again is P-21 (#186), which also looks into interpx's own part (the LU decomposition of every simplex, the cells). Against SciPy the construction is 4.4 times faster in the median. The grid is slower than SciPy in many cases even with `from_owned` (the cost of the checks. §3, #146).
    * **The construction of RBF over all points got slower by #170 (the estimate of the condition number)** (the ratio of P-17's run `c5fe29c` to the run before it): 1.7 to 2.3 times at 10 points, 1.35 to 1.63 times at 100, 1.25 to 1.33 times at 1000, 1.14 times at 3000. At the 10 838 weather stations 0.93 to 1.07 times, indistinguishable from the spread of a single measurement. No case is counted as slower by the threshold (1.10 times), but at 100 points it is level with SciPy (the median ratio interpx / SciPy 0.83 to 1.05; before, 0.6 to 0.7). Where interpx is slower in the median are `gaussian` (1.05) and `multiquadric` (1.05), and for both the fastest repeat is interpx's (93.5 µs and 94.9 µs, 73.8 µs and 74.9 µs), so by the way of counting too they are not counted as slower. The owner accepted the cost (§6.3).
    * **SciPy's construction of `thin_plate_spline` at 3000 points takes 1.75 to 2.3 seconds a call, straddling the line of long calls (2 seconds)** (the review of #171; in P-17's run). In this PR's first run (`eff8242`), the first call of round 2 took 1.99 seconds and was timed several times, and those timings were slow at 2.2 to 3.7 seconds. The rule for running one round only decided by the median, so round 3 did not run, and the value reported was 2.35 seconds, the median of 1.79 and 2.91 seconds, which read as 43% up from the previous run's 1.64 seconds. It is not that the first call is slower than the later ones (5 constructions in a row in one process put the first at 1.02 to 1.12 times the rest). Now that one round only is limited to pairs whose call was timed by one call, it is 1.94 seconds, the median of 3 rounds (1.94, 2.31, 1.75 seconds). Round 2 is slow again, and that this call moves by 20 to 30% from round to round remains.
    * **`LinearND` on real data** (10 838 weather stations, 1204 queries): construction interpx 7.96 ms, SciPy 30.6 ms, `spade` 5.04 ms, `scirs2-interpolate` 513 ms; evaluation interpx 81 µs, SciPy 1.18 ms, `spade` 1.36 ms, `scirs2-interpolate` 147 ms. The errors are the same for all (largest 22.1 degrees, root mean square 2.83 degrees, 1 query outside the hull with no value).
    * **`CloughTocher2D`** (M4-2): evaluation is the fastest in all 7 cases (1.9 to 9.6 times faster than SciPy; at 1000 points and 250 000 pixels interpx 13 ms, SciPy 24.6 ms, `scirs2-interpolate` 797 ms). Construction is 2.7 to 53 times faster than SciPy's (156 ms and 542 ms at 100 000 points). SciPy stops at its default `tol=1e-6`; interpx solves through (§6.3). `scirs2-interpolate`'s construction is faster than interpx's at 10 and 50 points (3.34 µs and 9.68 µs, 28.6 µs and 39.1 µs), but its values are those of another interpolant, 0.12 of the largest value from SciPy's, so these are not counted as losses. Its evaluation takes 23 seconds at 10 000 points and 296 seconds at 100 000; the latter is the longest run of this run (41% of the limit).
      * **Values**: in the 6 synthetic cases within $1.7 \times 10^{-8}$ of the largest value from SciPy's (what SciPy's `tol` leaves). On the weather data $1.3 \times 10^{-2}$ of the largest value apart. One edge of the triangulation differs (in 21 659 triangles): its 4 stations lie on one circle in their decimal coordinates, and rounded to `f64` the predicate is off by $4 \times 10^{-17}$. `convx` decides exactly on the rounded coordinates and returns the division that meets the Delaunay condition; Qhull returns the one that misses it by that much. The gradients are the solution of one system through all edges, so the one edge moves the values at queries far from it. SciPy's gradients are nearly the exact solution on Qhull's division (within $2.5 \times 10^{-9}$ in value of a direct solve; with `tol=1e-14` and `maxiter=10^6` SciPy warns that it has not converged, and the values do not change).
      * **The errors on the weather data are 699 degrees at most and 28 degrees in root mean square**, the same as SciPy's (`scirs2-interpolate`: 22.3 and 2.89 degrees). Pairs of stations 0.0003 to 0.0005 degrees apart differ by 3 to 4 degrees, and the gradient estimated there is about 5600 degrees Celsius a degree of longitude or latitude (the median magnitude is 11 a degree). That is the method, not a fault of a library.
      * **In the one case of 100 000 points, 3 of interpx's 9 processes took 2 to 3.5 times as long** (construction and evaluation, the fastest repeat too). The `from_owned` column's 296 ms and 41 ms come from this; the interpx column (156 ms, 21.4 ms) is not affected, being the median of the rounds. In the same rounds interpx was not slower on the other large cases (`LinearND` at 100 000 points, RBF at 3000), and SciPy on the same case was 1.2 to 1.4 times slower in round 2 only. 20 processes alone and 9 through the driver were all on the fast side (construction 145 to 171 ms, evaluation 20 to 23 ms): it does not reproduce. The cause is not known. It may be the phenomenon of P-19 (#181), a time that depends on the process, so the figures were added to #181.
    * **Against the previous run (M4-1's `58c8390`)**: the losses went 216 → 214 (against libraries with the same values, 206 → 202); 13 left and 11 entered. Of the 11, 2 are `CloughTocher2D`'s construction against `scirs2-interpolate` (other values). The rest are comparisons on the grid and in one variable whose ratio moved between 1.10 and 1.51 times (grid evaluation against `interpn`: 7 left and 6 entered; `Linear1D`'s construction at 10 000 points against `ndarray-interp`: 2 left and 1 entered; and so on). For `Linear1D`'s construction at 10 000 points interpx took 19.5 to 22.9 µs in both runs, `ndarray-interp`'s medians of the rounds moved between 14.8 and 22.3 µs, and its fastest repeat stayed at 13.2 to 13.5 µs. interpx's and SciPy's times are the same between the two runs (median per-case ratios 1.006 and 1.008).
    * **The difference from the run before that (M4-1's `959c7c5`, with the cells from the number of points)**: losses 209 → 216, with 11 out and 18 in. The 3 constructions of `LinearND` against `scirs2-interpolate` are there because it is compared now. The rest are grid and one-variable comparisons whose ratio moved between 1.10 and 1.85 times (grid evaluations against `interpn` 4 out and 7 in, constructions of `Linear1D` against `ndarray-interp` 1 out and 3 in, and so on).
    * **The difference from an earlier run (M4-1's `c519bf1`, with the barycentric coordinates from the inverse)**: losses 216 → 209, with 14 out and 7 in. The losses of `LinearND` stay 16 (construction against `spade`). The rest are grid and one-variable comparisons whose ratio moved between 1.10 and 1.53 times (grid evaluations against `interpn` 4 out and 2 in, constructions of `Linear1D` against `ndarray-interp` 3 out and 1 in, and so on).
    * **The difference from an earlier run (P-18, `79b6ced`)**: losses 197 → 216. 16 of the 19 more are the construction of `LinearND` (new cases). Besides, 6 went out and 10 came in. All are comparisons whose ratio moved between 1.05 and 1.39 times, most by the median of one side moving from round to round (`ndarray-interp`'s construction of `Linear1D` at 10 000 points went 15.3 → 22.1 µs with its fastest repeat the same at 13.48 µs; `scirs2-interpolate`'s construction of the grid's 3-D `linear` at `max200` went 1542 → 1331 µs, interpx's 1630 → 1628 µs). Where interpx's own fastest moved too is the construction of `CubicSpline` at 10 000 ordered points, 150 → 168 µs. The code of `CubicSpline` did not change between the two runs. It is read as one of the constructions at 10 000 points that move by round (P-19, #181), and is looked into in that row.
    * **The difference from the run before that (P-17, `c5fe29c`)** (P-18, #172; the warm-up before construction is timed was put in): losses 197 → 197, with 20 out and 20 in.
      * **12 of the 20 out are timings without a warm-up that measured the slow time.** All 8 constructions in one variable slower than SciPy are out: `Linear1D` at 10 000 points, 4 cases (interpx 93.4 → 19.6 µs, SciPy 50.5 → 44.0 µs), `Nearest1D` at 10 000 points, 2 cases (61.1 → 11.9 µs, SciPy 53.5 → 53.0 µs), `BSpline1D` of degree 1 at 10 000 points, 2 cases (66.0 → 10.8 µs, SciPy staying at 23.1 µs). So are the 2 cases of `CubicSpline` at 1000 points against `ndarray-interp` and `scirs2-interpolate` (38.6 → 16.9 µs), the 6.60 times of `Linear1D` at 10 000 points against `ndarray-interp`, and the 1.96 times of the grid's 3-D `linear` at `max40` against SciPy.
      * **10 of the 20 in are the other library's time going down.** The other libraries' runners measured the slow time too: `scirs2-interpolate`'s constructions at 10 000 points went `Akima1D` 142.6 → 48.4 µs, `Pchip` 173.9 → 42.7 µs, `CubicSpline` 334.8 → 145.5 µs. So `Akima1D` and `Pchip` at 10 000 points came in at 2.63 to 2.90 times `scirs2-interpolate` (at 1000 points they were already counted, at 1.7 to 1.9 times; the reason is §3's "Making coefficients"). The 3 cases of `CubicSpline` at 10 000 points (1.20 to 1.28 times) are against a library whose values differ from SciPy's (`scirs2-interpolate` has another end condition), and are not counted as losses. The 3 grid constructions against `scirs2-interpolate` (`max40` and `max200`, 1.94 to 2.27 times) are also because the other went down, from 2.7 to 15 µs to 1.4 to 3.3 µs.
      * **The rest are comparisons that sat around the threshold (1.10 times)** (grid evaluations against `interpn` 6 in and 4 out, 1.11 to 1.16 times, and so on).
      * Between the two runs SciPy's times hardly moved (the median of the per-case ratios: construction 0.982, evaluation 0.984). interpx's: construction 0.979 (a tenth of the cases at 0.63 times or less), evaluation 1.001.
      * In the run with the warm-up at 250 ms first (`7c7bbbd`), `Linear1D` at 10 000 points did not get past the slow time (medians of the rounds 44 to 63 µs, fastest 21 to 24 µs), and it was made 1 second (§8.3's "A warm-up before construction is timed").
    * **The difference from an earlier run (`056f117`)** (in P-17's run): losses 197 → 197, with 10 in and 10 out. 19 are comparisons that sat around the threshold (1.10 times) (1.10 to 1.24 times; the grid's `linear` evaluations against `interpn` 5 in and 3 out, and so on. Of those out, the 1.24 times of the elevation's `grid/nearest` construction against `scirs2-interpolate` is because the median of `scirs2-interpolate`'s one round moved from 12.4 to 16.2 ms, while its fastest repeat is the same, 12.0 and 12.4 ms). The remaining one (the construction of `CubicSpline` at 1000 ordered points, 2.22 times `ndarray-interp`, 38.6 µs and 17.4 µs) was 2.1 times in the median in the previous run too (37.5 µs and 17.5 µs). This construction of interpx takes either about 15 µs or about 36 to 38 µs depending on the process, and in the previous run one round's fastest repeat was 15.0 µs, below `ndarray-interp`'s fastest, so it was not counted as slower (below).
    * **Some constructions of interpx differed by more than 2 times from process to process** (P-18, #172). Over the three runs (`056f117`, `eff8242`, `c5fe29c`), 32 of the 287 constructions of interpx had a largest median of a round more than twice the fastest repeat of some round (`CubicSpline` at 1000 points, 15 µs and 38 µs; the grid's 3-D `linear` at `max40`, 2.7 µs and 15 µs; and others). The cause was in the measurement (§8.3's "A warm-up before construction is timed"). In the run with the warm-up (`79b6ced`), counted the same way within one run, 1 of 291 is above 2 times (6 in `c5fe29c`). That one is the grid's 3-D `linear` at `max40` with 100 queries, 9.3 µs in one round only (2.9 µs and 3.1 µs in the other two).
      * **After the warm-up too, some constructions at 10 000 points take 1.3 to 1.6 times as long in some rounds** (P-19). For `Akima1D` and `Pchip` at 10 000 points, there are rounds with a median of 83 to 91 µs and rounds with 100 to 138 µs, and the fastest repeat of every round is 81 to 91 µs (`scirs2-interpolate` on the same cases is 42 to 48 µs and does not do this). Measured in the runner in windows of 50 ms for 3 seconds, `Akima1D` at 10 000 points is about 90 µs, with one stretch of a few hundred ms at 1.3 to 1.6 times. When the timed loop of 50 ms falls in that stretch, the median of the round goes up. Why it goes up was not looked into.
    * **The difference from an earlier run (#167, `a23222e`)**: losses 205 → 197, with 12 out and 4 in. All were comparisons that sat around the threshold (1.10 times) (evaluations against `interpn` 5 out and 1 in, constructions against `scirs2-interpolate` 5 out and 2 in, and so on). Between the two runs SciPy's times did not change (the median of the per-case ratios: construction 0.996, evaluation 1.002; 80% of the cases within 0.93 to 1.03). interpx's are the same (0.997 and 0.999). **Pinning SciPy's process to the core again did not move its times**: SciPy's times so far, unpinned as they may have been, were the same as numbers.
    * **SciPy's memory is right from this run on** (P-16, #166). In `docs/performance.md` until then, SciPy was 4.3 MB in every case. What was measured was the `python.exe` of the `uv` environment (a launcher process of 4 MB), and the actual interpreter was its child. An error from the first run of the comparison (#120), found when SciPy read as 4 MB in the cases on real data (the queries alone are 400 MB). Now the SciPy runner says its process number, and the driver measures that process's memory, pins it to the core, and reads back that it is pinned (a run not pinned is a failure, not a number). There is a script that checks the driver itself (`compare/perf/check_driver.py`: it starts a process holding 300 MB the way SciPy is started, and checks that the driver reads 300 MB and that the process is on the core set). The corrected numbers (MB; when the data was read → highest):

      | Case (the run of `056f117`) | SciPy | interpx |
      | :--- | ---: | ---: |
      | Cubic spline in one variable, 10 000 points | 71 → 73 | 3.4 → 5.9 |
      | `linear` in 2 dimensions, 800 × 800 | 76 → 77 | 14.5 → 14.9 |
      | `cubic` in 2 dimensions, 400 × 400 | 72 → 161 | 4.2 → 8.4 |
      | `linear` in 4 dimensions (values of 200 MB) | 271 → 273 | 235 → 405 |
      | `NearestND`, 100 000 points | 77 → 92 | 14.0 → 18.2 |
      | RBF over all 3000 points | 71 → 144 | 3.3 → 80 |
      | RBF over 50 neighbors, 100 000 points | 73 → 99 | 6.7 → 11.9 |
      | The elevation's `linear` (8.65 million grid points, 25.91 million queries) | 555 → 1618 | 497 → 972 |
      | The elevation's `cubic` | 555 → 4087 | 510 → 1042 |
      | The weather's RBF `linear` over all 10 838 points | 71 → 1012 | 5 → 949 |

      * SciPy's process is about 70 MB doing nothing (the interpreter, NumPy, and SciPy). The difference in small cases is mostly that.
      * **On a `linear` grid with large values, interpx uses more** (405 MB and 273 MB): interpx keeps a copy of the values, and SciPy refers to the array passed. It is the cost of the decision of §3 (no borrowing form); `from_owned` does not copy. The tables until then showed SciPy at 4.3 MB, and this comparison could not be read.
      * **In the construction of splines, SciPy uses far more**: 4087 MB and 1042 MB for the elevation's `cubic`, 161 MB and 8.4 MB at 400 × 400. That is making and solving the sparse matrix over all the axes (measured per way of solving in P-13: solving axis by axis within the same SciPy adds 264 MB in the construction of the elevation's `cubic`, the default iterative solver 3362 MB. §6.4).
      * RBF over all points peaks at one matrix for both (about 940 MB at 10 838 points). `scirs2-interpolate` takes 3764 MB on the same case.
    * **The "runs that failed" got reasons** (P-15, #165). The 8 until then (`grid/linear/4d/max10`: 10 × 5 × 2 × 1 points; `scirs2-interpolate` and `interpn`) showed only "exit code 101". Run alone with a backtrace, **all 8 were panics of this repository's runners**: both runners read `axis[1]` to see whether an axis is ascending, and fall over on an axis of one point. The libraries had not even been called (for `scirs2-interpolate` the panic's location was inside `ndarray`, so it was first read as a panic inside the library; the backtrace had no function of the library). Fixed and passed to the libraries, both refuse: `scirs2-interpolate` with "Dimension 3 has less than 2 points", `interpn` with "All grids must have at least 2 entries". Axes of one point are accepted by SciPy and interpx, and not by these two. This was found by measuring only now. The Rust runners now tell the driver about panics and the libraries' refusals with the reason.
    * **No run was stopped at the time limit** (P-17, #168). In the previous run, `scirs2-interpolate` on the weather's RBF over all points (`thin_plate_spline`) was stopped. In the run before it, it finished (construction 273 s), and run alone it finished in 760 s (construction 252 s). The harness built it once for the values, once for the time, and once more for the time of the evaluation, and whether it finished or stopped was decided by the state of the machine. Now that a long construction is built once, the same run finishes in 310 s (construction 308 s; 285 s in this PR's first run). It is the longest run, at 34% of the limit. No run exceeded half the limit (450 s).
    * Going through all cases once takes, in the sum of the processes' times, 23 minutes for the 277 cases of synthetic data and 42 minutes for the 14 cases of real data (about 100 minutes before. Most of it is solving the RBF over all points on one thread: one construction takes SciPy 106 to 117 seconds, `scirs2-interpolate` 277 to 308 seconds, interpx 16 to 18 seconds).
  * **The results on real data** (M5-11, #160. The data and how it is divided are in §8.4.7. The times and errors in the tables below are from the first run measured (#167, `a23222e`). In the runs after the driver was fixed, too, SciPy's and interpx's times are the same within 0.93 to 1.06 times per case, and the errors the same. The current values are in `docs/performance.md`. Only one, `scirs2-interpolate`'s RBF over all points, was stopped at the time limit in a later run: P-17 above):
    * **Elevation** (a grid of 3601 × 2401 points, 25.91 million held-out queries):

      | Method | Library | Construction | Evaluation | RMS error at the held-out points |
      | :--- | :--- | ---: | ---: | ---: |
      | `nearest` | interpx | 16.3 ms | 346 ms | 0.178 m |
      | | SciPy | 7.0 µs | 1.32 s | 0.178 m |
      | | `interpn` | 0 | 716 ms | 0.178 m |
      | | `scirs2-interpolate` | 13.1 ms | 150 s | 0.178 m |
      | `linear` | interpx | 16.6 ms | 370 ms | 0.0601 m |
      | | SciPy | 6.9 µs | 1.07 s | 0.0601 m |
      | | `ndarray-interp` | 8.3 µs | 373 ms | 0.0601 m |
      | | `interpn` | 0 | 841 ms | 0.0601 m |
      | | `scirs2-interpolate` | 13.3 ms | 75.8 s | 0.0601 m |
      | `cubic` | interpx | 76.7 ms | 725 ms | 0.0507 m |
      | | SciPy | 14.6 s | 33.0 s | 0.0507 m |
      | | `interpn` (another interpolant) | 0 | 1.64 s | 0.0517 m |
      | `quintic` | interpx | 101 ms | 1.30 s | 0.0517 m |
      | | SciPy | 21.3 s | 37.8 s | 0.0517 m |
      | `pchip` | interpx | 79.2 ms | 1.13 s | 0.0490 m |
      | | SciPy | does not run | — | — |

      * In evaluation interpx is the fastest, or level with the fastest, for every method (`linear` the same as `ndarray-interp`). The parallel batch (8 threads) is 5 to 6 times faster than that (725 → 140 ms for `cubic`).
      * The construction of `nearest` and `linear` is slower in interpx (16 ms; SciPy and `ndarray-interp` a few µs, `interpn` 0). The others only borrow or refer to the data; interpx copies it and checks that it is finite. The same kind of loss as on synthetic data, with its reason in §3 (`from_owned` takes 5.2 ms, the rest being the checks).
      * `cubic` and `quintic` are 1 to 2 orders faster than SciPy in both construction and evaluation (77 ms against 14.6 s for construction). Most of this gap is solving all the axes together or axis by axis (measured in P-13: solving axis by axis inside SciPy, the construction of the elevation's `cubic` takes 569 ms, 1/24 of the default iterative solver's 13.5 s; the remaining 7.4 times is the implementation. Changing only the solver to the direct one ends for lack of memory. §6.4).
      * **SciPy's `pchip` does not run at this size**: it tries to allocate an array of 695 GiB (queries × the length of an axis) and gets `MemoryError`. SciPy's `pchip` (and `*_legacy`) builds a one-dimensional interpolant along the axis at every evaluation. interpx finds it from the values around the point alone, and finishes in 1.13 s.
      * The errors are the same for every library with the same method (they compute the same function). By method, `pchip` is the smallest, then `cubic`, and `quintic` is no better than `cubic`. Terrain is not a smooth function, so raising the degree does not help.
      * **What `interpn`'s `cubic` is** (read in the source of 0.11.2, `src/multicubic/rectilinear.rs`): a local cubic Hermite interpolation that solves no linear system. The slope at each grid point is decided by a second-order centered difference from the two neighbors (for uneven grids, the left and right secants weighted by the width of the interval on the opposite side; Veldman and Rinzema 1992), and a cubic is made per interval from the values and slopes at its two ends. On an even grid the slopes are those of Catmull–Rom. The end interval is a quadratic, taking "the third derivative is 0 at the end". It is continuous up to the first derivative and keeps no coefficients (which is why its construction is 0). It is a different function from the splines of SciPy and interpx (continuous up to the second derivative, coefficients solved from all the data), and its values differ by $5 \times 10^{-3}$ of the largest value. Adding this method as a choice in interpx was made a row after v0.1.0 (A-13, #163. The owner's proposal).
    * **Weather** (10 838 stations, 1 204 held out. Longitude and latitude in degrees, as they are, are the coordinates):

      | Case | Library | Construction | Evaluation | RMS error |
      | :--- | :--- | ---: | ---: | ---: |
      | Nearest | interpx | 0.68 ms | 0.16 ms | 3.42 °C |
      | | SciPy | 1.83 ms | 0.52 ms | 3.42 °C |
      | RBF `linear`, all points | interpx | 18.3 s | 19.0 ms | 2.74 °C |
      | | SciPy | 114 s | 55.7 ms | 2.74 °C |
      | | `scirs2-interpolate` | 314 s | 258 ms | 2.74 °C |
      | RBF `linear`, 50 neighbors | interpx | 0.83 ms | 19.2 ms | 2.75 °C |
      | | SciPy | 1.89 ms | 44.5 ms | 2.75 °C |

      * The construction of RBF over all points is about 6 times faster in interpx than in SciPy, and about 17 times faster than `scirs2-interpolate`. One construction of RBF over all points takes more than 2 seconds for every kernel, so it is 1 round (one measurement).
      * **The cubic and quintic kernels break down as interpolation on this data, in every library** (P-14, #164). The RMS error at the held-out points (the values of the first run):

        | Kernel | System | interpx | SciPy | `scirs2-interpolate` | Difference between interpx's and SciPy's values |
        | :--- | :--- | ---: | ---: | ---: | ---: |
        | `linear` | all points | 2.737 | 2.737 | 2.737 | $2 \times 10^{-9}$ |
        | `thin_plate_spline` | all points | 5.6696 | 5.6694 | 5.409 | $5 \times 10^{-4}$ |
        | `cubic` | all points | 27.05 | 24.36 | 8.48 | 0.31 |
        | `cubic` | 50 neighbors | 25.34 | 25.34 | — | $2 \times 10^{-7}$ |
        | `quintic` | all points | 13 586 | 4 916 | $2.8 \times 10^{15}$ | 3.7 |
        | `quintic` | 50 neighbors | 198.8 | 199.6 | — | 0.20 |

        (°C. The difference of values as a ratio to the largest value. The nearest station alone, as an interpolant, has an RMS error of 3.42 °C.) No error and no warning.
      * **The cause is stations very close together** (measured in P-14). The distance to the nearest other station is at least 0.0003 degrees, under 0.002 degrees at 1% of the stations, and 0.22 degrees in the median. The reciprocal condition number of the system over all points (interpx's `reciprocal_condition()`, of the balanced system. §6.3), and how far the interpolant misses its own data points (largest; without smoothing it should be 0):

        | Kernel | Reciprocal condition number | Missing the data points by (interpx) | The same (SciPy) | RMS at the held-out points (interpx) |
        | :--- | ---: | ---: | ---: | ---: |
        | `linear` | $6.9 \times 10^{-11}$ | $1 \times 10^{-7}$ | $2 \times 10^{-8}$ | 2.74 |
        | `thin_plate_spline` | $5.0 \times 10^{-17}$ | 0.06 | 0.007 | 5.67 |
        | `cubic` | $6.8 \times 10^{-21}$ | 33 | 11 | 27.1 |
        | `quintic` | $3.6 \times 10^{-24}$ | 56 000 | 5 200 | 13 600 |

        As far as double precision can solve it is $\varepsilon = 2.2 \times 10^{-16}$; `thin_plate_spline` falls slightly below, and `cubic` and `quintic` far below. On systems not solved, which of interpx and SciPy is nearer the data is a difference in the rounding of the solve, not a meaningful one (§6.3). The higher the power of the distance in the kernel, the more alike the rows of two close points, and the worse the condition. Checked by changing the data (with SciPy; RMS at the held-out points):

        | Kernel | As it is | Coordinates in the unit square | Stations closer than 0.1 degrees merged (9512 stations) | Both |
        | :--- | ---: | ---: | ---: | ---: |
        | `linear` | 2.74 | 2.76 | 2.75 | 2.77 |
        | `thin_plate_spline` | 5.67 | 4.35 | 3.09 | 3.19 |
        | `cubic` | 24.4 | 16.4 | 3.79 | 4.21 |
        | `quintic` | 4 922 | 10 350 | 11 940 | 1 670 |

        Merging close stations makes `cubic` and `thin_plate_spline` usable. Shrinking the coordinates alone is not enough. `quintic`, even merged and shrunk, has a reciprocal condition number (LAPACK's value of the unbalanced system) of $3 \times 10^{-24}$, and with this number and arrangement of points it cannot be used in double precision. `linear` is solved on every data set.
      * **What interpx does**: it does not refuse the system, and returns the estimate of the reciprocal condition number (`reciprocal_condition()`. §6.3). The user knows the system is not solved when it falls below $\varepsilon$. Synthetic data and the goldens did not show this: the points are scattered in a square of side 2 or less, and the goldens have at most 64 points.
    * **What real data showed, and what it did not**:
      * The order of speed did not change from synthetic data. No new loss in evaluation, and the losses in construction are of the same kind (copies and checks).
      * What synthetic data did not have: SciPy's `pchip` does not run on large queries, the RBF kernels of high powers break down on real coordinates (P-14), SciPy's memory numbers were wrong (P-16), and how to measure runs that take minutes (§8.4.7).
    * **A number from one round did not agree with its neighbors** (trial run): interpx and the `from_owned` variant, whose evaluation code is the same, read 1.06 s and 0.35 s for the evaluation of `nearest`, and 0.71 s and 1.63 s for `cubic`. In the run of 3 rounds, both are in the same range from round to round (`nearest`: 0.383, 0.345, 0.346 and 0.345, 0.347, 0.390; `cubic`: 0.779, 0.709, 0.725 and 0.804, 0.710, 0.748). Those two values of the trial run did not reproduce. Why they were read is not known (the trial run ran just after writing the files of 5 cases).
  * **The first results** (2026-10-13; 265 cases, 3 rounds; before the rows of P):
    * The 4 cases of evaluation in one variable (10 000 points in no order) were because `ndarray-interp` assumed even spacing and hit the interval in one step. The same form (guess and check) was put in, limited to evenly spaced grids, and 633 → 191 µs (the other 260 µs. P-2, #123). 12 cases of unevenly spaced data were added. There too interpx is the fastest (609 µs for linear on 10 000 points in no order, the other 1.92 ms).
    * In evaluation interpx is the fastest for one variable (96 of 100 cases), `PPoly`, `BSpline`, the curve, `NearestND` (all), and RBF (52 of 53). Only one case of evaluation is slower than SciPy. Evaluation on the grid is slower than `interpn` in 26 of 79 (1.1 to 1.6 times).
    * Construction loses on the grid (60 of 79) and in one variable (30 of 100) (after P-1: 59 of 79 and 23 of 100. The reasons for what remains are in §3's "The differences that remain, and why"). interpx copies the data and checks every number, and the others only borrow (on a grid of 200 MB, the copy 30 ms and the check of finiteness 10 ms; SciPy 20 µs).
    * The losses were made rows P-1 to P-5 of the roadmap (#122 to #126), to be done before `rayon`, `ndarray`, and `serde` (M5-1 to M5-3).
    * Not counted as losses (in the tables): the construction of the grid's `pchip` being 53 times slower than SciPy is because SciPy finds the slopes at every evaluation (the evaluation of the same case: SciPy 1.26 s, interpx 0.53 ms). `interpn`'s `cubic` is a local cubic through 4 points of the grid that solves no spline coefficients (the values differ). `scirs2-interpolate`'s `CubicSpline::new_not_a_knot` has values that differ from SciPy's not-a-knot ($10^{-2}$ of the largest value). SciPy's grid `cubic` differs from interpx's values by about $10^{-4}$ (SciPy solves the coefficients by an iterative method).
    * There were 8 cases interpx could not build (the 4-dimensional grid 10 × 5 × 2 × 1 of SciPy's benchmark; an axis of one point). P-4 made them accepted.

  Not taken: `peroxide` (its interpolation is part of a large numerical library), `splines` and `enterpolation` (for curves and animation, with nothing corresponding to SciPy's types), `interp` (linear only), `ninterp` (linear and nearest on grids; overlaps `interpn`), `rbf-interp`, `makima_spline`, and `bspline` (not updated since 2019 to 2022). For the nearest scattered point there is nothing to compare with, since this crate uses `kiddo`. `BSpline`, `PPoly`, the curve, and RBF's neighbors mode have nothing to compare with in Rust, and are compared with SciPy only.

### 8.4 Data sets for verification and benchmarks

The data the tests (§8.1) and benchmarks (§8.3) use is divided by purpose into 4 layers. No row (`docs/roadmap.md`) adds data not here silently. An addition is a design change first.

| Layer            | Purpose                          | Reference and scale                                | Run                                  |
| :--------------- | :------------------------------- | :------------------------------------------------- | :----------------------------------- |
| **A Correctness** | numerical agreement with SciPy 1.18 | small and deterministic. Committed in `compare/goldens/` | `cargo test`                     |
| **B Convergence and properties** | difference from analytic solutions, orders of convergence, invariants | functions whose values are known analytically. `proptest` | `cargo test` |
| **C Benchmarks** | speed and memory                 | the same parameters as SciPy's benchmarks + large Rust-only sizes | by hand (`just bench`, `compare/perf`) |
| **D Adversarial** | behavior in degenerate, close, and boundary cases | small and deterministic                        | `cargo test`                         |

#### 8.4.1 Layer A: correctness

**Franke's 2D functions (domain $[0,1]^2$).** The six surfaces used in Franke's (1979) comparative study of interpolating scattered data, originally McLain's (1974) functions. Let $r^2 = (x - \tfrac12)^2 + (y - \tfrac12)^2$.

| ID | Name                 | Definition                                                                                                                                                                   |
| :- | :------------------- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| F1 | Overall              | $0.75\,e^{-\frac{(9x-2)^2 + (9y-2)^2}{4}} + 0.75\,e^{-\frac{(9x+1)^2}{49} - \frac{(9y+1)^2}{10}} + 0.5\,e^{-\frac{(9x-7)^2 + (9y-3)^2}{4}} - 0.2\,e^{-(9x-4)^2 - (9y-7)^2}$ |
| F2 | Cliff                | $\dfrac{\tanh\bigl(9(y - x)\bigr) + 1}{9}$                                                                                                                                   |
| F3 | Saddle               | $\dfrac{1.25 + \cos(5.4\,y)}{6 + 6\,(3x - 1)^2}$                                                                                                                             |
| F4 | Gentle               | $\tfrac13\,e^{-\frac{81}{16} r^2}$                                                                                                                                           |
| F5 | Steep                | $\tfrac13\,e^{-\frac{81}{4} r^2}$                                                                                                                                            |
| F6 | Sphere               | $\dfrac{\sqrt{64 - 81\,r^2}}{9} - \tfrac12$ (the sphere of radius $\tfrac89$ centered at $(\tfrac12, \tfrac12, -\tfrac12)$. On the unit square the expression under the root is always positive) |

* **The record of the check**: F1 to F5 agree in two independent implementations (UQTestFuns, and Burkardt's `test_interp_2d` (MIT)) (the coefficient $\tfrac{81}{16}$ = 5.0625 of F4 and $\tfrac{81}{4}$ = 20.25 of F5, both with the factor $\tfrac13$). For F6, UQTestFuns lacks the $\tfrac19$, but the formula with $\tfrac19$ is the one consistent with the description of a sphere of radius $\tfrac89$, and Burkardt's implementation agrees with it. For the second term of F1 the square is not visible in Burkardt's excerpt, so UQTestFuns's formula is taken. Franke's (1979) original has not been obtained, so the first row that makes goldens checks again against the original or the public code of Renka and Brown, and adds the result to this section.
* **The arrangement of scattered points**: Franke's 100 and 33 points, and Lawson's 25 points (the 25 points are Lawson's arrangement, not Franke's). The coordinates are copied mechanically from Burkardt's `test_interp_2d` (MIT). Some of Franke's 100 points lie outside the unit square, so some queries on the unit square are outside the convex hull.
* **The scope of the goldens**: the goldens of layer A are in `compare/goldens/`, per family: `oned` / `cubic_spline` / `bspline1d` / `calculus` / `curve` / `grid` / `grid_pchip` / `scattered` / `rbf` / `rbf_system` (`cubic_spline` has every case of `CubicSpline`, and `calculus` the derivatives and integrals of `PPoly` and `BSpline`; both are split from `oned` to stay within the limit of one file. `bspline1d` has every case of `make_interp_spline`, split in M1-12 for the same reason. `grid` has the cases that compare six methods on the same grids, and `grid_pchip` the cases chosen for `pchip` (swapped axes, steps and flat values, axes of 2 and 3 points, and so on). `grid` is at 75% of the limit, so it is split in M2-3. The cases with axes of 2 and 3 points, which SciPy's grid does not accept, take as reference the values of SciPy's `PchipInterpolator` applied from the last axis on. `rbf` has the interpolated values and `rbf_system` the parts of RBF (matrix, right-hand side, change of coordinates, monomials, evaluation row), the output of SciPy's internal functions. The points are 12 or fewer (the matrix grows with the square of the number of points). `rbf` is 546 KB, so it is split in M3-2). There are also `datasets` with input only (classic data, arrangements of nodes, adversarial input), and `functions` with the values of analytic functions. The 1D cases have, besides the values, the first and second derivatives, the antiderivative, definite integrals, and the coefficients of the piecewise polynomial (or the knots and control points of the B-spline). `calculus` has the derivatives and antiderivatives themselves (the coefficients of `PPoly`, the knots and control points of `BSpline`), derivatives of order 3 and more, second antiderivatives, definite integrals that cross outside the domain, and definite integrals that cross periods. Behavior deliberately changed from SciPy (the "Differences from SciPy" of §6.1) is not put in the goldens and is checked by analytic solutions and property tests. M1-7 adds to `oned` infinite queries (`numpy.interp`), `nearest` / `nearest-up` with `fill_value="extrapolate"`, and values of 2 columns (`interp1d`). M2-1 adds to `grid` query points exactly halfway for `Nearest`, values of 2 columns, a 3-dimensional grid with only some axes descending, and `Extend` on uneven grids and in 3 dimensions. `Clamp` and `Periodic` have nothing corresponding in SciPy and are checked by property tests. M1-13 adds to `curve` each of the three ways of deciding parameters (with the parameters themselves), closed curves (periodic, degrees 2, 3, 5), points in 4 dimensions, point sequences of exactly $k + 1$ points, and point sequences with extremely uneven spacing. M1-12 adds to `bspline1d` even degrees (2 and 4) and degree 6 and more, data of exactly $k + 1$ points, periodic degrees 1 to 5 and periodic data with no more points than the degree, cubics with nonzero derivatives at the ends, cubics with not-a-knot at one end only (the reference SciPy's `CubicSpline`), values of 2 columns, and unevenly spaced data. The goldens of the banded solver are in `banded` (the solutions of LAPACK's `gbsv`, obtained from `scipy.linalg.lapack.dgbsv`, which `make_interp_spline` calls. `scipy.linalg.solve_banded` calls `gtsv`, dedicated to tridiagonal systems, when the band is one on each side, so it is not the reference. They also have the number of columns where LAPACK interchanged rows, so that even the choice of pivots is checked. Cyclic banded matrices are the solutions of the dense `scipy.linalg.solve`. Random banded matrices, matrices where row interchanges actually happen, ill-conditioned matrices, $n = 1$ to $3$, bands of full width (dense matrices), right-hand sides of several columns). M1-11 adds to `CubicHermiteSpline` 2 points, 3 points, unevenly spaced data, data with all slopes 0, data on a line given the line's slope, slopes that do not fit the data (an example that overshoots), values of 2 columns, and periodic extrapolation (data whose ends agree and data whose ends do not). M1-10 adds to both methods of `Akima1D` 2, 3, and 4 points, data on a line and flat stretches (points where the formula does not decide the slope), and values of 2 columns. M1-9 adds to `PchipInterpolator` 2 points, 3 points, data with flat stretches and changes of sign, monotone data, and values of 2 columns. M1-8 adds to `CubicSpline` 2 points (per boundary condition), 3 points (not-a-knot, periodic, mixed), periodic with 4 points, mixes of not-a-knot and derivatives, values of 2 columns, and data with extremely uneven spacing. The goldens of the tridiagonal solver are in `tridiagonal` (the solutions of SciPy's `solve_banded`: random diagonally dominant matrices, the matrices of the 4 kinds of boundary conditions made from the classic data, cyclic matrices, ill-conditioned matrices, $n = 1$ to $3$, right-hand sides of several columns). The reference of the grid's `cubic` / `quintic` is the values of `cubic_legacy` / `quintic_legacy`, as in §6.4. M2-2 adds to `grid` the coefficients themselves (`make_interp_spline` applied per axis), axes with exactly the minimum number of points, descending axes, and values of 2 columns.
* **Size limits**: the goldens are limited to 1 MB per file and 8 MB for the directory in all (both uncompressed). The goldens are included in the published package (so that `cargo test` passes on the crate as downloaded). The crates.io limit of 10 MB is the compressed size, and the goldens compress to about 1/5, so 8 MB is about 1.6 MB compressed. The limit per file is for diffs and review, and a family that exceeds it is split into files. A row that would bring the total over 90% of the limit decides the limit again first, without dropping cases. The history, and the alternatives not taken, are in ADR 0010.

**Classic 1D data.** The values are copied from primary sources (original papers and texts), and the source is recorded in the golden's `meta`.

| Data                     | Points | What it looks at                                                         | Source                     |
| :----------------------- | :----- | :----------------------------------------------------------------------- | :------------------------- |
| Fritsch–Carlson          | 11     | sharp turns and sparse points. Overshoot, and PCHIP's monotonicity       | Fritsch and Carlson (1980) |
| Titanium properties (titanium) | 49 | unevenly spaced nodes. The standard data of spline approximation        | de Boor and Rice           |
| $y = x^2$ + one close point | 12  | sensitivity to a point whose $y$ does not fit at a close $x$             | Burkardt `test_interp` p08 |
| Runge's function         | —      | the oscillation at the ends of $f(x) = \dfrac{1}{1 + 25x^2}$ ($x \in [-1, 1]$, evenly spaced nodes) | Runge |
| Step                     | —      | 0 for $x < 0.5$, 1 otherwise. Monotonicity, the overshoot of `CubicSpline` | —                        |
| Chirp                    | —      | $\sin(x^2)$ ($x \in [0, 10]$). Distortion of phase and high frequencies  | —                          |

#### 8.4.2 Layer B: convergence and properties

* **Orders of convergence**: on smooth functions the error is measured with $h$ halved, and $p = \log_2 (e_h / e_{h/2})$ found. The expected values (smooth data) are as follows.

  | Interpolant                                 | Expected order |
  | :------------------------------------------ | :------------- |
  | `Linear1D`, grid `linear`                   | 2              |
  | `CubicSpline` (NotAKnot), grid `cubic`      | 4              |
  | `PchipInterpolator`                         | 3              |
  | `BSpline1D` (degree $k$)                    | $k + 1$        |
  | grid `quintic`                              | 6              |

  Of the interpolants not above, `LinearND` (M4-1) and `CloughTocher2D` (M4-2) are in the table below.
* **The orders measured, and the types not in the table** (M5-4, #157; `tests/convergence.rs`). The function is $\sin(1.3x + 0.4) + 0.5\,e^{0.7x}$ (on $[0, 3]$; not symmetric, with an extremum). The points are both evenly spaced and not (the widest interval about 3 times the narrowest). The orders read over the last 2 halvings:

  | Interpolant | Order | Measured | Grounds |
  | :--- | :--- | :--- | :--- |
  | `Linear1D` | 2 | 1.98–1.99 | as in the table |
  | `Nearest1D` | 1 | 0.97–1.02 | the error of the value is slope × half an interval |
  | `CubicSpline` (not-a-knot) | 4 | 3.96–3.99 | as in the table |
  | the same, first derivative / second derivative | 3 / 2 | 2.96–2.99 / 1.96–2.0 | each differentiation lowers it by 1 |
  | the same, definite integral over the whole interval | 4 | 3.97–4.04 | the same as the values |
  | `CubicSpline` (natural) | 2 | 2.00 | it imposes "second derivative = 0" at the ends, which the function does not have, so it drops by 2 at the ends (de Boor, *A Practical Guide to Splines*, chapter IV) |
  | `PchipInterpolator` | 3 | 2.97–2.99 | as in the table (Fritsch and Carlson 1980; the estimate of the slopes is of order 2) |
  | `Akima1D` (Akima / makima) | 3 | 2.82–3.14 / 2.91–2.96 | the slope is a weighted mean of neighboring secants, of order 2, so the values are of order 3 (Akima 1970). The weights change from point to point, so Akima's order wavers from halving to halving |
  | `CubicHermiteSpline` (given the function's slopes) | 4 | 4.00 | cubic Hermite interpolation per interval (error $h^4 \max\lvert f^{(4)}\rvert / 384$) |
  | `BSpline1D` (degree $k = 1 \dots 5$) | $k + 1$ | within 0.1 of $k + 1$ for $k \le 4$; 5.72–5.96 for $k = 5$ | as in the table. Degree 5 is measured down to the spacing just before the error approaches the rounding of the solution ($10^{-11}$) |
  | Grid `nearest` / `linear` / `cubic` / `pchip` (2-D) | 1 / 2 / 4 / 3 | 0.94–1.09 / 1.87–1.96 / 3.60–3.93 / 2.85–3.03 | tensor products of the one-dimensional methods. In 2 dimensions the spacing cannot be made as fine, and it is somewhat far from the limit |
  | Grid `quintic` (2-D) | 6 | 5.42 → 5.77 → 5.87 (even), 4.55 → 5.05 → 5.70 (uneven) | still on its way up to the order. The test looks at the last order and at the rise |
  | Grid `linear` / `cubic` (3-D, Genz's oscillatory) | 2 / 4 | within widths 0.2 / 0.45 | the function of §8.4.4 |
  | `LinearNDInterpolator` (2-D, Franke's first function) | 2 | 1.55 → 2.01 → 2.01 (a lattice with its inner points moved by up to 0.15 h). On Halton points, 16 times the points bring the error down 12.8 times | piecewise linear interpolation on a triangulation whose angles stay away from zero is off by $C h^2$ (Ciarlet, *The Finite Element Method for Elliptic Problems*, 1978, theorem 3.1.6) |
  | `CloughTocher2DInterpolator` (2-D, Franke's first function, points moved the same way) | set from the measurement (right) | 2.38 on average over $h = 1/10 \to 1/160$ (error $2.06 \times 10^{-2} \to 2.79 \times 10^{-5}$). SciPy iterated to convergence gives the same errors on the same points within 15% | With exact values and gradients at the points the patch reproduces quadratics and is of order 3 (the reduced Hsieh–Clough–Tocher element, whose derivative across an edge is linear along it; Ciarlet, "Interpolation error estimates for the reduced Hsieh-Clough-Tocher triangle", *Math. Comp.* 32, 1978; the patch's reproduction of a quadratic is a test inside the module). The order with the estimated gradients has not been checked against a source (the abstract of Renka and Cline 1984 speaks of accuracy on examples, and no order was found; the paper itself has not been read). With the gradients replaced by the exact ones of the function, SciPy goes 3.25 → 3.66 → 3.69 → 3.05 ($h = 1/20 \to 1/160$): the estimate of the gradients is what sets the order. The test checks that the average is in 2.2 to 2.6 |

  * The widths of the tests (the tolerance from the expected order) were decided from the values measured and fixed in the tests. They are narrower than the distance to the next order (1).
  * **The points where the error is measured are at least 3 in the narrowest cell.** On a 3-dimensional grid, measuring with a lattice of 1 point or fewer per cell let the largest error (in a cell at the not-a-knot end) fall between the points measured, and the order read as 5.9 or 7.0 (correctly 4).
  * **RBF has no single order.** Scattered points have no spacing to halve, and the order depends on the kernel (and, for kernels with a shape parameter, on that too). What the tests look at is that, interpolating Franke's F1 on Halton points, 8 times the points (100 → 800) lowers the error by at least a decided factor and keeps it within a decided bound. The factors measured are from 12 times for `linear` to 20 000 times for `multiquadric` solving over all points, and 9 to 350 times over 30 neighbors. For `NearestNDInterpolator`, 16 times the points at least halves the error.
  * **The Gaussian gets worse past a point as points are added.** The error of the Gaussian over all points ($\varepsilon = 5$), with the points in the order given, was $2.9 \times 10^{-4}$ / $9.4 \times 10^{-5}$ / $2.0 \times 10^{-4}$ at 400 / 800 / 1600 points. SciPy too gets worse at 1600 points on the same data: $1.0 \times 10^{-4}$ / $1.8 \times 10^{-5}$ / $4.6 \times 10^{-5}$. The system is ill-conditioned (coefficients of $10^{9}$ to $10^{11}$). The error looks decided by rounding, but the condition number was not measured, and the mechanism has not been confirmed ("What is known" below). That adding points makes it worse happens in SciPy just the same.
    * **interpx's error with the points in the order given looked 3 to 5 times SciPy's, but that was one point within a range decided by the order of the points** (P-12, #159). Given the same points in another order, both errors move by 1 to 2 orders. The range and median of the largest error, with 12 orders given to interpx and 24 to SciPy:

      | Points | interpx's range | interpx's median | SciPy's range | SciPy's median |
      | ---: | :--- | ---: | :--- | ---: |
      | 400 | $4.3 \times 10^{-5}$ – $1.4 \times 10^{-3}$ | $1.3 \times 10^{-4}$ | $4.7 \times 10^{-5}$ – $4.9 \times 10^{-4}$ | $1.1 \times 10^{-4}$ |
      | 800 | $7.3 \times 10^{-6}$ – $9.4 \times 10^{-5}$ | $3.4 \times 10^{-5}$ | $7.0 \times 10^{-6}$ – $6.2 \times 10^{-4}$ | $2.4 \times 10^{-5}$ |
      | 1600 | $1.5 \times 10^{-5}$ – $2.3 \times 10^{-4}$ | $3.4 \times 10^{-5}$ | $2.2 \times 10^{-5}$ – $3.7 \times 10^{-4}$ | $5.6 \times 10^{-5}$ |

      The ranges overlap. The ratios of the medians (interpx / SciPy) are 1.2, 1.4, 0.6, within 1.7 times either way (at 1600 points SciPy's is 1.65 times larger). The numbers of orders differ (12 and 24), so the ends of the ranges cannot simply be compared: the more orders, the wider the range comes out (the $9.4 \times 10^{-5}$ at 800 points was the largest of interpx's 12 orders. That it first read as "outside the range" against only 6 orders of SciPy's was because 6 orders of SciPy's were too few; with 24 it goes up to $6.2 \times 10^{-4}$). The difference of the solves (interpx's LU, SciPy's symmetric factorization. §6.3) is not visible within this range.
    * The interpolant then also misses the data points by rounding (the residual at the data points is $7 \times 10^{-5}$ at 800 points and $2 \times 10^{-4}$ at 1600). `multiquadric` is well conditioned at the same numbers of points, its error does not move with the order of the points ($1.2$ to $1.4 \times 10^{-8}$ at 1600 points), its values are SciPy's, and its residual is at most $2 \times 10^{-9}$.
    * **The numbers of this section are from the Windows development machine.** With the same code and the same points they change by machine: on the Linux machine of the review (#161), the largest error at 800 points in the order given was $3.67 \times 10^{-5}$ ($9.4 \times 10^{-5}$ here; the same value in debug and release builds).
    * **What is known** (not a proof of the cause): in the same review, moving the values by random amounts of relative $10^{-16}$ to $10^{-12}$ left the largest error at $3$ to $5 \times 10^{-5}$. The error does not look like an amplification of the data's perturbation, but decided on the side of the matrix (the rounding of the factorization). That the values change by machine (rounding changes when the paths of fused multiply-add or SIMD differ) does not contradict this either.
    * **A larger $\varepsilon$ is no way out.** The condition improves, but the bumps narrow and do not reach between the points, and the error grows for another reason. The largest error on the same 800 points at $\varepsilon = 2$ / 5 / 10 / 20 / 40 was $2.8 \times 10^{-1}$ / $9.4 \times 10^{-5}$ / $2.6 \times 10^{-5}$ / $8.5 \times 10^{-3}$ / $3.9 \times 10^{-1}$ on this development machine, and $2.1 \times 10^{-1}$ / $3.7 \times 10^{-5}$ / $2.6 \times 10^{-5}$ / $8.5 \times 10^{-3}$ / $3.9 \times 10^{-1}$ on the review's Linux machine. From 10 up the two machines agree (well conditioned, not decided by rounding), and at 5 and below they differ (ill-conditioned). For this problem the bottom is around 10; 20 is about 300 times worse and 40 about 15 000 times. The first version of the rustdoc said "avoided with a large `epsilon`", which was wrong. What is safe is neighbors (the Gaussian over 30 neighbors gives $5.4 \times 10^{-5}$ at 800 points), or another kernel.
    * Ill-conditioned systems are not made goldens: their values move by 1 to 2 orders with the order of the points and rounding, so they cannot be compared with SciPy's values within a decided width. `tests/convergence.rs` only checks the Gaussian's error over all points to be at most $10^{-3}$ at 800 points (about 10 times the largest measured). The rustdoc of `RbfKernel::Gaussian` states this property.
  * **Changes that do not change the order are not found by this test.** Swapping the weights of PCHIP's inner slopes or changing makima's weights leaves the order at 3 (both are found by the tests that compare with SciPy's values). What the tests of order find are errors that pass through the data but lower the order, such as making the end slope a secant, or not extending the secants beyond the ends.
* **Analytic solutions**: derivatives, indefinite integrals, and definite integrals are measured against functions whose results are known analytically: $\sin$, polynomials, exponentials, and so on.
* **Properties**: the `proptest` items of §8.1.
* Functions: in 1D the classic functions above, $\sin$, and exponentials. In 2D Franke's F1 to F6.

#### 8.4.3 Layer C: benchmarks

The comparison uses **the same parameters as SciPy 1.18's benchmarks (`benchmarks/benchmarks/interpolate.py`)**. The parameters are checked again whenever the SciPy version is raised.

| Subject                                      | Parameters of SciPy's benchmark                                                                                                    |
| :------------------------------------------- | :--------------------------------------------------------------------------------------------------------------------------------- |
| 1D (`Interpolate1d`)                         | $n \in \{10, 50, 100, 1000, 10000\}$, methods linear / nearest / zero / slinear / quadratic / cubic. Exponentially decaying data, $4n$ evaluation points |
| `PPoly` / `BSpline`                          | cubic with 55 intervals. Evaluation at $n \in \{10, 100, 1000, 10000\}$ points                                                    |
| `RegularGridInterpolator`                    | dimensions $\{2, 3, 4\}$, largest points per axis $\{10, 40, 200\}$, evaluation points $\{10, \dots, 10000\}$, axis direction $\pm 1$ |
| Grid cubic / quintic (2D)                    | cubic: $\{10, 40, 100, 200, 400\}$ per axis, quintic: $\{10, 40\}$                                                                 |
| `griddata` (`nearest` / `linear` / `cubic`)  | 1000 random points in 2D, output grid $\{10, 100, 1000\}$ per axis                                                                 |
| `NearestNDInterpolator`                      | $n \in \{10, 100, 1000\}$, output grid $\{50, 100, 500\}$ per axis                                                                 |
| `CloughTocher2D`                             | $n \in \{10, 50, 100\}$, $50 \times 50$ evaluation                                                                                 |
| `RBFInterpolator`                            | `neighbors` `None` or 50, $N \in \{10, 100, 1000\}$, all 8 kernels, `epsilon = 5.0`, random points in $[-1, 1]^2$                  |

* **Large Rust-only sizes**: `NearestND`, `LinearND`, and `CloughTocher2D` with $N$ up to $10^5$ (`CloughTocher2D` from 1000 points, evaluated at the pixels of a 500 × 500 image), grids up to several hundred points per axis. RBF's global mode at $N \le 10^4$ (the memory is $O(N^2)$), and for larger $N$ only with `neighbors`.
* The comparison with the existing pure Rust crates is done in M5-5 (§8.3).
* Construction and evaluation are measured separately. For batch evaluation, `eval_batch_into` and `eval_batch_into_par` are compared with SciPy's vectorized calls.

#### 8.4.4 Functions of dimension $D \ge 3$

* **Genz's (1984) family of functions** ($x \in [0,1]^D$, $u_i = \tfrac12$). The parameters are fixed at $a_i = 2 + i/D$ ($i = 0, \dots, D-1$). This is a choice of implementation, not Genz's normalization of difficulty. The formulas are the standard forms but have not been checked against the original; the first row that uses this family to verify interpolation checks them and records the result.

  | Family                        | Formula                                               |
  | :---------------------------- | :---------------------------------------------------- |
  | Oscillatory                   | $\cos\bigl(2\pi u_1 + \sum_i a_i x_i\bigr)$           |
  | Product peak                  | $\prod_i \bigl(a_i^{-2} + (x_i - u_i)^2\bigr)^{-1}$   |
  | Corner peak                   | $\bigl(1 + \sum_i a_i x_i\bigr)^{-(D+1)}$             |
  | Gaussian peak                 | $\exp\bigl(-\sum_i a_i^2 (x_i - u_i)^2\bigr)$         |
  | Continuous                    | $\exp\bigl(-\sum_i a_i \lvert x_i - u_i \rvert\bigr)$ |

* **Rosenbrock** (a curved valley): $\sum_{i=1}^{D-1} \bigl[100\,(x_{i+1} - x_i^2)^2 + (1 - x_i)^2\bigr]$ is used as one case.
* The dimensions are $D \in \{3, 5, 8\}$ (`kiddo` recommends low dimensions, fewer than about 10).

#### 8.4.5 Sampling of scattered points

All four kinds are deterministic, with a fixed seed.

1. **Halton sequence** (quasi-random): irregular observation points that fill space.
2. **Uniform random**: clusters and sparse areas arise naturally.
3. **Grid with gaps**: a square grid with 30 to 70% of its points removed.
4. **Geometric distributions**: distributions corresponding to Qhull's `rbox` (inside a cube, on a sphere, on the surface of a cube, and so on). Only the kinds of distribution match; the random generator is NumPy's `Generator` (it does not match the generator of `convx`'s `benches/sets.txt`).

#### 8.4.6 Layer D: adversarial input

The expected behavior follows ADR 0005 and 0006.

| Case                   | Data                                          | Expected behavior                                                                          |
| :--------------------- | :-------------------------------------------- | :----------------------------------------------------------------------------------------- |
| Cocircular, cospherical | 4 or more points on a circle (2D), 5 or more on a sphere (3D) | succeeds, with a finite result. `LinearND` reproduces affine functions exactly |
| Collinear, coplanar    | point sets on one line or one plane           | `Err(DegenerateGeometry)`                                                                  |
| Very thin simplices    | flat simplices with one extremely short edge  | barycentric coordinates finite, within the tolerance                                       |
| Exact duplicate points | points with the same coordinates              | `Err(DuplicatePoints)` for `RBFInterpolator`. The others merge and accept them, using the value of the smallest index |
| Very close (not duplicate) | pairs of points about $10^{-12}$ apart    | `Ok` with finite values, or `Err(SolverFailed(Singular))`. No panic, no infinite loop      |
| Queries outside the convex hull | outside the point set                | `None` (`fill` for `eval_or`). No infinite loop                                            |
| NaN queries            | NaN in a query                                | `None` (`fill` for `eval_or`)                                                              |
| Input data not finite  | NaN / Inf in coordinates or values            | `Err(NonFinite)`                                                                           |

#### 8.4.7 Real-world data

* Included in v0.1.0. But not in `cargo test`; used only in the manual comparison of `compare/perf`.
* The fetching script (`compare/perf/real.py`, `just compare-fetch`) is in `compare/perf/`, and only data whose license was confirmed at its source is used. The raw data is not committed (it goes to `compare/perf/out/real/`, a directory that is ignored).
* **History**: at first it was to be chosen in M5-5, but M5-5 was closed without choosing, and this was found in the cross-check of M5-4 (#157). Including it in v0.1.0 was decided anew by the owner on 2026-10-13 (that "the results do not change" from a comparison of synthetic data alone is a guess until measured). The data sets were decided in the Grill of M5-11 (#160).
* **The data is not cut down** (the owner's decision). It is used whole, in the unit its source distributes it in. A size chosen by us would look chosen for convenience. What is used and how it is divided is decided before measuring and written here, and not changed after seeing the results.

  | | Elevation (grid) | Weather (scattered points) |
  | :--- | :--- | :--- |
  | Data | tile `BA32_10000_0401` of Auckland Part 1 LiDAR 1m DEM (2024) | the daily maximum temperature (`TMAX`) of 2025-01-01 from GHCN-Daily |
  | Source | the public bucket of Toitū Te Whenua LINZ, `s3://nz-elevation` (AWS Open Data) | NOAA's public bucket `s3://noaa-ghcn-pds` (AWS Open Data) |
  | License | CC BY 4.0 (the licensor is Regional Software Holdings Limited; confirmed in the survey's `collection.json` and at https://registry.opendata.aws/nz-elevation/) | CC0 1.0 (confirmed at https://registry.opendata.aws/noaa-ghcn/. NOAA asks for attribution and forbids implying NOAA's endorsement) |
  | Pinning the version | the file does not change. The SHA-256 is the value in the source's STAC (`06b2c8bd…0b92`), checked at every fetch | **cannot be pinned** (below) |
  | Size | 4800 × 7200 cells (a 1 m grid, 34.56 million points). No cell without a value. Elevations $-9.6$ to $195.7$ m | 12 185 rows for that day. 12 095 passed the source's quality checks, 12 042 stations with distinct coordinates |
  | Coordinates | the northing and easting of the cell centers (NZTM2000, m). The rows go from north to south, so the northing axis is reordered into ascending order before being passed | longitude and latitude (degrees) as two-dimensional coordinates as they are (not spherical distances. Using the data as it is means choosing no transformation) |
  | Data and held-out points | every second row and column of cells (and the last row and column) make the grid: 3601 × 2401 points. The other 25.91 million cells are held out as the query points. The last row and column go into the grid so that every held-out point is inside the grid | in the order of the file, the last of every 10 stations is held out: 1 204 stations. The other 10 838 stations are the data |
  | What runs | all 5 methods of the grid | `NearestNDInterpolator`, `LinearNDInterpolator` (M4-1), `CloughTocher2DInterpolator` (M4-2), and `RBFInterpolator` with the 4 kernels without a shape parameter (`linear`, `thin_plate_spline`, `cubic`, `quintic`), with one system over all points and with 50 neighbors |

  * **How the tile was chosen**: of the tiles containing the summit of Maungawhau (Mount Eden) (longitude 174.7644, latitude $-36.8773$), the one from the most recently acquired survey (the 2016–2018 survey has a tile of the same sheet too). **The mountain straddles the south edge of the tile** (the summit is about 110 m from the edge). Cutting out the mountain alone would mean deciding the extent ourselves, so one whole tile is used. The mountain's name is not put in the name.
  * **How the day and element were chosen**: January 1 of the latest completed year. The element with the most stations after precipitation (which is 0 at most stations). Of stations with the same coordinates, the first in the file's order is taken (53 rows drop; RBF does not accept duplicate points).
  * **Kernels with a shape parameter are not run**: their value would be ours to choose, and the result would follow the choice (the Gaussian example of §8.4.2).
  * **No upper bound on the number of points is decided beforehand**: RBF over all points has 10 838 points, a little above the guide of §8.4.3 ($10^4$). It is run rather than excluded by a bound beforehand. What is stopped at the time limit (15 minutes for one library on one case) is reported as "did not finish".
  * **Not taken**: R's `volcano` (the elevation of the same mountain; no license is stated for the data). The 2013 layer of the LINZ Data Service (downloading needs an account). Jura (no license statement found, and the source's personal site is gone). The USGS geochemical data release (DOI 10.5066/P1FXE9GK, CC0; the distribution site is behind a bot-protection check and cannot be fetched by a script. The check is not circumvented).
* **GHCN-Daily cannot have its version and checksums pinned.** The files of the public bucket are rewritten every day (the version changes daily, like `3.35-upd-2026101117`), and values of past days can be corrected. No dated archive is distributed. So the version string at the time of fetching and the SHA-256 of the 3 files are recorded in `out/real/weather.json` and put in the report. If a file fetched later differs from the record, the script stops (it does not silently measure on other data). There is no guarantee that the same files can be fetched later. A data set of scattered points with a pinned version was not found in the range searched ("Not taken" above).
* **What is compared**: besides the times of construction and evaluation and the memory (as for the other cases), the error at the held-out points (largest and RMS) and the number of queries that returned no value. The error is decided by the method of interpolation and the nature of the data, not by the merit of a library (with the same method, every library should give the same error). If the error differs, that library computes something else.
* **How long runs are measured**: a construction or evaluation that takes 2 seconds or more is timed by that one call (done 6 times like short calls, one case of RBF over all points took 11 minutes for SciPy alone. 3 measurements were 106 to 124 seconds, and repeating does not change the numbers). That one call is the construction that made the interpolator to write the values (for evaluation, the evaluation that found the values) itself. It is not rebuilt to be timed, and the time of evaluation is measured with that interpolator too (P-17, #168). Before, it was built once for the values, once for the time, and once more for the time of the evaluation. A construction of 250 seconds made 3 times put the run at the edge of the time limit.
  * **Runs near the time limit are named**: the time of the process is recorded per run (one library on one case, one process). The report says what fraction of the limit the longest finished run took, and names the runs that exceeded half the limit (450 seconds). The same run took 760 seconds alone and over 900 seconds inside a run of the whole comparison. The state of the machine moves it by 20% or more, and how far it can move is not known, so half the limit leaves a margin of 2 times (the owner's decision).
  * **A stopped run is run once more**: a run stopped at the time limit is recorded as such, and not started in the next round. After all the runs it is run once more alone, and reported as "stopped" only if it stops both times (the owner's decision). A run that stopped the first time and finished the second is not counted as a failure, and both are put in the table of runs near the time limit.
* **Long runs are measured for one round only**: on a case, a library whose construction or evaluation took 2 seconds or more in one call, and whose call was timed by one call, is measured in the first round only. What decides is whether the runner timed that call by one call, not the median reported (the review of #171. The owner's decision): SciPy's `thin_plate_spline` at 3000 points takes 1.75 to 2.3 seconds a call, straddling 2 seconds; in one round the first call took 1.99 seconds and was timed several times, but those timings were slow at 2.2 to 3.7 seconds, so the rule deciding by the median did not run the next round, and the value reported became the median of two, one normal round and one slow (1.64 seconds read as 2.35). The rounds are repeated to see the spread between processes, and for calls this long the spread is a few percent. Measuring all 3 rounds took 4 hours for the cases on real data alone, most of it one library solving the same system again. The rule is decided by time, not by the name of a library (SciPy's RBF over all points and interpx's alike are 1 round). Slow libraries are not removed from the comparison, and the time limit is not set to where they fall: "how many seconds the same problem takes" is a number this comparison produces.

#### 8.4.8 Measures

1. **Approximation error**: RMSE $= \sqrt{\tfrac1M \sum (f_{\mathrm{true}} - f_{\mathrm{interp}})^2}$, the largest error $\max \lvert f_{\mathrm{true}} - f_{\mathrm{interp}} \rvert$, the order of convergence (§8.4.2).
2. **Agreement with SciPy**: judged by the categories of §8.2 (no single tolerance across the board).
3. **Throughput**: the time of construction and the time of evaluation per evaluation point are reported separately. The APIs are `eval_batch`, `eval_batch_into`, `eval_batch_par`, `eval_batch_into_par`.
4. **Memory**: the peak allocation at construction and at evaluation are reported separately (§8.3).

---

## 9. Operations, Milestones, Roadmap

### 9.1 Operations
* **MSRV**: 1.89 / **edition**: 2021 / **license**: `MIT OR Apache-2.0`
* **CI (GitHub Actions)**: `fmt`, `clippy -D warnings`, tests on Linux / Windows / macOS, an MSRV job, a feature matrix (`--no-default-features` / `--all-features`).

### 9.2 Release policy
**No staged releases.** `v0.1.0` is published at once when all the main features are in (so that publication does not create constraints of backward compatibility early). The following are internal milestones, not published versions. Until `convx` is published, M1 to M3, which do not depend on `convx`, go first.

- [ ] **M1: Core + 1D + Curve**
  - [ ] The `Interpolator` trait, `InterpolateError`, `Extrapolate` / `BoundaryCondition`
  - [ ] `PPoly` / `BSpline` (derivatives, integrals, conversion, the Thomas algorithm)
  - [ ] `Linear1D`, `Nearest1D`, `CubicSpline` (NotAKnot / Natural / Clamped / Periodic / Derivative)
  - [ ] `PchipInterpolator`, `Akima1D` (makima), `CubicHermiteSpline`, `BSpline1D`
  - [ ] `BSplineCurve`
- [ ] **M2: Grid**
  - [ ] `RegularGridInterpolator` (`nearest` / `linear` / `slinear` / `cubic` / `quintic` / `pchip`)
- [ ] **M3: Scattered without `convx`**
  - [ ] `NearestNDInterpolator` (`kiddo`)
  - [ ] `RBFInterpolator` (`faer` LU, 8 kernels, polynomial terms, local neighbors)
- [ ] **M4: Scattered with `convx`** (started after `convx` is published)
  - [ ] `LinearNDInterpolator` (Delaunay + walk from cells, §6.2)
  - [ ] `CloughTocher2DInterpolator` (gradient estimate + Bernstein-Bézier)
- [ ] **M5: Integration and quality**
  - [ ] `rayon` parallel batch evaluation, `serde` (interoperability with `ndarray` after v0.1.0. §7)
  - [ ] SciPy 1.18 golden test suite, `proptest`
  - [ ] Comparison benchmarks (SciPy and existing Rust crates, speed and memory)
  - [ ] docs (recommended dimensions, RBF memory, handling of duplicate points, and so on), CI
  - [ ] **On completion, publish `v0.1.0`**

### 9.3 Roadmap (after v0.1.0)
* `roots` / `solve` (`PPoly` / `BSpline`), `NdPPoly`
* Derivatives on grids and of scattered data
* A wrapper with a dynamic dimension (grid, scattered)
* A pseudo-inverse fallback for RBF
* `BarycentricInterpolator`, `KroghInterpolator`
* The FITPACK family (`UnivariateSpline`, `splrep`, `splprep`, `RectBivariateSpline`, and others)
* Looking into `no_std` (`alloc`)
