# 0005: Scattered geometry: `convx` boundary, duplicates, walk, and serde

- Status: Accepted
- Date: 2026-10-11
- Issue: N/A (grill, pre-Issue)

## Context

`convx` merges duplicate points to the lowest index (`representative()`), reports a boundary as `u32::MAX` (not a public constant), splits cospherical groups into simplices, and does not guarantee `simplices` order across versions. `convx` is unpublished and has no license yet. A point-location walk can cycle in degenerate cases.

## Decision

- `NearestND`, `LinearND`, and `CloughTocher2D` accept duplicate points and use the value of the lowest-index point. `RBFInterpolator` returns `Err(DuplicatePoints)` (a duplicate makes the system singular).
- Point location: jump with `kiddo` to the nearest vertex, then walk. The step limit is `n_simplices`. On reaching the limit or detecting a cycle, scan all simplices with the broad tolerance. If none contains the point, return `None`.
- Tolerances follow SciPy: `eps = 100 * f64::EPSILON`, `eps_broad = f64::EPSILON.sqrt()`. A point on the hull boundary within tolerance is inside.
- `NO_NEIGHBOR = u32::MAX` is defined once in interpx until `convx` exports a constant (a request is made to `convx`).
- `serde` stores vertices, neighbors, values, and parameters with a format version, and validates lengths and index ranges on load. It does not re-triangulate, so saved data does not depend on `convx`'s ordering. It does not re-check the Delaunay property.
- `CloughTocher2D` reproduces SciPy's iterative global gradient estimate (`tol = 1e-6`, `maxiter = 400`, both parameters). Non-convergence is not an error; iteration count and `converged` are readable.
- `convx` must be published with a license before interpx is published; development uses a path dependency with a version.

## Consequences

Saved models survive a `convx` upgrade. Degenerate input terminates. Behavior at the hull boundary matches SciPy.
