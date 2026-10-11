"""Layer-C benchmark cases (docs/design.ja.md section 8.4.3).

The parameters are those of SciPy's own benchmark suite
(`benchmarks/benchmarks/interpolate.py`), so that interpx and SciPy are timed on
the same sizes, methods, and kernels. The runners that use this table are added
by roadmap row M5-5 (`compare/perf/`). `RUST_ONLY` lists the larger sizes that
are timed for interpx alone.

Re-check the table against SciPy's file whenever the SciPy baseline moves.
"""

SCIPY = {
    "interp1d": {
        "n_samples": [10, 50, 100, 1000, 10000],
        "methods": ["linear", "nearest", "zero", "slinear", "quadratic", "cubic"],
        "data": "exponential decay; evaluate at 4 * n_samples points",
    },
    "ppoly_bspline": {
        "intervals": 55,
        "k": 3,
        "n_eval": [10, 100, 1000, 10000],
    },
    "regular_grid": {
        "ndim": [2, 3, 4],
        "max_coord_size": [10, 40, 200],
        "n_samples": [10, 100, 1000, 10000],
        "flipped": [1, -1],
        "method": "linear",
    },
    "regular_grid_cubic_2d": {"n_per_axis": [10, 40, 100, 200, 400], "method": "cubic"},
    "regular_grid_quintic_2d": {"n_per_axis": [10, 40], "method": "quintic"},
    "griddata": {
        "n_points": 1000,
        "n_grids_per_axis": [10, 100, 1000],
        "methods": ["nearest", "linear", "cubic"],
        "dim": 2,
    },
    "nearest_nd": {"n_samples": [10, 100, 1000], "grid_per_axis": [50, 100, 500]},
    "clough_tocher": {"n_samples": [10, 50, 100], "eval_grid": [50, 50]},
    "rbf_interpolator": {
        "neighbors": [None, 50],
        "n_samples": [10, 100, 1000],
        "kernels": [
            "linear",
            "thin_plate_spline",
            "cubic",
            "quintic",
            "multiquadric",
            "inverse_multiquadric",
            "inverse_quadratic",
            "gaussian",
        ],
        "epsilon": 5.0,
        "domain": [-1.0, 1.0],
        "dim": 2,
    },
}

RUST_ONLY = {
    # Sizes SciPy's suite does not reach. RBF in global mode is O(N^2) in memory
    # (about 800 MB at N = 10^4), so larger N uses `neighbors` only.
    "nearest_nd": {"n_samples": [10_000, 100_000]},
    "linear_nd": {"n_samples": [10_000, 100_000]},
    "regular_grid": {"n_per_axis": [400, 800]},
    "rbf_global": {"n_samples": [3_000, 10_000]},
    "rbf_neighbors": {"n_samples": [10_000, 100_000], "neighbors": 50},
}
