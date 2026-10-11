"""Generate `goldens/rbf.json`: layer-A goldens for `RBFInterpolator`.

Reference: SciPy 1.18 `RBFInterpolator` with all eight kernels. For each kernel:
the defaults, `smoothing > 0`, a polynomial degree one above the kernel's
minimum, and the local `neighbors` mode. A case that SciPy rejects stores the
exception class as `expected_error`.

Run: `just gen-goldens`
"""

import numpy as np
from scipy.interpolate import RBFInterpolator

from common import classic_data, datasets, golden

GENERATOR = "generate_rbf.py"
SEED = 5

# Minimum polynomial degree per kernel in SciPy (-1: no polynomial needed).
MIN_DEGREE = {
    "linear": 0,
    "thin_plate_spline": 1,
    "cubic": 1,
    "quintic": 2,
    "multiquadric": 0,
    "inverse_multiquadric": -1,
    "inverse_quadratic": -1,
    "gaussian": -1,
}
SCALE_FREE = ("linear", "thin_plate_spline", "cubic", "quintic")
EPSILON = 3.0


def main() -> None:
    rng = np.random.default_rng(SEED)
    f33 = np.array(classic_data.FRANKE_33)
    h64 = datasets.halton(64, 2)
    h50 = datasets.halton(50, 3)
    data = [
        ("franke_33/f1", f33, datasets.franke(1, f33[:, 0], f33[:, 1]), "franke_f1"),
        ("halton_64/f4", h64, datasets.franke(4, h64[:, 0], h64[:, 1]), "franke_f4"),
        ("halton_50_d3/gaussian_peak", h50, datasets.genz("gaussian_peak", h50), "genz_gaussian_peak"),
    ]
    cases = []
    for label, pts, vals, function in data:
        d = pts.shape[1]
        q = np.vstack([-0.05 + 1.1 * rng.random((28, d)), pts[:2]])
        for kernel, min_degree in MIN_DEGREE.items():
            epsilon = None if kernel in SCALE_FREE else EPSILON
            variants = {
                "default": {},
                "smoothing": {"smoothing": 1e-3},
                "degree": {"degree": min_degree + 1},
                "neighbors": {"neighbors": 20},
            }
            for vname, extra in variants.items():
                kwargs = {"kernel": kernel, **extra}
                if epsilon is not None:
                    kwargs["epsilon"] = epsilon
                options = {
                    "kernel": kernel,
                    "epsilon": epsilon,
                    "smoothing": extra.get("smoothing", 0.0),
                    "degree": extra.get("degree"),
                    "neighbors": extra.get("neighbors"),
                    "min_degree": min_degree,
                    "function": function,
                }
                name = f"rbf/{kernel}/{vname}/{label}"
                common = {"points": pts, "values": vals, "queries": q}
                try:
                    f = RBFInterpolator(pts, vals, **kwargs)
                    cases.append(golden.case(name, "scipy.interpolate.RBFInterpolator", options, **common, expected=f(q)))
                except Exception as exc:  # recorded so that the difference is visible in review
                    print(f"  {name}: {type(exc).__name__}: {exc}")
                    cases.append(
                        golden.case(name, "scipy.interpolate.RBFInterpolator", options, **common, expected_error=type(exc).__name__)
                    )

    path = golden.write("rbf", GENERATOR, SEED, cases)
    print(f"wrote {path.name} ({len(cases)} cases)")


if __name__ == "__main__":
    main()
