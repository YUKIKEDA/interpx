"""Generate `goldens/rbf.json` and `goldens/rbf_system.json`: layer-A goldens
for `RBFInterpolator`.

Reference: SciPy 1.18 `RBFInterpolator` with all eight kernels. For each kernel:
the defaults, `smoothing > 0`, a polynomial degree one above the kernel's
minimum, and the local `neighbors` mode. A case that SciPy rejects stores the
exception class as `expected_error`.

`rbf_system.json` holds the parts the interpolant is made of, from SciPy's
internal functions (`system_cases`): the matrix and the right-hand side of its
system, the shift and scale of the polynomial terms, the powers of the
monomials, and the rows of some evaluations. The points are few: the matrix
has the square of their number of entries.

Run: `just gen-goldens`
"""

import numpy as np
from scipy.interpolate import RBFInterpolator, _rbfinterp_np, _rbfinterp_pythran

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


def system_case(name, kernel, epsilon, degree, smoothing, points, values, queries, note) -> dict:
    """The parts of one interpolant. `degree` is -1 for no polynomial terms;
    `epsilon` is 1 for a kernel that takes none."""
    n = len(points)
    powers = _rbfinterp_np._monomial_powers(points.shape[1], degree, np)
    lhs, rhs, shift, scale = _rbfinterp_pythran._build_system(
        points, values[:, None], np.full(n, smoothing), kernel, epsilon, powers
    )
    rows = _rbfinterp_pythran._build_evaluation_coefficients(queries, points, kernel, epsilon, powers, shift, scale)
    # The two functions that `_build_system` is made of agree with it.
    assert np.array_equal(_rbfinterp_pythran._kernel_matrix(points * epsilon, kernel), lhs[:n, :n] - smoothing * np.eye(n)) or smoothing != 0.0, name
    assert np.array_equal(_rbfinterp_pythran._polynomial_matrix((points - shift) / scale, powers), lhs[:n, n:]), name
    options = {"kernel": kernel, "epsilon": epsilon, "degree": degree, "smoothing": smoothing, "note": note}
    return golden.case(
        name,
        "scipy.interpolate._rbfinterp_pythran._build_system",
        options,
        points=points,
        values=values,
        queries=queries,
        powers=powers,
        matrix=lhs,
        rhs=rhs.ravel(),
        shift=shift,
        scale=scale,
        # The rows of the evaluations at the queries: what the case is asked for.
        expected=rows,
    )


def system_cases() -> list[dict]:
    rng = np.random.default_rng(SEED + 1)
    out = []
    pts = rng.uniform(-1.0, 3.0, (9, 2))
    vals = rng.uniform(-1.0, 1.0, 9)
    q = np.vstack([rng.uniform(-1.5, 3.5, (4, 2)), pts[:2]])
    # Every kernel with its default degree, an `epsilon` that is not 1 where
    # the kernel takes one, and smoothing.
    for kernel, min_degree in MIN_DEGREE.items():
        epsilon = 1.0 if kernel in SCALE_FREE else 0.7
        out.append(system_case(f"system/{kernel}/default", kernel, epsilon, max(min_degree, 0), 0.25, pts, vals, q, "default degree"))
    # Every degree from none to 3, on one kernel that needs no terms.
    for degree in (-1, 0, 1, 2, 3):
        pts3 = rng.uniform(0.0, 1.0, (12, 2))
        out.append(system_case(f"system/gaussian/degree_{degree}", "gaussian", 1.3, degree, 0.0, pts3, rng.uniform(-1.0, 1.0, 12), q[:3], "every degree"))
    # One and three dimensions.
    line = rng.uniform(-2.0, 2.0, (7, 1))
    out.append(system_case("system/quintic/d1", "quintic", 1.0, 2, 0.0, line, rng.uniform(-1.0, 1.0, 7), rng.uniform(-3.0, 3.0, (4, 1)), "one dimension"))
    space = rng.uniform(0.0, 5.0, (11, 3))
    out.append(system_case("system/cubic/d3", "cubic", 1.0, 2, 0.1, space, rng.uniform(-1.0, 1.0, 11), rng.uniform(0.0, 5.0, (4, 3)), "three dimensions, 10 monomials"))
    # An axis on which all points agree keeps its scale.
    flat = np.column_stack([rng.uniform(0.0, 4.0, 8), np.full(8, 2.5)])
    out.append(system_case("system/linear/flat_axis", "linear", 1.0, 0, 0.0, flat, rng.uniform(-1.0, 1.0, 8), np.array([[1.0, 2.5], [3.0, 4.0]]), "flat axis"))
    # Two points at the same place: the distance 0, where the logarithm of the
    # thin plate spline is not finite.
    twice = np.vstack([pts[:6], pts[:1]])
    out.append(system_case("system/thin_plate_spline/coincident", "thin_plate_spline", 1.0, 1, 0.5, twice, rng.uniform(-1.0, 1.0, 7), q[4:], "coincident points"))
    return out


def chosen_cases() -> list[dict]:
    """Cases chosen for the global solve (M3-3). Appended after the others,
    with a generator of their own, so that those stay as they are."""
    rng = np.random.default_rng(SEED + 2)
    reference = "scipy.interpolate.RBFInterpolator"
    out = []
    f33 = np.array(classic_data.FRANKE_33)
    vals = datasets.franke(1, f33[:, 0], f33[:, 1])
    q = np.vstack([-0.05 + 1.1 * rng.random((24, 2)), f33[:2]])

    def options(kernel, **more):
        base = {"kernel": kernel, "epsilon": None, "smoothing": 0.0, "degree": None, "neighbors": None, "min_degree": MIN_DEGREE[kernel]}
        base.update(more)
        return base

    # An `epsilon` for a kernel that does not depend on a scale, with
    # smoothing: the same interpolant as the smoothing over `epsilon^p`
    # without it, which is checked here and is what interpx is given.
    for kernel, power in (("linear", 1), ("thin_plate_spline", 2), ("cubic", 3), ("quintic", 5)):
        epsilon, smoothing = 2.5, 0.1
        want = RBFInterpolator(f33, vals, kernel=kernel, epsilon=epsilon, smoothing=smoothing)(q)
        same = RBFInterpolator(f33, vals, kernel=kernel, smoothing=smoothing / epsilon**power)(q)
        assert np.allclose(want, same, rtol=1e-9, atol=1e-9), kernel
        out.append(golden.case(
            f"rbf/{kernel}/epsilon_as_smoothing/franke_33/f1", reference,
            options(kernel, scipy_epsilon=epsilon, scipy_smoothing=smoothing, smoothing=smoothing / epsilon**power, function="franke_f1"),
            points=f33, values=vals, queries=q, expected=want,
        ))
    # No polynomial terms, for a kernel that needs none and for one that does.
    for kernel, epsilon in (("gaussian", 3.0), ("thin_plate_spline", None)):
        kwargs = {"kernel": kernel, "degree": -1}
        if epsilon is not None:
            kwargs["epsilon"] = epsilon
        out.append(golden.case(
            f"rbf/{kernel}/no_polynomial/franke_33/f1", reference,
            options(kernel, epsilon=epsilon, degree=-1, function="franke_f1"),
            points=f33, values=vals, queries=q, expected=RBFInterpolator(f33, vals, **kwargs)(q),
        ))
    # Heavy smoothing: close to the polynomial that fits best.
    out.append(golden.case(
        "rbf/thin_plate_spline/heavy_smoothing/franke_33/f1", reference,
        options("thin_plate_spline", smoothing=50.0, function="franke_f1"),
        points=f33, values=vals, queries=q, expected=RBFInterpolator(f33, vals, smoothing=50.0)(q),
    ))
    # One dimension, and four.
    line = np.sort(rng.uniform(-2.0, 2.0, 15))[:, None]
    lv = np.sin(2.0 * line[:, 0])
    ql = np.linspace(-2.5, 2.5, 21)[:, None]
    out.append(golden.case(
        "rbf/cubic/default/line_15/sine", reference, options("cubic", function="sine"),
        points=line, values=lv, queries=ql, expected=RBFInterpolator(line, lv, kernel="cubic")(ql),
    ))
    h4 = datasets.halton(40, 4)
    v4 = datasets.genz("product_peak", h4)
    q4 = rng.random((20, 4))
    out.append(golden.case(
        "rbf/multiquadric/default/halton_40_d4/product_peak", reference,
        options("multiquadric", epsilon=EPSILON, function="genz_product_peak"),
        points=h4, values=v4, queries=q4, expected=RBFInterpolator(h4, v4, kernel="multiquadric", epsilon=EPSILON)(q4),
    ))
    # Two values at every point. `values` and `expected` hold the first
    # column, as in every other case.
    both = np.column_stack([vals, datasets.franke(4, f33[:, 0], f33[:, 1])])
    case = golden.case(
        "rbf/thin_plate_spline/columns/franke_33/f1_f4", reference, options("thin_plate_spline", function="franke_f1_f4"),
        points=f33, values=vals, queries=q, expected=RBFInterpolator(f33, vals)(q),
    )
    case["values_columns"] = golden.encode(both)
    case["expected_columns"] = golden.encode(RBFInterpolator(f33, both)(q))
    out.append(case)
    # Points along an axis: a polynomial of degree 1 in two variables is not
    # determined by them, and SciPy says so.
    t = np.linspace(0.0, 1.0, 8)
    flat = np.column_stack([t, np.zeros(8)])
    try:
        RBFInterpolator(flat, np.sin(3.0 * t))
        raise SystemExit("points along an axis were expected to be singular")
    except np.linalg.LinAlgError as exc:
        out.append(golden.case(
            "rbf/thin_plate_spline/default/along_an_axis/sine", reference, options("thin_plate_spline", function="sine"),
            points=flat, values=np.sin(3.0 * t), queries=flat[:2], expected_error=type(exc).__name__,
        ))
    return out


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

    path = golden.write("rbf_system", GENERATOR, SEED, system_cases())
    print(f"wrote {path.name}")
    cases += chosen_cases()
    path = golden.write("rbf", GENERATOR, SEED, cases)
    print(f"wrote {path.name} ({len(cases)} cases)")


if __name__ == "__main__":
    main()
