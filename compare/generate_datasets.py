"""Generate `goldens/datasets.json` and `goldens/functions.json`.

`datasets.json` holds inputs only: the classic 1D data, the classic 2D node sets,
and the adversarial point sets of docs/design.ja.md section 8.4. `functions.json`
holds values of the analytic test functions (Franke, Genz, Rosenbrock) at fixed
points, so that the Rust copies in `tests/common/functions.rs` can be checked.

Run: `just gen-goldens`
"""

import numpy as np

from common import classic_data, datasets, golden

GENERATOR = "generate_datasets.py"
SEED = 0

BURKARDT = "John Burkardt, test_interp / test_interp_2d (MIT), https://people.sc.fsu.edu/~jburkardt/f_src/"


def classic_cases() -> list[dict]:
    one_d = [
        ("fritsch_carlson", classic_data.FRITSCH_CARLSON, "Fritsch and Carlson (1980), SIAM J. Numer. Anal. 17(2)"),
        ("titanium", classic_data.TITANIUM, "de Boor and Rice (1968); de Boor, A Practical Guide to Splines"),
        ("square_with_near_point", classic_data.SQUARE_WITH_NEAR_POINT, "no origin cited"),
    ]
    nodes = [
        ("franke_100", classic_data.FRANKE_100, "Franke's 100 node set"),
        ("franke_33", classic_data.FRANKE_33, "Franke's 33 node set"),
        ("lawson_25", classic_data.LAWSON_25, "Lawson's 25 node set"),
    ]
    cases = []
    for name, pairs, origin in one_d:
        a = np.array(pairs)
        cases.append(
            golden.case(name, "data", {"kind": "classic_1d", "origin": origin, "transcribed_from": BURKARDT}, x=a[:, 0], y=a[:, 1])
        )
    for name, pairs, origin in nodes:
        cases.append(
            golden.case(name, "data", {"kind": "nodes_2d", "origin": origin, "transcribed_from": BURKARDT}, points=np.array(pairs))
        )
    for adv in datasets.adversarial():
        name = adv.pop("name")
        tags = adv.pop("tags")
        cases.append(golden.case(name, "data", {"kind": "adversarial", "tags": tags}, **adv))
    return cases


def function_cases() -> list[dict]:
    rng = np.random.default_rng(SEED)
    g = np.linspace(0.0, 1.0, 5)
    xx, yy = np.meshgrid(g, g, indexing="ij")
    pts2 = np.column_stack([xx.ravel(), yy.ravel()])
    pts2 = np.vstack([pts2, rng.random((15, 2))])
    cases = []
    for k in range(1, 7):
        cases.append(
            golden.case(f"franke_f{k}", "datasets.franke", {"k": k}, queries=pts2, expected=datasets.franke(k, pts2[:, 0], pts2[:, 1]))
        )
    for d in (3, 5, 8):
        pts = rng.random((12, d))
        for kind in datasets.GENZ_KINDS:
            cases.append(
                golden.case(
                    f"genz_{kind}_d{d}",
                    "datasets.genz",
                    {"kind": kind, "d": d, "a": datasets.genz_a(d), "u": 0.5},
                    queries=pts,
                    expected=datasets.genz(kind, pts),
                )
            )
        cases.append(golden.case(f"rosenbrock_d{d}", "datasets.rosenbrock", {"d": d}, queries=pts, expected=datasets.rosenbrock(pts)))
    return cases


def main() -> None:
    counts = {
        "FRITSCH_CARLSON": 11,
        "TITANIUM": 49,
        "SQUARE_WITH_NEAR_POINT": 12,
        "FRANKE_100": 100,
        "FRANKE_33": 33,
        "LAWSON_25": 25,
    }
    for name, n in counts.items():
        assert len(getattr(classic_data, name)) == n, name
    for name, cases in (("datasets", classic_cases()), ("functions", function_cases())):
        path = golden.write(name, GENERATOR, SEED, cases)
        print(f"wrote {path.name} ({len(cases)} cases)")


if __name__ == "__main__":
    main()
