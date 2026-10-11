"""Generate `goldens/pipeline_check.json`.

A small case set that proves the golden pipeline end to end (writer, file
format, Rust loader, tolerance helper). The reference is `numpy.interp`; no
interpx code exists yet when this file is first used. Real families have their
own generators (`generate_<family>.py`).

Run: `just gen-goldens`
"""

import numpy as np

from common import golden

GENERATOR = "generate_pipeline_check.py"
SEED = 0


def case(name, x, y, queries, **options):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    q = np.asarray(queries, dtype=np.float64)
    return {
        "name": name,
        "reference": "numpy.interp",
        "options": options,
        "x": x,
        "y": y,
        "queries": q,
        "expected": np.interp(q, x, y, **options),
    }


def main() -> None:
    rng = np.random.default_rng(SEED)

    x_small = [0.0, 1.0, 2.0, 4.0]
    y_small = [0.0, 2.0, 1.0, 5.0]
    # Knots, interior points, both ends, outside on each side, and a NaN query.
    q_small = [0.0, 0.25, 1.0, 1.5, 3.0, 4.0, -1.0, 5.0, float("nan")]

    x_rand = np.sort(rng.uniform(-10.0, 10.0, size=16))
    y_rand = rng.normal(size=16)
    q_rand = rng.uniform(-12.0, 12.0, size=32)

    cases = [
        case("small_clamped", x_small, y_small, q_small),
        case("small_fill", x_small, y_small, q_small, left=-1.0, right=-2.0),
        case("random_clamped", x_rand, y_rand, q_rand),
    ]
    path = golden.write("pipeline_check", GENERATOR, SEED, cases)
    print(f"wrote {path} ({len(cases)} cases)")


if __name__ == "__main__":
    main()
