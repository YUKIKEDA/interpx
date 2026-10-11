"""Write golden JSON files for interpx.

File format (read by `tests/common/mod.rs`):

    {
      "meta": {"generator", "scipy", "numpy", "seed"},
      "cases": [{"name", "reference", "options", ...inputs..., "expected" | "expected_error"}]
    }

Real numbers are JSON numbers. A non-finite value is the string "NaN",
"Infinity", or "-Infinity", because JSON has no spelling for them. Floats are
written with Python's shortest round-trip repr, so the Rust side reads the exact
same f64.
"""

import json
import math
import sys
from pathlib import Path

import numpy as np
import scipy

SCIPY_BASELINE = "1.18"
GOLDENS = Path(__file__).resolve().parent.parent / "goldens"


def check_baseline() -> None:
    """Fail unless SciPy is the baseline named in docs/design.ja.md."""
    if not scipy.__version__.startswith(SCIPY_BASELINE + "."):
        sys.exit(f"SciPy {SCIPY_BASELINE}.x is the baseline, got {scipy.__version__}")


def encode(value):
    """Return a JSON-safe copy of `value` (arrays, scalars, nested lists)."""
    if isinstance(value, np.ndarray):
        return encode(value.tolist())
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    if isinstance(value, (np.floating, float)):
        v = float(value)
        if math.isnan(v):
            return "NaN"
        if math.isinf(v):
            return "Infinity" if v > 0 else "-Infinity"
        return v
    if isinstance(value, np.integer):
        return int(value)
    return value


def write(name: str, generator: str, seed: int, cases: list[dict]) -> Path:
    """Write `compare/goldens/<name>.json` and return its path."""
    check_baseline()
    doc = {
        "meta": {
            "generator": generator,
            "scipy": scipy.__version__,
            "numpy": np.__version__,
            "seed": seed,
        },
        "cases": encode(cases),
    }
    GOLDENS.mkdir(parents=True, exist_ok=True)
    path = GOLDENS / f"{name}.json"
    # Fixed key order and a trailing newline keep a regeneration byte-identical.
    text = json.dumps(doc, indent=2, allow_nan=False) + "\n"
    path.write_text(text, encoding="utf-8", newline="\n")
    return path
