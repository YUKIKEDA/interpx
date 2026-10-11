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

Layout: a list of scalars is one line, a matrix is one row per line. This keeps
the files small (the packaged `compare/goldens/` has a size budget, see
`MAX_FILE_BYTES`).
"""

import json
import math
import sys
from pathlib import Path

import numpy as np
import scipy

SCIPY_BASELINE = "1.18"
GOLDENS = Path(__file__).resolve().parent.parent / "goldens"

# docs/design.ja.md section 8.4 and the packaging budget: each file at most 1 MB,
# the whole directory at most 4 MB (crates.io accepts 10 MB).
MAX_FILE_BYTES = 1_000_000
MAX_TOTAL_BYTES = 4_000_000


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
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def _scalar(v) -> bool:
    return v is None or isinstance(v, (bool, int, float, str))


def dumps(obj, indent: int = 0) -> str:
    """Serialize `obj` with scalar lists on one line and matrices one row per line."""
    pad = " " * indent
    inner = " " * (indent + 2)
    if isinstance(obj, dict):
        if not obj:
            return "{}"
        items = [f"{inner}{json.dumps(k)}: {dumps(v, indent + 2)}" for k, v in obj.items()]
        return "{\n" + ",\n".join(items) + "\n" + pad + "}"
    if isinstance(obj, list):
        if all(_scalar(v) for v in obj):
            return "[" + ", ".join(json.dumps(v, allow_nan=False) for v in obj) + "]"
        if all(isinstance(v, list) and all(_scalar(w) for w in v) for v in obj):
            rows = [inner + dumps(v) for v in obj]
            return "[\n" + ",\n".join(rows) + "\n" + pad + "]"
        items = [inner + dumps(v, indent + 2) for v in obj]
        return "[\n" + ",\n".join(items) + "\n" + pad + "]"
    return json.dumps(obj, allow_nan=False)


def case(name: str, reference: str, options: dict | None = None, **fields) -> dict:
    """Build one case. `fields` are the inputs and `expected` or `expected_error`."""
    out = {"name": name, "reference": reference, "options": options or {}}
    out.update(fields)
    return out


def write(name: str, generator: str, seed: int, cases: list[dict]) -> Path:
    """Write `compare/goldens/<name>.json` and return its path."""
    check_baseline()
    names = [c["name"] for c in cases]
    if len(set(names)) != len(names):
        sys.exit(f"{name}: case names are not unique")
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
    text = dumps(doc) + "\n"
    if len(text.encode("utf-8")) > MAX_FILE_BYTES:
        sys.exit(f"{path.name} is {len(text)} bytes, over the {MAX_FILE_BYTES} budget")
    # A fixed key order and a trailing newline keep a regeneration byte-identical.
    path.write_text(text, encoding="utf-8", newline="\n")
    return path


def check_total_size() -> None:
    """Fail when the goldens directory is over the packaging budget."""
    total = sum(p.stat().st_size for p in GOLDENS.glob("*.json"))
    if total > MAX_TOTAL_BYTES:
        sys.exit(f"goldens total {total} bytes, over the {MAX_TOTAL_BYTES} budget")
    print(f"goldens total: {total} bytes")
