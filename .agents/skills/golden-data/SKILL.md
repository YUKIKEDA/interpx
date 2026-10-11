---
name: golden-data
description: >-
  Add or regenerate SciPy 1.18 golden data for interpx tests. Use when a new
  interpolator needs numeric-parity tests, when a golden file must change, or
  when the pinned SciPy version moves.
---

# Golden data

Golden data is the committed record of what SciPy 1.18 returns, so that `cargo test` can check interpx without Python. The rules are in `.cursor/rules/scipy-parity.mdc` and `.cursor/rules/numeric-tests.mdc`. The tolerances are in `docs/design.ja.md` §8.2.

## Files

- `compare/pyproject.toml` and `compare/uv.lock`: a `uv` project (`package = false`). SciPy is `1.18.*` and Python is 3.12 or newer (SciPy 1.18 does not support 3.11). The lock records the exact SciPy and NumPy.
- `compare/common/golden.py`: `golden.write(name, generator, seed, cases)` checks that SciPy is 1.18, encodes non-finite numbers, and writes `compare/goldens/<name>.json`. `golden.check_baseline()` is the version guard.
- `compare/generate_<family>.py`: one generator per family (`oned`, `curve`, `grid`, `scattered`, `rbf`, ...). `compare/generate_pipeline_check.py` is the first, a `numpy.interp` case set that proves the pipeline.
- `compare/goldens/<name>.json`: the committed output.
- `tests/common/mod.rs`: `load::<Case>(name)` (returns `Result`; tests `.expect` it), `Real` (a number or `"NaN"`, `"Infinity"`, `"-Infinity"`), `reals`, `Tol`, `DIRECT_F64`, `close`, `assert_close`.
- `just gen-goldens`: runs every generator. Add one line per new generator.

## Contract

- Run a generator with `uv run --directory compare python -X utf8 <script>` (the recipe does this).
- A golden file records: `meta` (generator, SciPy, NumPy, seed) and `cases`. A case has `name`, `reference` (the SciPy call), `options`, the inputs, and either `expected` or `expected_error` (the error class SciPy raises).
- Random inputs use a NumPy `Generator` with a fixed seed. The file does not record the Python or OS version, so regenerating on another machine gives the same bytes.
- The Rust side reads the JSON in `tests/` through `tests/common`. `cargo test` never runs Python and never writes the files. Floats are written with Python's shortest round-trip repr and read with `serde_json`'s `float_roundtrip`, so the Rust side sees the same `f64`.
- The package `include` list has `/compare/goldens/` and nothing else of `compare/`. The patterns are anchored with a leading `/`, because `tests/` alone also matches `compare/.venv/**/tests/`.
- Check regeneration: run `just gen-goldens` twice and compare file hashes. The second run must change nothing.

## Add a case

1. Find the SciPy call that is the reference. Record which class and options in the case.
2. Cover the families `.cursor/rules/spec-checks.mdc` asks for: boundary conditions, extrapolation modes, dtypes where the type is generic, dimensions, degenerate and near-degenerate inputs, and the error cases the design names.
3. Include cases where each term of the algorithm dominates (nearly coincident knots, large aspect ratio, values over many orders of magnitude).
4. Run the generator. Read the diff of the JSON: only the intended cases are added.
5. Write the Rust test with fixed tolerances from §8.2, and do the mutation check from `numeric-tests.mdc`: copy the JSON aside, change one expected value (and, separately, one NaN query), see the test fail, and restore the copy.

## Regenerate

Regeneration is its own change. The PR says why, lists the cases whose expected values moved, and for each says whether SciPy changed or the generator did. Never edit a JSON file by hand.

## Move the SciPy version

A `design` Issue first: update §5 of `docs/design.ja.md` and the differences listed against SciPy. Then regenerate everything, and review the moved cases one by one.
