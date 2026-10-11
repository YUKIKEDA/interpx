---
name: golden-data
description: >-
  Add or regenerate SciPy 1.18 golden data for interpx tests. Use when a new
  interpolator needs numeric-parity tests, when a golden file must change, or
  when the pinned SciPy version moves.
---

# Golden data

Golden data is the committed record of what SciPy 1.18 returns, so that `cargo test` can check interpx without Python. The rules are in `.cursor/rules/scipy-parity.mdc` and `.cursor/rules/numeric-tests.mdc`. The tolerances are in `docs/design.ja.md` §8.2.

> The pipeline itself (`compare/pyproject.toml`, `compare/generate_*.py`, `compare/goldens/`, `tests/common/`, and the `gen-goldens` recipe) is created by roadmap row M1-1. Until that row is `Done`, this skill states the contract. Do not create the files ahead of that row (`.cursor/rules/layout.mdc`).

## Contract

- Environment: `compare/` is a `uv` project. SciPy and NumPy are pinned in `pyproject.toml` and `uv.lock`. SciPy is 1.18. `uv run --directory compare python -X utf8 <script>`.
- One generator per family, `compare/generate_<family>.py` (for example `oned`, `curve`, `grid`, `scattered`, `rbf`). Each writes `compare/goldens/<family>.json`.
- A golden file records: the SciPy and NumPy versions, the seed, the generator name, and for each case the input (`x`, `y` or `values`, options), the query points, and the expected output. Cases that SciPy rejects record the error class as an expected error.
- Random inputs use a NumPy `Generator` with a fixed seed. No data depends on the machine.
- The Rust side reads the JSON in `tests/` through `tests/common`. `cargo test` never runs Python and never writes the files.
- The package `include` list has `compare/goldens/` and nothing else of `compare/`.

## Add a case

1. Find the SciPy call that is the reference. Record which class and options in the case.
2. Cover the families `.cursor/rules/spec-checks.mdc` asks for: boundary conditions, extrapolation modes, dtypes where the type is generic, dimensions, degenerate and near-degenerate inputs, and the error cases the design names.
3. Include cases where each term of the algorithm dominates (nearly coincident knots, large aspect ratio, values over many orders of magnitude).
4. Run the generator. Read the diff of the JSON: only the intended cases are added.
5. Write the Rust test with fixed tolerances from §8.2, and do the mutation check from `numeric-tests.mdc`.

## Regenerate

Regeneration is its own change. The PR says why, lists the cases whose expected values moved, and for each says whether SciPy changed or the generator did. Never edit a JSON file by hand.

## Move the SciPy version

A `design` Issue first: update §5 of `docs/design.ja.md` and the differences listed against SciPy. Then regenerate everything, and review the moved cases one by one.
