# 0010: The size budget of the golden data

- Status: Accepted
- Date: 2026-10-12
- Issue: #70

## Context

The golden data (`compare/goldens/*.json`, SciPy's results) is part of the published package, so that `cargo test` passes on the crate as downloaded. Its size has two limits, checked by the generator (`compare/common/golden.py`) and by `tests/goldens.rs`: 1 MB a file and 4 MB for the directory. The 4 MB was chosen against the 10 MB that crates.io accepts.

After M1-13 the directory holds 3.72 MB. Each row of M1 that added an interpolator added 40 to 120 KB, and nine rows of M2 to M4 are still to come. Twice a row reduced its cases to stay under the limit before it knew whether the limit was the right one.

The 10 MB of crates.io is the size of the compressed package, not of its files. Measured on the 3.72 MB: gzip brings the goldens to 0.77 MB, about a fifth. The sources, tests, and benches are 0.71 MB before compression. The 4 MB limit therefore kept the package at a tenth of what crates.io accepts, and pressed on the test data for no gain.

## Decision

- **The directory may hold 8 MB**, measured as the files are, not compressed. That is about 1.6 MB compressed, well inside the 10 MB.
- **A file may hold 1 MB**, as before. This limit is for diffs and reviews, not for the package: a family that outgrows its file gets one of its own (`docs/conventions.md`).
- **The goldens stay in the package.**
- **A limit stays.** A row whose cases would take the directory over 90% of it opens a Task to decide again before it adds them, and does not drop cases to fit.

Not taken:

- **Leaving the goldens out of the package.** The package would be smaller, and `cargo test` on the published crate would fail in every test that reads a golden file. Making those tests skip when the files are missing adds a second way for the suite to pass, one that checks nothing.
- **Measuring the compressed size.** It is the real constraint, but both the generator and the Rust test would have to compress the files, and the Rust side would need a dependency for it. The raw size with a factor of five in hand is enough.

## Consequences

The rows of M2 to M4 can add the cases their tests need. At the pace of M1 the directory reaches about 4.7 MB at v0.1.0.

The limit is no longer a reason to thin out a case. The curve cases that M1-13 had thinned to fit (queries at the midpoints between the data, and two of the three parametrizations in 3D) are restored in the change that raises the limit. The banded matrices that M1-12 had made smaller (80 rows for 150) stay as they are: a larger matrix of the same kind adds no case. A reason remains in the 1 MB per file and in the time `cargo test` takes to parse the files.

If the directory nears 8 MB, the question is asked again with the sizes of that day: this record is what it starts from.
