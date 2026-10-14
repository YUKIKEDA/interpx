# Gates named in CONTRIBUTING.md.
set windows-shell := ["powershell.exe", "-NoLogo", "-Command"]

lint:
    cargo fmt --check
    cargo clippy --all-targets -- -D warnings

test:
    cargo test

coverage:
    cargo llvm-cov --fail-under-lines 90

# Time benchmarks (criterion). Not a gate; CI does not run it. Extra arguments go
# to criterion, for example `just bench --save-baseline m1`.
bench *args:
    cargo bench -- {{args}}

# The comparison with SciPy and other Rust crates (not a gate; long: the skill says what sets the time;
# `.agents/skills/perf-compare/SKILL.md`). Extra arguments go to the driver,
# for example `just compare-perf --filter "^grid/"`.
compare-perf *args:
    uv run --directory compare python -X utf8 perf/run.py {{args}}

# Fetch the real-world data of the comparison into `compare/perf/out/real/`,
# and check it (docs/design.ja.md section 8.4.7). Its cases are in
# `compare-perf` once it is there.
compare-fetch:
    uv run --directory compare python -X utf8 perf/real.py

# `docs/performance.md` from the last run of `compare-perf`.
compare-report:
    uv run --directory compare python -X utf8 perf/report.py

# Regenerate golden data with the pinned SciPy. `cargo test` only reads the files.
# Each family adds a line here when its generator exists. The last line checks the set.
gen-goldens:
    uv run --directory compare python -X utf8 -W error generate_datasets.py
    uv run --directory compare python -X utf8 -W error generate_oned.py
    uv run --directory compare python -X utf8 -W error generate_curve.py
    uv run --directory compare python -X utf8 -W error generate_grid.py
    uv run --directory compare python -X utf8 -W error generate_scattered.py
    uv run --directory compare python -X utf8 -W error generate_rbf.py
    uv run --directory compare python -X utf8 -W error generate_tridiagonal.py
    uv run --directory compare python -X utf8 -W error generate_banded.py
    uv run --directory compare python -X utf8 check_goldens.py
