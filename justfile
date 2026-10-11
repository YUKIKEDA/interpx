# Gates named in CONTRIBUTING.md.
set windows-shell := ["powershell.exe", "-NoLogo", "-Command"]

lint:
    cargo fmt --check
    cargo clippy --all-targets -- -D warnings

test:
    cargo test

coverage:
    cargo llvm-cov --fail-under-lines 90

# Regenerate golden data with the pinned SciPy. `cargo test` only reads the files.
# Each family adds a line here when its generator exists.
gen-goldens:
    uv run --directory compare python -X utf8 generate_pipeline_check.py
