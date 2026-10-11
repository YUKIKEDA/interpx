# 0004: Thread safety, parallelism, allocation, and reproducibility

- Status: Accepted
- Date: 2026-10-11
- Issue: N/A (grill, pre-Issue)

## Context

The first design said interpolators have "no lifetime (`'static`)", which is not the same as `T: 'static`. A shared hint for a hinted search would break thread safety. The trait had no parallel batch method.

## Decision

- Interpolators are owned and have no lifetime parameter. `T: 'static` is not required. All public interpolators are `Send + Sync`, with no interior mutability. `&self` evaluation is pure.
- A search hint is a caller-owned value or a batch-local variable.
- `eval_batch_par` and `eval_batch_into_par` are separate methods behind the `rayon` feature, not default trait methods. Their results are bit for bit equal to the serial forms. `rayon` is used only at the outermost level.
- `eval`, `eval_or`, `eval_into`, and `eval_batch_into` do not allocate. `tests/alloc.rs` checks this with a counting allocator.
- Explicit SIMD is not in v0.1.0 unless measurement shows a bottleneck; if added, it is safe-crate based and bit-identical to the scalar path.
- `unsafe_code = "forbid"` for interpx's own code.

## Consequences

A caller can share an interpolator across threads. Tests can compare serial and parallel results exactly. Hot paths stay allocation-free and auditable.
