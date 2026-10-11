## Summary

<!-- What and why. -->

## Related

<!-- `Closes #N` on its own line when this PR finishes that Issue. `Part of #N` when it does not. N/A if none. -->

Closes #N

## Test plan

<!-- Commands and cases. -->

- [ ] `just lint`
- [ ] `just test`
- [ ] `just bench` (hot path only)

## Verification

<!-- What you actually ran. For a hot-path change: criterion numbers against the named baseline, or why speed does not change. N/A if not yet. -->

## Risk / Rollback

N/A

## Checklist

- [ ] PR title is Conventional Commits (`type(scope): subject`; type and scope in English)
- [ ] Matches `AGENTS.md`, `.cursor/rules/`, and `docs/roadmap.md`
- [ ] No `unwrap` / `expect` / `panic` / `dbg` on library paths
- [ ] No `unsafe` in interpx's own code
- [ ] Evaluation paths do not allocate (`tests/alloc.rs` passes)
- [ ] rustdoc in English for new public API
- [ ] A new numeric behavior matches SciPy 1.18, or the difference is written in `docs/design.ja.md` or an ADR
- [ ] This PR does not leave one path optimized or supported and its pair not
- [ ] Any deferral is a roadmap row that names when, or the item is already under Intentionally out of scope
- [ ] Docs updated when behavior or workflow changed, or N/A
