# 0011: SciPy is the verification reference only

- Status: Accepted
- Date: 2026-10-20
- Issue: #215

## Context

ADR 0001 made SciPy 1.18 the reference: the golden data is what it returns, and the tests hold interpx to it. The published pages of 0.1.0 went further. The README's first line, the `Cargo.toml` description, and the crate rustdoc said "covering the main classes of SciPy 1.18's `scipy.interpolate`"; `scipy` was a keyword; a "Coming from SciPy" section came before the first example; the rustdoc of most types began "SciPy's `X`", mapped SciPy's arguments onto interpx's, and listed the differences. The rules asked for all of it (`reader-docs.mdc`, `rust-docs.mdc`, `scipy-parity.mdc`, the `publish-review` skill), since the first design.

Read from the crate page, that is a port of SciPy, and a relation to the SciPy project that does not exist. The owner, 2026-10-20: SciPy is used to verify interpx, not for compatibility, and not to offer the same API.

## Decision

- **SciPy 1.18 is the reference interpx is verified against.** The golden data, the numeric tests, the tolerances of `docs/design.ja.md` §8.2, and the comparison (`docs/performance.md`) keep it. ADR 0001 stands.
- **It is not a compatibility target and not an API to reproduce.** Behavior is chosen for interpx; a difference from SciPy found in verification is fixed or written in the design, as before, but it is not a defect by being a difference.
- **The published documentation does not put SciPy forward.** The README, the crate rustdoc, and the description start from what interpx computes. The README names SciPy in one short section, "Verification": the values are checked against SciPy 1.18, and where interpx computes something else, what and why. Public rustdoc names neither SciPy nor NumPy. `scipy` is not a keyword.
- **No disclaimer of affiliation.** There is no relation, so there is nothing to disclaim (the owner).
- **Type names stay.** `RBFInterpolator`, `PchipInterpolator`, `CloughTocher2DInterpolator`, and the rest are names of methods, not of SciPy, and a rename breaks the API.
- **The records keep their references.** Code comments, private items, tests, goldens, the design text, the ADRs, and the comparison page name SciPy where it is the source of a method or the reference of a result. The published 0.1.0 section of `CHANGELOG.md` stays as it was published.

## Consequences

`scipy-parity.mdc` becomes `scipy-reference.mdc`, about verification only: it no longer asks a type to follow SciPy's defaults and method names, to start from SciPy's signature, or to state its SciPy counterpart. `reader-docs.mdc`, `rust-docs.mdc`, `scope.mdc`, and the `publish-review` and `golden-data` skills follow. The pages themselves change in V-2, and reach crates.io with the next published version.
