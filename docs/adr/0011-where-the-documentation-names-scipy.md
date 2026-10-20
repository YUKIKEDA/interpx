# 0011: Where the documentation names SciPy

- Status: Accepted
- Date: 2026-10-20
- Issue: #215

## Context

ADR 0001 made SciPy 1.18 the reference the results are verified against: the golden data is what it returns, and the tests hold interpx to it. The published pages of 0.1.0 named SciPy much more widely: in the README's first line, the `Cargo.toml` description, the crate rustdoc, the keywords, a section before the first example, and the rustdoc of most types. The rules asked for that, since the first design.

The owner, 2026-10-20: SciPy is used to verify interpx.

## Decision

- **SciPy 1.18 verifies the results.** The golden data, the numeric tests, the tolerances of `docs/design.ja.md` §8.2, and the comparison (`docs/performance.md`) keep it. ADR 0001 stands.
- **The published documentation names SciPy in two places only:** the README's Verification section (the values are checked against SciPy 1.18, and where interpx computes something else, what and why) and the README's link to the comparison page. The README's first line, the crate rustdoc, the rustdoc of the public API, the `description`, and the `keywords` say what interpx computes, in its own terms.
- **What SciPy is used for is written; what it is not is never written,** in the pages, the rules, the design, an Issue, or a pull request. A denial invites the suspicion it denies (the owner).
- **Type names stay.** They are names of methods, and a rename breaks the API.
- **The records keep their references.** Code comments, private items, tests, goldens, the design text, the ADRs, and the comparison page name SciPy where it is the source of a method or the reference of a result. The published 0.1.0 section of `CHANGELOG.md` stays as it was published.

## Consequences

`scipy-parity.mdc` becomes `scipy-reference.mdc`, about verification. `reader-docs.mdc`, `rust-docs.mdc`, `scope.mdc`, and the `publish-review` and `golden-data` skills follow. The pages change in V-2 and reach crates.io with the next published version.
