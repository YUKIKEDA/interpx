---
name: publish-review
description: >-
  Read interpx's README, comparison docs, and public rustdoc before crates.io
  publish. Use when preparing a release, running cargo package, running
  cargo publish --dry-run, or checking that a stranger can understand the
  crate from the docs alone.
---

# Publish review

Run this on the commit that will be published, before `cargo package` and before `cargo publish --dry-run`. A single failure stops those commands. Do not run `cargo publish` without `--dry-run` (`.cursor/rules/publish.mdc`).

## Read first

- `.cursor/rules/reader-docs.mdc`
- `.cursor/rules/rust-docs.mdc`
- `.cursor/rules/scipy-reference.mdc`
- Japanese pages in `.cursor/rules/reader-docs.mdc` are inspected by someone other than their author, with `.agents/skills/natural-japanese/SKILL.md` in score mode. The agent running this review does not score or rewrite them.

Do not restate those files here. Judge the pages against them.

## Pages

- `README.md` and `README.ja.md`
- `docs/comparison.md` and `docs/comparison.ja.md`, if they exist
- Crate-level `//!` and every `///` on an item an external crate can name with the default features
- `Cargo.toml` `description`

`docs/design.ja.md`, `docs/roadmap.md`, and `docs/adr/` are the decision record. Do not fail them for keeping Issue and milestone IDs.

## How to judge

Read as someone who has the crate page and has not opened this repository's Issue tracker.

Fail a page when any of these is true:

- A sentence is false against the code or `Cargo.toml`. The same wording in the README, the Japanese README, the crate rustdoc, and `description` does not make it true.
- A comparison number has no committed source and the sentence does not say where it came from (hardware, SciPy and crate versions, data, command).
- The README is not in the order `reader-docs` requires, or the English and Japanese sections diverge.
- The README still contains comparison tables or figures.
- A reader must know an Issue or PR number, a milestone or roadmap ID, `roadmap.md`, `AGENTS.md`, `.cursor/rules`, Grill, or a developer `just` recipe.
- Public rustdoc fails `.cursor/rules/rust-docs.mdc`, including Examples and Errors. `missing_docs` passing does not clear this.
- The first example is not pasteable after the documented dependency, or it uses `unwrap`.
- The README omits, in the place `reader-docs` gives it, any of: what the crate does (first); the boundary condition as an argument, the fill value at the call, `Extrapolate`, the input checks, owned types and `Send + Sync`, `f64` only for scattered types, C-order `values` for grids, NaN handling, the O(N²) RBF memory cost, the 0.x break policy, and the out-of-scope list (after the first example, under what a caller must know).
- The README names SciPy, NumPy, or their functions outside its Verification section and the link to the comparison page, or the Verification section leaves out a place where interpx computes something other than SciPy 1.18 (the differences in `docs/design.ja.md`).
- Public rustdoc, the crate `//!`, the `Cargo.toml` description, or `keywords` names SciPy, NumPy, or their functions (ADR 0011).
- The README links relatively at a file that is not in the package, or embeds an image.
- A Japanese page was scored or rewritten by the agent that edited it, or a reported gloss of the English remains.
- A public item that an external crate can name with the default features is missing from the README Usage section, has no usage example there, or has no example in its rustdoc.
- A public signature, error type, or rustdoc link exposes a `convx`, `kiddo`, or `faer` type (`.cursor/rules/api-contract.mdc`).
- `cargo package --list` differs from the allow list in `Cargo.toml`.

## Report

For each failure, name the file, the text, the rule it breaks, and why a new reader stops. Fix the failures in the publish Issue. When none remain, the package build and the dry-run are allowed.
