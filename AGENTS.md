# interpx: agent entry

Pure Rust interpolation covering the main classes of `scipy.interpolate` (SciPy 1.18). This file is the entry point. Steps for people are in `CONTRIBUTING.md`. Enforced rules are in `.cursor/rules/`.

## Source

| What                           | Where                                                         |
| ------------------------------ | ------------------------------------------------------------- |
| Design (canonical, Japanese)   | `docs/design.ja.md`                                           |
| Order and status of work       | `docs/roadmap.md`                                             |
| Modules and directory tree     | `docs/conventions.md`                                         |
| Why a decision was made        | `docs/adr/` (`0001-...`, one decision per file)               |
| Build, test, and gate commands | `justfile`                                                    |

The current row's ID lives only in `docs/roadmap.md`.

## Required rules

Read the rules that match the work. `alwaysApply` rules apply to every change. After every merge comes the retrospective of `workflow`: what the work taught goes into these files, not into an agent's memory.

| Rule                        | One line                                                                      |
| --------------------------- | ----------------------------------------------------------------------------- |
| `workflow`                  | Grill, Issue, branch, PR, human merge, retrospective; gates; language         |
| `consent`                   | A question is not a decision; options need a recommendation                   |
| `git`                       | Destructive git needs a named order; commit and push only when asked          |
| `conventional-commits`      | `type(scope): subject`; no `perf` type                                        |
| `pull-requests`             | Title and six-heading body from the template                                  |
| `github-closing`            | A closing verb before `#N` closes it; `Closes #N` alone on a line             |
| `defer`                     | Deferral only as a named roadmap row                                          |
| `scope`                     | Requirements decide scope; no half-shipped capability; one v0.1.0 release     |
| `layout`                    | Single crate; `compare/` is isolated; no empty files ahead of need            |
| `scratch`                   | Remove scratch files in the turn that made them                               |
| `record-corrections`        | A correction lands in the rule, design, ADR, or skill it belongs to           |
| `open-ends`                 | A PR leaves no path untimed and no figure unexplained                         |
| `similar-findings`          | Fix the same kind of problem where it also occurs                             |
| `powershell-shell`          | Shell syntax on this setup                                                    |
| `no-dead-code`              | Remove everything that only served a removed path                             |
| `rust`, `rust-api`, `types` | Safety, naming, invalid states unrepresentable                                |
| `rust-docs`                 | English rustdoc with `# Errors` and examples                                  |
| `numerics`                  | f32/f64, NaN, float comparison, stability                                     |
| `api-contract`              | Owned, `Send + Sync`, typed errors, no `convx` types in the API               |
| `scipy-parity`              | SciPy 1.18 is the reference; goldens come from the pinned script              |
| `numeric-tests`             | Independent reference, fixed tolerances, mutation checks                      |
| `spec-checks`               | A test implements the whole design statement                                  |
| `perf`                      | No allocation on evaluation; layout; dispatch; parallel and SIMD bit equality |
| `bench`                     | criterion, allocation test, memory, comparison with SciPy and Rust crates     |
| `dev-docs`                  | Mermaid and LaTeX in the design                                               |
| `reader-docs`               | README order; Japanese pages are checked by someone else                      |
| `publish`                   | What "ready" means; order of release steps                                    |

## Skills

Skills live in `.agents/skills/`. `skills-lock.json` pins the ones that come from outside.

- `grilling`: when the design branches (before an Issue).
- `github-closing`: before `gh pr create`, to check that the PR closes only what it finishes.
- `golden-data`: add or regenerate SciPy golden data.
- `perf-compare`: run the comparison with SciPy and other Rust crates.
- `cpu-profile`: sample a release binary on Windows.
- `publish-review`: read README, comparison docs, and public rustdoc before a release.
- `natural-japanese`: write and score Japanese pages (not by their author).

## External references

Rust API Guidelines and the rustdoc book. The reference for numerics is SciPy 1.18.
