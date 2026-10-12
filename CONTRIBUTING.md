# Contributing

## Order

```
Grill (when the design branches) → Issue → branch → work → PR → human review → squash merge → retrospective
```

1. **Grill** when the API, the meaning, or a milestone boundary branches. Skill: `.agents/skills/grilling/SKILL.md`.
2. **Issue** from a template in `.github/ISSUE_TEMPLATE/`. One Issue is one PR. Blank issues are off.

   | Template | Label           | Commit type            |
   | -------- | --------------- | ---------------------- |
   | bug      | `type:bug`      | `fix`                  |
   | feat     | `type:feat`     | `feat`                 |
   | refactor | `type:refactor` | `refactor`             |
   | test     | `type:test`     | `test`                 |
   | design   | `type:docs`     | `docs`                 |
   | task     | `type:chore`    | `chore`, `build`, `ci` |
   | spike    | `type:spike`    | `chore`, `test`        |

3. A change to what the spec means is a **design** Issue. It updates `docs/design.ja.md` and lands before the implementation Issue. The Japanese file is canonical.
4. **Acceptance text** goes on the Issue. `docs/roadmap.md` holds ID, title, Issue, and status only. A row that reads `Set after Grill on #n` is not started.
5. **Branch** from `main`: `type/<issue-number>-<slug>`, for example `feat/12-pchip`.
6. **PR** from `.github/pull_request_template.md`. Write `Closes #N` alone on a line under `## Related`. The title is a Conventional Commits line.
7. A **human reviews and squash-merges**. Agents do not merge. After the merge, delete the local and remote branch.
8. **Retrospective** on the merged PR, before the next feature branch (below).

## Retrospective

After a squash merge, write three points as a comment on the merged PR. Write "none" when a point is empty.

1. Whether a lesson belongs in `.cursor/rules/`, a skill, `docs/design.ja.md`, or `docs/adr/`. Something the work relied on that no file states is a lesson.
2. Which Issue will carry it, if any.
3. Whether to refactor. A yes opens a `refactor` Issue, and that Issue is next, before the next feature. Postponing it needs a roadmap row that names when.

The lessons of point 1 then go into the files they name, in a `chore` Issue and one PR, opened before the next feature branch. An agent does this without being asked. There is no retrospective log.

Decisions that need their reason kept become an ADR: `docs/adr/NNNN-slug.md` with Status, Date, Issue, Context, and Decision.

## Gates

| Command         | What                                                                |
| --------------- | ------------------------------------------------------------------- |
| `just lint`     | `cargo fmt --check` and `cargo clippy --all-targets -- -D warnings` |
| `just test`     | `cargo test` (reads golden JSON; does not run Python)               |
| `just coverage` | line coverage at least 90% (`cargo-llvm-cov`)                       |

`cargo bench` and the comparison runs in `compare/perf/` are not gates. Performance claims need numbers (`.cursor/rules/bench.mdc`).

CI also runs the tests on Linux, Windows, and macOS with `--no-default-features` and `--all-features`, and on the MSRV (1.89).

## Release

`v0.1.0` is published once, after M1 to M5 of `docs/roadmap.md` are done. There is no earlier release. See `.cursor/rules/publish.mdc`.

## Language

Identifiers, rustdoc, README, and the Cargo description are English. Issue and PR bodies and commit subjects may be Japanese. The canonical spec is `docs/design.ja.md`; the English `docs/design.md` is a roadmap row before v0.1.0.

## Local files

Remove what you made when you no longer need it (`.cursor/rules/scratch.mdc`). `.dev/` is ignored and for local notes.
