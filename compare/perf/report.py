"""Makes the summary of a comparison run: `docs/performance.md` and the
figures of memory over time under `docs/performance/`.

    uv run --directory compare python -X utf8 perf/report.py

Reads `out/results.jsonl`, `out/series/`, and `out/environment.json`, which
`run.py` wrote. The summary is what is committed; `out/` is not.
"""

import json
import statistics
from collections import defaultdict
from pathlib import Path

import cases as case_table

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
DOCS = HERE.parent.parent / "docs"
FIGURES = DOCS / "performance"

LIBRARIES = ("interpx", "interpx_owned", "scipy", "scirs2", "ndarray_interp", "interpn")
# The libraries interpx is compared with: not itself built another way.
OTHERS = ("scipy", "scirs2", "ndarray_interp", "interpn")
RUST = ("scirs2", "ndarray_interp", "interpn")
TITLES = {
    "interpx": "interpx",
    "interpx_owned": "interpx, `from_owned`",
    "scipy": "SciPy",
    "scirs2": "scirs2-interpolate",
    "ndarray_interp": "ndarray-interp",
    "interpn": "interpn",
}
# interpx counts as slower when its median is this much over another
# library's, and its fastest repeat is over that library's fastest too.
LOSS = 1.10
# Values further than this from SciPy's (as a share of the largest value) are
# those of another interpolant, not of another rounding.
SAME = 1e-6

FAMILIES = (
    ("oned", "One variable", "Data of `n` points, evaluated at `4 n` points: in order (`sorted`) and in no order (`shuffled`)."),
    ("ppoly", "`PPoly`", "A cubic of 55 intervals, evaluated at points in no order."),
    ("bspline", "`BSpline`", "A cubic B-spline of 55 intervals, evaluated at points in no order."),
    ("curve", "`BSplineCurve`", "A cubic curve through `n` points of the plane, evaluated at `4 n` parameters in order."),
    ("grid", "Regular grid", "Points in no order inside the grid. `maxN`: the first axis has `N` points and each axis half the one before (SciPy's suite); `eachN`: every axis has `N` points; `flipped`: every axis given from its high end."),
    ("nearest_nd", "`NearestNDInterpolator`", "`n` points of the unit square, evaluated at the pixels of a square image."),
    ("rbf", "`RBFInterpolator`", "`n` points of a square, evaluated at 1000 points in no order (`scattered`) or at the 10 000 pixels of an image (`image`). `global`: all points; `k50`: the 50 nearest."),
)
# The cases whose memory over time is reported and drawn.
DRAWN = case_table.MEMORY
COLORS = {"loaded": "#888888", "built": "#1f77b4", "evaluated": "#d62728", "timed": "#888888"}


def load() -> dict:
    """The runs, by case and library: the rows of all rounds."""
    runs = defaultdict(lambda: defaultdict(list))
    for line in (OUT / "results.jsonl").read_text(encoding="utf-8").splitlines():
        row = json.loads(line)
        runs[row["case"]][row["library"]].append(row)
    return runs


def measured(rows: list[dict]) -> list[dict]:
    return [row for row in rows if "build" in row]


def times(rows: list[dict], what: str) -> dict | None:
    rows = measured(rows)
    if not rows:
        return None
    return {
        "median": statistics.median(row[what]["median"] for row in rows),
        "min": min(row[what]["min"] for row in rows),
        "max": max(row[what]["max"] for row in rows),
    }


def seconds(value: float) -> str:
    for unit, scale in (("s", 1.0), ("ms", 1e-3), ("µs", 1e-6), ("ns", 1e-9)):
        if value >= scale:
            return f"{value / scale:.3g} {unit}"
    return "0"


def megabytes(value: float | None) -> str:
    return "" if value is None else f"{value / 1e6:.1f}"


def losses(runs: dict) -> list[dict]:
    found = []
    for case, by_library in runs.items():
        for what in ("build", "eval"):
            ours = times(by_library.get("interpx", []), what)
            if ours is None:
                continue
            for library in OTHERS:
                theirs = times(by_library.get(library, []), what)
                if theirs is None or theirs["median"] <= 0.0:
                    continue
                ratio = ours["median"] / theirs["median"]
                if ratio > LOSS and ours["min"] > theirs["min"]:
                    apart = max((row["apart"] or 0.0) for row in measured(by_library[library]))
                    found.append({"case": case, "what": what, "library": library, "ratio": ratio, "ours": ours, "theirs": theirs, "apart": apart})
    return sorted(found, key=lambda loss: -loss["ratio"])


def method_of(case: str) -> str:
    """The family of a case, with its method where the family has several."""
    parts = case.split("/")
    return "/".join(parts[:2]) if parts[0] in ("oned", "grid", "rbf") else parts[0]


def grouped(lost: list[dict]) -> list[str]:
    """The losses by method, operation, and library: how many, and how large."""
    groups = defaultdict(list)
    for loss in lost:
        same = loss["library"] == "scipy" or loss["apart"] < SAME
        groups[(method_of(loss["case"]), loss["what"], loss["library"], same)].append(loss)
    lines = [
        "| method | operation | against | same values | cases | slower by | the largest |",
        "| :--- | :--- | :--- | :--- | ---: | ---: | :--- |",
    ]
    for (method, what, library, same), found in sorted(groups.items(), key=lambda item: -max(loss["ratio"] for loss in item[1])):
        ratios = [loss["ratio"] for loss in found]
        worst = max(found, key=lambda loss: loss["ratio"])
        span = f"{min(ratios):.3g}x to {max(ratios):.3g}x" if len(found) > 1 else f"{ratios[0]:.3g}x"
        lines.append(
            f"| `{method}` | {'construction' if what == 'build' else 'evaluation'} | {TITLES[library]} | {'yes' if same else 'no'} | {len(found)} | {span} | `{worst['case']}`: {seconds(worst['ours']['median'])} against {seconds(worst['theirs']['median'])} |"
        )
    return lines


def tally(runs: dict) -> list[str]:
    """For each family and operation: in how many cases interpx is the fastest."""
    lines = ["| family | operation | cases | interpx fastest | slower than SciPy | slower than a Rust crate |", "| :--- | :--- | ---: | ---: | ---: | ---: |"]
    lost = {(loss["case"], loss["what"], loss["library"]) for loss in losses(runs)}
    for family, title, _ in FAMILIES:
        names = [case for case, by in runs.items() if case.split("/")[0] == family and times(by.get("interpx", []), "eval")]
        for what in ("eval", "build"):
            to_scipy = sum((case, what, "scipy") in lost for case in names)
            to_rust = sum(any((case, what, lib) in lost for lib in RUST) for case in names)
            fastest = sum(not any((case, what, lib) in lost for lib in OTHERS) for case in names)
            lines.append(f"| {title} | {'evaluation' if what == 'eval' else 'construction'} | {len(names)} | {fastest} | {to_scipy} | {to_rust} |")
    return lines


def table(runs: dict, family: str, what: str) -> list[str]:
    names = [case for case, by_library in runs.items() if next(iter(by_library.values()))[0]["family"] == family]
    present = [lib for lib in LIBRARIES if any(times(runs[case].get(lib, []), what) for case in names)]
    lines = ["| case | " + " | ".join(TITLES[lib] for lib in present) + " | interpx against the fastest other |", "| :--- | " + " | ".join("---:" for _ in present) + " | ---: |"]
    for case in names:
        cells, others = [], []
        for lib in present:
            got = times(runs[case].get(lib, []), what)
            if got is None:
                cells.append("")
                continue
            cells.append(seconds(got["median"]))
            if lib in OTHERS:
                others.append(got["median"])
        ours = times(runs[case].get("interpx", []), what)
        if ours and others and min(others) > 0.0:
            ratio = min(others) / ours["median"]
            verdict = f"{ratio:.2g}x faster" if ratio >= 1.0 else f"**{1.0 / ratio:.2g}x slower**"
        else:
            verdict = ""
        lines.append(f"| `{case.split('/', 1)[1]}` | " + " | ".join(cells) + f" | {verdict} |")
    return lines


def memory_table(runs: dict, cases: list[str]) -> list[str]:
    lines = [
        "| case | library | data read | while building | built | while evaluating | highest |",
        "| :--- | :--- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for case in cases:
        for lib in LIBRARIES:
            rows = measured(runs.get(case, {}).get(lib, []))
            if not rows:
                continue

            def middle(key: str) -> float | None:
                values = [row["rss"][key] for row in rows if row["rss"][key] is not None]
                return statistics.median(values) if values else None

            base = middle("loaded")

            def over(key: str) -> str:
                value = middle(key)
                return "" if value is None or base is None else f"+{max(value - base, 0.0) / 1e6:.1f}"

            lines.append(f"| `{case}` | {TITLES[lib]} | {megabytes(base)} | {over('build_peak')} | {over('built')} | {over('eval_peak')} | {megabytes(middle('overall_peak'))} |")
    return lines


def figure(runs: dict, case: str) -> str | None:
    """Draws the memory of each library over the time of its run, and returns the file's name."""
    libraries = [lib for lib in LIBRARIES if measured(runs.get(case, {}).get(lib, []))]
    if not libraries:
        return None
    width, panel, left, top = 760, 120, 70, 24
    parts = []
    for index, lib in enumerate(libraries):
        row = measured(runs[case][lib])[0]
        path = OUT / "series" / f"{case.replace('/', '__')}.{lib}.{row['round']}.csv"
        samples = [tuple(map(int, line.split(","))) for line in path.read_text(encoding="utf-8").splitlines()]
        end = row["marks_ns"]["timed"]
        samples = [(at, rss) for at, rss in samples if at <= end] or samples
        high = max(rss for _, rss in samples)
        y0 = top + index * (panel + 34)
        plot_w, plot_h = width - left - 20, panel

        def x(at: int) -> float:
            return left + plot_w * at / end

        def y(rss: int) -> float:
            return y0 + plot_h - plot_h * rss / high

        points = " ".join(f"{x(at):.1f},{y(rss):.1f}" for at, rss in samples)
        parts.append(f'<rect x="{left}" y="{y0}" width="{plot_w}" height="{plot_h}" fill="none" stroke="#cccccc"/>')
        parts.append(f'<polyline points="{points}" fill="none" stroke="#222222" stroke-width="1.2"/>')
        for name in ("loaded", "built", "evaluated"):
            at = row["marks_ns"][name]
            parts.append(f'<line x1="{x(at):.1f}" y1="{y0}" x2="{x(at):.1f}" y2="{y0 + plot_h}" stroke="{COLORS[name]}" stroke-dasharray="4 3"/>')
            parts.append(f'<text x="{x(at) + 3:.1f}" y="{y0 + 11}" font-size="10" fill="{COLORS[name]}">{name}</text>')
        parts.append(f'<text x="{left}" y="{y0 - 5}" font-size="12" font-weight="bold">{TITLES[lib]}</text>')
        parts.append(f'<text x="{left - 6}" y="{y0 + 10}" font-size="10" text-anchor="end">{high / 1e6:.0f} MB</text>')
        parts.append(f'<text x="{left - 6}" y="{y0 + plot_h}" font-size="10" text-anchor="end">0</text>')
        parts.append(f'<text x="{left + plot_w}" y="{y0 + plot_h + 12}" font-size="10" text-anchor="end">{end / 1e6:.0f} ms</text>')
    height = top + len(libraries) * (panel + 34)
    svg = f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" font-family="sans-serif">\n<rect width="100%" height="100%" fill="white"/>\n' + "\n".join(parts) + "\n</svg>\n"
    name = case.replace("/", "-") + ".svg"
    FIGURES.mkdir(parents=True, exist_ok=True)
    (FIGURES / name).write_text(svg, encoding="utf-8")
    return name


def main() -> None:
    runs = load()
    env = json.loads((OUT / "environment.json").read_text(encoding="utf-8"))
    rounds = max(len(rows) for by_library in runs.values() for rows in by_library.values())
    lost = losses(runs)
    failed = sorted({(case, lib, row["failed"]) for case, by in runs.items() for lib, rows in by.items() for row in rows if "failed" in row})

    doc = [
        "# Performance: interpx against SciPy and other Rust crates",
        "",
        "<!-- Written by compare/perf/report.py from a run of compare/perf/run.py. Do not edit by hand: run it again. -->",
        "",
        f"Measured on {env['date']}, at commit `{env['commit'][:10]}`{' with uncommitted changes' if env['dirty'] else ''}.",
        "",
        "## How it was measured",
        "",
        f"- Machine: {env['machine']}, {env['cores']} logical cores, {env['memory_bytes'] / 2**30:.0f} GB. {env['os']}.",
        f"- {env['rustc']}; Python {env['python']}, NumPy {env['numpy']}, SciPy {env['scipy']}.",
        "- The other Rust crates, at the versions in `compare/perf/runners/*/Cargo.toml`: `scirs2-interpolate` (the widest of the Rust crates: one variable, grids, RBF), `ndarray-interp` (the most used for one variable; also bilinear), `interpn` (grids).",
        f"- Every library ran each case in a process of its own, pinned to one core, on one thread, on the same data. The libraries took turns case by case, {rounds} rounds.",
        "- `interpx, from_owned` is interpx built from data that is handed over (`RegularGridInterpolator::from_owned`) and so not copied. Its column is there for the grids, next to the construction from borrowed data (`new`), which is the one the comparison counts. Each of its samples is given data made just before it, which is not in the nearest cache as the data of the other constructions is: on small grids it reads slower than `new` for that reason.",
        "- A time is that of one call: construction, or the evaluation of all the queries of the case (`eval_batch_into` for interpx, the call with an array for the others). It is the median of the repeats of a process, and then the median of the rounds.",
        f"- interpx is called slower where its median is more than {LOSS:.2f} times another library's and its fastest repeat is slower than that library's fastest.",
        "- Memory is the resident set size of the process, sampled every millisecond from outside it.",
        f"- Command: `{env['command']}`.",
        "",
        "## Where interpx is slower",
        "",
    ]
    doc += tally(runs) + [""]
    if lost:
        doc += [
            f"{len(lost)} comparisons (a case, an operation, and a library), by method. \"Same values\" says whether that library's values are SciPy's up to rounding: where they are not, it computes another interpolant under the same name, and its time is not that of the same work.",
            "",
        ]
        doc += grouped(lost)
        doc += [
            "",
            "<details><summary>Every comparison where interpx is slower</summary>",
            "",
            "\"Values apart\" is the distance of that library's values from SciPy's, as a share of the largest value.",
            "",
            "| case | operation | against | interpx | the other | slower by | values apart |",
            "| :--- | :--- | :--- | ---: | ---: | ---: | ---: |",
        ]
        for loss in lost:
            apart = "" if loss["library"] == "scipy" else f"{loss['apart']:.0e}"
            doc.append(f"| `{loss['case']}` | {'construction' if loss['what'] == 'build' else 'evaluation'} | {TITLES[loss['library']]} | {seconds(loss['ours']['median'])} | {seconds(loss['theirs']['median'])} | {loss['ratio']:.2f}x | {apart} |")
        doc += ["", "</details>"]
    else:
        doc.append("Nowhere.")
    doc.append("")
    if failed:
        doc += ["## Runs that failed", "", "| case | library | how |", "| :--- | :--- | :--- |"]
        doc += [f"| `{case}` | {TITLES[lib]} | {how} |" for case, lib, how in failed]
        doc.append("")

    for family, title, note in FAMILIES:
        doc += [f"## {title}", "", note, "", "### Evaluation", ""] + table(runs, family, "eval") + ["", "### Construction", ""] + table(runs, family, "build") + [""]

    doc += [
        "## Memory",
        "",
        "Megabytes. \"Data read\" is the process with the data of the case in memory, before anything is built. The columns after it are what was added to that: the highest while building, with the interpolator built, the highest while evaluating. \"Highest\" is the highest of the whole run, repeats included.",
        "",
    ]
    large = [case for case in runs if case in DRAWN]
    doc += memory_table(runs, large) + [""]
    doc += ["### Over time", "", "The resident set size of each library's process over its run, with the ends of the phases. The time after `evaluated` is the timed repeats.", ""]
    for case in large:
        name = figure(runs, case)
        if name:
            doc += [f"`{case}`", "", f"![memory over time, {case}](performance/{name})", ""]

    apart = defaultdict(float)
    for case, by_library in runs.items():
        for lib, rows in by_library.items():
            for row in measured(rows):
                if row["apart"] is not None:
                    apart[(row["family"], lib)] = max(apart[(row["family"], lib)], row["apart"])
    doc += [
        "## Do the libraries compute the same thing",
        "",
        f"The largest distance of a library's values from SciPy's, over the cases of a family, as a share of the largest value. Under {SAME:.0e} the two compute the same interpolant.",
        "",
        "| family | " + " | ".join(TITLES[lib] for lib in LIBRARIES if lib != "scipy") + " |",
        "| :--- | " + " | ".join("---:" for lib in LIBRARIES if lib != "scipy") + " |",
    ]
    for family, title, _ in FAMILIES:
        cells = [f"{apart[(family, lib)]:.0e}" if (family, lib) in apart else "" for lib in LIBRARIES if lib != "scipy"]
        doc.append(f"| {title} | " + " | ".join(cells) + " |")
    doc.append("")
    (DOCS / "performance.md").write_text("\n".join(doc), encoding="utf-8", newline="\n")
    print(f"{len(runs)} cases, {len(lost)} comparisons where interpx is slower, {len(failed)} failed runs")


if __name__ == "__main__":
    main()
