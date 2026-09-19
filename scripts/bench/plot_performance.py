#!/usr/bin/env python3
"""Build docs/PERFORMANCE.md and docs/figures/*.svg from docs/bench/*.csv."""

from __future__ import annotations

import csv
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BENCH_DIR = ROOT / "docs" / "bench"
FIG_DIR = ROOT / "docs" / "figures"
OUT_MD = ROOT / "docs" / "PERFORMANCE.md"

INT_FIELDS = (
    "nx",
    "n_elements",
    "n_nodes",
    "n_elements_out",
    "n_nodes_out",
    "mpi_ranks",
    "slurm_nodes",
    "omp_threads",
    "repeat",
    "bytes",
)
FLOAT_FIELDS = ("time_s", "melems_s", "gib_s")

KERNEL_ORDER = ["generate", "write", "read", "promote", "refine"]
ELEM_ORDER = ["HEX8", "TET4"]


def parse_row(raw: dict) -> dict | None:
    if not raw.get("kernel"):
        return None
    row = dict(raw)
    for key in INT_FIELDS:
        row[key] = int(float(row[key])) if row.get(key) not in (None, "") else 0
    for key in FLOAT_FIELDS:
        row[key] = float(row[key]) if row.get(key) not in (None, "") else 0.0
    return row


def load_rows(extra: list[Path]) -> list[dict]:
    paths: list[Path] = []
    if BENCH_DIR.exists():
        paths.extend(sorted(BENCH_DIR.glob("*.csv")))
    for path in extra:
        path = path.resolve()
        if path not in paths:
            paths.append(path)
    rows: list[dict] = []
    for path in paths:
        if not path.is_file():
            continue
        with path.open(newline="") as fh:
            for raw in csv.DictReader(fh):
                row = parse_row(raw)
                if row:
                    rows.append(row)
    return rows


def group_median(rows: list[dict]) -> list[dict]:
    buckets: dict[tuple, list[dict]] = defaultdict(list)
    for row in rows:
        key = (row["kernel"], row["elem_type"], row["nx"], row["mpi_ranks"])
        buckets[key].append(row)
    out = []
    for _, items in sorted(buckets.items()):
        base = dict(items[0])
        base["time_s"] = statistics.median(item["time_s"] for item in items)
        base["melems_s"] = statistics.median(item["melems_s"] for item in items)
        base["gib_s"] = statistics.median(item["gib_s"] for item in items)
        base["n_repeats"] = len(items)
        out.append(base)
    return out


def efficiency(series: list[dict]) -> dict[int, float]:
    if not series:
        return {}
    by_rank = {row["mpi_ranks"]: row["time_s"] for row in series}
    p0 = min(by_rank)
    t0 = by_rank[p0]
    out = {}
    for p, t in by_rank.items():
        denom = (p / p0) * t
        out[p] = t0 / denom if denom > 0 else 0.0
    return out


def _log_span(values: list[float], pad: float = 0.08) -> tuple[float, float]:
    positives = [v for v in values if v > 0]
    if not positives:
        return -1.0, 1.0
    lo = math.log10(min(positives))
    hi = math.log10(max(positives))
    if abs(hi - lo) < 1e-12:
        lo -= 0.5
        hi += 0.5
    span = hi - lo
    return lo - pad * span, hi + pad * span


def _lin_span(values: list[float], pad: float = 0.08) -> tuple[float, float]:
    if not values:
        return 0.0, 1.0
    lo = min(values)
    hi = max(values)
    if abs(hi - lo) < 1e-12:
        if hi == 0:
            return -1.0, 1.0
        return lo * 0.9 if lo > 0 else lo * 1.1, hi * 1.1 if hi > 0 else hi * 0.9
    span = hi - lo
    lo2 = lo - pad * span
    if min(values) >= 0:
        lo2 = max(0.0, lo2)
    return lo2, hi + pad * span


def _ticks_log(lo: float, hi: float) -> list[float]:
    start = int(math.floor(lo))
    stop = int(math.ceil(hi))
    ticks = [10.0**e for e in range(start, stop + 1) if lo - 1e-9 <= e <= hi + 1e-9]
    return ticks or [10.0 ** ((lo + hi) / 2.0)]


def _ticks_lin(lo: float, hi: float, n: int = 5) -> list[float]:
    if hi <= lo:
        return [lo]
    return [lo + (hi - lo) * i / (n - 1) for i in range(n)]


def _fmt(v: float) -> str:
    av = abs(v)
    if av == 0:
        return "0"
    if av >= 100 or av < 0.01:
        return f"{v:.2g}"
    return f"{v:.3g}"


def write_svg(
    series: list[dict],
    ykey: str,
    ylabel: str,
    path: Path,
    logx: bool,
    logy: bool,
    ideal: bool,
) -> None:
    series = sorted(series, key=lambda r: r["mpi_ranks"])
    xs = [float(r["mpi_ranks"]) for r in series]
    ys = [float(r[ykey]) for r in series]
    w, h = 640.0, 420.0
    l, r, t, b = 72.0, 24.0, 24.0, 56.0
    pw, ph = w - l - r, h - t - b

    x_log = [math.log2(x) if x > 0 else 0.0 for x in xs] if logx else xs
    if logx:
        x_lo, x_hi = min(x_log) - 0.4, max(x_log) + 0.4
        if x_hi <= x_lo:
            x_lo, x_hi = x_lo - 0.5, x_hi + 0.5
    else:
        x_lo, x_hi = _lin_span(xs)
    if logy:
        y_lo, y_hi = _log_span(ys)
    else:
        y_lo, y_hi = _lin_span(ys)

    def xpix(x: float) -> float:
        u = math.log2(x) if logx and x > 0 else x
        src_lo, src_hi = (x_lo, x_hi)
        return l + (u - src_lo) / (src_hi - src_lo) * pw

    def ypix(y: float) -> float:
        u = math.log10(y) if logy and y > 0 else y
        src_lo, src_hi = (y_lo, y_hi)
        return t + ph - (u - src_lo) / (src_hi - src_lo) * ph

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{int(w)}" height="{int(h)}" viewBox="0 0 {int(w)} {int(h)}">',
        "<style>text{font-family:Helvetica,Arial,sans-serif;font-size:12px;fill:#222} .muted{fill:#666;font-size:11px}</style>",
        f'<rect x="0" y="0" width="{int(w)}" height="{int(h)}" fill="#fff"/>',
        f'<rect x="{l}" y="{t}" width="{pw}" height="{ph}" fill="#fafafa" stroke="#ccc"/>',
    ]

    if logx:
        xticks = sorted({2**e for e in range(int(math.floor(min(x_log))), int(math.ceil(max(x_log))) + 1)})
        xticks = [x for x in xticks if x >= 1]
        if not xticks:
            xticks = xs
    else:
        xticks = _ticks_lin(x_lo, x_hi)
    yticks = _ticks_log(y_lo, y_hi) if logy else _ticks_lin(y_lo, y_hi)

    for xv in xticks:
        if xv <= 0 and logx:
            continue
        xp = xpix(xv)
        parts.append(f'<line x1="{xp:.1f}" y1="{t}" x2="{xp:.1f}" y2="{t+ph}" stroke="#eee"/>')
        parts.append(f'<text class="muted" x="{xp:.1f}" y="{t+ph+16}" text-anchor="middle">{_fmt(xv)}</text>')
    for yv in yticks:
        if yv <= 0 and logy:
            continue
        yp = ypix(yv)
        parts.append(f'<line x1="{l}" y1="{yp:.1f}" x2="{l+pw}" y2="{yp:.1f}" stroke="#eee"/>')
        parts.append(f'<text class="muted" x="{l-8}" y="{yp+4:.1f}" text-anchor="end">{_fmt(yv)}</text>')

    if ideal and ykey == "time_s" and xs and ys and xs[0] > 0 and ys[0] > 0:
        iys = [ys[0] * xs[0] / x for x in xs]
        d = " ".join(f"{xpix(x):.1f},{ypix(y):.1f}" for x, y in zip(xs, iys) if y > 0)
        if d:
            parts.append(f'<polyline fill="none" stroke="#999" stroke-dasharray="6 4" stroke-width="1.5" points="{d}"/>')
            parts.append(f'<text class="muted" x="{l+8}" y="{t+16}">dashed: ideal</text>')

    pts = [(x, y) for x, y in zip(xs, ys) if y > 0 or not logy]
    if len(pts) >= 2:
        d = " ".join(f"{xpix(x):.1f},{ypix(y):.1f}" for x, y in pts)
        parts.append(f'<polyline fill="none" stroke="#0b6e99" stroke-width="2" points="{d}"/>')
    for x, y in pts:
        parts.append(f'<circle cx="{xpix(x):.1f}" cy="{ypix(y):.1f}" r="4" fill="#0b6e99"/>')

    parts.append(f'<text x="{l+pw/2}" y="{h-12}" text-anchor="middle">MPI ranks (1 rank / socket, 72 OMP threads)</text>')
    parts.append(
        f'<text x="16" y="{t+ph/2}" text-anchor="middle" transform="rotate(-90 16 {t+ph/2})">{ylabel}</text>'
    )
    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n")


def md_table(series: list[dict], eff: dict[int, float]) -> str:
    if not series:
        return "_No rows yet._\n"
    lines = [
        "| ranks | nodes | N | elems | time [s] | Melem/s | GiB/s | efficiency |",
        "|------:|------:|--:|------:|---------:|--------:|------:|-----------:|",
    ]
    for row in sorted(series, key=row_key):
        p = row["mpi_ranks"]
        e = eff.get(p, float("nan"))
        lines.append(
            "| {ranks} | {nodes} | {nx} | {ne} | {t:.4g} | {m:.4g} | {g:.4g} | {e:.3f} |".format(
                ranks=p,
                nodes=row["slurm_nodes"] or "—",
                nx=row["nx"],
                ne=row["n_elements"],
                t=row["time_s"],
                m=row["melems_s"],
                g=row["gib_s"],
                e=e,
            )
        )
    return "\n".join(lines) + "\n"


def row_key(row: dict) -> tuple:
    return (row["mpi_ranks"], row["nx"])


def write_md(groups: dict[tuple[str, str, int], list[dict]], figures: dict[str, Path]) -> None:
    lines = [
        "# PERFORMANCE.md",
        "",
        "Hybrid strong scaling on CSCS Alps (GH200): **1 MPI rank per socket**, 72 OpenMP threads.",
        "Points: 1 socket (1×72), 1 node (4×72), 2 / 4 / 8 / 16 nodes (8 / 16 / 32 / 64 ranks).",
        "",
        "## Method",
        "",
        "- Driver: `smesh_bench` (`src/benchmark/drivers/smesh_bench.exe.cpp`). Wall time is `max` over ranks after `MPI_Barrier`.",
        "- Each kernel is repeated `SMESH_BENCH_REPEAT` times (default 3). Tables use the **median**.",
        "- Default sizes: HEX8 `N=512` (`N³` hexes); TET4 `SMESH_BENCH_TET_N=256` (`6 N³` tets). Override with `SMESH_BENCH_N`.",
        "- Throughput: `n_elements / time` (Melem/s). IO also reports GiB/s from `n_nodes * sdim * sizeof(geom_t) + n_elements * nxe * sizeof(idx_t)`.",
        "- Efficiency: `T_1 / ((P / P_1) * T_P)` relative to the smallest rank count in the CSV (the 1-socket point when present).",
        "- Alps build: `-DSMESH_ENABLE_MPI=ON -DSMESH_ENABLE_OPENMP=ON -DCMAKE_BUILD_TYPE=Release`.",
        "- Submit: `scripts/slurm/alps/submit_strong_scaling.sh`. CSV files land in `docs/bench/`. Regenerate this page with `scripts/bench/plot_performance.py`.",
        "- Rows with `nodes = —` (`slurm_nodes=0`) are not Alps jobs (local smoke or interactive runs).",
        "",
        "`create_cube` uses serial fill when `comm_size == 1` and distributed create when `comm_size > 1`. "
        "The 1-socket generate point is a different kernel than 4+ ranks. Refine and promote use the production MPI paths.",
        "",
        "OpenMP helps SS refine graphs. Distributed generate is MPI-decomposed and is not itself an OpenMP loop.",
        "",
    ]

    if not groups:
        lines += [
            "## Results",
            "",
            "No CSV rows in `docs/bench/` yet. Run the Alps jobs, copy `docs/bench/*.csv` into the repo, and rerun `scripts/bench/plot_performance.py`.",
            "",
        ]
        OUT_MD.write_text("\n".join(lines))
        return

    ordered: list[tuple[str, str, int]] = []
    for kernel in KERNEL_ORDER:
        for elem in ELEM_ORDER:
            nxs = sorted({k[2] for k in groups if k[0] == kernel and k[1] == elem})
            for nx in nxs:
                ordered.append((kernel, elem, nx))
    seen = set(ordered)
    for key in ordered:
        series = groups.get(key)
        if not series:
            continue
        kernel, elem, nx = key
        title = f"{kernel} {elem} N={nx}"
        note = " (local / non-Slurm rows present)" if any(row["slurm_nodes"] == 0 for row in series) else ""
        lines += [f"## {title}{note}", ""]
        lines.append(md_table(series, efficiency(series)))
        fig_time = figures.get(f"{kernel}_{elem}_{nx}_time")
        fig_thr = figures.get(f"{kernel}_{elem}_{nx}_thr")
        fig_io = figures.get(f"{kernel}_{elem}_{nx}_io")
        if fig_time:
            lines += ["", f"![{title} time]({fig_time.relative_to(OUT_MD.parent)})"]
        if fig_thr:
            lines += ["", f"![{title} throughput]({fig_thr.relative_to(OUT_MD.parent)})"]
        if fig_io:
            lines += ["", f"![{title} GiB/s]({fig_io.relative_to(OUT_MD.parent)})"]
        lines.append("")

    for key in groups:
        if key in seen:
            continue
        kernel, elem, nx = key
        lines += [f"## {kernel} {elem} N={nx}", "", md_table(groups[key], efficiency(groups[key])), ""]

    OUT_MD.write_text("\n".join(lines).rstrip() + "\n")


def main(argv: list[str]) -> int:
    extra = [Path(a) for a in argv[1:]]
    rows = load_rows(extra)
    med = group_median(rows)
    groups: dict[tuple[str, str, int], list[dict]] = defaultdict(list)
    for row in med:
        groups[(row["kernel"], row["elem_type"], row["nx"])].append(row)

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    BENCH_DIR.mkdir(parents=True, exist_ok=True)
    figures: dict[str, Path] = {}
    for (kernel, elem, nx), series in groups.items():
        tpath = FIG_DIR / f"{kernel}_{elem}_{nx}_time.svg"
        mpath = FIG_DIR / f"{kernel}_{elem}_{nx}_thr.svg"
        write_svg(series, "time_s", "time [s]", tpath, True, True, True)
        write_svg(series, "melems_s", "Melem/s", mpath, True, False, False)
        figures[f"{kernel}_{elem}_{nx}_time"] = tpath
        figures[f"{kernel}_{elem}_{nx}_thr"] = mpath
        if kernel in ("read", "write"):
            ipath = FIG_DIR / f"{kernel}_{elem}_{nx}_io.svg"
            write_svg(series, "gib_s", "GiB/s", ipath, True, False, False)
            figures[f"{kernel}_{elem}_{nx}_io"] = ipath

    write_md(groups, figures)
    print(f"wrote {OUT_MD} ({sum(len(v) for v in groups.values())} median series)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
