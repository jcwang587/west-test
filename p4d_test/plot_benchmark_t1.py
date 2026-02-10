#!/usr/bin/env python3
"""Plot p4d_test benchmarks for the 1-thread case (t1_g*)."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator


WALL_RE_TEMPLATE = r"{label}\s*:\s*.*?CPU\s+(.*?)\s+WALL"
DURATION_RE = re.compile(r"([0-9]+(?:\.[0-9]+)?)\s*([hms])")


def parse_duration_seconds(text: str) -> float:
    """Convert a duration like '1m 7.25s' or '44.42s' to seconds."""
    total = 0.0
    for value, unit in DURATION_RE.findall(text):
        value_f = float(value)
        if unit == "h":
            total += value_f * 3600.0
        elif unit == "m":
            total += value_f * 60.0
        elif unit == "s":
            total += value_f
    if total <= 0.0:
        raise ValueError(f"Could not parse duration: '{text}'")
    return total


def extract_wall_seconds(log_path: Path, label: str) -> float:
    """Extract wall time (seconds) from the last timing line for a label."""
    content = log_path.read_text(encoding="utf-8", errors="ignore")
    pattern = re.compile(WALL_RE_TEMPLATE.format(label=re.escape(label)))
    matches = pattern.findall(content)
    if not matches:
        raise ValueError(f"Could not find '{label}' timing in {log_path}")
    return parse_duration_seconds(matches[-1].strip())


def collect_t1_benchmark_data(
    base_dir: Path,
) -> tuple[list[int], list[float], list[float]]:
    """Return g-values, PWSCF wall times, and WSTAT wall times for t1_g* tests."""
    entries: list[tuple[int, float, float]] = []
    for test_dir in sorted(base_dir.glob("t1_g*")):
        if not test_dir.is_dir():
            continue
        g_text = test_dir.name.split("_g")[-1]
        if not g_text.isdigit():
            continue
        g_value = int(g_text)

        pw_log = test_dir / "pw.out"
        wstat_log = test_dir / "wstat.out"
        if not pw_log.exists() or not wstat_log.exists():
            print(
                f"Skipping {test_dir.name}: missing pw.out or wstat.out",
                file=sys.stderr,
            )
            continue

        pw_wall = extract_wall_seconds(pw_log, "PWSCF")
        wstat_wall = extract_wall_seconds(wstat_log, "WSTAT")
        entries.append((g_value, pw_wall, wstat_wall))

    if not entries:
        raise RuntimeError("No valid t1_g* benchmark outputs found under p4d_test.")

    entries.sort(key=lambda item: item[0])
    g_vals = [item[0] for item in entries]
    pw_vals = [item[1] for item in entries]
    wstat_vals = [item[2] for item in entries]
    return g_vals, pw_vals, wstat_vals


def make_plot(
    g_vals: list[int], pw_vals: list[float], wstat_vals: list[float], out_path: Path
) -> None:
    """Create and save the benchmark plot."""
    x_vals = g_vals
    all_vals = pw_vals + wstat_vals
    g_ref = g_vals[0]
    ideal_ref = wstat_vals[0]
    ideal_vals = [ideal_ref * (g_ref / g) for g in g_vals]
    pw_speedup = [pw_vals[0] / v for v in pw_vals]
    wstat_speedup = [wstat_vals[0] / v for v in wstat_vals]

    fig, ax = plt.subplots(figsize=(5.8, 5.8), facecolor="white")
    ax.set_facecolor("#f8fafc")

    ax.plot(
        x_vals,
        pw_vals,
        color="#1f77b4",
        marker="o",
        markersize=6.5,
        markeredgecolor="white",
        markeredgewidth=1.1,
        linewidth=2.2,
        label="PWSCF",
        zorder=4,
    )
    ax.plot(
        x_vals,
        wstat_vals,
        color="#e4572e",
        marker="s",
        markersize=6.8,
        markeredgecolor="white",
        markeredgewidth=1.1,
        linewidth=2.2,
        label="WSTAT",
        zorder=4,
    )
    ax.plot(
        x_vals,
        ideal_vals,
        linestyle="--",
        color="#222222",
        linewidth=1.5,
        dashes=(5, 3),
        alpha=0.9,
        label="_nolegend_",
        zorder=3,
    )

    for x, y, s in zip(x_vals[1:], pw_vals[1:], pw_speedup[1:]):
        ax.annotate(
            f"{s:.2f}x",
            xy=(x, y),
            xytext=(0, 9),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=10.5,
            color="#1f77b4",
            fontweight="semibold",
            zorder=5,
        )

    for x, y, s in zip(x_vals[1:], wstat_vals[1:], wstat_speedup[1:]):
        ax.annotate(
            f"{s:.2f}x",
            xy=(x, y),
            xytext=(0, 9),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=10.5,
            color="#e4572e",
            fontweight="semibold",
            zorder=5,
        )

    ax.set_xlabel("Number of GPUs", fontsize=14)
    ax.set_ylabel("Wall time (s)", fontsize=14)
    ax.set_xscale("log", base=2)
    ax.set_xlim(1, 8)
    ax.set_xticks(g_vals, labels=[str(v) for v in g_vals])
    ax.set_yscale("log")
    ax.set_yticks([10, 100, 1000], labels=["10", "100", "1000"])
    ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=(2, 3, 4, 5, 6, 7, 8, 9)))
    ax.set_ylim(10, max(1000, max(all_vals) * 1.1))
    ax.set_box_aspect(1)
    ax.grid(which="major", linestyle="-", linewidth=0.8, color="#c6ccd6", alpha=0.65)
    ax.grid(which="minor", linestyle="-", linewidth=0.55, color="#dfe4ec", alpha=0.55)

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.15)
    ax.tick_params(axis="both", colors="black", width=0.95, labelsize=12)

    leg = ax.legend(
        loc="upper right",
        fontsize=12,
        frameon=True,
        framealpha=0.95,
        borderpad=0.45,
        handlelength=1.8,
    )
    leg.get_frame().set_edgecolor("#c6ccd6")
    leg.get_frame().set_linewidth(0.9)
    leg.get_frame().set_facecolor("white")

    fig.tight_layout()
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot p4d_test benchmark results for the 1-thread tests (t1_g*)."
    )
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory containing t1_g* folders (default: script directory).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "benchmark_t1.png",
        help="Output PNG path (default: p4d_test/benchmark_t1.png).",
    )
    args = parser.parse_args()

    g_vals, pw_vals, wstat_vals = collect_t1_benchmark_data(args.base_dir)
    make_plot(g_vals, pw_vals, wstat_vals, args.output)

    print("Plotted 1-thread benchmark:")
    print(f"  g values:     {g_vals}")
    print(f"  PWSCF (s):    {[round(v, 2) for v in pw_vals]}")
    print(f"  WSTAT (s):    {[round(v, 2) for v in wstat_vals]}")
    print(f"  Output image: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
