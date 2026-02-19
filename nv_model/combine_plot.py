#!/usr/bin/env python3
"""
Combine GPU utilization and memory plots from p4d_nv_1g_1o and p4d_nv_8g_8o
into one figure. Run from nv_model/ or pass paths to gpu.log files.
"""
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import pandas as pd
from pathlib import Path

COLS = [
    "timestamp",
    "index",
    "util_gpu",
    "util_mem",
    "mem_used",
    "mem_total",
    "power",
    "temp",
]

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
LOG_1G = ROOT / "p4d_nv_1g_1o" / "gpu.log"
LOG_8G = ROOT / "p4d_nv_8g_8o" / "gpu.log"


def load_log(path):
    df = pd.read_csv(path, header=None, names=COLS, skipinitialspace=True)
    df["timestamp"] = pd.to_datetime(
        df["timestamp"],
        format="%Y/%m/%d %H:%M:%S.%f",
        errors="coerce",
    )
    df["index"] = pd.to_numeric(df["index"], errors="coerce")
    for c in ["util_gpu", "util_mem", "mem_used", "mem_total", "power", "temp"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.dropna(subset=["timestamp", "index"]).copy()


def plot_panel(ax_util, ax_mem, df, title):
    for gpu_id, g in df.groupby("index"):
        g = g.sort_values("timestamp")
        ax_util.plot(g["timestamp"], g["util_gpu"], label=f"GPU {int(gpu_id)}")
        ax_mem.plot(g["timestamp"], g["mem_used"], label=f"GPU {int(gpu_id)}")
    ax_util.set_ylabel("GPU Utilization (%)")
    ax_mem.set_ylabel("Memory Used (MB)")
    ax_util.set_title(title)
    ax_util.legend(ncol=2, fontsize=7, loc="upper right")


def main():
    df1 = load_log(LOG_1G)
    df8 = load_log(LOG_8G)

    # 2 columns: left = 1 GPU run, right = 8 GPU run; each column has util (top) + mem (bottom)
    fig, ((ax1_util, ax8_util), (ax1_mem, ax8_mem)) = plt.subplots(
        nrows=2,
        ncols=2,
        sharex="col",
        figsize=(15, 6),
        constrained_layout=True,
    )

    plot_panel(ax1_util, ax1_mem, df1, "1 GPU")
    plot_panel(ax8_util, ax8_mem, df8, "8 GPUs")

    # Time format on bottom row
    for ax in (ax1_mem, ax8_mem):
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
        ax.xaxis.get_offset_text().set_visible(False)

    def _k_formatter(x, pos):
        if abs(x) >= 1000:
            v = x / 1000.0
            if abs(v - round(v)) < 1e-9:
                return f"{int(round(v))}k"
            return f"{v:.1f}k"
        return f"{int(x)}"

    for ax in (ax1_mem, ax8_mem):
        ax.yaxis.set_major_formatter(FuncFormatter(_k_formatter))

    out_png = SCRIPT_DIR / "gpu_util_and_mem_combined.png"
    fig.savefig(out_png, dpi=300)
    print(f"Saved: {out_png}")


if __name__ == "__main__":
    main()
