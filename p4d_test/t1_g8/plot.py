import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import pandas as pd

# Adjust names if your nvidia-smi version outputs slightly different headers
cols = [
    "timestamp",
    "index",
    "util_gpu",
    "util_mem",
    "mem_used",
    "mem_total",
    "power",
    "temp",
]


df = pd.read_csv("gpu.log", header=None, names=cols, skipinitialspace=True)
df["timestamp"] = pd.to_datetime(
    df["timestamp"],
    format="%Y/%m/%d %H:%M:%S.%f",
    errors="coerce",
)
df["index"] = pd.to_numeric(df["index"], errors="coerce")
for c in ["util_gpu", "util_mem", "mem_used", "mem_total", "power", "temp"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

df = df.dropna(subset=["timestamp", "index"]).copy()

# Two stacked subplots: utilization (top), memory used (bottom)
fig, (ax1, ax2) = plt.subplots(
    nrows=2,
    ncols=1,
    sharex=True,
    figsize=(8, 5),
    constrained_layout=True,
)

for gpu_id, g in df.groupby("index"):
    g = g.sort_values("timestamp")
    ax1.plot(g["timestamp"], g["util_gpu"], label=f"GPU {int(gpu_id)}")
    ax2.plot(g["timestamp"], g["mem_used"], label=f"GPU {int(gpu_id)}")

# Format x-axis as time only and remove the extra date offset text (e.g. "12")
ax2.xaxis.set_major_locator(mdates.AutoDateLocator())
ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
ax2.xaxis.get_offset_text().set_visible(False)

ax1.set_ylabel("GPU Utilization (%)")
ax2.set_ylabel("Memory Used (MB)")
ax2.set_xlabel("Time")

# Format memory ticks as 1k, 4k, etc (k = 1000 MiB)
def _k_formatter(x, pos):
    if abs(x) >= 1000:
        v = x / 1000.0
        # show no decimals if integer, else one decimal
        if abs(v - round(v)) < 1e-9:
            return f"{int(round(v))}k"
        return f"{v:.1f}k"
    return f"{int(x)}"

ax2.yaxis.set_major_formatter(FuncFormatter(_k_formatter))

# Legend: 2 columns (=> 4 rows for 8 GPUs), inside top plot (upper-left)
ax1.legend(ncol=2, fontsize=8, loc="upper left")

out_png = "gpu_util_and_mem.png"
fig.savefig(out_png, dpi=300)

print(f"Saved: {out_png}")
