#!/usr/bin/env python
"""All-slides CopyKAT spatial grid, mirroring puree_perspot_spatial_grid.png layout.
Each spot colored by CopyKAT call: aneuploid (tumor) / diploid / not.defined.
Reuses puree_viz_common for coords, spot sizing, tissue order/colors."""
import os, glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from puree_viz_common import coords, spot_diameter, TISSUE_ORDER, TISSUE_COLOR as tcol

ROOT = "/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore"
CK = f"{ROOT}/data/processed/copykat"
FIG = "/scratch/rprest2/Spatial-MetScore/figures/copykat"
os.makedirs(FIG, exist_ok=True)

CALL_COLOR = {"aneuploid": "#C44E52", "diploid": "#4C72B0", "not.defined": "#CFCFCF"}
CALL_ORDER = ["aneuploid", "diploid", "not.defined"]

meta = pd.read_csv(f"{ROOT}/data/raw/patient_metadata.csv")
meta["sample"] = meta["Sample_ID"].str.upper()
tt = dict(zip(meta["sample"], meta["Tissue_Type"]))

# discover samples with a prediction file
samples = []
for d in sorted(os.listdir(CK)):
    if os.path.exists(f"{CK}/{d}/{d}_copykat_prediction.csv"):
        samples.append(d)
samples = sorted(samples, key=lambda s: (TISSUE_ORDER.index(tt.get(s, "Primary PDAC")), s))
print("n samples", len(samples), flush=True)

def cat_scatter(ax, x, y, colors, diameter_fullres, fudge=1.05, pad_frac=0.04):
    x0, x1 = float(x.min()), float(x.max()); y0, y1 = float(y.min()), float(y.max())
    px = (x1 - x0) * pad_frac or diameter_fullres
    py = (y1 - y0) * pad_frac or diameter_fullres
    ax.set_xlim(x0 - px, x1 + px); ax.set_ylim(y0 - py, y1 + py)
    ax.set_aspect("equal"); ax.figure.canvas.draw()
    bbox = ax.get_window_extent(); dpi = ax.figure.dpi
    data_w = (x1 + px) - (x0 - px)
    px_per_data = bbox.width / data_w if data_w > 0 else 1.0
    diameter_pts = diameter_fullres * px_per_data * 72.0 / dpi
    s = (diameter_pts * fudge) ** 2
    return ax.scatter(x, y, c=colors, s=s, linewidths=0)

n = len(samples); ncol = 7; nrow = int(np.ceil(n / ncol))
fig = plt.figure(figsize=(ncol * 2.3, nrow * 2.4))
gs = GridSpec(nrow, ncol, figure=fig, hspace=0.35, wspace=0.15)
for i, s in enumerate(samples):
    ax = fig.add_subplot(gs[i // ncol, i % ncol])
    pred = pd.read_csv(f"{CK}/{s}/{s}_copykat_prediction.csv")
    pred.columns = ["bc", "call"]
    d = pred.merge(coords(s), on="bc", how="left").dropna(subset=["pxl_row", "pxl_col"])
    colors = d["call"].map(CALL_COLOR).fillna("#CFCFCF").values
    cat_scatter(ax, d["pxl_col"].values, -d["pxl_row"].values, colors, spot_diameter(s))
    frac_an = (d["call"] == "aneuploid").mean()
    ax.set_title(f"{s}\n{tt.get(s,'?')[:12]} ({frac_an*100:.0f}% an)", fontsize=6)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_edgecolor(tcol.get(tt.get(s), "k")); sp.set_linewidth(1.5)

handles = [Patch(facecolor=CALL_COLOR[c], label=c) for c in CALL_ORDER]
fig.legend(handles=handles, title="CopyKAT call", loc="lower center",
           ncol=3, frameon=False, fontsize=10, title_fontsize=10, bbox_to_anchor=(0.5, -0.01))
fig.suptitle("CopyKAT per-spot CNV call — all %d slides" % n, y=0.995, fontsize=13)
fig.savefig(f"{FIG}/copykat_spatial_grid.png", dpi=140, bbox_inches="tight")
plt.close(fig); print("grid done", flush=True)
