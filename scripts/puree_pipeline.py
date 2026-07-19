#!/usr/bin/env python
"""
Unified PUREE pipeline for Spatial-MetScore: tumor-purity scoring, spatial
visualization, and GMM-based tumor/normal classification, all in one script.

Subcommands
-----------
  validate                        Reproduce PUREE's shipped test case (sanity check).
  pseudobulk   --output PATH       Sum raw counts per sample -> one purity per slide (54 rows).
  perspot      --output PATH       Score every spot individually across all 54 slides.
  viz-cohort                       Cohort-level spatial grid + boxplot + per-spot-vs-pseudobulk figures.
  viz-persample                    One PNG per sample + one combined PNG per patient.
  gmm          --k {2,3}           Gaussian-mixture tumor/normal split on pooled per-spot purity.

Heavy modes (pseudobulk/perspot/viz-*) are meant to run via SBATCH
(see scripts/slurm/puree.sh), never on the login node.
"""
import argparse
import json
import os
import sys
from collections import defaultdict

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from joblib import load

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
ROOT = "/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore"
DATA_ROOT = f"{ROOT}/data/raw/Samples"
PUREE_OUT = f"{ROOT}/data/processed/puree"
PATIENT_META = f"{ROOT}/data/raw/patient_metadata.csv"

PUREE_DIR = "/scratch/rprest2/Spatial-MetScore/tools/PUREE"
CONV = os.path.join(PUREE_DIR, "data/gene_id_conversion_table.csv")
RANK_UNIV = os.path.join(PUREE_DIR, "data/ranking_universe.csv")
PRED_GENES = os.path.join(PUREE_DIR, "data/predictive_genes.csv")
MODEL = os.path.join(PUREE_DIR, "models/model.joblib")
IMPUTER = os.path.join(PUREE_DIR, "models/median_imputer.joblib")

SAMPLES_MANIFEST = "/scratch/rprest2/Spatial-MetScore/scripts/copykat_samples.txt"
FIG_DIR = "/scratch/rprest2/Spatial-MetScore/figures/puree"

TISSUE_ORDER = ["Primary PDAC", "Liver metastasis", "Lung metastasis", "Peritoneal metastasis"]
TISSUE_COLOR = {"Primary PDAC": "#4C72B0", "Liver metastasis": "#DD8452",
                 "Lung metastasis": "#55A868", "Peritoneal metastasis": "#C44E52"}


def samples():
    return [x.strip() for x in open(SAMPLES_MANIFEST) if x.strip()]


# ============================================================================
# Section 1: scoring (validate / pseudobulk / perspot)
# ============================================================================

def read_h5(sample):
    """Return (spots x genes) CSR, spot barcodes, ensembl gene ids. QC nCount>=100."""
    h5 = os.path.join(DATA_ROOT, sample, "filtered_feature_bc_matrix.h5")
    with h5py.File(h5, 'r') as f:
        g = f['matrix']
        data = g['data'][:]; indices = g['indices'][:]; indptr = g['indptr'][:]; shape = g['shape'][:]
        M = sp.csc_matrix((data, indices, indptr), shape=tuple(shape))  # genes x spots (10x CSC)
        barcodes = [b.decode() for b in g['barcodes'][:]]
        ens = [x.decode() for x in g['features']['id'][:]]
    ens = [e.split('.')[0] for e in ens]  # strip ensembl version suffix
    counts = np.asarray(M.sum(axis=0)).ravel()
    keep = counts >= 100
    M = M[:, keep]; barcodes = [b for b, k in zip(barcodes, keep) if k]
    return M.T.tocsr(), barcodes, ens  # spots x genes


def score(predictor, gene_ids="ENSEMBL"):
    selected = pd.read_csv(PRED_GENES)['ENSEMBL_ID']
    universe = pd.read_csv(RANK_UNIV)['ENSEMBL_ID']
    predictor.clean_data(CONV, universe, input_ids=gene_ids)
    predictor.rank_normalize_data(rank_method='min')
    predictor.filter_genes(selected)
    imp = load(IMPUTER)
    vals = imp.transform(predictor.data)
    predictor.data = pd.DataFrame(vals, index=predictor.data.index, columns=predictor.data.columns)
    predictor.predict_purities(MODEL)
    return predictor.purities


def _predictor():
    sys.path.insert(0, PUREE_DIR)
    from scripts.purity_predictor import PurityPredictor
    return PurityPredictor


def run_validate():
    PurityPredictor = _predictor()
    p = PurityPredictor(); p.read_data(os.path.join(PUREE_DIR, "tests/expression.tsv"))
    got = score(p, "ENSEMBL")
    exp = pd.read_csv(os.path.join(PUREE_DIR, "tests/purities_expected.tsv"), sep='\t', index_col=0)
    m = got.join(exp, lsuffix='_got', rsuffix='_exp'); m['abs_diff'] = (m.iloc[:, 0] - m.iloc[:, 1]).abs()
    print(m.to_string()); print("MAX ABS DIFF:", m['abs_diff'].max())
    assert m['abs_diff'].max() < 1e-6, "VALIDATION FAILED"
    print("VALIDATION PASSED")


def run_pseudobulk(outp):
    PurityPredictor = _predictor()
    rows = {}
    for s in samples():
        M, bc, ens = read_h5(s)                    # spots x genes
        pb = np.asarray(M.sum(axis=0)).ravel()      # gene sums
        ser = pd.Series(pb, index=ens).groupby(level=0).sum()  # collapse dup ensembl
        rows[s] = ser
        print("pseudobulk", s, M.shape, flush=True)
    df = pd.DataFrame(rows).T.fillna(0.0)           # samples x genes (full union)
    print("pseudobulk matrix", df.shape, flush=True)
    p = PurityPredictor(); p.data = df
    got = score(p, "ENSEMBL"); got.to_csv(outp)
    print("written", outp, got.shape); print(got.describe().to_string())


def run_perspot(outp):
    PurityPredictor = _predictor()
    universe = pd.read_csv(RANK_UNIV)['ENSEMBL_ID'].tolist()
    mats = []; index = []; ens_ref = None
    for s in samples():
        M, bc, ens = read_h5(s)                     # spots x genes CSR
        if ens_ref is None:
            ens_ref = ens
        else:
            assert ens == ens_ref, "gene order differs at %s" % s
        mats.append(M); index += ["%s|%s" % (s, b) for b in bc]
        print("loaded", s, M.shape, flush=True)
    combined = sp.vstack(mats, format='csr'); del mats
    print("combined", combined.shape, flush=True)
    # projection matrix: full gene set -> universe cols (sum duplicate ensembl)
    col_of = defaultdict(list)
    for j, e in enumerate(ens_ref):
        col_of[e].append(j)
    rows = []; cols = []; present = np.zeros(len(universe), dtype=bool)
    for k, e in enumerate(universe):
        for j in col_of.get(e, []):
            rows.append(j); cols.append(k); present[k] = True
    P = sp.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)),
                       shape=(len(ens_ref), len(universe)))
    proj = (combined @ P)                            # spots x universe sparse
    del combined
    arr = np.asarray(proj.todense(), dtype=np.float32); del proj
    arr[:, ~present] = np.nan                         # genes absent from input -> NaN (PUREE convention)
    df = pd.DataFrame(arr, index=index, columns=universe); del arr
    print("universe matrix", df.shape, flush=True)
    p = PurityPredictor(); p.data = df
    got = score(p, "ENSEMBL"); got.index.name = 'barcode'
    got['sample'] = [b.split('|')[0] for b in got.index]
    got.to_csv(outp)
    print("written", outp, got.shape, flush=True)
    print(got.groupby('sample')['purity'].agg(['count', 'mean', 'std']).to_string())


# ============================================================================
# Section 2: spatial-viz helpers
# ============================================================================

def coords(sample, samples_root=DATA_ROOT):
    """Spot pixel coordinates; handles both SpaceRanger v1 and v2 tissue_positions formats."""
    base = f"{samples_root}/{sample}/spatial"
    f1, f2 = f"{base}/tissue_positions_list.csv", f"{base}/tissue_positions.csv"
    if os.path.exists(f1):
        c = pd.read_csv(f1, header=None,
                         names=["bc", "in_tissue", "row", "col", "pxl_row", "pxl_col"])
    else:
        c = pd.read_csv(f2)
        c.columns = [x.strip() for x in c.columns]
        c = c.rename(columns={"barcode": "bc", "array_row": "row", "array_col": "col",
                               "pxl_row_in_fullres": "pxl_row", "pxl_col_in_fullres": "pxl_col"})
    return c[["bc", "pxl_row", "pxl_col"]]


def spot_diameter(sample, samples_root=DATA_ROOT):
    with open(f"{samples_root}/{sample}/spatial/scalefactors_json.json") as fh:
        return json.load(fh)["spot_diameter_fullres"]


def spot_scatter(ax, x, y, c, diameter_fullres, cmap="viridis", vmin=None, vmax=None,
                  fudge=1.05, pad_frac=0.04, **kwargs):
    """Scatter spots at their TRUE relative size: marker diameter in points is derived
    from the axes' actual rendered pixel size and the real Visium spot diameter (data
    units), not a hardcoded constant — so spots tile edge-to-edge regardless of panel
    size or sample. Call AFTER final gridspec/subplot layout, before savefig."""
    x0, x1 = float(x.min()), float(x.max())
    y0, y1 = float(y.min()), float(y.max())
    px = (x1 - x0) * pad_frac or diameter_fullres
    py = (y1 - y0) * pad_frac or diameter_fullres
    ax.set_xlim(x0 - px, x1 + px)
    ax.set_ylim(y0 - py, y1 + py)
    ax.set_aspect("equal")
    ax.figure.canvas.draw()
    bbox = ax.get_window_extent()
    dpi = ax.figure.dpi
    data_w = (x1 + px) - (x0 - px)
    px_per_data = bbox.width / data_w if data_w > 0 else 1.0
    diameter_px = diameter_fullres * px_per_data
    diameter_pts = diameter_px * 72.0 / dpi
    s = (diameter_pts * fudge) ** 2
    return ax.scatter(x, y, c=c, cmap=cmap, vmin=vmin, vmax=vmax, s=s, linewidths=0, **kwargs)


def _load_meta():
    meta = pd.read_csv(PATIENT_META)
    meta["sample"] = meta["Sample_ID"].str.upper()
    tt = dict(zip(meta["sample"], meta["Tissue_Type"]))
    pat = dict(zip(meta["sample"], meta["patient"]))
    return tt, pat


def run_viz_cohort():
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from matplotlib.patches import Patch
    from scipy.stats import pearsonr

    os.makedirs(FIG_DIR, exist_ok=True)
    ps = pd.read_csv(f"{PUREE_OUT}/puree_perspot_allslides.csv")
    ps["bc"] = ps["barcode"].str.split("|").str[1]
    pb = pd.read_csv(f"{PUREE_OUT}/puree_pseudobulk_54.csv"); pb.columns = ["sample", "pb_purity"]
    tt, _ = _load_meta()
    samps = sorted(ps["sample"].unique(), key=lambda s: (TISSUE_ORDER.index(tt.get(s, "Primary PDAC")), s))

    # Figure 1: spatial purity maps, grid
    n = len(samps); ncol = 7; nrow = int(np.ceil(n / ncol))
    fig = plt.figure(figsize=(ncol * 2.3, nrow * 2.4))
    gs = GridSpec(nrow, ncol, figure=fig, hspace=0.35, wspace=0.15)
    sc = None
    for i, s in enumerate(samps):
        ax = fig.add_subplot(gs[i // ncol, i % ncol])
        d = ps[ps["sample"] == s].merge(coords(s), on="bc", how="left").dropna(subset=["pxl_row", "pxl_col"])
        sc = spot_scatter(ax, d["pxl_col"].values, -d["pxl_row"].values, d["purity"].values,
                           spot_diameter(s), cmap="viridis", vmin=0.2, vmax=1.0)
        ax.set_title(f"{s}\n{tt.get(s,'?')[:12]}", fontsize=6)
        ax.set_xticks([]); ax.set_yticks([])
        for spn in ax.spines.values():
            spn.set_edgecolor(TISSUE_COLOR.get(tt.get(s), "k")); spn.set_linewidth(1.5)
    cbar = fig.colorbar(sc, ax=fig.axes, shrink=0.4, pad=0.01, aspect=30)
    cbar.set_label("PUREE purity")
    fig.suptitle("PUREE per-spot tumor purity — all 54 slides", y=0.995, fontsize=13)
    fig.savefig(f"{FIG_DIR}/puree_perspot_spatial_grid.png", dpi=140, bbox_inches="tight")
    plt.close(fig); print("fig1 done", flush=True)

    # Figure 2: per-sample distributions by tissue
    ps["tissue"] = ps["sample"].map(tt)
    order = sorted(samps, key=lambda s: (TISSUE_ORDER.index(tt.get(s, "Primary PDAC")),
                   ps[ps["sample"] == s]["purity"].median()))
    fig, ax = plt.subplots(figsize=(15, 5))
    data = [ps[ps["sample"] == s]["purity"].values for s in order]
    bp = ax.boxplot(data, positions=range(len(order)), widths=0.7, showfliers=False,
                     patch_artist=True, medianprops=dict(color="k"))
    for patch, s in zip(bp["boxes"], order):
        patch.set_facecolor(TISSUE_COLOR.get(tt.get(s), "grey")); patch.set_alpha(0.8)
    ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=90, fontsize=7)
    ax.set_ylabel("PUREE per-spot purity")
    ax.set_title("Per-spot purity distribution by slide (grouped/colored by tissue)")
    ax.legend(handles=[Patch(facecolor=TISSUE_COLOR[t], label=t) for t in TISSUE_ORDER],
              loc="upper left", fontsize=8)
    fig.tight_layout(); fig.savefig(f"{FIG_DIR}/puree_perspot_boxplot.png", dpi=150, bbox_inches="tight")
    plt.close(fig); print("fig2 done", flush=True)

    # Figure 3: per-spot mean vs pseudobulk
    agg = ps.groupby("sample")["purity"].agg(["mean", "median", "std", "count"]).reset_index().merge(pb, on="sample")
    agg["tissue"] = agg["sample"].map(tt)
    r, _ = pearsonr(agg["mean"], agg["pb_purity"])
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(agg["pb_purity"], agg["mean"], c=agg["tissue"].map(TISSUE_COLOR), s=55, edgecolor="k", linewidth=0.4)
    lims = [0.2, 0.9]; ax.plot(lims, lims, "--", c="grey", lw=1)
    ax.set_xlabel("Pseudobulk purity (54)"); ax.set_ylabel("Per-spot mean purity")
    ax.set_title(f"Per-spot mean vs pseudobulk   r={r:.2f}")
    ax.set_xlim(lims); ax.set_ylim(lims)
    fig.tight_layout(); fig.savefig(f"{FIG_DIR}/puree_perspot_vs_pseudobulk.png", dpi=150, bbox_inches="tight")
    plt.close(fig); print("fig3 done", flush=True)

    agg.to_csv(f"{PUREE_OUT}/puree_perspot_sample_summary.csv", index=False)
    print("all done", agg.shape, "r_perspot_vs_pb=%.3f" % r, flush=True)
    print(agg[["sample", "tissue", "mean", "median", "std", "count", "pb_purity"]].to_string(index=False))


def run_viz_persample():
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(FIG_DIR, exist_ok=True)
    ps = pd.read_csv(f"{PUREE_OUT}/puree_perspot_allslides.csv")
    ps["bc"] = ps["barcode"].str.split("|").str[1]
    tt, pat = _load_meta()
    samps = sorted(ps["sample"].unique())
    VMIN, VMAX = 0.2, 1.0

    # Part A: one PNG per sample, mirroring copykat/<sample>/<sample>_copykat_spatial.png
    for s in samps:
        d = ps[ps["sample"] == s].merge(coords(s), on="bc", how="left").dropna(subset=["pxl_row", "pxl_col"])
        outdir = f"{FIG_DIR}/{s}"; os.makedirs(outdir, exist_ok=True)
        fig, ax = plt.subplots(figsize=(5, 5))
        sc = spot_scatter(ax, d["pxl_col"].values, -d["pxl_row"].values, d["purity"].values,
                           spot_diameter(s), cmap="viridis", vmin=VMIN, vmax=VMAX)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{s}  —  {tt.get(s,'?')}\nPUREE purity", fontsize=11)
        cb = fig.colorbar(sc, ax=ax, shrink=0.75); cb.set_label("PUREE purity")
        fig.tight_layout()
        fig.savefig(f"{outdir}/{s}_puree_spatial.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
    print("part A done: per-sample pngs written", len(samps), flush=True)

    # Part B: one combined figure per patient, all their slides side by side
    os.makedirs(f"{FIG_DIR}/by_patient", exist_ok=True)
    patients = sorted(set(pat.get(s) for s in samps if pat.get(s) is not None),
                       key=lambda x: int(''.join(ch for ch in x if ch.isdigit())))
    for p in patients:
        psamp = sorted([s for s in samps if pat.get(s) == p],
                        key=lambda s: TISSUE_ORDER.index(tt.get(s, "Primary PDAC")) if tt.get(s) in TISSUE_ORDER else 99)
        n = len(psamp)
        if n == 0:
            continue
        fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 4.6))
        if n == 1:
            axes = [axes]
        sc = None
        for ax, s in zip(axes, psamp):
            d = ps[ps["sample"] == s].merge(coords(s), on="bc", how="left").dropna(subset=["pxl_row", "pxl_col"])
            sc = spot_scatter(ax, d["pxl_col"].values, -d["pxl_row"].values, d["purity"].values,
                               spot_diameter(s), cmap="viridis", vmin=VMIN, vmax=VMAX)
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_title(f"{s}\n{tt.get(s,'?')}", fontsize=10, color=TISSUE_COLOR.get(tt.get(s), "k"))
            for spn in ax.spines.values():
                spn.set_edgecolor(TISSUE_COLOR.get(tt.get(s), "k")); spn.set_linewidth(2)
        cb = fig.colorbar(sc, ax=axes, shrink=0.7, pad=0.01); cb.set_label("PUREE purity")
        fig.suptitle(f"{p} — PUREE per-spot tumor purity across sites", fontsize=13, y=1.02)
        fig.savefig(f"{FIG_DIR}/by_patient/{p}_puree_combined.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
    print("part B done: per-patient pngs written", len(patients), flush=True)


# ============================================================================
# Section 3: GMM tumor/normal classification (k=2 or k=3)
# ============================================================================

def run_gmm(k):
    from sklearn.mixture import GaussianMixture

    csv = f"{PUREE_OUT}/puree_perspot_allslides.csv"
    df = pd.read_csv(csv)
    x = df['purity'].values.reshape(-1, 1)
    print("n spots", len(x), "range", float(x.min()), float(x.max()), flush=True)

    fits = {}
    for kk in (1, 2, 3):
        g = GaussianMixture(n_components=kk, n_init=10, random_state=0, covariance_type='full').fit(x)
        fits[kk] = dict(bic=g.bic(x), aic=g.aic(x), model=g)
        print("k=%d bic=%.1f aic=%.1f" % (kk, g.bic(x), g.aic(x)), flush=True)

    if k == 2:
        g2 = fits[2]['model']
        means = g2.means_.ravel(); stds = np.sqrt(g2.covariances_.ravel()); weights = g2.weights_.ravel()
        order = np.argsort(means); normal_i, tumor_i = int(order[0]), int(order[1])

        p_tumor = g2.predict_proba(x)[:, tumor_i]
        label = np.where(p_tumor >= 0.5, 'tumor', 'normal')

        grid = np.linspace(means[normal_i], means[tumor_i], 20001).reshape(-1, 1)
        pg = g2.predict_proba(grid)[:, tumor_i]
        threshold = float(grid[np.argmin(np.abs(pg - 0.5)), 0])
        print("threshold purity=%.4f" % threshold, flush=True)

        df['p_tumor'] = p_tumor; df['label'] = label
        df.to_csv(f"{PUREE_OUT}/puree_perspot_gmm_labels.csv", index=False)

        per = df.groupby('sample').agg(n_spots=('label', 'size'),
            n_tumor=('label', lambda s: (s == 'tumor').sum())).reset_index()
        per['frac_tumor'] = per['n_tumor'] / per['n_spots']
        per.to_csv(f"{PUREE_OUT}/puree_perslide_tumor_fraction_gmm.csv", index=False)
        print(per.to_string(), flush=True)

        counts, edges = np.histogram(x.ravel(), bins=80, range=(0, 1), density=True)
        params = dict(n_spots=int(len(x)),
            bic={str(kk): float(fits[kk]['bic']) for kk in fits},
            aic={str(kk): float(fits[kk]['aic']) for kk in fits},
            means=means.tolist(), stds=stds.tolist(), weights=weights.tolist(),
            normal_i=normal_i, tumor_i=tumor_i, threshold=threshold,
            n_tumor=int((label == 'tumor').sum()), n_normal=int((label == 'normal').sum()),
            frac_tumor_overall=float((label == 'tumor').mean()),
            hist_counts=counts.tolist(), hist_edges=edges.tolist())
        json.dump(params, open(f"{PUREE_OUT}/gmm_purity_params.json", "w"), indent=2)

    elif k == 3:
        g = fits[3]['model']
        means = g.means_.ravel(); stds = np.sqrt(g.covariances_.ravel()); weights = g.weights_.ravel()
        order = np.argsort(means)  # low -> high
        print("components (sorted low->high):", flush=True)
        for r, i in enumerate(order):
            print("  comp%d  mean=%.3f  std=%.3f  weight=%.3f" % (r, means[i], stds[i], weights[i]), flush=True)

        hard = g.predict(x)  # component index per spot
        gA = np.where(np.isin(hard, [order[1], order[2]]), 'tumor', 'normal')  # normal = lowest only
        gB = np.where(hard == order[2], 'tumor', 'normal')                    # normal = lowest 2

        df['comp'] = pd.Series(hard).map({order[0]: 'low', order[1]: 'mid', order[2]: 'high'}).values
        df['label_A'] = gA; df['label_B'] = gB
        df.to_csv(f"{PUREE_OUT}/puree_perspot_gmm3_labels.csv", index=False)

        def perslide(lab):
            t = pd.Series(lab)
            d = df.assign(_l=t.values).groupby('sample')['_l'].agg(
                n_spots='size', n_tumor=lambda s: (s == 'tumor').sum())
            d = d.reset_index(); d['frac_tumor'] = d['n_tumor'] / d['n_spots']; return d

        pA = perslide(gA); pB = perslide(gB)
        pA.to_csv(f"{PUREE_OUT}/puree_perslide_gmm3_groupA.csv", index=False)
        pB.to_csv(f"{PUREE_OUT}/puree_perslide_gmm3_groupB.csv", index=False)
        print("\nGrouping A (normal=low only)   overall frac_tumor=%.3f, per-slide range %.3f-%.3f" % (
            (gA == 'tumor').mean(), pA['frac_tumor'].min(), pA['frac_tumor'].max()), flush=True)
        print("Grouping B (normal=low+mid)    overall frac_tumor=%.3f, per-slide range %.3f-%.3f" % (
            (gB == 'tumor').mean(), pB['frac_tumor'].min(), pB['frac_tumor'].max()), flush=True)

        counts, edges = np.histogram(x.ravel(), bins=80, range=(0, 1.35), density=True)
        params = dict(means=means.tolist(), stds=stds.tolist(), weights=weights.tolist(),
            order=order.tolist(), bic=float(g.bic(x)), aic=float(g.aic(x)),
            hist_counts=counts.tolist(), hist_edges=edges.tolist())
        json.dump(params, open(f"{PUREE_OUT}/gmm3_purity_params.json", "w"), indent=2)

    print("DONE", flush=True)


# ============================================================================
# CLI
# ============================================================================

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    sub.add_parser("validate")

    p_pb = sub.add_parser("pseudobulk"); p_pb.add_argument("--output", required=True)
    p_ps = sub.add_parser("perspot"); p_ps.add_argument("--output", required=True)

    sub.add_parser("viz-cohort")
    sub.add_parser("viz-persample")

    p_gmm = sub.add_parser("gmm"); p_gmm.add_argument("--k", type=int, choices=[2, 3], required=True)

    a = ap.parse_args()

    if a.cmd in ("validate", "pseudobulk", "perspot"):
        os.chdir(PUREE_DIR)
        if a.cmd == "validate":
            run_validate()
        elif a.cmd == "pseudobulk":
            run_pseudobulk(a.output)
        elif a.cmd == "perspot":
            run_perspot(a.output)
    elif a.cmd == "viz-cohort":
        run_viz_cohort()
    elif a.cmd == "viz-persample":
        run_viz_persample()
    elif a.cmd == "gmm":
        run_gmm(a.k)
