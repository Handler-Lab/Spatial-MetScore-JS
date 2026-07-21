#!/usr/bin/env python
"""Phase 1c + Phase 2 builder.
Consolidates the audit findings into:
  - reference_audit.csv         (one row per atlas)
  - label_harmonization.csv     (native label -> coarse type, per atlas, with n_cells)
  - harmonized_composition.csv   (cells per coarse type per atlas + total)
  - harmonized_composition.png   (grouped/stacked bar chart)
Native label counts are the observed value_counts from the Phase 1a/1b audits.
"""
import csv, os, json
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUTDATA = "/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/data/processed/rctd"
AUDITDIR = os.path.join(OUTDATA, "audit")
HARMDIR  = os.path.join(OUTDATA, "harmonization")
FIGDIR   = "/scratch/rprest2/Spatial-MetScore/figures/deconvolution"
for d in (AUDITDIR, HARMDIR, FIGDIR): os.makedirs(d, exist_ok=True)

# ---------------------------------------------------------------------------
# 1) reference_audit.csv
# ---------------------------------------------------------------------------
audit_rows = [
 dict(atlas="peng_2019", tissue="Pancreas (PDAC primary)", object_class="AnnData (h5ad)",
      n_cells=57530, raw_counts_location=".raw.h5ad :: X (float32, integer-valued)",
      gene_id_type="symbol", symbol_column="var/_index (SYMBOL col is numeric codes-> use _index)",
      chosen_label_column="Cell_type", n_native_labels=10,
      notes="Labels live in .annotated.h5ad (57423 cells); barcode set/order DIFFERS from .raw.h5ad (57530) -> join counts<-labels by barcode; annotated raw/X is NOT integer."),
 dict(atlas="hlca", tissue="Lung (normal/healthy reference)", object_class="AnnData (h5ad)",
      n_cells=584944, raw_counts_location="raw/X (float32, integer-valued); top-level X is lognorm",
      gene_id_type="Ensembl (var/_index); symbol in var/feature_name",
      symbol_column="var/feature_name (& raw/var/feature_name)",
      chosen_label_column="ann_level_3", n_native_labels=25,
      notes="HLCA core. cell_type=50 (CL ontology), ann_level_1=4, ann_level_2=11, ann_level_3=25, ann_finest=61. ann_level_3 chosen as coarse-but-informative."),
 dict(atlas="lambrechts_2018", tissue="Lung (NSCLC tumor)", object_class="loom (HDF5)",
      n_cells=52698, raw_counts_location="matrix (float64, integer-valued)",
      gene_id_type="symbol", symbol_column="row_attrs/Gene",
      chosen_label_column="ClusterName", n_native_labels=40,
      notes="8 compartment looms present; Thienpont_Tumors_52k is the full atlas (52698 cells). 'cancer cells pt N' = NSCLC malignant cells (patient-specific)."),
 dict(atlas="loveless_2025", tissue="Pancreas (PDAC atlas incl. mets)", object_class="Seurat (RNA assay)",
      n_cells=726107, raw_counts_location="RNA@counts (integer)",
      gene_id_type="symbol", symbol_column="rownames",
      chosen_label_column="Clusters", n_native_labels=14,
      notes="FILE DOUBLE-GZIPPED: read via readRDS(gzcon(gzfile(fp,'rb'))). 36601 genes. Has 'If.metastatic..location' + DiseaseState metadata. PDAC/pancreas, not liver."),
]
audit_df = pd.DataFrame(audit_rows)
audit_df.to_csv(os.path.join(AUDITDIR, "reference_audit.csv"), index=False)
print("=== reference_audit.csv ===")
print(audit_df[["atlas","tissue","n_cells","chosen_label_column","n_native_labels"]].to_string(index=False))

# ---------------------------------------------------------------------------
# 2) label_harmonization.csv  --  native label -> coarse type (DRAFT)
# ---------------------------------------------------------------------------
# Coarse vocabulary (draft): Epithelial, Acinar, Endocrine, Fibroblast, Endothelial,
#   Mural(Pericyte/SMC), T_NK, B_Plasma, Myeloid, Mast, Alveolar, (Hepatocyte=none), drop
# Each entry: (atlas, native_label, n_cells, coarse_type)

AUDIT_SRC = os.environ.get("AUDIT_SRC",
    "/scratch/rprest2/.claude-science/jobs/b1d97dcf-ce0b-403b-87b2-65b19ffbfd5d/audit_out")
def parse_labels(fname, section):
    """Parse a '## labels::SECTION' block of '<count>\\t<label>' lines from an audit dump."""
    path=os.path.join(AUDIT_SRC, fname); out={}; grab=False
    for line in open(path):
        s=line.rstrip("\n")
        if s.startswith("## "):
            grab = ("::"+section) in s or s.strip().endswith(section)
            continue
        if grab and "\t" in s:
            cnt,lab = s.split("\t",1)
            try: out[lab.strip()] = int(cnt.strip())
            except ValueError: pass
    return out

peng = parse_labels("peng_2019_labels.txt", "Cell_type")
peng_map = {
 "Ductal cell type 2":"Epithelial","Ductal cell type 1":"Epithelial","Endothelial cell":"Endothelial",
 "Fibroblast cell":"Fibroblast","Stellate cell":"Fibroblast","Macrophage cell":"Myeloid","T cell":"T_NK",
 "B cell":"B_Plasma","Acinar cell":"Acinar","Endocrine cell":"Endocrine"}

loveless = {  # Clusters
 "DUCTAL":277301,"TNK":102067,"MYELOID":76253,"FIBROBLASTS":66192,"ENDOTHELIAL":52388,
 "PERICYTES":35774,"ACINAR":35109,"CYCLING DUCTAL":23276,"ENDOCRINE":17833,"CYCLING TNK":14607,
 "B CELLS":12880,"PLASMA":6117,"MAST":3688,"CYCLING. MYELOID":2622}
loveless_map = {
 "DUCTAL":"Epithelial","CYCLING DUCTAL":"Epithelial","TNK":"T_NK","CYCLING TNK":"T_NK",
 "MYELOID":"Myeloid","CYCLING. MYELOID":"Myeloid","FIBROBLASTS":"Fibroblast","ENDOTHELIAL":"Endothelial",
 "PERICYTES":"Mural","ACINAR":"Acinar","ENDOCRINE":"Endocrine","B CELLS":"B_Plasma","PLASMA":"B_Plasma","MAST":"Mast"}

hlca = parse_labels("hlca_labels.txt", "ann_level_3")
# map keys are matched case-insensitively / flexibly below via _mapkey()
hlca_map_raw = {
 "Macrophages":"Myeloid","Monocytes":"Myeloid","Dendritic cells":"Myeloid",
 "Basal":"Epithelial","Secretory":"Epithelial","Multiciliated lineage":"Epithelial",
 "Submucosal Secretory":"Epithelial","Rare":"Epithelial",
 "AT2":"Alveolar","AT1":"Alveolar",
 "T cell lineage":"T_NK","Innate lymphoid cell NK":"T_NK","B cell lineage":"B_Plasma",
 "EC capillary":"Endothelial","EC venous":"Endothelial","EC arterial":"Endothelial",
 "Lymphatic EC mature":"Endothelial","Lymphatic EC differentiating":"Endothelial",
 "Lymphatic EC proliferating":"Endothelial",
 "Fibroblasts":"Fibroblast","Myofibroblasts":"Fibroblast",
 "SM activated stress response":"Mural","Smooth muscle FAM83D+":"Mural",
 "Mast cells":"Mast","None":"drop"}

lam = {  # ClusterName (Tumors_52k) - trailing spaces stripped
 "CD8+ T cells":12040,"CD4+ T cells":9617,"macrophages":8074,"cancer cells pt 5":2799,
 "follicular B cells":2675,"cancer cells pt 4":2045,"natural killer cells":1741,
 "regulatory T cells":1513,"cancer cells pt 2":1365,"plasma B cells":1239,"cancer cells pt 3":1078,
 "MALT B cells":892,"Langerhans cells":714,"tumour endothelial cell":618,"mast cells":613,
 "normal endothelial cell":569,"cuboidal alveolar type 2 (AT2) cells":520,"granulocytes":478,
 "lower quality endothelial cell":320,"COL12A1-expressing fibroblasts":315,
 "lower quality alveolar cell":310,"flat alveolar type 1 (AT1) cells":303,
 "monocyte-derived dendritic cells":298,"COL4A2-expressing fibroblasts":266,
 "lower quality fibroblasts":240,"PLA2G2A-expressing fibroblasts":219,"COPD-injured alveolar cells":205,
 "GABARAP-expressing fibroblasts":195,"cross-presenting dendritic cells":192,
 "normal  lung fibroblasts":175,"respiratory epithelial cells":168,"epithelial cell":164,
 "cancer cells pt 1":160,"secretory club cells":136,"plasmacytoid dendritic cells":96,
 "erythroblasts":88,"lymphatic EC":85,"basal cells":68,"TFPI2-expressing fibroblasts":55,
 "lower quality epithelial cell":50}
def lam_coarse(lab):
    l=lab.lower()
    if "cancer cells" in l: return "drop"   # NSCLC malignant - wrong tumor identity for PDAC lung mets (user decision)
    if "t cells" in l or "natural killer" in l: return "T_NK"
    if "b cells" in l: return "B_Plasma"
    if "macrophage" in l or "dendritic" in l or "langerhans" in l or "granulocyte" in l or "monocyte" in l: return "Myeloid"
    if "mast" in l: return "Mast"
    if "fibroblast" in l: return "Fibroblast"
    if "endothelial" in l or "lymphatic ec" in l: return "Endothelial"
    if "alveolar" in l or "at2" in l or "at1" in l: return "Alveolar"
    if "epithelial" in l or "club" in l or "basal" in l: return "Epithelial"
    if "erythroblast" in l: return "drop"
    return "UNMAPPED"
lam_map = {k: lam_coarse(k) for k in lam}

hlca_map = hlca_map_raw

def _norm(s):
    return " ".join(str(s).strip().split()).lower()

def resolve(mp, lab):
    """Look up lab in mp with exact then whitespace/case-normalized fallback."""
    if lab in mp: return mp[lab]
    nmap = {_norm(k): v for k, v in mp.items()}
    return nmap.get(_norm(lab), "UNMAPPED")

harm_rows = []
for atlas, counts, mp in [("peng_2019",peng,peng_map),("hlca",hlca,hlca_map),
                          ("lambrechts_2018",lam,lam_map),("loveless_2025",loveless,loveless_map)]:
    assert counts, f"{atlas}: parsed 0 native labels — check AUDIT_SRC path"
    for lab, n in counts.items():
        harm_rows.append(dict(atlas=atlas, native_label=lab, n_cells=n, coarse_type=resolve(mp, lab)))
harm_df = pd.DataFrame(harm_rows)
harm_df.to_csv(os.path.join(HARMDIR, "label_harmonization.csv"), index=False)
print("\n=== label_harmonization.csv (n rows = %d) ===" % len(harm_df))
print("UNMAPPED labels:", harm_df[harm_df.coarse_type=="UNMAPPED"]["native_label"].tolist() or "none")

# ---------------------------------------------------------------------------
# 3) composition summary (cells per coarse type per atlas)
# ---------------------------------------------------------------------------
comp = (harm_df[harm_df.coarse_type!="drop"]
        .groupby(["coarse_type","atlas"])["n_cells"].sum().unstack(fill_value=0))
comp["TOTAL"] = comp.sum(axis=1)
comp = comp.sort_values("TOTAL", ascending=False)
comp.to_csv(os.path.join(HARMDIR, "harmonized_composition.csv"))
print("\n=== harmonized_composition.csv ===")
print(comp.to_string())

# flag coarse types below RCTD ~25-cell floor within an atlas (where present but tiny)
floor = 25
tiny = []
for ct in comp.index:
    for atlas in ["peng_2019","hlca","lambrechts_2018","loveless_2025"]:
        v = comp.loc[ct, atlas] if atlas in comp.columns else 0
        if 0 < v < floor: tiny.append((ct, atlas, int(v)))
print("\nCoarse types below 25-cell floor in an atlas:", tiny or "none")

# ---------------------------------------------------------------------------
# 4) bar chart
# ---------------------------------------------------------------------------
atlases = ["peng_2019","loveless_2025","hlca","lambrechts_2018"]
palette = {"peng_2019":"#4C72B0","loveless_2025":"#55A868","hlca":"#C44E52","lambrechts_2018":"#8172B3"}
plot_df = comp.drop(columns=["TOTAL"])
fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))
# panel A: grouped absolute (log)
ax=axes[0]; x=np.arange(len(plot_df.index)); w=0.2
for i,a in enumerate(atlases):
    vals = plot_df[a].values if a in plot_df.columns else np.zeros(len(x))
    ax.bar(x+(i-1.5)*w, vals, w, label=a, color=palette[a])
ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(plot_df.index, rotation=45, ha="right")
ax.set_ylabel("cells (log scale)"); ax.set_title("A. Harmonized cells per coarse type, by atlas")
ax.legend(frameon=False, fontsize=9); ax.spines[['top','right']].set_visible(False)
# panel B: composition proportions within atlas (stacked)
ax=axes[1]
prop = plot_df.div(plot_df.sum(axis=0), axis=1).fillna(0)
bottom=np.zeros(len(atlases))
import matplotlib.cm as cm
colors = cm.tab20(np.linspace(0,1,len(prop.index)))
for j,ct in enumerate(prop.index):
    vals=[prop.loc[ct,a] if a in prop.columns else 0 for a in atlases]
    ax.bar(atlases, vals, bottom=bottom, label=ct, color=colors[j]); bottom+=np.array(vals)
ax.set_ylabel("proportion of atlas cells"); ax.set_title("B. Coarse-type composition within each atlas")
ax.set_xticklabels(atlases, rotation=20, ha="right")
ax.legend(bbox_to_anchor=(1.01,1), loc="upper left", frameon=False, fontsize=8)
ax.spines[['top','right']].set_visible(False)
fig.tight_layout()
figpath=os.path.join(FIGDIR, "harmonized_composition.png")
fig.savefig(figpath, dpi=150, bbox_inches="tight")
print("\nwrote figure:", figpath)
print("wrote:", os.path.join(AUDITDIR,"reference_audit.csv"))
print("wrote:", os.path.join(HARMDIR,"label_harmonization.csv"))
print("wrote:", os.path.join(HARMDIR,"harmonized_composition.csv"))
