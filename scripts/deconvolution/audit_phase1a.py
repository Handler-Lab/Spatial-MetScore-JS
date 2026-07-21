#!/usr/bin/env python
"""Phase 1a audit: peng_2019, hlca, lambrechts_2018.
Reads obs/var structure with h5py only (no scanpy dep). h5ad and loom are both HDF5.
Reports: dims, raw-counts location, gene-ID namespace, candidate label cols + value_counts.
"""
import h5py, numpy as np, os, gzip, shutil, json, sys

BASE = "/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/data/raw/scRNAseq"
OUT  = "audit_out"
os.makedirs(OUT, exist_ok=True)
summary = {}

def h5ad_obs_cols(f):
    """Return list of obs column names (AnnData >=0.7 layout)."""
    obs = f["obs"]
    cols = []
    # column-order attr if present
    for k in obs.keys():
        if k in ("__categories",): continue
        cols.append(k)
    return cols

def read_h5ad_col_valuecounts(f, col, top=60):
    """Return dict label->count for an obs column, handling categorical group or flat dataset."""
    obs = f["obs"]
    node = obs[col]
    try:
        if isinstance(node, h5py.Group) and "categories" in node and "codes" in node:
            cats = node["categories"][:]
            codes = node["codes"][:]
        elif isinstance(node, h5py.Dataset):
            # could be categorical via __categories (besca/old anndata) or plain
            arr = node[:]
            if "__categories" in obs and col in obs["__categories"]:
                cats = obs["__categories"][col][:]
                codes = arr
            else:
                # plain (string or numeric) -> treat values directly
                vals = arr
                if vals.dtype.kind in ("S","O"):
                    vals = np.array([v.decode() if isinstance(v,bytes) else v for v in vals])
                uniq, cnt = np.unique(vals, return_counts=True)
                order = np.argsort(-cnt)
                return {str(uniq[i]): int(cnt[i]) for i in order[:top]}, len(uniq)
        else:
            return {}, 0
        cats = np.array([c.decode() if isinstance(c,bytes) else c for c in cats])
        codes = np.asarray(codes)
        valid = codes[codes >= 0]
        cnt = np.bincount(valid, minlength=len(cats))
        order = np.argsort(-cnt)
        d = {str(cats[i]): int(cnt[i]) for i in order if cnt[i] > 0}
        return d, int((cnt > 0).sum())
    except Exception as e:
        return {"ERROR": str(e)}, -1

def counts_location(f):
    """Heuristic: report shape/dtype of X and raw/X to identify integer counts."""
    info = {}
    for path in ["X", "raw/X", "layers/counts"]:
        if path in f:
            node = f[path]
            if isinstance(node, h5py.Group):  # sparse CSR/CSC
                dt = node["data"].dtype
                shape = tuple(node.attrs.get("shape", [None,None]))
                # sample a few data values to test integrality
                sample = node["data"][:1000]
                is_int = bool(np.allclose(sample, np.round(sample)))
                info[path] = {"sparse": True, "dtype": str(dt), "shape": [int(x) for x in shape] if shape[0] is not None else None, "sample_integer": is_int}
            else:
                dt = node.dtype; shape = node.shape
                sample = node[0, :min(1000, shape[1])] if len(shape)==2 else node[:1000]
                is_int = bool(np.allclose(sample, np.round(sample)))
                info[path] = {"sparse": False, "dtype": str(dt), "shape": list(shape), "sample_integer": is_int}
    return info

def var_geneids(f, key="var"):
    info = {}
    if key in f:
        v = f[key]
        info["var_keys"] = list(v.keys())
        # try index / feature_name / SYMBOL / ENSEMBL samples
        for cand in ["_index","index","feature_name","SYMBOL","ENSEMBL","gene_symbols","gene_ids"]:
            if cand in v:
                node = v[cand]
                try:
                    if isinstance(node, h5py.Group) and "categories" in node:
                        s = node["categories"][:5]
                    else:
                        s = node[:5]
                    s = [x.decode() if isinstance(x,bytes) else str(x) for x in s]
                    info[f"sample::{cand}"] = s
                except Exception as e:
                    info[f"sample::{cand}"] = f"ERR {e}"
    return info

# ---------- peng_2019 ----------
print("="*70); print("PENG_2019"); print("="*70)
peng = {}
raw_fp = f"{BASE}/peng_2019/StdWf1_PRJCA001063_CRC_besca2.raw.h5ad"
ann_fp = f"{BASE}/peng_2019/StdWf1_PRJCA001063_CRC_besca2.annotated.h5ad"
with h5py.File(raw_fp,"r") as f:
    peng["raw_file_counts"] = counts_location(f)
    peng["raw_file_var"] = var_geneids(f)
    peng["raw_n_obs"] = int(f["obs"][list(f["obs"].keys())[0]].shape[0])
with h5py.File(ann_fp,"r") as f:
    peng["ann_obs_cols"] = h5ad_obs_cols(f)
    peng["ann_counts"] = counts_location(f)
    peng["ann_var"] = var_geneids(f)
    peng["ann_n_obs"] = int(f["obs"][list(f["obs"].keys())[0]].shape[0])
    for col in ["Cell_type","celltype0","celltype1","celltype2","celltype3","dblabel"]:
        if col in f["obs"]:
            vc, n = read_h5ad_col_valuecounts(f, col)
            peng[f"labels::{col}"] = {"n_unique": n, "counts": vc}
    # barcode order check
    def barcodes(ff):
        for bk in ["_index","index","CELL"]:
            if bk in ff["obs"]:
                node = ff["obs"][bk]
                arr = node[:] if isinstance(node,h5py.Dataset) else node["categories"][:]
                return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in arr])
        return None
    ann_bc = barcodes(f)
with h5py.File(raw_fp,"r") as f:
    raw_bc = None
    for bk in ["_index","index","CELL"]:
        if bk in f["obs"]:
            node = f["obs"][bk]
            arr = node[:] if isinstance(node,h5py.Dataset) else node["categories"][:]
            raw_bc = np.array([x.decode() if isinstance(x,bytes) else str(x) for x in arr]); break
if ann_bc is not None and raw_bc is not None:
    peng["barcode_order_match"] = bool(len(ann_bc)==len(raw_bc) and np.array_equal(ann_bc, raw_bc))
    peng["barcode_setequal"] = bool(set(ann_bc)==set(raw_bc))
summary["peng_2019"] = peng
print(json.dumps(peng, indent=2)[:4000])

# ---------- hlca ----------
print("="*70); print("HLCA"); print("="*70)
hlca = {}
with h5py.File(f"{BASE}/hlca/hlca_core.h5ad","r") as f:
    hlca["obs_cols"] = h5ad_obs_cols(f)
    hlca["counts"] = counts_location(f)
    hlca["var_X"] = var_geneids(f, "var")
    if "raw" in f and "var" in f["raw"]:
        hlca["raw_var_keys"] = list(f["raw"]["var"].keys())
        # sample gene id from raw/var index + feature_name if present
        rv = f["raw"]["var"]
        for cand in ["_index","index","feature_name"]:
            if cand in rv:
                node = rv[cand]
                try:
                    s = (node["categories"][:5] if isinstance(node,h5py.Group) and "categories" in node else node[:5])
                    hlca[f"raw_var_sample::{cand}"] = [x.decode() if isinstance(x,bytes) else str(x) for x in s]
                except Exception as e:
                    hlca[f"raw_var_sample::{cand}"] = f"ERR {e}"
        hlca["raw_counts"] = counts_location(f["raw"])
    hlca["n_obs"] = int(f["obs"][list(f["obs"].keys())[0]].shape[0])
    for col in ["cell_type","ann_level_1","ann_level_2","ann_level_3","ann_finest_level","scanvi_label"]:
        if col in f["obs"]:
            vc, n = read_h5ad_col_valuecounts(f, col)
            hlca[f"labels::{col}"] = {"n_unique": n, "counts": vc}
summary["hlca"] = hlca
print(json.dumps({k:v for k,v in hlca.items() if not k.startswith("labels::")}, indent=2)[:3000])
for k in hlca:
    if k.startswith("labels::"):
        print(f"\n-- {k} (n_unique={hlca[k]['n_unique']}) --")
        print(json.dumps(hlca[k]["counts"], indent=1)[:2500])

# ---------- lambrechts_2018 (loom.gz) ----------
print("="*70); print("LAMBRECHTS_2018"); print("="*70)
lam = {}
loom_gz = None
for root,_,files in os.walk(f"{BASE}/lambrechts_2018"):
    for fn in files:
        if fn.endswith(".loom.gz"):
            loom_gz = os.path.join(root, fn)
lam["loom_gz_path"] = loom_gz
if loom_gz:
    tmp_loom = os.path.join(os.environ.get("TMPDIR","/tmp"), "lambrechts.loom")
    print("gunzip ->", tmp_loom)
    with gzip.open(loom_gz,"rb") as fi, open(tmp_loom,"wb") as fo:
        shutil.copyfileobj(fi, fo, length=16*1024*1024)
    with h5py.File(tmp_loom,"r") as f:
        lam["root_keys"] = list(f.keys())
        if "matrix" in f:
            m = f["matrix"]
            lam["matrix_shape"] = list(m.shape); lam["matrix_dtype"] = str(m.dtype)
            s = m[:5,:1000] if len(m.shape)==2 else m[:1000]
            lam["matrix_sample_integer"] = bool(np.allclose(s, np.round(s)))
        lam["row_attrs"] = list(f["row_attrs"].keys()) if "row_attrs" in f else []
        lam["col_attrs"] = list(f["col_attrs"].keys()) if "col_attrs" in f else []
        # gene id sample
        for g in ["Gene","GeneName","gene","var_names","Accession"]:
            if "row_attrs" in f and g in f["row_attrs"]:
                s = f["row_attrs"][g][:5]
                lam[f"gene_sample::{g}"] = [x.decode() if isinstance(x,bytes) else str(x) for x in s]
        # candidate cell-type cols in col_attrs
        if "col_attrs" in f:
            for c in f["col_attrs"].keys():
                node = f["col_attrs"][c]
                try:
                    arr = node[:]
                    if arr.dtype.kind in ("S","O"):
                        arr = np.array([x.decode() if isinstance(x,bytes) else x for x in arr])
                        uniq, cnt = np.unique(arr, return_counts=True)
                        if 2 <= len(uniq) <= 200:  # looks categorical
                            order = np.argsort(-cnt)
                            lam[f"colattr::{c}"] = {"n_unique": int(len(uniq)),
                                "counts": {str(uniq[i]): int(cnt[i]) for i in order[:60]}}
                except Exception as e:
                    lam[f"colattr::{c}"] = f"ERR {e}"
    os.remove(tmp_loom)
summary["lambrechts_2018"] = lam
print(json.dumps({k:v for k,v in lam.items() if not k.startswith("colattr::")}, indent=2)[:3000])
for k in lam:
    if k.startswith("colattr::"):
        print(f"\n-- {k} (n_unique={lam[k]['n_unique']}) --")
        print(json.dumps(lam[k]["counts"], indent=1)[:2000])

# ---------- dump native label vocab files ----------
for atlas, d in summary.items():
    with open(f"{OUT}/{atlas}_labels.txt","w") as fo:
        fo.write(f"# {atlas} native label vocabularies\n")
        for k,v in d.items():
            if k.startswith("labels::") or k.startswith("colattr::"):
                fo.write(f"\n## {k}  (n_unique={v.get('n_unique')})\n")
                for lab,ct in v["counts"].items():
                    fo.write(f"{ct}\t{lab}\n")

with open(f"{OUT}/phase1a_summary.json","w") as fo:
    json.dump(summary, fo, indent=2, default=str)
print("\n=== WROTE", OUT, "===")
print(os.listdir(OUT))
