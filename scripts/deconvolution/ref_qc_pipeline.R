#!/usr/bin/env Rscript
# =====================================================================
# Composite RCTD reference QC pipeline (all 4 atlases, full cell counts)
#   1) load raw counts + chosen native label from each atlas -> Seurat
#   2) apply label_harmonization.csv (native -> coarse), drop 'drop' rows
#   3) merge on shared genes; PCA
#   4) UMAP uncorrected  (pre)   + UMAP Harmony(batch=atlas)  (post)
#      colored by atlas and by coarse type  -> 4 UMAP panels
#   5) cross-atlas coarse-type mean-profile correlation heatmap
#   6) FindAllMarkers (presto) on coarse labels -> dotplot + csv
# =====================================================================
suppressWarnings(suppressMessages({
  library(Seurat); library(Matrix); library(dplyr); library(ggplot2)
  library(harmony); library(schard); library(pheatmap); library(RColorBrewer)
}))
set.seed(1)
RAW  <- "/group/jshandl-g00/Spatial-MetScore/data/raw/scRNAseq"
PROC <- "/group/jshandl-g00/Spatial-MetScore/data/processed/rctd"
FIG  <- "/scratch/rprest2/Spatial-MetScore/figures/deconvolution"
dir.create(FIG, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(PROC,"harmonization"), recursive=TRUE, showWarnings=FALSE)
HARM <- read.csv(file.path(PROC,"harmonization","label_harmonization.csv"),
                 stringsAsFactors=FALSE, check.names=FALSE)
norm_lab <- function(x) tolower(trimws(gsub("\\s+"," ",x)))
map_for <- function(atlas){
  d <- HARM[HARM$atlas==atlas,]
  setNames(d$coarse_type, norm_lab(d$native_label))
}
apply_map <- function(labels, atlas){
  m <- map_for(atlas); out <- m[norm_lab(labels)]
  out[is.na(out)] <- "UNMAPPED"; unname(out)
}
tt <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), msg))

# ---- 1. loaders: return Seurat obj with meta col 'coarse' + 'atlas' -------
make_obj <- function(counts, labels, atlas){
  if(!inherits(counts,"CsparseMatrix")) counts <- as(counts, "CsparseMatrix")
  stopifnot(ncol(counts)==length(labels))
  coarse <- apply_map(labels, atlas)
  keep <- coarse!="drop" & coarse!="UNMAPPED"
  tt(sprintf("  %s: %d cells, dropping %d (drop/unmapped) -> %d kept",
             atlas, length(coarse), sum(!keep), sum(keep)))
  counts <- counts[, keep, drop=FALSE]; coarse <- coarse[keep]
  colnames(counts) <- paste0(atlas,"|",seq_len(ncol(counts)))
  o <- CreateSeuratObject(counts=counts, project=atlas)
  o$atlas <- atlas; o$coarse <- coarse; o
}

tt("loading peng_2019 (raw counts .raw.h5ad, labels .annotated by barcode)")
peng_raw <- schard::h5ad2seurat(file.path(RAW,"peng_2019","StdWf1_PRJCA001063_CRC_besca2.raw.h5ad"), use.raw=FALSE)
peng_ann <- schard::h5ad2seurat(file.path(RAW,"peng_2019","StdWf1_PRJCA001063_CRC_besca2.annotated.h5ad"), use.raw=FALSE)
# join labels (annotated) onto raw counts by barcode intersection
bc <- intersect(colnames(peng_raw), colnames(peng_ann))
tt(sprintf("  peng barcodes: raw=%d ann=%d shared=%d", ncol(peng_raw), ncol(peng_ann), length(bc)))
pc <- GetAssayData(peng_raw, layer="counts")[, bc]
plabs <- peng_ann$Cell_type[match(bc, colnames(peng_ann))]
peng <- make_obj(pc, as.character(plabs), "peng_2019")

tt("loading hlca (raw counts via use.raw=TRUE; remap Ensembl rownames -> symbols)")
hlca_f <- file.path(RAW,"hlca","hlca_core.h5ad")
hlca_s <- schard::h5ad2seurat(hlca_f, use.raw=TRUE)   # raw/X integer counts, Ensembl rownames
hl_labs <- as.character(hlca_s$ann_level_3)
hlca_cts <- as(GetAssayData(hlca_s, layer="counts"), "CsparseMatrix")
# Ensembl -> symbol using raw/var/feature_name (categorical codes)
suppressMessages(library(rhdf5))
cats  <- as.character(h5read(hlca_f,"raw/var/feature_name/categories"))
codes <- as.integer(h5read(hlca_f,"raw/var/feature_name/codes"))
ens   <- as.character(h5read(hlca_f,"raw/var/_index"))
sym   <- cats[codes+1]              # plain character vector
names(sym) <- ens
stopifnot(all(rownames(hlca_cts) %in% names(sym)))
rownames(hlca_cts) <- as.character(unname(sym[rownames(hlca_cts)]))
keep_g <- which(!duplicated(rownames(hlca_cts)))   # integer index (avoids array-class dispatch)
hlca_cts <- hlca_cts[keep_g, , drop=FALSE]          # drop 8 dup symbols
hlca <- make_obj(hlca_cts, hl_labs, "hlca")
rm(hlca_s, hlca_cts); gc()

tt("loading lambrechts_2018 Tumors_52k loom")
library(hdf5r)
LP <- file.path(RAW,"lambrechts_2018","3bd84635-f4c8-4d45-9cdc-ea9a03a32fa0","Thienpont_Tumors_52k_v4_R_fixed.loom.gz")
tmp_loom <- tempfile(fileext=".loom"); system2("gzip",c("-dc",shQuote(LP)),stdout=tmp_loom)
lf <- H5File$new(tmp_loom,"r")
lm <- lf[["matrix"]][,]                       # cells x genes
genes <- lf[["row_attrs/Gene"]][]; cn <- lf[["col_attrs/ClusterName"]][]
lf$close_all()
lm <- t(lm); rownames(lm) <- genes            # genes x cells
lamb <- make_obj(as(lm,"CsparseMatrix"), as.character(cn), "lambrechts_2018")

tt("loading loveless_2025 Seurat (double-gzipped rds)")
con <- gzcon(gzfile(file.path(RAW,"loveless_2025","scAtlas.rds.gz"),"rb"))
lv <- readRDS(con); close(con)
lov <- make_obj(GetAssayData(lv, assay="RNA", layer="counts"),
                as.character(lv$Clusters), "loveless_2025")
rm(lv, peng_raw, peng_ann); gc()

# ---- 2/3. merge on shared genes (KEEP layers SPLIT) ----------------------
# NOTE: 4 atlases x ~1.4M cells joined into ONE sparse dgCMatrix overflows the
# 32-bit column-pointer limit (nnz > 2^31-1). Seurat v5 defers layer-joining
# for exactly this reason: NormalizeData/FindVariableFeatures run per-layer, so
# we normalize on FULL genes (correct library sizes), pick HVGs, THEN subset to
# a ~2000-gene panel and only then JoinLayers -- the panel-joined matrix is far
# under 2^31 nnz. No cells are dropped; all 1.4M go through PCA/UMAP/Harmony.
tt("merging on shared genes (split layers)")
objs <- list(peng, hlca, lamb, lov)
shared <- Reduce(intersect, lapply(objs, rownames))
tt(sprintf("  shared genes across 4 atlases: %d", length(shared)))
objs <- lapply(objs, function(o) o[shared,])
merged <- merge(objs[[1]], objs[-1])          # layers stay split: counts.<atlas>
rm(objs, peng, hlca, lamb, lov); gc()
tt(sprintf("  merged (split): %d genes x %d cells", nrow(merged), ncol(merged)))
saveRDS(table(merged$atlas, merged$coarse), file.path(PROC,"harmonization","merged_atlas_x_coarse.rds"))

# ---- normalize (full genes, per-layer) + HVG -----------------------------
merged <- NormalizeData(merged, verbose=FALSE)          # per-layer, full-gene library sizes
merged <- FindVariableFeatures(merged, nfeatures=2000, verbose=FALSE)
hvg <- VariableFeatures(merged)

# canonical markers (known biology) -- fold into the panel so DotPlot has them
canon <- c("EPCAM","KRT19","KRT8","MKI67",         # epithelial/cycling
           "PRSS1","CPA1","CTRB1",                 # acinar
           "INS","CHGA","CHGB",                    # endocrine
           "COL1A1","DCN","PDGFRB","ACTA2","RGS5", # fibro/mural
           "PECAM1","VWF",                         # endothelial
           "CD3D","CD8A","NKG7","IL7R",            # T/NK
           "MS4A1","CD79A","MZB1",                 # B/plasma
           "LYZ","CD68","CD14","FCGR3A",           # myeloid
           "TPSAB1","CPA3",                        # mast
           "SFTPC","SFTPB","AGER")                 # alveolar
canon <- canon[canon %in% rownames(merged)]
panel <- union(hvg, canon)
tt(sprintf("  panel genes (HVG u canonical): %d", length(panel)))

# ---- subset to panel, THEN join (panel-joined sparse is under 2^31 nnz) ---
tt("subsetting to panel + JoinLayers")
sub <- subset(merged, features=panel)
sub <- JoinLayers(sub)
VariableFeatures(sub) <- hvg
rm(merged); gc()

# ---- ScaleData(HVG) + PCA + UMAP pre + Harmony + UMAP post ---------------
sub <- ScaleData(sub, features=hvg, verbose=FALSE)
sub <- RunPCA(sub, features=hvg, npcs=30, verbose=FALSE)
tt("UMAP uncorrected (pre)")
sub <- RunUMAP(sub, dims=1:30, reduction="pca", reduction.name="umap_pre", verbose=FALSE)
tt("Harmony integration (batch=atlas)")
sub <- RunHarmony(sub, group.by.vars="atlas", verbose=FALSE)
tt("UMAP harmony (post)")
sub <- RunUMAP(sub, dims=1:30, reduction="harmony", reduction.name="umap_post", verbose=FALSE)

# CHECKPOINT: all expensive compute (merge/PCA/2xUMAP/Harmony) is done -- persist
# now so a downstream plotting failure never forces a full recompute.
tt("checkpoint: saveRDS computed object before plotting")
dir.create(FIG, recursive=TRUE, showWarnings=TRUE)
saveRDS(sub, "/scratch/rprest2/composite_ref_qc_merged.rds")
tt("  checkpoint written")

pal_atlas <- setNames(brewer.pal(4,"Set1"), sort(unique(sub$atlas)))
ncoarse <- length(unique(sub$coarse))
pal_coarse <- setNames(colorRampPalette(brewer.pal(12,"Paired"))(ncoarse), sort(unique(sub$coarse)))
mkumap <- function(red, grp, pal, ttl){
  # raster=TRUE renders all 1.4M cells onto a high-dpi canvas (no cells dropped)
  DimPlot(sub, reduction=red, group.by=grp, cols=pal, raster=TRUE,
          raster.dpi=c(1200,1200), pt.size=1.2, shuffle=TRUE) +
    ggtitle(ttl) + theme(legend.position="right")
}
save_png <- function(p, name, w=8, h=6){
  ggsave(file.path(FIG,name), p, width=w, height=h, dpi=300, bg="white", create.dir=TRUE)
  tt(paste("wrote", name))
}
save_png(mkumap("umap_pre","atlas", pal_atlas,"Pre-integration (merged, uncorrected) - by reference"),  "umap_pre_by_atlas.png")
save_png(mkumap("umap_pre","coarse",pal_coarse,"Pre-integration (merged, uncorrected) - by cell type"), "umap_pre_by_celltype.png", 9)
save_png(mkumap("umap_post","atlas", pal_atlas,"Post-integration (Harmony) - by reference"),  "umap_post_by_atlas.png")
save_png(mkumap("umap_post","coarse",pal_coarse,"Post-integration (Harmony) - by cell type"), "umap_post_by_celltype.png", 9)

# ---- 5. cross-atlas coarse-type mean-profile correlation -----------------
tt("cross-atlas mean-profile correlation heatmap")
sub$atlas_coarse <- paste(sub$atlas, sub$coarse, sep=" | ")
expr <- GetAssayData(sub, layer="data")[hvg,]
grp  <- sub$atlas_coarse
means <- sapply(split(seq_len(ncol(expr)), grp), function(ix) Matrix::rowMeans(expr[,ix,drop=FALSE]))
# require >=25 cells per group
n_per <- table(grp); means <- means[, names(n_per)[n_per>=25], drop=FALSE]
cc <- cor(means, method="pearson")
ann <- data.frame(coarse=sub("^.*\\| ","",colnames(cc)), atlas=sub(" \\|.*$","",colnames(cc)))
rownames(ann) <- colnames(cc)
png(file.path(FIG,"crossatlas_profile_correlation.png"), width=2400, height=2200, res=200)
pheatmap(cc, annotation_row=ann, annotation_col=ann["coarse"],
         show_rownames=TRUE, show_colnames=FALSE, fontsize_row=6,
         main="Cross-atlas coarse-type mean-profile correlation (HVG, Pearson)")
dev.off(); tt("wrote crossatlas_profile_correlation.png")

# ---- 6. markers (presto-accelerated FindAllMarkers) on coarse ------------
tt("FindAllMarkers on coarse labels (presto)")
Idents(sub) <- "coarse"
mk <- FindAllMarkers(sub, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.5, verbose=FALSE)
write.csv(mk, file.path(PROC,"harmonization","celltype_markers.csv"), row.names=FALSE)
tt(sprintf("  markers: %d rows across %d types", nrow(mk), length(unique(mk$cluster))))
topn <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n=5) %>% pull(gene) %>% unique()
dp <- DotPlot(sub, features=topn, group.by="coarse") + RotatedAxis() +
  ggtitle("Top coarse-type markers (data-derived)") +
  theme(axis.text.x=element_text(size=7))
save_png(dp, "celltype_marker_dotplot.png", 14, 7)

# canonical-marker dotplot (known biology)
dpc <- DotPlot(sub, features=canon, group.by="coarse") + RotatedAxis() +
  ggtitle("Canonical markers by coarse type") + theme(axis.text.x=element_text(size=8))
save_png(dpc, "celltype_canonical_dotplot.png", 14, 7)

saveRDS(sub, "/scratch/rprest2/composite_ref_qc_merged.rds")
tt("DONE")
