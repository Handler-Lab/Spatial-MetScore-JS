#!/usr/bin/env Rscript
# BANKSY clustering — single-slide prototype on PT-1A
# Runs BOTH regimes: lambda=0.2 (cell-typing) and lambda=0.8 (spatial domains)
# Env: spatial_R_scratch (Banksy 1.9.1, Seurat 5.5.1, SeuratWrappers 0.4.0)

suppressPackageStartupMessages({
  library(Banksy)
  library(Seurat)
  library(SeuratWrappers)
  library(ggplot2)
  library(patchwork)
})

set.seed(1000)
options(future.globals.maxSize = 8 * 1024^3)

DATA_DIR <- "/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/data/raw/Samples/PT-1A/"
SAMPLE   <- "PT-1A"
OUT      <- "./out"; dir.create(OUT, showWarnings = FALSE)

# ---- parameters ----
K_GEOM     <- 18      # spatial neighbours (Visium: 6=1st ring, 18=1st+2nd)
NPCS       <- 20
RES        <- 0.5     # Leiden/Louvain resolution
LAMBDAS    <- c(cell_typing = 0.2, domains = 0.8)

# ============================================================
# 1. Load + QC  (match python pipeline: min_counts = 100)
# ============================================================
seu <- Load10X_Spatial(data.dir = DATA_DIR, filename = "filtered_feature_bc_matrix.h5")
cat(sprintf("Loaded %s: %d spots x %d genes (pre-QC)\n", SAMPLE, ncol(seu), nrow(seu)))

seu <- subset(seu, subset = nCount_Spatial >= 100)
cat(sprintf("After min_counts>=100 filter: %d spots\n", ncol(seu)))

# ============================================================
# 2. Normalise + HVGs
# ============================================================
seu <- NormalizeData(seu, normalization.method = "LogNormalize",
                     scale.factor = 1e4, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
HVG <- VariableFeatures(seu)                       # fixed gene set, reused by both regimes
cat(sprintf("HVGs selected: %d\n", length(HVG)))

# ============================================================
# 3. BANKSY per lambda regime
# ============================================================
run_banksy <- function(seu, lambda, tag) {
  cat(sprintf("\n=== BANKSY %s (lambda=%.1f) ===\n", tag, lambda))
  aname <- paste0("BANKSY_", tag)
  rname <- paste0("pca_", tag)
  uname <- paste0("umap_", tag)
  cname <- paste0("banksy_", tag)

  # Reset to the base expression assay; the previous regime left its own
  # BANKSY_* assay as default, which breaks variable-feature lookup here.
  DefaultAssay(seu) <- "Spatial"
  # use_agf=FALSE: the AGF (complex-matrix) path hits a Banksy+data.table+Matrix
  # incompatibility in this env. Core BANKSY (mean-neighbour features) is used.
  # Resetting DefaultAssay above restores Spatial's 2000 VariableFeatures,
  # so features="variable" resolves correctly for each regime.
  seu <- RunBanksy(seu, lambda = lambda, assay = "Spatial", slot = "data",
                   features = "variable", k_geom = K_GEOM, use_agf = FALSE,
                   assay_name = aname, verbose = FALSE)
  seu <- RunPCA(seu, assay = aname, features = rownames(seu[[aname]]),
                npcs = NPCS, reduction.name = rname, verbose = FALSE)
  seu <- RunUMAP(seu, reduction = rname, dims = 1:NPCS,
                 reduction.name = uname, verbose = FALSE)
  seu <- FindNeighbors(seu, reduction = rname, dims = 1:NPCS, verbose = FALSE)
  seu <- FindClusters(seu, resolution = RES, cluster.name = cname, verbose = FALSE)

  n_clust <- length(unique(seu[[cname]][, 1]))
  cat(sprintf("%s: %d clusters\n", tag, n_clust))
  print(table(seu[[cname]][, 1]))
  seu
}

for (i in seq_along(LAMBDAS)) {
  seu <- run_banksy(seu, LAMBDAS[i], names(LAMBDAS)[i])
}

# ============================================================
# 4. Figures
# ============================================================
plot_one <- function(seu, tag, lambda) {
  cname <- paste0("banksy_", tag)
  uname <- paste0("umap_", tag)
  Idents(seu) <- cname
  sp <- SpatialDimPlot(seu, group.by = cname, pt.size.factor = 1.6, label = TRUE,
                       label.size = 3) +
        ggtitle(sprintf("%s — spatial (lambda=%.1f, %s)", SAMPLE, lambda, tag))
  um <- DimPlot(seu, reduction = uname, group.by = cname, label = TRUE) +
        ggtitle(sprintf("%s — UMAP (lambda=%.1f)", SAMPLE, lambda))
  list(spatial = sp, umap = um)
}

p02 <- plot_one(seu, "cell_typing", 0.2)
p08 <- plot_one(seu, "domains",     0.8)

# H&E reference (tissue image with spots faint)
he <- tryCatch(
  SpatialFeaturePlot(seu, features = "nCount_Spatial", alpha = c(0.1, 0.1)) +
    ggtitle(paste0(SAMPLE, " — H&E")),
  error = function(e) p08$spatial + ggtitle(paste0(SAMPLE, " — (H&E unavailable)")))

ggsave(file.path(OUT, "PT-1A_banksy_spatial_l02_celltyping.png"),
       p02$spatial, width = 8, height = 7, dpi = 200)
ggsave(file.path(OUT, "PT-1A_banksy_spatial_l08_domains.png"),
       p08$spatial, width = 8, height = 7, dpi = 200)
ggsave(file.path(OUT, "PT-1A_banksy_umap_l02.png"), p02$umap, width = 7, height = 6, dpi = 200)
ggsave(file.path(OUT, "PT-1A_banksy_umap_l08.png"), p08$umap, width = 7, height = 6, dpi = 200)

# combined overview: H&E | domains-spatial | celltyping-spatial
combo <- he + p08$spatial + p02$spatial + plot_layout(ncol = 3)
ggsave(file.path(OUT, "PT-1A_banksy_overview.png"), combo, width = 22, height = 7, dpi = 200)

# ============================================================
# 5. Save outputs
# ============================================================
meta <- seu@meta.data[, c("nCount_Spatial", "nFeature_Spatial",
                          "banksy_cell_typing", "banksy_domains")]
meta$barcode <- rownames(meta)
write.csv(meta, file.path(OUT, "PT-1A_banksy_clusters.csv"), row.names = FALSE)

saveRDS(seu, file.path(OUT, "PT-1A_banksy.rds"))

# cross-tab of the two regimes
cat("\n=== cross-tab: cell_typing (rows) x domains (cols) ===\n")
print(table(seu$banksy_cell_typing, seu$banksy_domains))

cat("\nDONE. Outputs in", OUT, "\n")
