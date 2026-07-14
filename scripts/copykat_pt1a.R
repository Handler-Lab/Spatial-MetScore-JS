#!/usr/bin/env Rscript
# CopyKAT on PT-1A -- infer per-spot CNA and classify aneuploid (tumor) vs diploid (normal)
# Env: spatial_R_scratch (copykat 1.2.5, Seurat 5.5.1)
# Run from repo root, e.g.: Rscript scripts/copykat_pt1a.R

suppressPackageStartupMessages({
  library(Seurat)
  library(copykat)
  library(ggplot2)
  library(patchwork)
})

set.seed(1000)

DATA_DIR <- "/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/data/raw/Samples/PT-1A/"
SAMPLE   <- "PT-1A"
FIG_DIR  <- "figures/copykat/PT-1A"
DATA_OUT <- "/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/data/processed/copykat/PT-1A"
OUT      <- "./copykat_work"   # copykat writes many intermediate files here
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_OUT, showWarnings = FALSE, recursive = TRUE)
N_CORES  <- 8

# ============================================================
# 1. Load + QC (match BANKSY: min_counts >= 100)
# ============================================================
seu <- Load10X_Spatial(data.dir = DATA_DIR, filename = "filtered_feature_bc_matrix.h5")
cat(sprintf("Loaded %s: %d spots x %d genes (pre-QC)\n", SAMPLE, ncol(seu), nrow(seu)))
seu <- subset(seu, subset = nCount_Spatial >= 100)
cat(sprintf("After min_counts>=100 filter: %d spots\n", ncol(seu)))

# raw UMI counts, genes x spots (CopyKAT needs raw, not normalized)
raw <- GetAssayData(seu, assay = "Spatial", layer = "counts")
raw <- as.matrix(raw)
cat(sprintf("Raw matrix for CopyKAT: %d genes x %d spots\n", nrow(raw), ncol(raw)))

# ============================================================
# 2. Run CopyKAT
#    id.type="S" (gene symbols), genome hg20 (human)
#    no known normal cells supplied -> CopyKAT defines baseline
# ============================================================
oldwd <- getwd()
setwd(OUT)   # copykat writes its files into cwd
ck <- copykat(rawmat = raw, id.type = "S", sam.name = SAMPLE,
              ngene.chr = 5, win.size = 25, KS.cut = 0.1,
              distance = "euclidean", genome = "hg20",
              n.cores = N_CORES, plot.genes = "TRUE", output.seg = "FALSE")
setwd(oldwd)

pred <- data.frame(ck$prediction)
cat("\n=== CopyKAT prediction table ===\n")
print(table(pred$copykat.pred))
saveRDS(ck, file.path(DATA_OUT, paste0(SAMPLE, "_copykat.rds")))
write.csv(pred, file.path(DATA_OUT, paste0(SAMPLE, "_copykat_prediction.csv")), row.names = FALSE)

# copy copykat's own diagnostic files (CNA heatmap pdf, etc.) into figures/
ck_files <- list.files(OUT, pattern = "\\.(pdf|png|txt)$", full.names = TRUE)
file.copy(ck_files, FIG_DIR, overwrite = TRUE)

# ============================================================
# 3. Spatial overlay of aneuploid / diploid on tissue
# ============================================================
seu$copykat <- "not.defined"
m <- match(colnames(seu), pred$cell.names)
seu$copykat[!is.na(m)] <- pred$copykat.pred[m[!is.na(m)]]

cols <- c(aneuploid = "#D7263D", diploid = "#1B98E0", "not.defined" = "grey80")
sp <- SpatialDimPlot(seu, group.by = "copykat", pt.size.factor = 1.6) +
      scale_fill_manual(values = cols, name = "CopyKAT") +
      ggtitle(sprintf("%s -- CopyKAT aneuploid vs diploid", SAMPLE))
ggsave(file.path(FIG_DIR, paste0(SAMPLE, "_copykat_spatial.png")),
       sp, width = 8, height = 7, dpi = 200)

cat("\nDONE. Figures in", FIG_DIR, " | data in", DATA_OUT, "\n")
