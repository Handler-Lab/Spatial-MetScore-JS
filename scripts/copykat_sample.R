#!/usr/bin/env Rscript
# CopyKAT on a single Visium sample.
# Usage (from repo root): Rscript scripts/copykat_sample.R <SAMPLE_ID>
# Keeps only deliverables (spatial PNG, CNA heatmap PDF, prediction CSV/RDS,
# one-line summary); deletes CopyKAT's large intermediate bin/gene-by-cell matrices.

suppressPackageStartupMessages({
  library(Seurat)
  library(copykat)
  library(ggplot2)
})
set.seed(1000)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: copykat_sample.R <SAMPLE_ID>")
SAMPLE <- args[1]

DATA_ROOT <- "/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/data/raw/Samples"
PROC_ROOT <- "/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/data/processed/copykat_pearson"
DATA_DIR  <- file.path(DATA_ROOT, SAMPLE)
FIG_DIR   <- file.path("figures/copykat_pearson", SAMPLE)   # relative to repo root
DATA_OUT  <- file.path(PROC_ROOT, SAMPLE)
WORK      <- file.path(tempdir(), paste0("copykat_", SAMPLE))
for (d in c(FIG_DIR, DATA_OUT, WORK)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Load + QC (match BANKSY / metscore pipeline: min_counts >= 100) ----
seu <- Load10X_Spatial(data.dir = DATA_DIR, filename = "filtered_feature_bc_matrix.h5")
cat(sprintf("[%s] Loaded %d spots x %d genes (pre-QC)\n", SAMPLE, ncol(seu), nrow(seu)))
seu <- subset(seu, subset = nCount_Spatial >= 100)
n_postqc <- ncol(seu)
cat(sprintf("[%s] After min_counts>=100 filter: %d spots\n", SAMPLE, n_postqc))
raw <- as.matrix(GetAssayData(seu, assay = "Spatial", layer = "counts"))

# ---- 2. Run CopyKAT inside a scratch dir so intermediates are easy to drop ----
oldwd <- getwd(); setwd(WORK)
ck <- tryCatch(
  copykat(rawmat = raw, id.type = "S", sam.name = SAMPLE,
          ngene.chr = 5, win.size = 25, KS.cut = 0.1,
          distance = "euclidean", genome = "hg20",
          n.cores = 8, plot.genes = "TRUE", output.seg = "FALSE"),
  error = function(e) { message(sprintf("[%s] CopyKAT ERROR: %s", SAMPLE, conditionMessage(e))); NULL }
)
setwd(oldwd)

if (is.null(ck) || is.null(ck$prediction)) {
  summ <- data.frame(sample = SAMPLE, n_spots_postQC = n_postqc, n_classified = NA,
                     n_aneuploid = NA, n_diploid = NA, n_not_defined = NA, status = "failed")
  write.csv(summ, file.path(DATA_OUT, paste0(SAMPLE, "_copykat_summary.csv")), row.names = FALSE)
  unlink(WORK, recursive = TRUE)
  stop(sprintf("[%s] CopyKAT produced no prediction.", SAMPLE))
}

pred <- data.frame(ck$prediction)
write.csv(pred, file.path(DATA_OUT, paste0(SAMPLE, "_copykat_prediction.csv")), row.names = FALSE)
saveRDS(ck, file.path(DATA_OUT, paste0(SAMPLE, "_copykat.rds")))

# keep only the CNA heatmap (pdf/jpeg); drop giant *_results_*_by_cell.txt matrices
hm <- list.files(WORK, pattern = "heatmap\\.(pdf|jpeg|jpg|png)$", full.names = TRUE)
if (length(hm)) file.copy(hm, FIG_DIR, overwrite = TRUE)

# ---- 3. Spatial overlay of aneuploid / diploid on tissue ----
seu$copykat <- "not.defined"
m <- match(colnames(seu), pred$cell.names)
seu$copykat[!is.na(m)] <- pred$copykat.pred[m[!is.na(m)]]
cols <- c(aneuploid = "#D7263D", diploid = "#1B98E0", "not.defined" = "grey80")
tryCatch({
  sp <- SpatialDimPlot(seu, group.by = "copykat", pt.size.factor = 1.6) +
        scale_fill_manual(values = cols, name = "CopyKAT") +
        ggtitle(sprintf("%s - CopyKAT aneuploid vs diploid", SAMPLE))
  ggsave(file.path(FIG_DIR, paste0(SAMPLE, "_copykat_spatial.png")), sp,
         width = 8, height = 7, dpi = 200)
}, error = function(e) message(sprintf("[%s] spatial plot failed: %s", SAMPLE, conditionMessage(e))))

# ---- 4. One-line per-sample summary ----
# CopyKAT sometimes falls back to low-confidence cluster labels
# (e.g. "c1:diploid:low.conf", "c2:aneuploid:low.conf") instead of the
# plain "aneuploid"/"diploid" strings -- classify robustly by substring
# rather than exact match so these aren't silently dropped as NA.
call_class <- ifelse(grepl("aneuploid", pred$copykat.pred), "aneuploid",
              ifelse(grepl("diploid", pred$copykat.pred), "diploid",
                     "not.defined"))
tb <- table(factor(call_class, levels = c("aneuploid", "diploid", "not.defined")))
summ <- data.frame(sample = SAMPLE, n_spots_postQC = n_postqc, n_classified = nrow(pred),
                   n_aneuploid = as.integer(tb["aneuploid"]),
                   n_diploid   = as.integer(tb["diploid"]),
                   n_not_defined = as.integer(tb["not.defined"]),
                   status = "ok")
write.csv(summ, file.path(DATA_OUT, paste0(SAMPLE, "_copykat_summary.csv")), row.names = FALSE)

# ---- 5. Drop all intermediates ----
unlink(WORK, recursive = TRUE)
cat(sprintf("[%s] DONE: %d spots -> %d aneuploid / %d diploid / %d not.defined\n",
            SAMPLE, n_postqc, tb["aneuploid"], tb["diploid"], tb["not.defined"]))
