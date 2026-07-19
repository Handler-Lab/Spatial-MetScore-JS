#!/usr/bin/env Rscript
# Unified CopyKAT pipeline for Spatial-MetScore: per-sample CNA calling,
# pooled/combined multi-slide calling, and cohort-summary aggregation,
# all in one script.
#
# Subcommands
# -----------
#   sample     <SAMPLE_ID>   Run CopyKAT on one Visium slide (used by the
#                            array job, one task per sample).
#   combined                 Pool ALL slides into one Seurat/CopyKAT run so
#                            low-purity slides supply a shared diploid
#                            baseline (mitigates the per-slide reference-free
#                            baseline failure -- see copykat_puree_inversion_diag).
#   summarize                Aggregate every *_copykat_prediction.csv under
#                            --outdir into one cohort summary CSV, merged
#                            with patient_metadata.csv.
#
# Heavy modes (sample/combined) are meant to run via SBATCH
# (see scripts/slurm/copykat.sh), never on the login node.
#
# Usage (from repo root):
#   Rscript scripts/copykat_pipeline.R sample <SAMPLE_ID> [--outdir DIR] [--distance TYPE] [--ncores N]
#   Rscript scripts/copykat_pipeline.R combined [--outdir DIR] [--ncores N]
#   Rscript scripts/copykat_pipeline.R summarize [--outdir DIR] [--out CSV]

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(copykat)
  library(ggplot2)
})
set.seed(1000)

# ----------------------------------------------------------------------------
# Paths / constants (shared across modes)
# ----------------------------------------------------------------------------
DATA_ROOT       <- "/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/data/raw/Samples"
PATIENT_META    <- "/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/data/raw/patient_metadata.csv"
SAMPLES_MANIFEST <- "scripts/copykat_samples.txt"     # relative to repo root
PROC_ROOT_DEFAULT <- "/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/data/processed/copykat"
FIG_ROOT_DEFAULT  <- "figures/copykat"                # relative to repo root
CALL_COLORS <- c(aneuploid = "#D7263D", diploid = "#1B98E0", "not.defined" = "grey80")

samples_list <- function() {
  s <- readLines(SAMPLES_MANIFEST)
  s[nzchar(s)]
}

# robust classification -- handles CopyKAT's low-confidence label variants
# (e.g. "c1:diploid:low.conf", "c2:aneuploid:low.conf") by substring match
# rather than exact match, so these aren't silently dropped as NA.
classify_calls <- function(calls) {
  ifelse(grepl("aneuploid", calls), "aneuploid",
  ifelse(grepl("diploid",   calls), "diploid", "not.defined"))
}

# ----------------------------------------------------------------------------
# arg parsing helper: pull "--flag value" pairs out of a trailing args vector
# ----------------------------------------------------------------------------
get_opt <- function(args, flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) && i[1] < length(args)) return(args[i[1] + 1])
  default
}

# ============================================================================
# Mode 1: sample -- CopyKAT on a single Visium slide
# ============================================================================
run_sample <- function(args) {
  if (length(args) < 1) stop("Usage: copykat_pipeline.R sample <SAMPLE_ID> [--outdir DIR] [--distance TYPE] [--ncores N]")
  SAMPLE   <- args[1]
  PROC_ROOT <- get_opt(args, "--outdir", PROC_ROOT_DEFAULT)
  DISTANCE  <- get_opt(args, "--distance", "euclidean")
  NCORES    <- as.integer(get_opt(args, "--ncores", "8"))

  DATA_DIR <- file.path(DATA_ROOT, SAMPLE)
  FIG_DIR  <- file.path(FIG_ROOT_DEFAULT, SAMPLE)   # relative to repo root
  DATA_OUT <- file.path(PROC_ROOT, SAMPLE)
  WORK     <- file.path(tempdir(), paste0("copykat_", SAMPLE))
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
            distance = DISTANCE, genome = "hg20",
            n.cores = NCORES, plot.genes = "TRUE", output.seg = "FALSE"),
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
  tryCatch({
    sp <- SpatialDimPlot(seu, group.by = "copykat", pt.size.factor = 1.6) +
          scale_fill_manual(values = CALL_COLORS, name = "CopyKAT") +
          ggtitle(sprintf("%s - CopyKAT aneuploid vs diploid", SAMPLE))
    ggsave(file.path(FIG_DIR, paste0(SAMPLE, "_copykat_spatial.png")), sp,
           width = 8, height = 7, dpi = 200)
  }, error = function(e) message(sprintf("[%s] spatial plot failed: %s", SAMPLE, conditionMessage(e))))

  # ---- 4. One-line per-sample summary ----
  call_class <- classify_calls(pred$copykat.pred)
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
}

# ============================================================================
# Mode 2: combined -- pooled CopyKAT across ALL slides in one object
# ============================================================================
# Rationale: reference-free CopyKAT inverts calls on high-purity slides that
# lack an internal low-CNA "normal" cluster. Pooling every slide makes the
# abundant stromal/immune spots from lower-purity slides available as a shared
# diploid baseline for the whole cohort.
#
# Memory design: we never hold all 54 per-sample Seurat objects (with images,
# metadata, etc.) in memory simultaneously. Step 1 reads only the raw sparse
# counts matrix per sample (via Read10X_h5), cbinds them into one big sparse
# matrix, and discards each per-sample matrix as it's folded in. Step 3
# re-loads each sample fresh from disk (cheap, single-slide) only when needed
# for its spatial plot -- so at no point do we hold more than one full Seurat
# object plus the single pooled matrix/prediction in memory.
run_combined <- function(args) {
  PROC_ROOT <- get_opt(args, "--outdir", "/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/data/processed/copykat_combined")
  FIG_ROOT  <- "figures/copykat_combined"    # relative to repo root
  NCORES    <- as.integer(get_opt(args, "--ncores", "2"))  # forced low by default: bounds per-fork memory duplication in copykat step 3
  dir.create(PROC_ROOT, recursive = TRUE, showWarnings = FALSE)

  SAMPLES <- samples_list()
  cat(sprintf("Loading %d samples (sparse counts only, QC nCount>=100)...\n", length(SAMPLES)))

  # ---- 1. Read raw sparse counts per slide, QC, cbind into ONE sparse matrix ----
  mat_list <- vector("list", length(SAMPLES)); names(mat_list) <- SAMPLES
  for (s in SAMPLES) {
    m <- Read10X_h5(file.path(DATA_ROOT, s, "filtered_feature_bc_matrix.h5"))
    keep <- Matrix::colSums(m) >= 100
    m <- m[, keep, drop = FALSE]
    colnames(m) <- paste0(s, "_", colnames(m))
    mat_list[[s]] <- m
    cat(sprintf("  %-12s %d spots (post-QC)\n", s, ncol(m)))
  }

  # sanity check: all slides must share the same gene set (standard 10x human ref)
  gene_sets <- lapply(mat_list, rownames)
  stopifnot(all(sapply(gene_sets, identical, gene_sets[[1]])))

  raw <- do.call(cbind, mat_list)
  rm(mat_list); gc()
  cat(sprintf("Combined sparse matrix: %d genes x %d spots\n", nrow(raw), ncol(raw)))
  raw <- as.matrix(raw)   # copykat requires a dense matrix internally anyway
  gc()

  # ---- 2. Run CopyKAT once on the pooled matrix ----
  WORK <- file.path(tempdir(), "copykat_combined")
  dir.create(WORK, recursive = TRUE, showWarnings = FALSE)
  oldwd <- getwd(); setwd(WORK)
  ck <- copykat(rawmat = raw, id.type = "S", sam.name = "combined",
                ngene.chr = 5, win.size = 25, KS.cut = 0.1,
                distance = "euclidean", genome = "hg20",
                n.cores = NCORES, plot.genes = "FALSE", output.seg = "FALSE")
  setwd(oldwd)
  rm(raw); gc()

  pred <- data.frame(ck$prediction)
  write.csv(pred, file.path(PROC_ROOT, "combined_copykat_prediction.csv"), row.names = FALSE)
  saveRDS(ck, file.path(PROC_ROOT, "combined_copykat.rds"))
  rm(ck); gc()

  pred$class <- classify_calls(pred$copykat.pred)

  # ---- 3. Per-slide visualization + summary from the SHARED prediction ----
  summ_rows <- list()
  for (s in SAMPLES) {
    seu <- Load10X_Spatial(data.dir = file.path(DATA_ROOT, s), filename = "filtered_feature_bc_matrix.h5")
    seu <- subset(seu, subset = nCount_Spatial >= 100)
    pref <- paste0(s, "_", colnames(seu))
    m    <- match(pref, pred$cell.names)
    seu$copykat <- "not.defined"
    seu$copykat[!is.na(m)] <- pred$class[m[!is.na(m)]]

    fig_dir <- file.path(FIG_ROOT, s)
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    tryCatch({
      sp <- SpatialDimPlot(seu, group.by = "copykat", pt.size.factor = 1.6) +
            scale_fill_manual(values = CALL_COLORS, name = "CopyKAT") +
            ggtitle(sprintf("%s - pooled CopyKAT", s))
      ggsave(file.path(fig_dir, paste0(s, "_copykat_spatial.png")), sp, width = 8, height = 7, dpi = 200)
    }, error = function(e) message(sprintf("[%s] spatial plot failed: %s", s, conditionMessage(e))))

    cls <- factor(seu$copykat, levels = c("aneuploid", "diploid", "not.defined"))
    tb  <- table(cls)
    n_classified <- sum(!is.na(m))
    summ_rows[[s]] <- data.frame(
      sample = s, n_spots_postQC = ncol(seu), n_classified = n_classified,
      n_aneuploid = as.integer(tb["aneuploid"]),
      n_diploid   = as.integer(tb["diploid"]),
      n_not_defined = as.integer(tb["not.defined"]))
    rm(seu); gc()
  }

  summ <- do.call(rbind, summ_rows)
  summ$frac_aneuploid <- round(summ$n_aneuploid / summ$n_classified, 4)
  write.csv(summ, file.path(PROC_ROOT, "copykat_combined_summary_all_samples.csv"), row.names = FALSE)
  cat("DONE. Per-slide summary:\n"); print(summ)
}

# ============================================================================
# Mode 3: summarize -- aggregate per-sample predictions into cohort CSV
# ============================================================================
run_summarize <- function(args) {
  PROC_ROOT <- get_opt(args, "--outdir", PROC_ROOT_DEFAULT)
  OUT       <- get_opt(args, "--out", file.path(PROC_ROOT, "copykat_summary_all_samples.csv"))

  pred_files <- list.files(PROC_ROOT, pattern = "_copykat_prediction\\.csv$", recursive = TRUE, full.names = TRUE)
  if (!length(pred_files)) stop("No per-sample prediction CSVs found under ", PROC_ROOT)

  summ_files <- list.files(PROC_ROOT, pattern = "_copykat_summary\\.csv$", recursive = TRUE, full.names = TRUE)
  postqc_by_sample <- list()
  for (f in summ_files) {
    s <- read.csv(f, stringsAsFactors = FALSE)
    postqc_by_sample[[s$sample[1]]] <- s$n_spots_postQC[1]
  }

  rows <- do.call(rbind, lapply(pred_files, function(f) {
    sample <- sub("_copykat_prediction\\.csv$", "", basename(f))
    pred <- read.csv(f, stringsAsFactors = FALSE)
    cls <- classify_calls(pred$copykat.pred)
    tb <- table(factor(cls, levels = c("aneuploid", "diploid", "not.defined")))
    data.frame(sample = sample,
               n_spots_postQC = ifelse(is.null(postqc_by_sample[[sample]]), NA, postqc_by_sample[[sample]]),
               n_classified = nrow(pred),
               n_aneuploid = as.integer(tb["aneuploid"]),
               n_diploid = as.integer(tb["diploid"]),
               n_not_defined = as.integer(tb["not.defined"]),
               status = "ok")
  }))

  meta <- tryCatch(read.csv(PATIENT_META, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(meta)) {
    key <- intersect(c("Sample_ID", "sample", "sample_id"), names(meta))
    if (length(key)) {
      meta$.k <- toupper(gsub("[^A-Za-z0-9]", "", meta[[key[1]]]))
      rows$.k <- toupper(gsub("[^A-Za-z0-9]", "", rows$sample))
      keep <- intersect(c("Tissue_Type", "patient", "tissue"), names(meta))
      rows <- merge(rows, meta[, c(".k", keep)], by = ".k", all.x = TRUE)
      rows$.k <- NULL
    }
  }

  rows <- rows[order(rows$sample), ]
  rows$frac_aneuploid <- round(rows$n_aneuploid / rows$n_classified, 4)
  stopifnot(all(rows$n_aneuploid + rows$n_diploid + rows$n_not_defined == rows$n_classified))
  write.csv(rows, OUT, row.names = FALSE)
  cat(sprintf("Wrote %s with %d samples\n", OUT, nrow(rows)))
  print(rows)
}

# ============================================================================
# CLI
# ============================================================================
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) < 1) {
  stop("Usage: copykat_pipeline.R {sample <SAMPLE_ID>|combined|summarize} [options]")
}
mode <- argv[1]
rest <- if (length(argv) > 1) argv[-1] else character(0)

if (mode == "sample") {
  run_sample(rest)
} else if (mode == "combined") {
  run_combined(rest)
} else if (mode == "summarize") {
  run_summarize(rest)
} else {
  stop(sprintf("Unknown mode '%s'. Expected sample|combined|summarize", mode))
}
