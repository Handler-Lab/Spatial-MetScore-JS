#!/usr/bin/env Rscript
# Aggregate per-sample CopyKAT predictions into one summary CSV.
# Classifies directly from each sample's *_copykat_prediction.csv (robust to
# CopyKAT's low-confidence label variants, e.g. "c1:diploid:low.conf"),
# rather than trusting the per-sample summary's prior (buggy) exact-match tally.
# Usage (from repo root): Rscript scripts/copykat_summarize.R
PROC_ROOT <- "/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/data/processed/copykat"
META      <- "/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/data/raw/patient_metadata.csv"
OUT       <- file.path(PROC_ROOT, "copykat_summary_all_samples.csv")

pred_files <- list.files(PROC_ROOT, pattern = "_copykat_prediction\\.csv$",
                         recursive = TRUE, full.names = TRUE)
if (!length(pred_files)) stop("No per-sample prediction CSVs found under ", PROC_ROOT)

summ_files <- list.files(PROC_ROOT, pattern = "_copykat_summary\\.csv$",
                         recursive = TRUE, full.names = TRUE)
postqc <- setNames(rep(NA_integer_, length(summ_files)), summ_files)
postqc_by_sample <- list()
for (f in summ_files) {
  s <- read.csv(f, stringsAsFactors = FALSE)
  postqc_by_sample[[s$sample[1]]] <- s$n_spots_postQC[1]
}

rows <- do.call(rbind, lapply(pred_files, function(f) {
  sample <- sub("_copykat_prediction\\.csv$", "", basename(f))
  pred <- read.csv(f, stringsAsFactors = FALSE)
  cls <- ifelse(grepl("aneuploid", pred$copykat.pred), "aneuploid",
         ifelse(grepl("diploid", pred$copykat.pred), "diploid", "not.defined"))
  tb <- table(factor(cls, levels = c("aneuploid", "diploid", "not.defined")))
  data.frame(sample = sample,
             n_spots_postQC = ifelse(is.null(postqc_by_sample[[sample]]), NA, postqc_by_sample[[sample]]),
             n_classified = nrow(pred),
             n_aneuploid = as.integer(tb["aneuploid"]),
             n_diploid = as.integer(tb["diploid"]),
             n_not_defined = as.integer(tb["not.defined"]),
             status = "ok")
}))

meta <- tryCatch(read.csv(META, stringsAsFactors = FALSE), error = function(e) NULL)
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
