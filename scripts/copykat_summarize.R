#!/usr/bin/env Rscript
# Aggregate per-sample CopyKAT summaries into one CSV.
# Usage (from repo root): Rscript scripts/copykat_summarize.R
PROC_ROOT <- "/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/data/processed/copykat"
META      <- "/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/data/raw/patient_metadata.csv"
OUT       <- file.path(PROC_ROOT, "copykat_summary_all_samples.csv")

files <- list.files(PROC_ROOT, pattern = "_copykat_summary\\.csv$",
                    recursive = TRUE, full.names = TRUE)
if (!length(files)) stop("No per-sample summary CSVs found under ", PROC_ROOT)
rows <- do.call(rbind, lapply(files, function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  if (!"status" %in% names(df)) df$status <- "ok"
  df
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
write.csv(rows, OUT, row.names = FALSE)
cat(sprintf("Wrote %s with %d samples\n", OUT, nrow(rows)))
print(rows)
