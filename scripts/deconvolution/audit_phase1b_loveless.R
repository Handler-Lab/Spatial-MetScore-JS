#!/usr/bin/env Rscript
## Phase 1b: discover loveless_2025 scAtlas.rds.gz (32G). Report structure only.
suppressWarnings(suppressMessages({library(methods)}))
`%||%` <- function(a,b) if (is.null(a)) b else a
FP <- "/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/data/raw/scRNAseq/loveless_2025/scAtlas.rds.gz"
sink("loveless_discovery.txt", split=TRUE)
cat("=== loading", FP, "===\n"); cat("start:", format(Sys.time()), "\n")
cat("NOTE: file is double-gzipped (outer gzip wraps an inner gzipped RDS)\n")
t0 <- Sys.time()
## outer gzfile() decompresses layer 1 -> inner gzip stream; gzcon() decompresses layer 2 -> RDS bytes
obj <- tryCatch({
  con <- gzcon(gzfile(FP, "rb"))
  on.exit(close(con), add=TRUE)
  readRDS(con)
}, error=function(e) {
  cat("nested-connection read failed:", conditionMessage(e), "\n -> falling back to decompress-to-scratch\n")
  tmp <- file.path("/scratch/rprest2", paste0("loveless_inner_", Sys.getpid(), ".rds"))
  system2("gzip", c("-dc", shQuote(FP)), stdout=tmp)
  o <- readRDS(tmp); unlink(tmp); o
})
cat("load seconds:", round(as.numeric(Sys.time()-t0, units="secs")), "\n")
cat("class:", paste(class(obj), collapse=","), "\n\n")

int_check <- function(x) {
  s <- tryCatch(as.numeric(x[1:min(2000, length(x))]), error=function(e) NA)
  if (all(is.na(s))) return(NA)
  isTRUE(all.equal(s, round(s)))
}

if (inherits(obj, "Seurat")) {
  suppressWarnings(suppressMessages(library(Seurat)))
  cat("Seurat version of object; dims (genes x cells):", paste(dim(obj), collapse=" x "), "\n")
  cat("assays:", paste(names(obj@assays), collapse=", "), "\n")
  cat("DefaultAssay:", DefaultAssay(obj), "\n")
  for (a in names(obj@assays)) {
    as_obj <- obj@assays[[a]]
    slots <- tryCatch(slotNames(as_obj), error=function(e) character(0))
    cat(sprintf("  assay '%s' class=%s slots=%s\n", a, class(as_obj)[1], paste(slots, collapse=",")))
    cm <- tryCatch(SeuratObject::GetAssayData(obj, assay=a, layer="counts"), error=function(e) NULL)
    if (!is.null(cm) && length(cm)>0) {
      cat(sprintf("    counts layer dim=%s integer=%s\n", paste(dim(cm),collapse="x"), int_check(cm@x %||% cm[,1])))
    }
  }
  md <- obj@meta.data
  cat("\nmeta.data columns (", ncol(md), "):\n", sep=""); print(colnames(md))
  cat("\ngene id samples (rownames):\n"); print(head(rownames(obj), 8))
} else if (inherits(obj, "SingleCellExperiment") || inherits(obj, "SummarizedExperiment")) {
  suppressWarnings(suppressMessages(library(SingleCellExperiment)))
  cat("SCE dims (genes x cells):", paste(dim(obj), collapse=" x "), "\n")
  cat("assayNames:", paste(SummarizedExperiment::assayNames(obj), collapse=", "), "\n")
  for (a in SummarizedExperiment::assayNames(obj)) {
    m <- SummarizedExperiment::assay(obj, a)
    cat(sprintf("  assay '%s' dim=%s integer=%s\n", a, paste(dim(m),collapse="x"),
        int_check(if (inherits(m,"dgCMatrix")) m@x else m[,1])))
  }
  md <- as.data.frame(SummarizedExperiment::colData(obj))
  cat("\ncolData columns (", ncol(md), "):\n", sep=""); print(colnames(md))
  cat("\ngene id samples (rownames):\n"); print(head(rownames(obj), 8))
} else if (is.list(obj)) {
  cat("Object is a list. names:\n"); print(names(obj))
  cat("element classes:\n"); print(sapply(obj, function(e) class(e)[1]))
  md <- NULL
} else {
  cat("Unrecognized structure. str (max 2 levels):\n"); str(obj, max.level=2, list.len=40)
  md <- NULL
}

## value_counts on candidate categorical meta columns
if (exists("md") && !is.null(md)) {
  cat("\n=== candidate label columns value_counts ===\n")
  cand <- grep("cell|type|annot|ident|cluster|lineage|celltype|CellType|tissue|organ|site|compartment|major|minor|subset",
               colnames(md), ignore.case=TRUE, value=TRUE)
  cat("candidate columns:\n"); print(cand)
  for (col in cand) {
    v <- md[[col]]
    if (is.factor(v)) v <- as.character(v)
    if (is.character(v) || is.factor(v) || (is.numeric(v) && length(unique(v))<=200)) {
      nu <- length(unique(v))
      if (nu >= 1 && nu <= 300) {
        cat(sprintf("\n-- %s (n_unique=%d) --\n", col, nu))
        tb <- sort(table(v), decreasing=TRUE)
        print(head(tb, 80))
      }
    }
  }
}
cat("\n=== DONE", format(Sys.time()), "===\n")
sink()
