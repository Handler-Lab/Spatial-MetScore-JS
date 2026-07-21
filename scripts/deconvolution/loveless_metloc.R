#!/usr/bin/env Rscript
# Probe loveless_2025 for liver-met content + any hepatocyte-like label.
suppressWarnings(suppressMessages({library(Seurat)}))
FP <- "/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/data/raw/scRNAseq/loveless_2025/scAtlas.rds.gz"
cat("=== loading (double-gzipped) ===\n"); t0<-Sys.time()
con <- gzcon(gzfile(FP,"rb")); obj <- readRDS(con); close(con)
cat("load seconds:", round(as.numeric(Sys.time()-t0,units="secs")),"\n\n")
md <- obj@meta.data

show <- function(col){
  if(!col %in% colnames(md)){cat("--",col,"NOT PRESENT --\n\n");return(invisible())}
  cat("== ",col," (n_unique=",length(unique(md[[col]])),") ==\n",sep="")
  print(sort(table(md[[col]], useNA="ifany"), decreasing=TRUE)); cat("\n")
}
for(col in c("If.metastatic..location","DiseaseState","Treatment","Name","GSE.SRA..Study.","Study..Citation..PMID."))
  show(col)

# cross-tab: which Clusters appear in liver-met samples
loc <- md[["If.metastatic..location"]]
cat("=== cross-tab Clusters x If.metastatic..location ===\n")
print(table(md[["Clusters"]], loc, useNA="ifany"))

# any hepatocyte marker expression? check ALB, APOA1, TTR, TF, HP mean by cluster
mk <- c("ALB","APOA1","APOA2","TTR","TF","HP","SERPINA1","CYP2E1","APOB")
mk <- mk[mk %in% rownames(obj)]
cat("\n=== hepatocyte-marker mean expression (counts) by Clusters ===\n")
if(length(mk)){
  m <- GetAssayData(obj, assay="RNA", layer="counts")[mk,,drop=FALSE]
  agg <- t(apply(m, 1, function(r) tapply(r, md[["Clusters"]], mean)))
  print(round(agg,3))
} else cat("none of the hepatocyte markers found in rownames\n")
cat("\n=== DONE",format(Sys.time()),"===\n")
