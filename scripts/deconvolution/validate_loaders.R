#!/usr/bin/env Rscript
suppressWarnings(suppressMessages({library(Seurat); library(schard); library(Matrix)}))
RAW <- "/group/jshandl-g00/Spatial-MetScore/data/raw/scRNAseq"
chk <- function(tag, expr) { cat("\n#### ",tag," ####\n"); tryCatch(expr, error=function(e) cat("ERROR:",conditionMessage(e),"\n")) }
isint <- function(m){ s<-m@x[1:min(5000,length(m@x))]; all(s==round(s)) }

chk("peng RAW via schard", {
  o <- schard::h5ad2seurat(file.path(RAW,"peng_2019/StdWf1_PRJCA001063_CRC_besca2.raw.h5ad"))
  cat("dims:",paste(dim(o),collapse=" x "),"\n")
  cat("rownames[1:4]:",paste(head(rownames(o),4),collapse=", "),"\n")
  cat("colnames[1:2]:",paste(head(colnames(o),2),collapse=", "),"\n")
  cat("counts integer:",isint(GetAssayData(o,layer="counts")),"\n")
})
chk("peng ANNOTATED via schard (old anndata format - Cell_type decode?)", {
  o <- schard::h5ad2seurat(file.path(RAW,"peng_2019/StdWf1_PRJCA001063_CRC_besca2.annotated.h5ad"))
  cat("dims:",paste(dim(o),collapse=" x "),"\n")
  cat("meta cols:",paste(colnames(o@meta.data),collapse=", "),"\n")
  cat("colnames[1:2]:",paste(head(colnames(o),2),collapse=", "),"\n")
  if("Cell_type" %in% colnames(o@meta.data)){
    cat("Cell_type class:",class(o$Cell_type),"\n")
    print(table(o$Cell_type, useNA="ifany"))
  } else cat("!! Cell_type NOT in meta\n")
})
chk("hlca via schard use.raw=TRUE (expect Ensembl rownames, need remap)", {
  o <- schard::h5ad2seurat(file.path(RAW,"hlca/hlca_core.h5ad"), use.raw=TRUE)
  cat("dims:",paste(dim(o),collapse=" x "),"\n")
  cat("rownames[1:4]:",paste(head(rownames(o),4),collapse=", "),"\n")
  cat("counts integer:",isint(GetAssayData(o,layer="counts")),"\n")
  cat("ann_level_3 present:", "ann_level_3" %in% colnames(o@meta.data),"\n")
  if("ann_level_3" %in% colnames(o@meta.data)) cat("n ann_level_3:",length(unique(o$ann_level_3)),"\n")
})
# hlca symbol remap source: read feature_name (categorical) via rhdf5
chk("hlca feature_name remap (rhdf5)", {
  library(rhdf5)
  cats <- h5read(file.path(RAW,"hlca/hlca_core.h5ad"),"raw/var/feature_name/categories")
  codes<- h5read(file.path(RAW,"hlca/hlca_core.h5ad"),"raw/var/feature_name/codes")
  ens  <- h5read(file.path(RAW,"hlca/hlca_core.h5ad"),"raw/var/_index")
  sym <- cats[codes+1]
  cat("n genes:",length(sym)," ens[1:2]:",paste(head(ens,2),collapse=","),
      " sym[1:2]:",paste(head(sym,2),collapse=","),"\n")
  cat("dup symbols:",sum(duplicated(sym)),"\n")
})
cat("\nDONE\n")
