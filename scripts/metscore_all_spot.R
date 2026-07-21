suppressPackageStartupMessages({
  library(Matrix)
  library(singscore)
  library(tidyverse)
  library(data.table)
})

# Loading the data
data_dir <- "/group/jshandl-g00/Spatial-MetScore/data/processed"

info_df <- read.csv(file.path(data_dir, "Visium_AllSpots_InfoTable_noMS.csv"), stringsAsFactors = FALSE)
expr_df <- read.csv(file.path(data_dir, "Visium_Expression_Matrix_for_NMF.csv"), row.names = 1, check.names = FALSE)



met_high_df <- read.csv("/group/jshandl-g00/Spatial-MetScore/data/raw/met-high-genes-human.csv")
met_low_df  <- read.csv("/group/jshandl-g00/Spatial-MetScore/data/raw/met-low-genes-human.csv")
met_high_genes <- na.omit(met_high_df$SYMBOL)
met_low_genes  <- na.omit(met_low_df$SYMBOL)


# 1. 检查基因是否有交集
intersect_high <- length(intersect(rownames(expr_df), met_high_genes))
intersect_low <- length(intersect(rownames(expr_df), met_low_genes))

cat("============== 诊断结果 1：基因对齐情况 ==============\n")
cat("高表达基因在表达矩阵中找到了: ", intersect_high, " 个\n")
cat("低表达基因在表达矩阵中找到了: ", intersect_low, " 个\n")
if(intersect_high == 0 & intersect_low == 0) {
    cat("🚨 警告：基因完全没有对上！请检查大小写，或者是否一个是 Symbol，一个是 Ensembl ID！\n")
    cat("矩阵的基因名长这样(前5个): ", head(rownames(expr_df), 5), "\n")
    cat("你提供的基因名长这样(前5个): ", head(met_high_genes, 5), "\n")
}

# 2. 检查 Spot_ID 是否有交集
intersect_spots <- length(intersect(info_df$Spot_ID, colnames(expr_df)))

cat("\n============== 诊断结果 2：Spot 对齐情况 ==============\n")
cat("Info表 和 表达矩阵 成功匹配了: ", intersect_spots, " 个 Spot\n")
if(intersect_spots == 0) {
    cat("🚨 警告：Spot 名字完全没有对上！通常是 R 语言把横杠 '-' 变成了点 '.' 导致的。\n")
    cat("Info表的 Spot 长这样(前3个): ", head(info_df$Spot_ID, 3), "\n")
    cat("表达矩阵的 Spot 长这样(前3个): ", head(colnames(expr_df), 3), "\n")
}


expression_matrix <- as.matrix(expr_df)
visium_ranked <- rankGenes(expression_matrix)

met_scores <- simpleScore(visium_ranked, upSet = met_high_genes, downSet = met_low_genes)
rownames(met_scores) <- sub("-[0-9]+$", "", rownames(met_scores))


info_df$MET_Score <- met_scores$TotalScore[match(info_df$Spot_ID, rownames(met_scores))]

final_columns <- c(
  "Patient_Source", 
  "Sample_Name", 
  "Organ", 
  "Reads", 
  "MET_Score", 
  "ESTIMATE_TumorPurity", 
  "ESTIMATE_ImmuneScore", 
  "ESTIMATE_StromaScore"
)


rownames(info_df) <- info_df$Spot_ID
final_table <- info_df[, final_columns]


output_path <- file.path(data_dir, "Visium_Final_InfoTable_WithMetScore.csv")
write.csv(final_table, file = output_path, row.names = TRUE)

cat("Done! stored in:\n", output_path, "\n")

