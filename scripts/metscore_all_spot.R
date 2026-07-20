# Loading the data
data_dir <- "/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/data/processed"

info_df <- read.csv(file.path(data_dir, "Visium_AllSpots_InfoTable_noMS.csv"), stringsAsFactors = FALSE)
expr_df <- read.csv(file.path(data_dir, "Visium_Expression_Matrix_for_NMF.csv"), row.names = 1, check.names = FALSE)



met_high_df <- read.csv("/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/data/raw/met-high-genes-human.csv")
met_low_df  <- read.csv("/group/jshandl-g00/Spatial-MetScore/Spatial-MetScore/data/raw/met-low-genes-human.csv")
met_high_genes <- na.omit(met_high_df$SYMBOL)
met_low_genes  <- na.omit(met_low_df$SYMBOL)



expression_matrix <- as.matrix(expr_df)
visium_ranked <- rankGenes(expression_matrix)

met_scores <- simpleScore(visium_ranked, upSet = met_high_genes, downSet = met_low_genes)


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