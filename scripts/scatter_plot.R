library(ggplot2)
library(patchwork) 

# ==========================================
# 1. 读取大表 CSV 数据
# ==========================================
spatial_data <- read.csv("/group/jshandl-g00/Spatial-MetScore/data/processed/Visium_Final_InfoTable_WithMetScore.csv", header = TRUE)

# ==========================================
# 2. 计算 Pearson 相关系数 (R) 和 R平方 (R²)
# ==========================================
# 计算 Reads 和 MET_Score 的相关性 
cor_reads <- cor(spatial_data$Reads, spatial_data$MET_Score, method = "pearson", use = "complete.obs")
label_reads <- sprintf("Pearson R = %.3f | R-squared = %.3f", cor_reads, cor_reads^2)

# 计算 ESTIMATE_TumorPurity 和 MET_Score 的相关性
cor_purity <- cor(spatial_data$ESTIMATE_TumorPurity, spatial_data$MET_Score, method = "pearson", use = "complete.obs")
label_purity <- sprintf("Pearson R = %.3f | R-squared = %.3f", cor_purity, cor_purity^2)

# ==========================================
# 3. 准备画图并直接保存为上下排列的 PNG
# ==========================================
png(filename = "/scratch/hhuan40/Spatial-MetScore-JS/nmf/scatter_plots.png", 
    width = 12, height = 10, units = "in", res = 300, type = "cairo")

# 图1: x 为 Reads
plot1 <- ggplot(spatial_data, aes(x = Reads, y = MET_Score)) +
  geom_point(alpha = 0.2, size = 0.5, color = "#1f77b4") +  
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) + 
  theme_classic() +                            
  labs(title = "MetScore vs. Reads",
       subtitle = label_reads,
       x = "Reads",
       y = "MetScore")

# 图2: x 为 ESTIMATE_TumorPurity
plot2 <- ggplot(spatial_data, aes(x = ESTIMATE_TumorPurity, y = MET_Score)) +
  geom_point(alpha = 0.2, size = 0.5, color = "#d62728") +
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) + 
  theme_classic() +
  labs(title = "MetScore vs. Tumor Purity",
       subtitle = label_purity, 
       x = "Tumor Purity",
       y = "MetScore")

# 把两张图上下拼在一起，并写入 PNG
print(plot1 / plot2)

# 关闭画板，完成文件保存
dev.off()









# 1. 直接在现有数据中，对 Reads 进行 Log10 转换压缩 (加 1 是为了防止出现 log(0))
spatial_data$Log10_Reads <- log10(spatial_data$Reads + 1)

# 2. 重新计算 Pearson 相关系数 (R) 和 R平方
cor_reads <- cor(spatial_data$Log10_Reads, spatial_data$MET_Score, method = "pearson", use = "complete.obs")
label_reads <- sprintf("Pearson R = %.3f | R-squared = %.3f", cor_reads, cor_reads^2)

cor_purity <- cor(spatial_data$ESTIMATE_TumorPurity, spatial_data$MET_Score, method = "pearson", use = "complete.obs")
label_purity <- sprintf("Pearson R = %.3f | R-squared = %.3f", cor_purity, cor_purity^2)

# 3. 准备画图并直接保存为新的 PNG
png(filename = "/scratch/hhuan40/Spatial-MetScore-JS/nmf/scaled_scatter_plots.png", 
    width = 12, height = 10, units = "in", res = 300, type = "cairo")

# 图1: x 轴变为压缩后的 Log10_Reads
plot3 <- ggplot(spatial_data, aes(x = Log10_Reads, y = MET_Score)) +
  geom_point(alpha = 0.2, size = 0.5, color = "#1f77b4") +  
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) + 
  theme_classic() +                            
  labs(title = "MetScore vs. Log10(Reads)",
       subtitle = label_reads,
       x = expression(Log[10]("Reads")), 
       y = "MetScore")

# 图2: x 轴依然是 ESTIMATE_TumorPurity
plot4 <- ggplot(spatial_data, aes(x = ESTIMATE_TumorPurity, y = MET_Score)) +
  geom_point(alpha = 0.2, size = 0.5, color = "#d62728") +
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) + 
  theme_classic() +
  labs(title = "MetScore vs. Tumor Purity",
       subtitle = label_purity, 
       x = "Tumor Purity",
       y = "MetScore")

# 把两张图上下拼在一起，并写入 PNG
print(plot3 / plot4)

# 关闭画板，完成文件保存
dev.off()