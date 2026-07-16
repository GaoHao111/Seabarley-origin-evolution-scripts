# 0.加载所有需要的包
library(readxl)     # 读取Excel环境数据
library(vegan)      # RDA分析核心包
library(ggplot2)    # 绘图
library(ggrepel)    # 智能标签
library(dplyr)      # 数据操作
library(car)        # 用于VIF共线性诊断
library(corrplot)   # 可视化环境变量相关性
library(data.table) # 用于高效读取.raw文件

# 1. 运行前数据准备
# 读取等位基因数据
snp_raw <- fread("genotype_matrix.raw") 
snp_ids <- as.character(snp_raw$IID)

# 读取环境数据（Excel文件）
env_data <- read_excel("G2andM2.2SelectedDistribution_ClimateDataAll.xlsx")
env_ids <- as.character(env_data$Output)

# 快速验证样本匹配
common_ids <- intersect(snp_ids, env_ids)
cat("✅ 共同样本数:", length(common_ids), "\n")

# 创建已匹配的SNP矩阵（跳过质控，直接转换）
snp_matrix <- as.matrix(snp_raw[, 7:ncol(snp_raw), with = FALSE])
rownames(snp_matrix) <- snp_ids

# 创建已匹配的环境变量矩阵
# 确保env_data的Output为字符型
env_data$Output <- as.character(env_data$Output)
# 提取所有环境变量（从Bio1到Bio19）
env_matrix <- as.matrix(env_data[, c("Bio1", "Bio2", "Bio3", "Bio4", "Bio5", "Bio6", 
                                     "Bio7", "Bio8", "Bio9", "Bio10", "Bio11", "Bio12", 
                                     "Bio13", "Bio14", "Bio15", "Bio16", "Bio17", "Bio18", 
                                     "Bio19")])
rownames(env_matrix) <- env_data$Output

# 按共同ID排序并取子集（确保顺序完全一致）
snp_matrix <- snp_matrix[common_ids, ]
env_matrix <- env_matrix[common_ids, ]

# 2.去除原始气候变量间的高度相关性
# 标准化原始环境变量
env_scaled <- scale(env_matrix) # env_matrix是你的原始19变量矩阵

# 执行PCA
pca_result <- prcomp(env_scaled, center = FALSE, scale. = FALSE)
summary_pca <- summary(pca_result)

# 确定保留的主成分数（通常取特征值 > 1 的成分）
# 查看碎石图和累计贡献率
screeplot(pca_result, type = "lines", main = "环境变量PCA碎石图")
abline(h = 1, col = "red", lty = 2) # Kaiser准则：特征值>1

print("各主成分累计方差解释率：")
print(summary_pca$importance[3, ]) # 第三行是累计贡献率

# 假设我们选择前4个主成分（请根据你的累计贡献率图调整，通常3-5个）
n_pc <- 4
cat(sprintf("选择前 %d 个主成分，它们解释了 %.1f%% 的环境变异。\n", 
            n_pc, 100 * summary_pca$importance[3, n_pc]))

# 提取主成分得分作为新的环境变量
env_pca_scores <- pca_result$x[, 1:n_pc]
colnames(env_pca_scores) <- paste0("EnvPC", 1:n_pc)

# （关键）保存载荷矩阵以解读主成分的生态意义
pc_loadings <- pca_result$rotation[, 1:n_pc]
write.csv(pc_loadings, file = "环境变量PCA_载荷矩阵.csv", row.names = TRUE)
cat("主成分载荷矩阵已保存至 '环境变量PCA_载荷矩阵.csv'，请打开查看每个PC由哪些原始变量驱动。\n")

# 3.使用环境主成分运行RDA
snp_scaled <- scale(snp_matrix) # 你的SNP矩阵

# 检查 snp_scaled 矩阵（已标准化的SNP矩阵）
cat("检查 snp_scaled 矩阵:\n")
cat("  是否有 NA?", any(is.na(snp_scaled)), "\n")
cat("  是否有 NaN?", any(is.nan(snp_scaled)), "\n")
cat("  是否有 Inf?", any(is.infinite(snp_scaled)), "\n")

# 如果有问题，定位到具体行列
if(any(is.na(snp_scaled))) {
  na_indices <- which(is.na(snp_scaled), arr.ind = TRUE)
  cat("  发现", nrow(na_indices), "个NA值。前几个位置（行，列）:\n")
  print(head(na_indices))
}

# 检查 env_pca_scores 矩阵（环境主成分得分）
cat("\n检查 env_pca_scores 矩阵:\n")
cat("  是否有 NA?", any(is.na(env_pca_scores)), "\n")
cat("  是否有 NaN?", any(is.nan(env_pca_scores)), "\n")
cat("  是否有 Inf?", any(is.infinite(env_pca_scores)), "\n")

# 处理原始 snp_matrix 中的缺失值（关键步骤！）
# 使用该SNP在所有样本中的均值进行填充（对0,1,2编码而言，均值代表平均等位基因剂量）
impute_snp_missing <- function(x) {
  if(any(is.na(x))) {
    x[is.na(x)] <- mean(x, na.rm = TRUE) # 用列均值填充NA
  }
  return(x)
}

# 应用填充函数，对每一列（每个SNP）进行操作
snp_matrix_imputed <- apply(snp_matrix, 2, impute_snp_missing)
rownames(snp_matrix_imputed) <- rownames(snp_matrix) # 保持行名

# 检查填充结果
cat("填充后缺失值检查:\n")
cat("  snp_matrix_imputed 中是否有 NA?", any(is.na(snp_matrix_imputed)), "\n")
cat("  原始缺失值总数:", sum(is.na(snp_matrix)), "\n")
cat("  填充后缺失值总数:", sum(is.na(snp_matrix_imputed)), "\n")

# 对填充后的矩阵进行标准化（生成RDA需要的 snp_scaled）
snp_scaled <- scale(snp_matrix_imputed) # 现在这个矩阵里应该没有NA了

# 快速验证
cat("\n标准化后矩阵检查:\n")
cat("  snp_scaled 中是否有 NA?", any(is.na(snp_scaled)), "\n")
cat("  snp_scaled 中是否有 Inf?", any(is.infinite(snp_scaled)), "\n")

rda_result <- rda(snp_scaled ~ ., data = as.data.frame(env_pca_scores), scale = FALSE)

# 4.模型整体显著性检验（999次置换）
cat("=== RDA模型整体显著性检验 (999次置换) ===\n")
set.seed(123) # 确保结果可重复
rda_anova_terms <- anova(rda_result, permutations = 999, by = "terms")
print(rda_anova_terms)
global_r2 <- RsquareAdj(rda_result)$r.squared
global_p <- rda_anova_terms$`Pr(>F)`[1]

# 计算调整后R²和模型概况
rda_summary <- summary(rda_result)
prop_explained <- rda_summary$constr.chi / rda_summary$tot.chi
adj_r2 <- RsquareAdj(rda_result)$adj.r.squared

cat("\n*** 核心统计摘要 ***\n")
cat(sprintf("约束轴解释的遗传变异: %.2f%%\n", prop_explained * 100))
cat(sprintf("模型整体R²: %.3f\n", global_r2))
cat(sprintf("调整后R²: %.3f\n", adj_r2))
cat(sprintf("模型显著性 P值: %.4f\n", global_p))
if(global_p < 0.001) {
  cat("模型极显著 (P < 0.001)\n")
} else if(global_p < 0.01) {
  cat("模型高度显著 (P < 0.01)\n")
} else if(global_p < 0.05) {
  cat("模型显著 (P < 0.05)\n")
} else {
  cat("模型不显著 (P ≥ 0.05)\n")
}

# 5.准备绘图数据
site_scores <- scores(rda_result, display = "sites")
env_arrows <- scores(rda_result, display = "bp")

# 明确分组：前21个为M2.2，后7个为G2，直接在标签中包含样本数
plot_sites <- data.frame(
  RDA1 = site_scores[, 1],
  RDA2 = site_scores[, 2],
  Sample = rownames(site_scores),
  Group = rep(c("M2.2 (n=21)", "G2 (n=7)"), times = c(21, 7))
)

plot_env <- data.frame(
  RDA1 = env_arrows[, 1],
  RDA2 = env_arrows[, 2],
  Variable = rownames(env_arrows)
)

# 计算箭头缩放因子
multiplier <- 0.75 * max(abs(plot_sites[, 1:2])) / max(abs(plot_env[, 1:2]))
plot_env[, 1:2] <- plot_env[, 1:2] * multiplier

# 定义美观的颜色和形状（与新Group名称匹配）
group_colors <- c("M2.2 (n=21)" = "#e87070", "G2 (n=7)" = "#299d8f")
group_shapes <- c("M2.2 (n=21)" = 19, "G2 (n=7)" = 17)  # 19=实心圆, 17=实心三角

# 创建精美RDA双序图
library(ggplot2)
library(ggrepel)

rda_plot <- ggplot() +
  # 背景网格（非常淡的灰色，仅作参考）
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey88", linewidth = 0.25) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey88", linewidth = 0.25) +
  
  # 样本点 - 增大并添加描边
  geom_point(data = plot_sites, 
             aes(x = RDA1, y = RDA2, color = Group, shape = Group, fill = Group),
             size = 4.5,  # 增大点的大小
             stroke = 0.8,  # 描边粗细
             alpha = 0.85) +  # 略微透明
  
  # 环境变量箭头
  geom_segment(data = plot_env,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               color = "#555555",  # 深灰色箭头
               linewidth = 0.6,
               alpha = 0.8) +
  
  # 环境变量标签
  geom_text_repel(data = plot_env,
                  aes(x = RDA1 * 1.08, y = RDA2 * 1.08, label = Variable),
                  color = "#333333",
                  size = 3.8,
                  fontface = "bold",
                  segment.size = 0.2,
                  segment.color = "grey70",
                  box.padding = 0.35,
                  point.padding = 0.3) +
  
  # 坐标轴和标题
  labs(
    x = sprintf("RDA1 (%.1f%%)", 100 * rda_summary$cont$importance[2, 1]),
    y = sprintf("RDA2 (%.1f%%)", 100 * rda_summary$cont$importance[2, 2]),
    title = "RDA: Genetic Variation vs Environmental Gradients",
    subtitle = sprintf("EnvPC1 significant (P=%.3f) | Adj-R²=%.3f", global_p, adj_r2),
    color = "Population", shape = "Population", fill = "Population"
  ) +
  
  # 颜色、形状和填充的映射（简化版）
  scale_color_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  scale_fill_manual(values = group_colors) +
  
  # 精美主题设置
  theme_bw(base_size = 12) +
  theme(
    # 面板和背景
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey30", linewidth = 0.6),
    
    # 标题和文本
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14,
                              margin = margin(b = 8)),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40",
                                 margin = margin(b = 12)),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 10, color = "grey30"),
    
    # 图例
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.5, "cm"),
    legend.box.spacing = unit(0.3, "cm"),
    
    # 整体布局
    aspect.ratio = 1,
    plot.margin = margin(12, 12, 12, 12)
  ) +
  
  # 固定纵横比
  coord_fixed(ratio = 1, clip = "off")

# 显示图形
print(rda_plot)

# 保存高质量图像
ggsave("RDA_Biplot_Final.pdf", rda_plot, 
       width = 8.5, height = 7, device = cairo_pdf)
