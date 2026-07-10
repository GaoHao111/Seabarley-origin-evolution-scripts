suppressMessages(library(admixtools))
suppressMessages(library(tidyverse))
suppressMessages(library(igraph))
suppressMessages(library(future))

# 0. 全局设置与数据加载
set.seed(20260707) # 固定随机数种子保证顶刊要求的绝对可重复性
my_f2_dir <- "D:/ProfileAdress/93SB_f2_dir/"

cat("[INFO] 加载 f2 统计量数据...\n")
f2_data <- f2_from_precomp(my_f2_dir)

# 1. 定义八大竞争模型 (Hypothesis Models Database)

# Mod1: 核心假说 (SubA 来自 M2subgr1; SubB 来自 G2) - 无杂交
mod1 <- matrix(c(
  'R', 'G_anc',        'R', 'M_anc',
  'G_anc', 'G2',       'G_anc', 'G4SubB',
  'M_anc', 'M2subgr2', 'M_anc', 'M2subgr1_anc',
  'M2subgr1_anc', 'M2subgr1',
  'M2subgr1_anc', 'G4SubA'
), ncol = 2, byrow = TRUE)

# Mod2: 反证底盘 (SubA 来自 M2subgr2) - 无杂交
mod2 <- matrix(c(
  'R', 'G_anc',        'R', 'M_anc',
  'G_anc', 'G2',       'G_anc', 'G4SubB',
  'M_anc', 'M2subgr1', 'M_anc', 'M2subgr2_anc',
  'M2subgr2_anc', 'M2subgr2',
  'M2subgr2_anc', 'G4SubA'
), ncol = 2, byrow = TRUE)

# Mod3: 基部分化假说 (SubA 早于 M1/M2 分化)
mod3 <- matrix(c(
  'R', 'G_anc',        'R', 'M_anc',
  'G_anc', 'G2',       'G_anc', 'G4SubB',
  'M_anc', 'G4SubA',   'M_anc', 'M2_core',
  'M2_core', 'M2subgr1',
  'M2_core', 'M2subgr2'
), ncol = 2, byrow = TRUE)

# Mod4: M1底盘 + M2渗入
mod4 <- matrix(c(
  'R', 'G_anc',        'R', 'M_anc',
  'G_anc', 'G2',       'G_anc', 'G4SubB',
  'M_anc', 'M2subgr2_anc', 'M_anc', 'M2subgr1_anc',
  'M2subgr1_anc', 'M2subgr1',
  'M2subgr2_anc', 'M2subgr2',
  'M2subgr1_anc', 'G4SubA_anc',
  'M2subgr2_anc', 'G4SubA_anc',
  'G4SubA_anc', 'G4SubA'
), ncol = 2, byrow = TRUE)

# Mod5: M1底盘 + G2渗入 (最佳单杂交模型)
mod5 <- matrix(c(
  'R', 'G_anc',        'R', 'M_anc',
  'G_anc', 'G2_anc',   
  'G2_anc', 'G2',      'G_anc', 'G4SubB',
  'M_anc', 'M2subgr2', 'M_anc', 'M2subgr1_anc',
  'M2subgr1_anc', 'M2subgr1',
  'M2subgr1_anc', 'G4SubA_anc',
  'G2_anc', 'G4SubA_anc',
  'G4SubA_anc', 'G4SubA'
), ncol = 2, byrow = TRUE)

# Mod6: M1底盘 + SubB渗入 (多倍化后亚基因组串扰)
mod6 <- matrix(c(
  'R', 'G_anc',          'R', 'M_anc',
  'G_anc', 'G2',         
  'G_anc', 'G4SubB_anc', 'G4SubB_anc', 'G4SubB',
  'M_anc', 'M2subgr2',   'M_anc', 'M2subgr1_anc', 
  'M2subgr1_anc', 'M2subgr1',
  'M2subgr1_anc', 'G4SubA_anc',                     
  'G4SubB_anc', 'G4SubA_anc',                       
  'G4SubA_anc', 'G4SubA'
), ncol = 2, byrow = TRUE)

# Mod7: M2底盘 + G2渗入 (最严苛的控制变量反证)
mod7 <- matrix(c(
  'R', 'G_anc',        'R', 'M_anc',
  'G_anc', 'G2_anc',   
  'G2_anc', 'G2',      'G_anc', 'G4SubB',
  'M_anc', 'M2subgr1', 'M_anc', 'M2subgr2_anc',
  'M2subgr2_anc', 'M2subgr2',
  'M2subgr2_anc', 'G4SubA_anc',
  'G2_anc', 'G4SubA_anc',
  'G4SubA_anc', 'G4SubA'
), ncol = 2, byrow = TRUE)

# Mod8: 多重共祖基部模型
mod8 <- matrix(c(
  'R', 'G_anc',        'R', 'M_anc',
  'G_anc', 'G2',       'G_anc', 'G4SubB',
  'M_anc', 'M2subgr2_anc', 'M2subgr2_anc', 'M2subgr1',   
  'M2subgr2_anc', 'M2subgr2',
  'M2subgr2_anc', 'G4SubA'                        
), ncol = 2, byrow = TRUE)

# 整合至列表，键名将直接用于图表标题
models_list <- list(
  "Mod1_Target_SubA_from_M1"         = mod1,
  "Mod2_Alt_SubA_from_M2"            = mod2,
  "Mod3_Alt_SubA_from_BasalM"        = mod3,
  "Mod4_Admix_M1_Intro_M2"           = mod4,
  "Mod5_Admix_M1_Intro_G2"           = mod5,
  "Mod6_Admix_Subgenome_Crosstalk"   = mod6,
  "Mod7_Control_M2_Intro_G2"         = mod7,
  "Mod8_Alt_M2_Basal_M1"             = mod8
)

# 2. 批量评估 8 个模型并生成汇总表
cat("[INFO] 开始批量运行并评估 8 个竞争模型...\n")

results_df <- imap_dfr(models_list, function(edges, model_name) {
  fit <- qpgraph(f2_data, edges, return_fstats = TRUE)
  
  fstats <- fit$f4
  if (is.null(fstats)) fstats <- fit$f3
  
  if (!is.null(fstats) && "z" %in% colnames(fstats)) {
    worst_row <- fstats %>% arrange(desc(abs(z))) %>% slice(1)
    pop_cols <- grep("^pop", colnames(worst_row), value = TRUE)
    worst_pair <- paste(worst_row[1, pop_cols], collapse = ",")
  } else {
    worst_pair <- "N/A"
  }
  
  tibble(
    Model = model_name,
    Score = round(fit$score, 3),
    Worst_Z = round(fit$worst_residual, 3),
    Worst_f4_Pair = worst_pair
  )
})

# 按拟合得分排名
final_table <- results_df %>% arrange(Score)

cat("\n================ 8大模型评估结果排名 =================\n")
print(as.data.frame(final_table))
write_csv(final_table, "Table_S1_qpGraph_8Models_Evaluation.csv")

# 3. 输出 8 个模型的拓扑图合集
cat("[INFO] 正在生成 8个模型的拓扑比对图 (Fig_S1_All_8_Models.pdf)...\n")
pdf("Fig_S1_All_8_Models.pdf", width = 8, height = 8)

for (model_name in names(models_list)) {
  # 提取对应指标用于标题
  stats <- results_df %>% filter(Model == model_name)
  title_str <- sprintf("%s\nScore: %.1f | Worst Z: %.2f", 
                       model_name, stats$Score, stats$Worst_Z)
  
  p <- plot_graph(models_list[[model_name]], title = title_str)
  print(p)
}
dev.off()

# 4. 基于 Mod5 的有约束二次杂交搜索
cat("\n[INFO] 根据评估，Mod5 为最优底盘。启动基于 Mod5 的 2-Admixture 探索...\n")

# 必须将矩阵转换为 igraph 对象以通过类型检查
mod5_igraph <- edges_to_igraph(mod5)

plan(multisession, workers = 16) # 启动多线程

search_2admix_res <- find_graphs(
  data = f2_data,
  outpop = NULL,                 
  initgraph = mod5_igraph,       # 强约束：以 Mod5 为强制起点
  numadmix = 2,                  # 允许第二条杂交边出现
  stop_gen = 2000,               
  stop_gen2 = 30
)

plan(sequential) # 释放计算资源
cat("[INFO] 二次约束搜索完成。\n")

# 5. 二次探索结果评估与可视化
# 提取最优 2-Admix 模型
sorted_res <- search_2admix_res %>% arrange(score)
best_2admix_edges <- sorted_res$edges[[1]]

fit_2admix <- qpgraph(f2_data, best_2admix_edges, return_fstats = TRUE)

cat("\n=============== Mod5 升级版 (2-Admix) 评估 ===============\n")
cat(sprintf("Score 从 Mod5 的 78.113 变化为: %.3f\n", fit_2admix$score))
cat(sprintf("Worst Z-score 从 8.591 变化为: %.3f\n", fit_2admix$worst_residual))
cat("==========================================================\n")

# 保存终极网络图
graph_title <- sprintf("Constrained 2-Admix Model (Init: Mod5)\nScore: %.2f | Worst Z: %.2f", 
                       fit_2admix$score, fit_2admix$worst_residual)

pdf("Fig_Main_Constrained_2Admix_Model.pdf", width = 8, height = 8)
print(plot_graph(best_2admix_edges, title = graph_title))
dev.off()

# 保存计算对象
saveRDS(search_2admix_res, "Data_Constrained_2Admix_FullSearch.rds")
cat("[SUCCESS] 所有分析流程已结束。文件已成功归档入库。\n")
