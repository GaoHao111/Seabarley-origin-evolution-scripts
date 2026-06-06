setwd("E:/ProfileAdress1/")
# ==============================================================================
# 严格 10-fold Bootstrap 交叉验证：防内存溢出版 (CBI + 双向 MAE)
# ==============================================================================

if(!require(terra)) install.packages("terra")
if(!require(ecospat)) install.packages("ecospat")

library(terra)
library(ecospat)

sp_name <- "Seabarley"

boyce_results    <- numeric(10)
mae_mixed_results <- numeric(10) 
mae_pres_results  <- numeric(10) 

# 设置 terra 处理大文件时的临时目录规则，防止撑爆 C 盘
terraOptions(memfrac = 0.6, tempdir = tempdir()) 

for (i in 0:9) {
  cat(sprintf("\n正在处理第 %d 次 Bootstrap 运行...\n", i + 1))
  
  asc_file <- paste0(sp_name, "_", i, ".asc")
  csv_file <- paste0(sp_name, "_", i, "_samplePredictions.csv")
  
  # 读入栅格
  r_map <- rast(asc_file)
  pts_data <- read.csv(csv_file, check.names = FALSE)
  
  # 严格筛选测试集
  test_data <- pts_data[pts_data[["Test or train"]] == "test", ]
  test_coords <- test_data[, c("X", "Y")]
  
  # 1. 计算 Continuous Boyce Index (CBI)
  boyce_eval <- ecospat.boyce(fit = r_map, obs = test_coords, nclass = 0, PEplot = FALSE)
  boyce_results[i + 1] <- boyce_eval$cor
  
  # 提取测试存在点的模型预测值
  pred_pres <- terra::extract(r_map, test_coords)[, 2]
  pred_pres <- na.omit(pred_pres)
  
  # 2. 计算【仅存在点 MAE】 
  mae_pres_results[i + 1] <- mean(1 - pred_pres)
  
  # 3. 计算【存在-背景 混合 MAE】 
  set.seed(42 + i)
  pred_bg <- spatSample(r_map, size = 10000, method = "random", na.rm = TRUE, xy = FALSE)[, 1]
  
  pred_all   <- c(pred_pres, pred_bg)
  actual_all <- c(rep(1, length(pred_pres)), rep(0, length(pred_bg)))
  mae_mixed_results[i + 1] <- mean(abs(actual_all - pred_all))
  
  # -------------------------------------------------------------
  # 【核心修复区域：强制内存清理】
  # -------------------------------------------------------------
  # 1. 从 R 环境中删除当次循环占用内存的巨大变量
  rm(r_map, pts_data, test_data, test_coords, pred_pres, pred_bg, pred_all, actual_all, boyce_eval)
  
  # 2. 强制 R 语言立即执行垃圾回收 (释放 RAM)
  gc(verbose = FALSE)
  
  # 3. 清理 terra 包在底层 C++ 产生的临时缓存文件
  terra::tmpFiles(remove = TRUE)
}
