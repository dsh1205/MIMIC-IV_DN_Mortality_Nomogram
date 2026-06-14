# XGBoost 分析：糖尿病肾病危重患者死亡风险预测
# 目标：量化14个预测因子的权重，绘制SHAP图
# 数据来源：MIdata3652.csv
# ================================
#检查xgboost包是否己经安装
if (!require("kernelshap",quietly = TRUE)){
  #如果未安装，则使用install.packages函数安装
  install.packages("kernelshap")}
# 加载包
library(xgboost)
library(pROC)
library(tidyr)
library(reshape2)

library(survival)
library(kernelshap)
library(shapviz)
library(ggplot2)
library(dplyr)

# 设置路径并读取数据
setwd("E:/rtest")
data <- read.csv("MIdata3652.csv", header = TRUE)

# 定义特征和生存变量
features <- c("Severe_liver_disease", "Malignant_cancer", "Cardiovascular_disease", 
"Congestive_heart_failure", "ICU_stay", "Renal_replacement_therapy", 
"International_normalized_ratio", "Magnesium","Chloride", "Potassium",
"Blood_urea_nitrogen", "Albumin", "Hemoglobin", "Age")
target_time <- "time"
target_status <- "status"

# 数据清洗
data_model <- data[, c(features, target_time, target_status)]
data_model <- na.omit(data_model)

# 构建 COX 模型
cox_formula <- as.formula(paste("Surv(", target_time, ",", target_status, ") ~", 
                                paste(features, collapse = " + ")))
cox_model <- coxph(cox_formula, data = data_model, x = TRUE)

# 定义预测函数（返回线性预测值）
pred_fun <- function(object, newdata) {
  predict(object, newdata = newdata, type = "lp")
}

# 使用 kernelshap 计算 SHAP 值（背景数据集取前 100 行以加快速度，可调整）
set.seed(123)
bg <- data_model[sample(nrow(data_model), 100), features]  # 背景样本
ks <- kernelshap(cox_model, X = data_model[, features], bg_X = bg, 
                 pred_fun = pred_fun, feature_names = features)

# 转换为 shapviz 对象
shp <- shapviz(ks)

# 特征重要性条形图（平均绝对 SHAP 值）
p_bar <- sv_importance(shp, kind = "bar") +
  labs(title = "Feature Importance (mean|SHAP value|)") +
  theme_minimal()
print(p_bar)

# SHAP summary 图（蜂群图）
p_bee <- sv_importance(shp, kind = "bee") +
  labs(title = "SHAP Summary (Feature Contribution)", x = "SHAP value") +
  theme_minimal()
print(p_bee)

# 若要手动绘制与您图片风格完全一致的散点图，可提取 SHAP 矩阵
shap_mat <- as.matrix(shp$S)
colnames(shap_mat) <- features
shap_long <- reshape2::melt(shap_mat, varnames = c("row", "Feature"), value.name = "SHAP")
X_long <- reshape2::melt(data_model[, features], varnames = c("row", "Feature"), value.name = "Value")
shap_long$Value <- X_long$Value

# 计算平均绝对 SHAP 排序
imp <- sort(colMeans(abs(shap_mat)), decreasing = TRUE)
shap_long$Feature <- factor(shap_long$Feature, levels = names(imp))

p_manual <- ggplot(shap_long, aes(x = SHAP, y = Feature, color = Value)) +
  geom_jitter(width = 0.01, height = 0.2, size = 1, alpha = 0.6) +
  scale_color_gradient2(low = "blue", mid = "green", high = "red") +
  labs(title = "SHAP Summary (Feature Contribution)", x = "SHAP value", color = "Feature value") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p_manual)

# 保存图像
ggsave("Feature_Importance.tif", p_bar, width = 8, height = 6, dpi = 300)
ggsave("SHAP_Summary.tif", p_manual, width = 10, height = 6, dpi = 300)



# 13. 分段回归分析（增强版：支持调试 + 备选拐点检测）
# ==========================================================================

# 安装并加载 strucchange 包
if (!require("strucchange", quietly = TRUE)) {
  install.packages("strucchange", repos = "https://cloud.r-project.org/")
}
library(strucchange)

# ===== 1. 调试：检查 shap_long 结构 =====
cat("\n=== shap_long 数据框结构 ===\n")
str(shap_long)
cat("\n列名:", paste(colnames(shap_long), collapse = ", "), "\n")
cat("\n前 3 行数据:\n")
print(head(shap_long, 3))

# 确保列名存在
if (!all(c("Feature", "SHAP", "Value") %in% colnames(shap_long))) {
  stop("shap_long 中缺少必要的列 (Feature / SHAP / Value)，请检查之前的 reshape2::melt 步骤")
}

# ===== 2. 定义分段回归函数（带备选方法）=====
find_breakpoint_robust <- function(data, feature_name) {
  # 提取该特征的数据
  df <- data %>%
    filter(Feature == feature_name) %>%
    select(value = Value, shap = SHAP) %>%
    na.omit()
  
  if (nrow(df) < 20) {
    cat(sprintf("特征 %s 样本量不足 (%d)，跳过\n", feature_name, nrow(df)))
    return(NULL)
  }
  
  # 检查值的变异程度
  if (length(unique(df$value)) < 5) {
    cat(sprintf("特征 %s 的值变异太小（仅有 %d 个不同值），跳过\n", 
                feature_name, length(unique(df$value))))
    return(NULL)
  }
  
  # 方法1：使用 breakpoints 函数
  lm_fit <- lm(shap ~ value, data = df)
  bp_obj <- tryCatch({
    breakpoints(lm_fit, breaks = 1, data = df)
  }, error = function(e) {
    cat(sprintf("breakpoints 函数出错 (%s): %s\n", feature_name, e$message))
    return(NULL)
  })
  
  # 如果 breakpoints 成功找到有效断点
  if (!is.null(bp_obj) && length(bp_obj$breakpoints) > 0) {
    bp_index <- bp_obj$breakpoints[1]
    if (bp_index >= 1 && bp_index <= nrow(df)) {
      df_sorted <- df %>% arrange(value)
      bp_value <- df_sorted$value[bp_index]
      # 计算两侧斜率
      left <- df_sorted[1:bp_index, ]
      right <- df_sorted[(bp_index+1):nrow(df_sorted), ]
      slope_left <- if (nrow(left) > 1) coef(lm(shap ~ value, data = left))["value"] else NA
      slope_right <- if (nrow(right) > 1) coef(lm(shap ~ value, data = right))["value"] else NA
      
      cat("\n========================================\n")
      cat("特征:", feature_name, "\n")
      cat("方法: breakpoints\n")
      cat("断点估计值:", round(bp_value, 3), "\n")
      cat("断点左侧斜率:", round(slope_left, 4), "\n")
      cat("断点右侧斜率:", round(slope_right, 4), "\n")
      cat("========================================\n")
      
      return(list(
        feature = feature_name,
        breakpoint = bp_value,
        slope_left = slope_left,
        slope_right = slope_right,
        data = df,
        method = "breakpoints"
      ))
    }
  }
  
  # 方法2：备选方法（基于残差平方和最小化）
  cat(sprintf("特征 %s 未通过 breakpoints 找到有效断点，尝试备选方法...\n", feature_name))
  df_sorted <- df %>% arrange(value)
  n <- nrow(df_sorted)
  best_rss <- Inf
  best_i <- NA
  
  for (i in 2:(n-1)) {
    left <- df_sorted[1:i, ]
    right <- df_sorted[(i+1):n, ]
    if (nrow(left) > 1 && nrow(right) > 1) {
      fit_left <- lm(shap ~ value, data = left)
      fit_right <- lm(shap ~ value, data = right)
      rss <- sum(resid(fit_left)^2) + sum(resid(fit_right)^2)
      if (rss < best_rss) {
        best_rss <- rss
        best_i <- i
      }
    }
  }
  
  if (!is.na(best_i) && best_i >= 2 && best_i <= n-1 && best_rss < Inf) {
    bp_value <- df_sorted$value[best_i]
    left <- df_sorted[1:best_i, ]
    right <- df_sorted[(best_i+1):n, ]
    slope_left <- coef(lm(shap ~ value, data = left))["value"]
    slope_right <- coef(lm(shap ~ value, data = right))["value"]
    
    cat("\n========================================\n")
    cat("特征:", feature_name, "\n")
    cat("方法: alternative (RSS minimization)\n")
    cat("断点估计值:", round(bp_value, 3), "\n")
    cat("断点左侧斜率:", round(slope_left, 4), "\n")
    cat("断点右侧斜率:", round(slope_right, 4), "\n")
    cat("========================================\n")
    
    return(list(
      feature = feature_name,
      breakpoint = bp_value,
      slope_left = slope_left,
      slope_right = slope_right,
      data = df,
      method = "alternative"
    ))
  } else {
    cat(sprintf("特征 %s 备选方法也未能找到合适断点\n", feature_name))
    return(NULL)
  }
}

# ===== 3. 获取可用的特征并筛选 =====
available_features <- unique(shap_long$Feature)
cat("\n可用的特征列表:\n")
print(available_features)

# 定义需要分析的连续型特征（与 features 中的名称严格一致）
features_of_interest <- c("Albumin", "Potassium", "Magnesium", "International_normalized_ratio",
                          "Chloride", "Blood_urea_nitrogen", "Hemoglobin", "Age")

# 取交集
features_to_analyze <- intersect(features_of_interest, available_features)
if (length(features_to_analyze) == 0) {
  cat("警告：没有找到需要分析的特征，请检查特征名称是否匹配。\n")
  cat("当前 shap_long 中的特征有:", paste(available_features, collapse = ", "), "\n")
} else {
  cat("将分析以下特征:\n")
  print(features_to_analyze)
}

# ===== 4. 执行分段回归 =====
breakpoints_alt <- list()
for (feat in features_to_analyze) {
  res <- find_breakpoint_robust(shap_long, feat)
  if (!is.null(res)) {
    breakpoints_alt[[feat]] <- res
  }
}

# ===== 5. 绘制断点依赖图 =====
nice_names <- c(
  "Albumin" = "Albumin (g/dL)",
  "Potassium" = "Potassium (mEq/L)",
  "Magnesium" = "Magnesium (mEq/L)",
  "International_normalized_ratio" = "INR",
  "Chloride" = "Chloride (mEq/L)",
  "Blood_urea_nitrogen" = "BUN (mg/dL)",
  "Hemoglobin" = "Hemoglobin (g/dL)",
  "Age" = "Age (years)"
)

if (length(breakpoints_alt) > 0) {
  for (feat in names(breakpoints_alt)) {
    res <- breakpoints_alt[[feat]]
    df <- res$data
    bp <- res$breakpoint
    x_label <- ifelse(feat %in% names(nice_names), nice_names[feat], feat)
    
    p_break <- ggplot(df, aes(x = value, y = shap)) +
      geom_point(alpha = 0.4, size = 1, color = "darkgray") +
      geom_smooth(method = "lm", se = FALSE, color = "red", size = 1.2) +
      geom_vline(xintercept = bp, linetype = "dashed", color = "blue", size = 0.8) +
      annotate("text", x = bp, y = max(df$shap) * 0.9,
               label = paste("Breakpoint =", round(bp, 2)),
               hjust = -0.1, color = "blue", size = 4) +
      labs(title = paste("SHAP Dependence with Breakpoint –", feat),
           x = x_label, y = "SHAP value") +
      theme_minimal()
    print(p_break)
    ggsave(paste0("SHAP_breakpoint_alt_", feat, ".png"), p_break, width = 6, height = 4, dpi = 300)
  }
  
  # ===== 6. 输出结果汇总表 =====
  results_alt_df <- data.frame(
    Feature = character(),
    Breakpoint = numeric(),
    Slope_left = numeric(),
    Slope_right = numeric(),
    Method = character(),
    stringsAsFactors = FALSE
  )
  for (feat in names(breakpoints_alt)) {
    res <- breakpoints_alt[[feat]]
    results_alt_df <- rbind(results_alt_df, data.frame(
      Feature = feat,
      Breakpoint = round(res$breakpoint, 3),
      Slope_left = round(res$slope_left, 4),
      Slope_right = round(res$slope_right, 4),
      Method = res$method,
      stringsAsFactors = FALSE
    ))
  }
  cat("\n\n分段回归结果汇总表：\n")
  print(results_alt_df)
  write.csv(results_alt_df, "SHAP_breakpoint_alternative_results.csv", row.names = FALSE)
} else {
  cat("\n未找到任何断点。可能原因：\n")
  cat("1. 特征值分布无明显的分段线性模式\n")
  cat("2. 样本量过小或值变异不足\n")
  cat("3. shap_long 中的特征名称与 features_of_interest 不一致\n")
  cat("请检查控制台输出的调试信息。\n")
}

cat("\n所有分析完成。\n")
