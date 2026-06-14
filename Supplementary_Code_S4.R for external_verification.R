# Logistic回归模型外部验证：性能评估、ROC曲线（含AUC 95%CI）、
# 校准曲线（含截距/斜率95%CI、Hosmer-Lemeshow检验p值）、决策曲线
# ============================================================

# 安装必要包（若未安装）
if (!require("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("boot", quietly = TRUE)) install.packages("boot")
if (!require("dcurves", quietly = TRUE)) install.packages("dcurves")
if (!require("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!require("pROC", quietly = TRUE)) install.packages("pROC")
if (!require("ResourceSelection", quietly = TRUE)) install.packages("ResourceSelection")

# 加载包
library(dplyr)
library(ggplot2)
library(boot)
library(dcurves)
library(patchwork)
library(pROC)
library(ResourceSelection)

# 设置工作路径（请修改为您的实际路径）
setwd("E:/rtest")

# 1. 读取数据（请根据实际文件名调整）
data <- read.csv("IMP.data.all.long10.csv", header = TRUE)

# 2. 定义预测变量和目标变量
features <- c("Age", "Hemoglobin", "Potassium", "Chloride", "Magnesium",
              "Albumin", "International_normalized_ratio", "Blood_urea_nitrogen",
              "Renal_replacement_therapy", "liver_disease", "Malignant_cancer",
              "Cardiovascular_disease", "Congestive_heart_failure")
# 1. 读取数据（请根据实际文件名调整）
#data <- read.csv("MIdata3652.csv", header = TRUE)

# 2. 定义预测变量和目标变量#features <- c("Age", "Hemoglobin", "Potassium", "Chloride", "Magnesium",
#            "Albumin", "International_normalized_ratio", "Blood_urea_nitrogen",
#              "Renal_replacement_therapy", "Severe_liver_disease", "Malignant_cancer",
#              "Cardiovascular_disease", "Congestive_heart_failure", "ICU_stay")
target <- "status"   # 住院死亡标志 (1=死亡, 0=存活)

# 检查变量存在性
if (!all(features %in% colnames(data))) stop("部分预测变量不在数据中")
if (!target %in% colnames(data)) stop("目标变量不在数据中")

# 提取完整数据并删除缺失值（也可用多重插补，此处为演示简便）
data_complete <- data[, c(features, target)] %>% na.omit()
cat("完整样本量:", nrow(data_complete), "\n")

# 将分类变量转为因子
data_complete$gender <- as.factor(data_complete$gender)

# 3. 划分训练集和测试集（70%训练，30%测试）
set.seed(123)
train_idx <- sample(1:nrow(data_complete), size = 0.7 * nrow(data_complete), replace = FALSE)
train_data <- data_complete[train_idx, ]
test_data <- data_complete[-train_idx, ]

# 4. 构建Logistic回归模型
formula <- as.formula(paste(target, "~", paste(features, collapse = " + ")))
model <- glm(formula, data = train_data, family = binomial)

# 5. 在测试集上预测死亡概率
pred_prob <- predict(model, newdata = test_data, type = "response")
obs <- test_data[[target]]

# ========== 6. 计算ROC、AUC及其95%置信区间（DeLong法） ==========
roc_obj <- roc(obs, pred_prob, ci = TRUE, boot.n = 2000, ci.method = "delong")
auc_value <- round(auc(roc_obj), 3)
auc_ci <- round(ci(roc_obj), 3)
cat("\n========== ROC 指标 ==========\n")
cat("AUC:", auc_value, " (95% CI:", auc_ci[1], "-", auc_ci[3], ")\n")

# 绘制ROC曲线（使用ggplot2）
roc_df <- data.frame(sensitivity = roc_obj$sensitivities,
                     specificity = roc_obj$specificities)
roc_df$fpr <- 1 - roc_df$specificity

p_roc <- ggplot(roc_df, aes(x = fpr, y = sensitivity)) +
  geom_line(color = "#2E86C1", size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40", size = 0.8) +
  labs(title = "ROC Curve (External Validation)",
       x = "1 - Specificity (False Positive Rate)",
       y = "Sensitivity (True Positive Rate)") +
  theme_minimal() +
  annotate("text", x = 0.75, y = 0.10, hjust = 0, size = 4,
           label = paste0("AUC = ", auc_value, "\n95% CI: ", auc_ci[1], "–", auc_ci[3])) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme(plot.title = element_text(hjust = 0.5))

print(p_roc)

# ========== 7. 计算校准指标及置信区间（Bootstrapping） ==========
logit_pred <- log(pred_prob / (1 - pred_prob))
calib_data <- data.frame(obs = obs, logit_pred = logit_pred)

set.seed(123)
boot_results <- boot(calib_data, statistic = function(data, indices) {
  d <- data[indices, ]
  fit <- glm(obs ~ logit_pred, family = binomial, data = d)
  c(coef(fit)[1], coef(fit)[2])
}, R = 500)

intercept_ci <- boot.ci(boot_results, type = "perc", index = 1)$percent[4:5]
slope_ci <- boot.ci(boot_results, type = "perc", index = 2)$percent[4:5]

brier_point <- mean((pred_prob - obs)^2)

cat("\n========== 校准指标 ==========\n")
cat("Brier Score:", round(brier_point, 4), "\n")
cat("Calibration Intercept:", round(boot_results$t0[1], 4), 
    " (95% CI:", round(intercept_ci[1], 4), "-", round(intercept_ci[2], 4), ")\n")
cat("Calibration Slope:", round(boot_results$t0[2], 4), 
    " (95% CI:", round(slope_ci[1], 4), "-", round(slope_ci[2], 4), ")\n")

# Brier评分的95%置信区间（Bootstrap）
set.seed(123)
brier_boot <- boot(test_data, statistic = function(data, indices, model, target) {
  boot_sample <- data[indices, ]
  pred <- predict(model, newdata = boot_sample, type = "response")
  obs_boot <- boot_sample[[target]]
  return(mean((pred - obs_boot)^2))
}, R = 500, model = model, target = target)

brier_ci <- boot.ci(brier_boot, type = "perc")
brier_lower <- brier_ci$percent[4]
brier_upper <- brier_ci$percent[5]
cat("Brier Score 95% CI:", round(brier_lower, 4), "-", round(brier_upper, 4), "\n")

# ========== 新增：Hosmer-Lemeshow 校准检验（p值） ==========
hl_test <- hoslem.test(obs, pred_prob, g = 10)
hl_pvalue <- hl_test$p.value
cat("\nHosmer-Lemeshow检验 p值:", round(hl_pvalue, 5), "\n")

# ========== 8. 绘制校准曲线（X/Y轴上限为0.5，含p值） ==========
n_groups <- 10
cal_data <- data.frame(pred = pred_prob, obs = obs)
cal_data$group <- cut(cal_data$pred,
                      breaks = quantile(cal_data$pred, probs = seq(0, 1, length.out = n_groups + 1), na.rm = TRUE),
                      include.lowest = TRUE)

cal_summary <- cal_data %>%
  group_by(group) %>%
  summarise(
    mean_pred = mean(pred),
    mean_obs = mean(obs),
    n = n(),
    lower = (mean_obs + qnorm(0.975)^2/(2*n) - qnorm(0.975)*sqrt((mean_obs*(1-mean_obs) + qnorm(0.975)^2/(4*n))/n)) / (1 + qnorm(0.975)^2/n),
    upper = (mean_obs + qnorm(0.975)^2/(2*n) + qnorm(0.975)*sqrt((mean_obs*(1-mean_obs) + qnorm(0.975)^2/(4*n))/n)) / (1 + qnorm(0.975)^2/n)
  ) %>%
  filter(n >= 5)

# 创建p值标签文本
hl_label <- paste0("Hosmer-Lemeshow p = ", round(hl_pvalue, 4))

p_cal <- ggplot(cal_summary, aes(x = mean_pred, y = mean_obs)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.02, alpha = 0.7, color = "gray30") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 1) +
  geom_smooth(method = "loess", se = FALSE, color = "blue", size = 1) +
  labs(title = "Calibration Curve with 95% CI",
       x = "Predicted Probability", y = "Observed Probability") +
  theme_minimal() +
  annotate("text", x = 0.05, y = 0.45, hjust = 0, size = 3.5,
           label = paste0("Intercept: ", round(boot_results$t0[1], 3),
                          " (95% CI: ", round(intercept_ci[1], 3), "-", round(intercept_ci[2], 3), ")\n",
                          "Slope: ", round(boot_results$t0[2], 3),
                          " (95% CI: ", round(slope_ci[1], 3), "-", round(slope_ci[2], 3), ")\n",
                          "Brier: ", round(brier_point, 3),
                          " (95% CI: ", round(brier_lower, 3), "-", round(brier_upper, 3), ")\n",
                          hl_label)) +
  coord_cartesian(xlim = c(0, 0.5), ylim = c(0, 0.5))   # 修改刻度上限为0.5
print(p_cal)

# ========== 9. 决策曲线分析 (DCA) ==========
dca_data <- data.frame(status = obs, pred = pred_prob)
dca_result <- dca(status ~ pred, data = dca_data, thresholds = seq(0, 0.5, by = 0.01))
p_dca <- plot(dca_result) +
  labs(title = "Decision Curve Analysis (External Validation)",
       x = "Threshold Probability", y = "Net Benefit") +
  theme_minimal()
print(p_dca)

# ========== 10. 组合三个图形（可选） ==========
combined_plot <- (p_roc | p_cal) / p_dca + 
  plot_annotation(title = "External Validation Performance",
                  theme = theme(plot.title = element_text(hjust = 0.5)))
print(combined_plot)

# 保存图形
ggsave("ROC_Curve.png", p_roc, width = 6, height = 5)
ggsave("Calibration_Curve_with_CI.png", p_cal, width = 6, height = 5)
ggsave("DCA_Curve.png", p_dca, width = 6, height = 5)
ggsave("Combined_Performance.png", combined_plot, width = 12, height = 10)

cat("\n所有图形已保存至工作目录。\n")

