
#下载用到的R包

#加载上述依赖的R包以后，再手动加载VRPM包才不会报错
install.packages("D:/R+Rsudio/RStudio/o1d/packages/VRPM_1.2.tar.gz",repos=N
#if (!require("regplot", quietly = TRUE)) { install.packages("regplot") }
#if (!require("DynNom", quietly = TRUE)) { install.packages("DynNom") }                

if (!require("rms", quietly = TRUE)) install.packages("rms")
if (!require("regplot", quietly = TRUE)) install.packages("regplot")

if(!require(devtools)) install.packages("devtools")
devtools::install_version("regplot", version = "1.1", repos = "http://cran.us.r-project.org")
1if (!require("survival", quietly = TRUE)) install.packages("survival")
if (!require("DynNom", quietly = TRUE)) install.packages("DynNom")
if (!require("compare", quietly = TRUE)) install.packages("compare")
if (!require("stargazer", quietly = TRUE)) install.packages("stargazer")
if (!require("rsconnect", quietly = TRUE)) install.packages("rsconnect")
#Cox回归模型，普通列线图+添加个案列线图+动态列线图+彩色条带式静态列线图
setwd("E:/rtest")                 
#加载需要的R包
library(rms)  #拟合棋型
library(regplot) #添加个案的列线图
library(DynNom) #动态列线图
library(shiny)
#library(VRPM) #彩色条带式静态列线图
library(survival) #生存函数
library(compare)
library(stargazer)
# 1. 重新加载和准备数据
mydata <- read.csv("MIdata3652.csv", sep = ",", header = TRUE)

# 2. 变量处理（保持原样）
# 2. 变量处理（保持原样）
mydata$status <- factor(mydata$status, labels = c("0", "1"))
mydata$ICU_stay <- factor(mydata$ICU_stay, labels = c("No", "Yes"))
mydata$Renal_replacement_therapy <- factor(mydata$Renal_replacement_therapy, labels = c("No","Yes"))
mydata$Cardiovascular_disease <- factor(mydata$Cardiovascular_disease ,labels=c("No","Yes"))
mydata$Congestive_heart_failure<- factor(mydata$Congestive_heart_failure,labels=c("No","Yes"))
mydata$Malignant_cancer<- factor(mydata$Malignant_cancer,labels=c("No","Yes"))
mydata$Severe_liver_disease<- factor(mydata$Severe_liver_disease,labels=c("No","Yes"))
mydata$BGR_cut2 <- factor(mydata$BGR_cut2,labels = c("1", "2", "3"))


# 3. 创建正确的生存对象
surv_object <- with(mydata, Surv(time, status == 1))

# 4. 设置datadist
dd <- datadist(mydata)
options(datadist = "dd")

# 5. 拟合简化模型
#age+HGB+alb+BUN+potassium+chloride+magnesium+inr+rrt+icustay+CHF+CVD+MC+SLD,
model <- cph(surv_object ~ Age + Hemoglobin + Albumin + Blood_urea_nitrogen + Potassium + Chloride 
             + Magnesium + International_normalized_ratio + Renal_replacement_therapy + ICU_stay 
             + Congestive_heart_failure  + Cardiovascular_disease + Malignant_cancer + Severe_liver_disease,
              data = mydata, x = TRUE, y = TRUE, surv = TRUE,
                   time.inc = 30, iter.max = 50)

# 6. 检查模型
print(model)

# 7. 定义生存函数
if(!all(coef(model) == 0)) { # 确保模型拟合成功
  #设置不同节点的生存函数---
  surv<-Survival(model)#拟合生存函数
  surv_30d <- function(x) surv(30,x)#30天生存函数
  surv_90d <- function(x) surv(90,x)#90天生存函数
  surv_180d <- function(x) surv(180,x)#180天生存函数
  surv_365d <- function(x) surv(365,x)#365天生存函数


  # 8. 绘制列线图#初始版本的列线图
  Nomogram_1 <- nomogram(model,fun =list(surv_30d,surv_90d,surv_180d,surv_365d),lp=FALSE,#模型
                         funlabel =c('30 day survival rate', '90 day survival rate', 
                                     '180 day survival rate', '365 day survival rate'), #风险预测轴名称
                         maxscale=100,
                         fun.at= c(0.1,seq(0.1,0.9,by=0.1),0.90)) # 生存概率范围 0-1，#风险预测轴概率取值
  
  
  # 9. 绘制结果#调整一些参数
  plot(Nomogram_1,
       xfrac=.35,#变量与图形的占比（调整变量与坐标抽距离）
       cex.var=1.2,#变量字体加祖
       cex.axis=1.0,#数轴：字体的大小
       tc1=-0.5,#数轴：刻度的长度
       lmgp=0.3,#数轴：文字与刻度的距离
       label.every=1,#数轴：划度下的文字，1=连续显示，2=隔一个显示一个
       col.grid =gray(c(0.8,0.95)))

  #10. 添加个案的列线图
  library(regplot) #添加个案的列线图
 # Nomogram_2<-cph(Surv(time,status)~icustay + age + alb + stage + CHF,x =TRUE ,y =TRUE,surv =TRUE,data =mydata)
  Nomogram_2<- cph(Surv(time, as.numeric(status) - 1) ~ Age + Hemoglobin + Albumin + Blood_urea_nitrogen + Potassium + Chloride 
                   + Magnesium + International_normalized_ratio + Renal_replacement_therapy + ICU_stay 
                   + Congestive_heart_failure  + Cardiovascular_disease + Malignant_cancer + Severe_liver_disease,
         data = mydata, x = TRUE, y = TRUE, surv = TRUE)
  
    # 简化版 regplot
  regplot(Nomogram_2,
          title = nullfile(),
          points = TRUE,
          # 生存概率设置
          failtime = c(365, 180, 90, 30),  # 预测时间点
          prfail = TRUE,  # 显示失败概率（死亡概率）
          droplines = TRUE,#droplines参数指定是否绘制垂直于x轴的线
          center = TRUE, #center参数指定是否将图形的中心点设置为零，对
          odds = FALSE,#odds参数指定是否显示比值或者概率
          interval="confidence", #interval参数指定使用置信区间还是预
          plots = c("density", "boxes"), #p1ots参数指定要绘制的图形
          dencol = "skyblue",#dencol参数指定密度图的颜色
          boxcol = "lightgreen", # 箱线图颜色
          clickable = TRUE, #是否启用交互式模式，这个大家一定要换成T, 
          observation=mydata[46,],#指定在图形中展示的样本，选择第几
          # 其他选项
          plot = TRUE,        # 立即绘图
          verbose = TRUE)     # 显示详细信息 
  
  
  #11. # 使用 DynNom 创建动态列线图
  #install.packages("fastmap")
  #library(fastmap)
  library(DynNom)
  mydata$status_numeric <- as.numeric(mydata$status) - 1  # 因为因子的第一个水平是"0"，第二个是"1"，所以减1后得到0和1
  
  # 确保模型公式正确
  Nomogram_3 <- cph(Surv(time, as.numeric(status) - 1) ~ Age + Hemoglobin + Albumin + Blood_urea_nitrogen + Potassium + Chloride 
                    + Magnesium + International_normalized_ratio + Renal_replacement_therapy + ICU_stay 
                    + Congestive_heart_failure  + Cardiovascular_disease + Malignant_cancer + Severe_liver_disease,
                   data = mydata, x = TRUE, y = TRUE, surv = TRUE)
  
  summary(Nomogram_3)  #确保所有变量都是正确的。
  # 启动动态列线图
  DynNom(Nomogram_3, 
         data = mydata,
         clevel = 0.95, #clevel参数指定要使用的置信水平
         m.summary = "formatted",
         covariate = "slider",
         DNtitle="Nomogram", #DNtitle参数指定图形的标题
         DNxlab="probability", #DNxlab参数指定x轴的标签
         DNylab=NULL, #DNylab参数指定y轴的标签
         DNlimits= NULL) #DNlimts参数指定x轴的限制范围
  
  ####建立动态列线图
  library (DynNom)
  library(rsconnect)
  # 重新生成应用
  
  DNbuilder(Nomogram_3,covariate="numeric")#从滑块改为输入
  
 
   # 安装和配置部署工具
  #install.packages('rsconnect')
  
  # 设置shinyapps.io账户
  library(rsconnect)
  rsconnect::setAccountInfo(
    name = "dsh1205",  # 您的账户名
    token = "60DE4560F9E5B3B4D3AAF9ADFB4483BF",        # "your-token"从shinyapps.io获取
    secret = "n7zZeKJ+Rc696DoYsH91qo6OOVmnIs2//otRib6t"#"your-secret",从shinyapps.io获取
  )
  ####确定路径
  dir <- getwd()  # 使用英文括号，没有空格
  print(dir)  # 检查当前工作目录
  path <-paste0(dir,"/DynNomapp")
  print(path)#检查完整路径
  # 或者使用 file.path() 函数（更推荐）
  #path <- file.path(dir, "DynNomapp")
 # print(path)
  
  # 检查目录是否存在
  if(dir.exists(path)) {
    print("目录存在，准备部署...")
    
    # 部署应用
   deployApp(appDir = path, 
#    deployApp(appDir = "DynNomapp"           
              appName = "DKD_Nomogram_App",
              account = "dsh1205")
  } else {
    print("目录不存在，请先创建应用")
  }###布置网络

  