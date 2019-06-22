#######################
## Qualitaive traits ##
#######################

# read single snps files
setwd("E:/GitHub/Bioinformatics-Project/Project1_R语言完成人群队列数据分析/part3_Association_study and Meta_analysis/")
df <- read.csv("Single.csv")
head(df)
#rs367896724   height disease
#1           1 161.4140       1
#2           1 160.3762       1
#3           1 179.7105       1
#4           1 143.7456       0
#5           0 177.9926       1
#6           1 174.3760       0

# select glm model to fit 
fit <- glm(disease ~ rs367896724, data=df, family = binomial())
s <- summary(fit)
# 查看结果
s
# 取出rs367896724的系数
s$coefficient[2, ]
# 计算OR值,离散变量中, OR=e^beta
OR <- exp(s$coefficient[2, 1])
OR

#########################
## Quantitative traits ##
#########################
df <- read.csv("Single.csv")
# 使用广义线性模型做拟合, 此时就不用binomial了
fit <- glm(height ~ rs367896724, data=df)
s <- summary(fit)
s$coefficient[2, ]
# 查看Beta值
s$coefficient[2, 1]
# [1] 0.2396126

