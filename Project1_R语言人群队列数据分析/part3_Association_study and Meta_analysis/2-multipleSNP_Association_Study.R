#######################
## Qualitaive traits ##
#######################
df <- read.csv("multiple.csv",stringsAsFactors = F)
head(df)
# rs367896724 rs540431307 rs555500075   height disease
# 1           1           0           1 161.4140       1
# 2           1           0           1 160.3762       1
# 3           1           0           1 179.7105       1
# 4           1           0           1 143.7456       0
# 5           0           0           1 177.9926       1
# 6           1           0           1 174.3760       0
fit <- glm(disease ~ rs367896724 + rs540431307 + rs555500075, 
           data=df, family = binomial())
# 使用summary函数获取拟合函数结果
s <- summary(fit)
# 查看3个SNP计算出的系数
s$coefficients[c(2:4),]
# 计算OR值
OR　<- exp(s$coefficients[c(2:4),1])
OR
# rs367896724 rs540431307 rs555500075 
# 1.195924    4.944932    4.179826 

#########################
## Quantitative traits ##
#########################
df <- read.csv("Multiple.csv")
fit <- glm(height ~ rs367896724 + rs540431307 + rs555500075, data=df)
s <- summary(fit)
s$coefficient[c(2:4), ]



