## 【目录】

#### 1. PCA 校验数据的群体结构
#### 2. 独立效应的SNP信号检测
+ 2.1 练习-检测区间内的SNP独立信号
+ 2.2 Conditional analysis 条件分析
+ 2.3 SNP效应
#### 3. 表型预测
+ 3.1 SNP Selection 选择与表型相关的SNPs
+ 3.2 建立预测模型
+ 3.3 评估预测准确度
#### 4. PRS 多基因风险评分
+ 4.1 PRS流程
+ 4.2 一些Notes以及影响因素
+ 4.3 PRS的常用软件

---

## 1. PCA 校验数据的群体结构

做GWAS之前，要先做PCA(主成分分析) 过滤差的样本，矫正群体结构。所以我们用R来实际对数据做一次PCA，看如何进行群体结构的校验。

+ PCA实验的数据： 2 代表纯合，1.5 杂合 

+ PCA PPT那个第二张代码（也加了颜色） 

是把自带样本和sample样本 放在同一个坐标系下，看PCA情况，可以看到我们的自带样本和东亚人群的sample重合

---

## 2. 独立效应的SNP信号检测

![1](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part4_Prediction/pic/1.png)

通过GWAS研究，我们可以得到曼哈顿图。从图中我们可以看出有多少位点与outcome相关（过P阈值）。但是并不是所有过阈值的点都是独立与outcome相关。如何寻找真正的causal SNPs 是接下来需要解决的问题。


#### 2.1 练习-检测区间内的SNP独立信号

![2](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part4_Prediction/pic/2.png)

大样本得出的GWAS，会有很多信号，很多snp， 在这个区域到底哪些snp真正有独立效应-对表型有关呢？ 这后续就需要筛选
+ 首先使用R语言来练习绘制SNP局部图(上图)，并标明最显著的SNP
```r
# reading data
df <- read.csv("dependent_effect.csv") # 数据3000行(人),500列(snp)
# choosing susceptibility SNPs
pos <- c(25, 300)
# generating simulated phenotype
sphe <- as.matrix(df[, pos]) %*% rep(1,2) + rnorm(nrow(df), sd = 5)
# boxplot(sphe~df[,25]) #看一下数据分布
df$y <- sphe

# function for association analysis of continuous traits
getAssoc_c <- function(data, yname, xname){
  # formula
  fo <- as.formula(paste(yname, "~", xname))
  # linear model
  fit <- lm(fo, data = data)
  s <- summary(fit)
  # association result
  assocRe <- coef(s)[2, c(1, 2, 4)]
  return(assocRe)
}

# association analysis for continuous traits
snplist <- names(df)[1:500] #对500个snp 调用函数进行分析
yname <- "y"
assocResult <- matrix(NA, 500, 4)
assocResult <- as.data.frame(assocResult)
names(assocResult) <- c("SNP", "Beta", "SE", "P")
for(i in 1:500){
  xname <- snplist[i]
  assocResult[i, ] <- c(xname, getAssoc_c(df, yname, xname))
}

##### 查看结果 assocResult 谁更显著
#which.max(-log10(p)) #用p的log来判断

#################

## plotting -log(P) vs. SNPs
# highlighting the true susceptibility SNPs by point color and size
color <- rep("black", 500)
color[pos] <- "red"
size <- rep(1, 500)
size[pos] <- 2
y <- -log10(as.numeric(assocResult$P)) # p取log后就变为 >1 的数，画图更显著些
plot(y, col = color, cex=size, xlab = "SNP", ylab = "-Log10(P)", pch=19) # 横轴代表snp, 纵轴为显著性
## plotting LD Pattern
library(pheatmap)
LD <- cor(df[, 1:500])^2
LD[upper.tri(LD)] <- 0
diag(LD) <- 0
pheatmap(LD, color = gray(100:0/100), cluster_rows = F, cluster_cols = F, labels_row = "", labels_col = "")
```

#### 2.2 Conditional analysis 条件分析
> 参考文献：Conditional and joint multiple-SNP analysis of GWAS summary statistics identifies additional variants influencing complex traits

**理解条件分析的原理和步骤**
> + http://felixfan.github.io/condition-VS-interaction/
> + https://www.cnblogs.com/chenwenyan/p/10278893.html

![conditionalAnalysis](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part4_Prediction/pic/3-conditionAnalysis.png)

GWAS结果中会有许多显著的信号，假设某个区域内，SNP位点 rsA 和 rsB都与某表型显著相关，已经确定rsA对表型是有影响效应的了，那么，rsB位点与表型的显著关联是它真的有独立效应，还是说由于它和rsA 高度连锁不平衡 而产生的假阳性的显著关联呢？
此时，我们就需要用 **`条件分析`** 来检测rsB位点是否对表型有独立的影响效应

+ **`条件分析的原理是这样的`**
假设现在经过第一轮关联分析后，我们发现rs002 和rs001 与表型的关联性都很高，为了排除连锁不平衡所产生的假阳性，想测试 rs002 SNP 对表型产生的效应是否独立于 rs001 SNP，那么首先构建以下这个模型：
**Y = b0 + b1.rs002 + b2.rs001**
+ 我们想要测试的 **rs002作为自变量** ，**rs001作为协变量** 共同构建一个回归模型，b0,b1,b2代表各SNP的效应剂量(0,1,2， 例如假设B是效应allele, 那么 基因型BB的效应剂量就是2，基因型AB的效应剂量就是１，AA则是０)　
+ 若回归模型最后计算出的结果显示 rs002的 **`P值`** 仍然小于0.05，那么可以认为rs002确实对表型有独立效应

```r
## read data
df <- read.csv("conditional_analysis.csv")
tag <- c()
remaining <- names(df)[1:500] # select 500 SNPs from df

#################
## first round ##
#################
# get p value for each SNP of remaining set
pvalue <- c()
for(i in 1:length(remaining)){
  fo <- formula(paste("y ~ ", remaining[i], sep = ""))
  fit <- lm(fo, data = df)
  s <- summary(fit)
  pvalue <- c(pvalue, coef(s)[2, 4])
}

# check if there is SNP with pvalue less than 0.05
sum(pvalue < 0.05)
## [1] 351
# update tag SNP set (the most significant SNP is chosen and added to tag SNP set)
tag <- c(tag,remaining[which.min(pvalue)])
# updata remaining SNP set (removing the chosen SNP)
remaining <- remaining[-which.min(pvalue)]

##################
## second round ##
##################

# get p value for each SNP of remaining set conditioning on the chosen SNP of tag set
pvalue <- c()
for(i in 1:length(remaining)){
  fo <- formula(paste("y ~ ", paste(tag, collapse = "+"), "+", remaining[i], sep = ""))
  fit <- lm(fo, data = df)
  s <- summary(fit)
  pvalue <- c(pvalue, coef(s)[3, 4])
}
# check if there is SNP with pvalue less than 0.05
sum(pvalue < 0.05)
```

#### 2.3 SNP效应
genotype - 理解为1个snp (AA,AB,BB) ，每增加一个B，allele 对表型增加/降低 一个效应

---

## 3. 表型预测
得到GWAS的显著关联的SNPs结果后，我们可以用这些结果来做多基因疾病的风险评估，以及其它复杂性状的预测

![predictPipeline](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part4_Prediction/pic/4-pipelineforpredict.png)

**`注意`** 
+ 做预测的时候我们的四个队列数据都是独立不同的（用于GWAS 的样本是不能用于预测分析的 ）
+ 预测模型的训练集/测试集是最好都是要独立的（但很多时候测试集可能做不到独立，可以通过交叉验证从训练集中取，但是 **训练集一定要独立** ）
+ 表型预测流程图

![pipeline](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part4_Prediction/pic/5-pipeline.png)

#### 3.1 SNP Selection  选择与表型相关的SNPs

![snpSelect](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part4_Prediction/pic/6-snpSelection.png)

+ 选择SNP的方法

![snpSelectMethods](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part4_Prediction/pic/7-methodforsnpSelect.png)


#### 3.2 建立预测模型
`Model`
+ PRS方法-线性叠加-直截了当-
+ 多分类表型-可以用多分类逻辑回归
+ 等等有多种方法，基本都可以

建模的时候，用前人做过的GWAS报道的snp Beta/ P/ SE值来用就行 ，也可以用自己的研究中的结果来作为参数输入模型

![model](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part4_Prediction/pic/8-model.png)


+ **`代码示例-连续表型-线性回归做预测`**

```r
# read data and selected predictors

df <- read.csv("prediction.csv") # 数据最后y1 y2是表型

tag <- read.csv("tag.csv") # 有独立效应的snp的列表

tag <- as.character(tag$x)

# divide all data into training set (80%) and testing set (20%)

train_index <- sample(1:nrow(df), round(nrow(df)*0.8))

test_index <- setdiff(1:nrow(df), train_index)

train <- df[train_index, ]

test <- df[test_index, ]

## for continuous trait

# build model in training set

fo1 <- formula(paste("y1 ~ ", paste(tag, collapse = "+"), sep = "")) # 针对y1 做一个线性回归

fit1 <- lm(fo1, data = train) # 要注意代码分清train和test

# predict in testing set and get the predicted value

pred1 <- predict(fit1, newdata = test) #predict是函数

head(pred1)

cor(pred1,test$y1)^2  # 查看预测表型和实际表型直接的相关性（加了平方的这个是一个表型方差）解释表型方差 来看预测是否准确。如果是100% 说明模型很好

plot(pred1,test$y1)

```

+ **`代码示例-离散表型(binary)-逻辑回归做预测`**
```r
## for binary trait

# build model in training set

fo2 <- formula(paste("y2 ~ ", paste(tag, collapse = "+"), sep = ""))

fit2 <- glm(fo2, data = train, family = binomial()) # 由于离散表型，所以用logit回归预测

# 要注意代码分清train和test

# predict in testing set and get the probability

pred2 <- predict(fit2, newdata = test, type = "response")

head(pred2)

# 之后还可以用AUC /ROC来看看模型情况

```

#### 3.3 评估预测准确度
对于不同类型的表型，我们可以使用相应的方法去评估表型预测的准确度
常见的有 **`MAD`** , **`ROC`** ,   **`AUC`** ， **`Sensitivity`** 和 **`Specificity`**

![predict_accuracy](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part4_Prediction/pic/9-prediction_accuracy.png)

+ 以binary表型为例，各评估指标的计算方法如下

![a1](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part4_Prediction/pic/10-assessment1.png)


![a2](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part4_Prediction/pic/10-assessment1.png)

+ PPV和NPV是横向维度

![ROC](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part4_Prediction/pic/11-ROC.png)


+ **示例- 计算sensitivity 和 specificity **
```r
# reading predicted result --- predicted probability

df <- read.csv("predicted_results.csv")

df1 <- df

df1$probability[df1$probability >= 0.5] <- 1

df1$probability[df1$probability < 0.5] <- 0

# get confusion matrix of df1

confm <- table(df1)

# Sensitivity

sens1 <- confm[2, 2]/sum(confm[, 2])

# Specificity

spec1 <- confm[1, 1]/sum(confm[, 1])

```

+ **示例- 计算ROC 和 AUC 曲线**
```r
# reading predicted result --- predicted probability

df <- read.csv("predicted_results.csv")

df1 <- df

df1$probability[df1$probability >= 0.5] <- 1 # 1认为他是患病的

df1$probability[df1$probability < 0.5] <- 0 # 0认为他是不患病的

# get confusion matrix of df1

confm <- table(df1) ## table做一个混淆矩阵

confm

# Sensitivity

sens1 <- confm[2, 2]/sum(confm[, 2])

# Specificity

spec1 <- confm[1, 1]/sum(confm[, 1])

```

![p1](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part4_Prediction/pic/12-p1.png)

![p2](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part4_Prediction/pic/13-p2.png)


**`注`** ： **用proc包可以直接计算ROC,AUC**

---

## 4. PRS 多基因风险评分
PRS的结果情况和GWAS结果的参数紧密相连，它的输入可以是已有的，其它GWAS研究所得到的参数值，不需要我们自己去做个GWAS 。

#### 4.1 PRS流程

![PRSpipeline](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part4_Prediction/pic/15-PRSpipeline.png)

#### 4.2 一些Notes以及影响因素

![note1](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part4_Prediction/pic/16-note1.png)

![note2](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part4_Prediction/pic/16-note2.png)

![factors](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part4_Prediction/pic/17-factors.png)

#### 4.3 PRS的常用软件
+ **`PRSice`** 是最常用的软件包

![PRSsoft](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part4_Prediction/pic/18-PRSsoft.png)

![prs](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part4_Prediction/pic/19-prs.png)
