[toc]

# Meta analysis 荟萃分析
> 参考阅读文献：https://link.springer.com/book/10.1007/978-3-319-21416-0

## 1. 简介
  单个群体的GWAS研究由于样本量限制，在解释遗传效应的时候可能不是那么的完善，因此我们就需要整合其它GWAS研究的结果，来补充完善我们自己的研究结论。
  整合其它GWAS研究结果有两种方式：
+ （1）直接从原始研究数据开始整合。我们把其它GWAS研究的原始SNP，表型等数据拿到，和我们自己的数据结合，做全基因组关联分析，来看结果情况，但是这样做的缺点是太麻烦了 
+ （2）直接从统计结果开始整合。使用一定的统计策略，将相关的多个GWAS研究的统计结果整合，来说明某个位点和表型的疾病相关性，这种直接从统计结果开始整合的研究方法就是Meta分析（不用再做一次GWAS了）
＋Meta分析流程

## 2. Meta分析常用方法

#### 2.1 固定效应模型的Meta分析
> https://link.springer.com/chapter/10.1007%2F978-3-319-21416-0_2

+ 固定效应模型的Meta分析是比较流行的
+ 常用的计算模型有：（1）Inverse variance based approach【 **`最常使用`** 】（2）Cochran–Mantel–Haenszel approach
+ **`Inverse variance based approach`** 计算模型的详细公式



#### 2.2 Meta分析示例练习：使用Inverse variance based approach完成Meta分析
**`Inverse variance based approach`** 这个方法的计算公式已经列出来了：

+ 第一步， 收集不同GWAS研究中的 **`Beta值`** (snp在样本中的效应) 和 **`se(stand error)`**
+ 第二步，根据以上表格中第二行的公式计算得到 - **`整合多个研究的新的β值`** 和 **`整合多个研究的新的SE值`**
+ 第三步，由第二步得到的β值和SE值，计算 **`Z-score`**
+ 第四步，由第三步得到的 Z-score, 计算出最终的整合结果的 **`P值`**

示例数据练习 **Association_result1.csv**， **Association_result1.csv**， **Association_result1.csv** 为三份GWAS研究的统计结果（含Beta, SE, P，rsid） ，现在要将这三份研究的统计结果做Meta分析，得到整合后的统计值

```r
## function for meta analysis
meta_test <- function(b,se){
  v <- se^2 
  w <- 1/v
  w_to <- sum(w)
  b_meta <- sum(w*b)/w_to
  se_meta <- sqrt(1/w_to)
  z_meta <- b_meta/se_meta
  p_meta <- pnorm(abs(z_meta),lower.tail=F)*2
  return(c(b_meta,se_meta,z_meta,p_meta))
}
## conducting meta analysis for three studies
# reading data
study1 <- read.csv("Association_result1.csv")
study2 <- read.csv("Association_result2.csv")
study3 <- read.csv("Association_result3.csv")

dim(study1)
dim(study2)
dim(study3)
# checking if one-to-one correspondence
sum(study1$SNP == study2$SNP)
sum(study2$SNP == study3$SNP)
# checking if the effect allele is consistent among all studies
sum(study1$A1 == study2$A1)
sum(study2$A1 == study3$A1)
# checking if the effect allele is consistent among all studies
sum(study1$A1 == study2$A1)
sum(study2$A1 == study3$A1)


# rearranging data for easy calculation
b <- cbind(study1$BETA, study2$BETA, study3$BETA)
se <- cbind(study1$SE, study2$SE, study3$SE)

# meta-analysis
meta <- matrix(NA, 1000, ncol=4)
meta <- as.data.frame(meta)
names(meta) <- c("beta","se","stat","p")
for(i in 1:1000){
  tempb <- b[i,]
  tempse <- se[i,]
  meta[i,] <- meta_test(tempb,tempse)
}

head(meta)
```

#### 2.4 随机效应模型的Meta分析



#### 2.5 Meta分析的软件包比较


