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
#assocResult[which.max(-log10(as.numeric(assocResult$P))),]

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
pheatmap(LD, color = gray(100:0/100), 
         cluster_rows = F, cluster_cols = F, 
         labels_row = "", labels_col = "")

