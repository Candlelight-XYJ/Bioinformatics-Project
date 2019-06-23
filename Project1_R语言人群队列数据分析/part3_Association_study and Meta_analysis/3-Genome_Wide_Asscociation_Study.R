## continuous traits
#就是把计算关联这个步骤封装成一个函数了
options(stringsAsFactors = F)
getAssoc_c <- function(data, yname, xname, covname=NULL){
  # yname 为表型名，xname为snp ，covname 为协变量
  if(is.null(covname)){
    fo <- as.formula(paste(yname,"~",xname))
  }else{
  # 如果使用协变量做矫正，那么公式就要加上协变量的名
    fo <- as.formula(paste(yname, "~", xname, "+", paste(covname,collapse = "+")))
  }
  # linear model  
  fit <- lm(fo, data = data)
  # association result
  s <- summary(fit)
  # 取出需要进行分析的量，p-value,beta,SE
  res <- coef(s)[2,c(1,2,4)]
  return(res)
}

## conducting association analysis
# reading data
df <- read.csv("Data_for_Association.csv")
head(df)[,1:5] 
head(df)[,99:103]
snplist <- names(df)[2:101]
# association analysis for continuous traits
yname <- "y1"
assocResult <- as.data.frame(matrix(NA,nrow(df),4))
colnames(assocResult) <- c("SNP", "Beta", "SE", "P")

for(i in 1:nrow(df)){
  print(i)
  xname <- snplist[i]
  assocResult[i,] <- c(xname, getAssoc_c(df, yname, xname))
}
head(assocResult)











