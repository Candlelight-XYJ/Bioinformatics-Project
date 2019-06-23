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

