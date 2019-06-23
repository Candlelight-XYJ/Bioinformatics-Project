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
