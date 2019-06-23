setwd("E:/GitHub/Bioinformatics-Project/Project1_R语言人群队列数据分析/part4_Prediction/")

#############
##         ##
#############
input<-read.table("pca_input.txt")
dim(input)
input[1:5,1:5]
a<-input[,-1]
rownames(a)<-input$V1
a[1:5,1:5]
b<-t(a)
b[1:5,1:5]
pca<-princomp(b,scores = T,cor = T)
s<-summary(pca)
dim(pca$loadings)
pc1<-pca$loadings[,1]
pc2<-pca$loadings[,2]
pc3<-pca$loadings[,3]
pc<-data.frame(pc1,pc2,pc3)
rownames(pc)<-rownames(a)
head(pc)
library(ggplot2)
ggplot(pc,aes(x=pc1,y=pc2))+geom_point()

###################


sample<-read.table("sample_group.txt",head=T)
index<-which(sample$sample %in% rownames(pc))
pc_sample<-data.frame(sample[index,],pc)
head(pc_sample)
ggplot(pc_sample,aes(x=pc1,y=pc2,shape=Group))+geom_point()
cluster<-kmeans(pc_sample[3:5],3,nstart = 20)
pc_sample$cluster<-as.factor(cluster$cluster)
head(pc_sample)
ggplot(pc_sample,aes(x=pc1,y=pc2,shape=Group,col=cluster))+geom_point()+theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
