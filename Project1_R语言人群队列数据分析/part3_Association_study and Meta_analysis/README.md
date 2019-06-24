## 【目录】

#### 1. 人群大队列的基本概念
+ 1.1 队列研究
+ 1.2 队列研究中常见的参数定义
+ 1.3 队列分类

#### 2. 基因遗传关联性研究 （Genetic association studies）
+ 2.1 单个SNP的关联研究
  + 2.1.1 针对离散表型
  + 2.1.2 针对数量性状
+ 2.2 多个SNP的关联研究
  + 2.2.1 离散表型
  + 2.2.2 数量性状

#### 3. Genome-wide association study 全基因组关联分析
+ 3.1 对离散表型/连续性状的关联分析
+ 3.2 GWAS数据的质量控制
+ 3.3 研究设计
+ 3.4 GWAS结果解释
  + 3.4.1 Manhattan plot 曼哈顿图
  + 3.4.2 QQ(Quantile-quantile plot )分位图
  + 3.4.3 Regional Manhattan Plot & Haploview 局部曼哈顿图
  + 3.4.4 对遗传关联显著性的解读
+ 3.5 GWAS中存在的问题


---

## 1. 人群大队列的基本概念
#### 1.1 队列研究
+ 队列研究是纵向研究的一种形式，对具有共同属性的人进行抽样组合，在一段时间内进行交叉分析。
+ 队列研究的目的主要有四点：检验病因假说，疾病预防，研究疾病的自然史，新药监测

+ 关于连锁分析（LD）和队列研究(Population)的区别
  + 连锁分析 主要集中在 单基因/罕见疾病 做研究
  + cohort/population based ， 用相关性检测，检测哪些遗传位点和表型有相关性
  + 人群based的study 关注点在频率普遍高的遗传变异或者是表型；而LD 的主要是罕见病，先天病
  + 人群based的数据收集，价格花费是比较高的，当然好处也很多

#### 1.2 队列研究中常见的参数定义
+ **Incidence(发病率)** = 发病率指在一定期间内，一定人群中某病新发生的病例出现的频率。是反映疾病对人群健康影响和描述疾病分布状态的一项测量指标。
发病率 =（某时期内某人群中某病新病例人数/同时期内暴露人口数）


+ **Prevalence(患病率)** = 指某特定时间内总人口中某病新旧病例所占比例。可按观察时间的不同分为：
  + 期间患病率（period prevalence） = 某观察期间一定人群中现患某病的新旧病例数/同期的平均人口数(或被观察人数)

  + 时点患病率（point prevalence） = 某一时点一定人群中现患某病新旧病例数/该时点人口数(或被观察人数)

  + 注意，期间患病率实际上等于某一特定期间开始时患病率加上该期间内的发病率（incidence）。


**`注意区分患病率和发病率`**

+ **Risk**,**Risk difference**, **Risk ratio**, **Odds ratio**
**`注意`** 回顾性的**病例对照研究不能计算RR**，只能计算OR，前瞻性的队列研究可以计算OR和RR
> 关于OR，RR 值的具体解释，可参考博客：http://blog.sina.com.cn/s/blog_44befaf60102uza6.html 的解释

Risk ratio(RR, 相对危险度):  暴露组某病发病率（或死亡率）/ 非暴露组该病发病率（或死亡率）之比，是反映暴露因素和疾病关联强度的一个指标。
意义：RR值说明了 **暴露组发病的危险性是非暴露组的多少倍**

Odds ratio(比值比):  **在病例对照研究中**，比值比（OR）指  病例组中暴露与非暴露人数的比值和对照组中暴露与非暴露人数的比值的比。**在队列研究中**，指的是暴露组中患病与非患病者的比值和非暴露组中患病与非患病者的比值的比。

**对于罕见疾病risk ratio 和 odd ratio 是近似相等的**

+ **例子**

![RS_OS](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part3_Association_study%20and%20Meta_analysis/picAssoc/1-RS_OS.png)


#### 1.3 队列分类
+ 根据选择标准分类
  + 出生队列
  + 暴露队列

</br>

+ 根据队列中研究对象入选的时间分类
  + 固定队列
  + 动态队列

</br>

+ 根据实验设计分类
  + 病例-对照 队列
  + 前瞻性队列
  + 回顾性队列：回顾性队列研究的研究对象是根据其在过去某时点的特征或暴露情况而入选并分组的，然后从已有的记录中追溯从那时开始到其后某一时点或直到研究当时为止这一期间内，每一成员的死亡或发病情况

</br>

+ 队列研究中常见的偏差
  + Selection Bias
  + Loss to follow-up Bias
  + Confounding Bias
     + External factors obscure or exaggerate the link between risk 
     + factors and diseases

  + Information Bias



---

## 2. 基因遗传关联性研究 （Genetic association studies）

![associationStudy](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part3_Association_study%20and%20Meta_analysis/picAssoc/2-associationStudy.png)

+ 基因关联研究

  + 检验表型或疾病状况与遗传变异之间的相关性

  + 识别与特定表型或疾病有关的基因或基因组区域

+ 复杂表型与疾病:由遗传和环境因素共同决定

+ 单核苷酸多态性(SNP)

  + 单个核苷酸的变异

  + 包括插入、删除、转换和转换

SNPs是关联研究中广泛用于检测的标记物，关联研究中使用字母字符(a,b,c ...)或数字(0,1,2)来表示snp展现的基因型

![snp_coding](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part3_Association_study%20and%20Meta_analysis/picAssoc/3-snp_coding.png)


#### 2.1 单个SNP的关联研究
##### 2.1.1 **针对离散表型**
对于离散的表型，我们使用logistic回归来计算SNP与表型性状的关联性

![logistic](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part3_Association_study%20and%20Meta_analysis/picAssoc/4-logistic.png)

+ 示例,计算SNP `rs367896724` 与disease的关联性
```r
#######################
##       Qualitaive traits      ##
#######################

# read single snps files
setwd("E:/GitHub/Bioinformatics-Project/Project1_R语言完成人群队列数据分析/part3_Association_study and Meta_analysis/")
df <- read.csv("Single.csv")
head(df)
#rs367896724 height disease
#1 1 161.4140 1
#2 1 160.3762 1
#3 1 179.7105 1
#4 1 143.7456 0
#5 0 177.9926 1
#6 1 174.3760 0

# select glm model to fit 
fit <- glm(disease ~ rs367896724, data=df, family = binomial())
s <- summary(fit)
# 查看结果
s
# 取出rs367896724的系数
s$coefficient[2, ]
# 计算OR值,离散变量中, OR=e^beta
OR <- exp(s$coefficient[2, 1])
# > OR
# [1] 1.427409
```
 **注意** 数据中的snp 0-隐性，1-杂合，2-显性， -1 - missing
+ 数据由于是非连续变量，所以使用逻辑回归（如果是连续变量 例如身高，则用lm/或glm, 线性回归）
+ **`对于逻辑回归，离散变量beta的意思是 每增加一个效应基因，对odd ratio的效应是多少`**
+ 疾病更关注的是OR， OR=e^β，因此本例中计算出OR值
+ **`logistic回归中，OR值=1，表示该因素对疾病的发生不起作用 OR值大于1，表示该因素是一个危险因素，OR值小于1，表示该因素是一个保护因素`** , 本例中的OR值为1.43>1，说明SNP rs367896724是一个危险因素


##### 2.1.2 **针对数量性状**
数量性状(例如身高，血压)使用线性回归，来做单个SNP的关联分析

![continuous](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part3_Association_study%20and%20Meta_analysis/picAssoc/5-continuous.png)

β1代表SNP对数量性状的影响

+ 示例

```r
###################
## Quantitative traits ##
###################
df <- read.csv("Single.csv")
# 使用广义线性模型做拟合, 此时就不用binomial了
fit <- glm(height ~ rs367896724, data=df)
s <- summary(fit)
s$coefficient[2, ]
# 查看Beta值
s$coefficient[2, 1]
# [1] 0.2396126
```
+ 对于线性回归，连续变量： **`beta值显示snp对疾病的影响大小。每增加一个allele ,效应是多少`** （例如假设B是效应等位基因，表型为y， 则每增加一个B,对y表型产生多大程度的影响）
+ 线性回归结束后，重点关注值：Estimate,std,p-value

#### 2.2 多个SNP的关联研究
##### 2.2.1 离散表型
针对离散表型，使用多元逻辑回归来做关联分析

![multiple_logistic](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part3_Association_study%20and%20Meta_analysis/picAssoc/6-multiple_logistic.png)

β0，β1，β2 ... βm 代表每个SNP对风险易感因素的影响效应 （也可以视作权重）

+ 示例，三个SNP和疾病(0,1)的相关性分析 
```r
##################
## Qualitaive traits ###
##################
df <- read.csv("multiple.csv",stringsAsFactors = F)
head(df)
# rs367896724 rs540431307 rs555500075 height disease
# 1 1 0 1 161.4140 1
# 2 1 0 1 160.3762 1
# 3 1 0 1 179.7105 1
# 4 1 0 1 143.7456 0
# 5 0 0 1 177.9926 1
# 6 1 0 1 174.3760 0
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
# 1.195924 4.944932 4.179826 
```
本例中，三个SNP的OR值都大于了1，说明这三个SNP对疾病而言都是危险因素

##### 2.2.2 数量性状
针对数量性状，使用多元线性回归来做关联分析

![multiple_continuous](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part3_Association_study%20and%20Meta_analysis/picAssoc/7-multiple_continuous.png)

+ 示例，3个SNP和身高的相关性分析
```r
df <- read.csv("Multiple.csv")

fit <- glm(height ~ rs367896724 + rs540431307 + rs555500075, data=df)

s <- summary(fit)

s$coefficient[c(2:4), ]
#               Estimate Std. Error     t value  Pr(>|t|)
# rs367896724  0.23749504  0.4750555  0.49993117 0.6171675
# rs540431307 -3.08198778  5.8607848 -0.52586606 0.5990279
# rs555500075  0.03321229  0.7011484  0.04736842 0.9622234
```

---

## 3. Genome-wide association study 全基因组关联分析 
**`GWAS 简单来说其实就是 关联分析（最简单的情况是用回归来做关联）/假设检验 在全基因组水平上的应用（对每个snp做一次）`**
+ 例如，假设某物种全基因组水平上有500万个SNP，想探究全基因水平上有哪些变异位点和我们研究的表型有关联，那么GWAS就是指我们对这500万个SNP分别做一次（共500万次）关联分析，最后按显著性排序，选最显著的。
+ GWAS对不同个体的全基因组遗传变异的研究，目的是看看是否有任何变异位点与某个性状特征相关联，当然，我们是比较关注snp和人类复杂特征之间的联系哒，GWAS就是很好的工具

#### 3.1 对离散表型/连续性状的关联分析
+ 示例，目前数据 `Data_for_Association.csv` 中有2504个样本，100个SNP，现在我们要对这100个SNP计算它和连续表型y1的关联性 （其实就是对100个snp循环做100次相关性分析，最后得出结果 组成一个大的结果框）
```r
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
```

+ 如果是离散的表型，就把上述代码中回归的线性回归模型换成适合离散数据的logistic回归，其它基本一样，此处不再演示


#### 3.2 GWAS数据的质量控制
+ **样品制备**

  + 高质量DNA

  + 减少批次效应的影响

+ **基因型分型后的检查**

  + 过滤掉性能差的芯片

  + 检查基因型和表型性别匹配情况

  + 过滤掉差的snp

  + 过滤掉差的样本

+ **统计检查**

  + 删除有太多缺失基因型数据的个体

  + 根据基因型缺失率来删除SNP
  + 根据MAF值来删除SNP

  + 删除那些没有通过哈迪-温伯格检验的snp

  + 根据孟德尔错误率排除个体或SNP

#### 3.3 研究设计
**单阶段研究**
+ 首先选择足够的样本对所有选择的snp进行一次基因分型
+ 然后分析各SNP与表型的关联性，计算效应大小

**多阶段研究**
+ 第一阶段(发现阶段):对所有SNP进行关联分析，找到阳性SNP(`通常p值<5e-8`)
+ 第二阶段(重复阶段):使用更大的样本对第一阶段发现的阳性snp进行基因型分型和重现
最好使用功能实验对结果进行进一步的解释和重现

#### 3.4 GWAS结果解释
 > 参考：https://www.jianshu.com/p/987859ae503c  如何理解GWAS中Manhattan plot和QQ plot所传递的信息
##### 3.4.1 Manhattan plot 曼哈顿图

![manhattan](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part3_Association_study%20and%20Meta_analysis/picAssoc/8-manhattan.jpg)

+ 根据每个SNP的染色体位置和 -log10p值，绘制与其关联的图
+ 红线表示全基因组层面的显著性阈值(P=5e-8)
+ 蓝线表示有意义的阈值(P=1e-5)

##### 3.4.2 QQ(Quantile-quantile plot )分位图

![qqnorm](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part3_Association_study%20and%20Meta_analysis/picAssoc/9-qqnorm.jpg)

+ QQ图应该前面是个直线，后面再飘起来，这样比较好，说明研究的SNP和表型有关（如上图所示）；如果QQ图从前面就开始飘，那么数据群体结构可能有问题，要再看看数据分布情况
+ **`例如下面这张图就是有问题的，需要进行数据结构的矫正`**




+ QQ图用于估计数量性状观测值与预测值之间的差异，如果直线出现了偏离，那么就认为这个SNP位点的偏离是由SNP突变所产生的遗传作用造成的，我们所研究的表型和基因型之间是存在着显著相关的自然选择作用的

+ QQ图的纵轴，是SNP位点的p-value值（加了Log,变为 -1og10(p-value)）（这是实际得到的结果，observed），
横轴是则是均匀分布的概率值（这是Expecte的结果），同样也是换算为-log10。 


##### 3.4.3 Regional Manhattan Plot & Haploview 局部曼哈顿图

##### 3.4.4 对遗传关联显著性的解读
GWAS结果中，我们会得到一些显著性关联的SNP位点，那么怎么去解释它们的关联性呢，他们是假阳性的吗？ 可以分成以下三类：
+ 直接关联：显著性的SNP确实是引起疾病易感性的真实原因
+ 间接关联：得到的显著性的SNP与真正引起疾病的基因之间是连锁不平衡的关系
+ 假阳性结果：例如由于人群结构的问题，造成了SNP与表型显著关联的假阳性现象

**`得到阳性结果后，一定要考虑以下几点`**
+ Whether the choice of phenotype is strict
+ Whether the research design is reasonable
+ Whether the statistical analysis is appropriate
+ Whether this result is false-positive
+ Whether this result is consistent with the credibility of biology
+ Can you repeat the result in another cohort


#### 3.5 GWAS中存在的问题
+ 群体分层与矫正
+ 多重检验与矫正
+ Permutation Test  - 可以用来设置p-value的阈值
+ Statistical power 统计效能







