[toc]

## 1. 基因－基因交互关系
#### 1.1 上位效应
+ 定义：上位性是指一个基因的作用依赖于一个或多个“修饰基因”的存在，即遗传背景。最初这个术语的意思是一个基因的表型效应被另一个基因所掩盖。
+ 上位性的分类：

![type](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part5_Interaction%20Analysis/pic/1.jpg)

  + positive epistasis
  + negative epistasis
  + sign epistasis
  + reciprocal sign epistasis
遗传－环境因素交互作用，基因之间的交互作用/上位效应 是多种多样的
上位性的定义中，positive-negative 这个定义不是绝对正确的，还是从synergistic－antagonistic角度看问题吧，只要偏离了，当两个效应同时存在时，产生的效果大于两者本身的叠加， 就算是上位性？

#### 1.2 解释表型方差
遗传因素之间的交互作用，用于解释表型方差
![2](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part5_Interaction%20Analysis/pic/2.png)


## 2. 检测遗传因素交互作用的方法(snp-snp之间互作)
+ 统计方法： 回归模型(线性回归/逻辑回归)
+ 其他机器学习/数据挖掘算法：随即森林，神经网络，等等

#### 2.1 回归模型检测
![3](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part5_Interaction%20Analysis/pic/3.png)

![4](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part5_Interaction%20Analysis/pic/4.png)

![5](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part5_Interaction%20Analysis/pic/5.png)


需要注意，统计学上有意义不代表生物学上有意义
**`根据x1x2的系数 来判断x1- x2 之间是否有交互作用`**
+ **示例练习题**
fit1 <- lm(height ~ rs1+rs2 + r1*r2,df)
数据和PPT 结果不一样，老师数据出来r1xr2得出的参数是NA， 原因是54043那个snp 的基因型有问题，只有0，1，没有2的，所以统计层面不能支撑住这个交互关系的检测
+ 另外一个例子 连续变量,自主模拟数据，说明交互效应对表型方差的解释性 ／ 离散变量一样的，只不过用的是逻辑回归而已
```r
x1<- rnorm(1000)
x2<- rnorm(1000)
x1x2 <-x1*x2
y <- x1*x2 +rnorm(1000) # 
hist（x1x2）# 看一下直方图 📊
sd(x1x2) ＃ 看一下标准差
cor(y,x1x2) 

fit1 <- lm(y~x1+x2) ＃ 无交互效应的模型
fit2 <- lm(y~x1+x2+x1*x2) 和 lm(y~x1*x2) 语法意义一样 ＃ 有交互效应的模型
summary(fit1)
summary(fit2)
```

+ 使用上面的fit1,fit2运行出来的检测模型， 预测一下df1和df2
```r
df1 <- as.data.frame(cbind(x1,x2,y))
df2 <- as.data.frame(cbind(newx1,newx2.newy))
phat1 <- predict(fit1,newdata=df2)
phar2 <- predict(fit2,newdata=df2)
head(phat1)
hist(phat1)
plot(phat1,df2$y) # 第一个模型 预测效果较差
plot(phat2,df2$y) # 第二个模型 明显效果更好，图类似于一个斜线
```
#### 2.2 回归模型的优点、缺点
+ 数据量不是很大的情况下，用统计模型基本就可以了/ 机器学习那些。。。数据量大的时候用

![6](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part5_Interaction%20Analysis/pic/6.png)


#### 2.3 其它模型
+ MDR

![7](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part5_Interaction%20Analysis/pic/7.png)


+ Single classification tree

![8](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part5_Interaction%20Analysis/pic/8.png)


+ 贪婪搜索

![9](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part5_Interaction%20Analysis/pic/9.png)


## 3. 遗传多效性

#### 3.1 混合杂合性

![10](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part5_Interaction%20Analysis/pic/10.png)


#### 3.2 基因多效性
**`一个基因可以对多个表型有影响`**

![11](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part5_Interaction%20Analysis/pic/11.png)

#### 3.3 检测基因多效性的方法

![12](https://github.com/Candlelight-XYJ/Bioinformatics-Project/blob/master/Project1_R%E8%AF%AD%E8%A8%80%E4%BA%BA%E7%BE%A4%E9%98%9F%E5%88%97%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90/part5_Interaction%20Analysis/pic/12.png)

MTAG方法不错，但是还是有差距
combined test 刘老师开发的新方法，下一步会在GPU 上优化


