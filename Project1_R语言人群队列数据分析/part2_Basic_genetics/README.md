#### 1. 群体遗传学的基础知识
+ 群体遗传学定义
+ 基础名词解释
+ Hardy-Weinberg Equilibrium 哈迪温伯格平衡
+ Linkage disequilibrium 连锁不平衡
+ Heritability 遗传力

#### 2. 练习题
+ 使用R语言计算等位基因和基因型的频率
+ 使用R语言绘制HWE图

---


## 1. 群体遗传学的基础知识
#### 群体遗传学定义
+ 遗传学的分支学科，也是进化生物学的一部分，研究种群内部和种群之间的差异
+ **目标**： 描述群体的 **`遗传结构`** 和 对群体的进化压力进行理论研究

#### 基础名词解释
+ **Locus(loci)** (位点): the place on a chromosome where a gene resides( 基因组上的一段区域DNA/一个位点/一个snp)
+ **Allele** (等位基因): two different forms of a gene located at a specific position on a specific chromosome(通常说的snp就是biallele)
+ **Multiple allele** (多等位基因): three or more alternative or allelic forms of a gene, only two of which can exist in any normal, diploid individual (例如重复片断的插入)
+ **Genotype** (基因型)： the genetic constitution of an organism, which comprises the entire complex of genes inherited from both parents
+ **Allele frequency**  (等位基因频率): the ratio of the number of a specified allele in a population to the total of all alleles at its genetic locus.
+ **Genotype frequency** （基因型频率） : the proportion of a given genotype in a population.
+ **Evolution** （进化）: A change in the frequency of an allele.

#### Hardy-Weinberg Equilibrium 哈迪温伯格平衡
+ 定义：“哈迪-温伯格定律”是指在理想状态下，各等位基因的频率和等位基因的基因型频率在遗传中是稳定不变的，即保持着基因平衡
</br>
+ 实现HWE定律需要 **`满足几个条件`** ：
  + Large population size（种群足够大）
  + Random mating（种群的个体间随机交配）
  + No mutations（没有突变）
  + No natural selection（没有自然选择）
  + No migration（没有迁移）
</br>
+ HWE中各基因频率和各基因型存在如下等式关系，且保持不变：
  + 设基因型Aa中，A的基因频率是p, a的基因频率是q 
  + **`p+q=1, 且 p^2 +2pq + q^2 = 1`**
</br>
+ 几个影响哈迪温伯格平衡定律的因素
  + Genetic drift (遗传漂变)
  + Mutation (突变)
  + Natural selection （自然选择）
  + Non-random mating (非随机交配)
  + Gene flow (基因流动)
</br>

**`需要注意的是：`** 尽管这些点会影响哈迪文伯格定律，但是大多数情况下，哈迪文伯格平衡是不受影响的一个不遵循哈迪巍峨伯格平衡的群体经过一代自由交配后，第二代会变成平衡的(这个是可证明的，练习题中有R 脚本证明)

+ 哈迪温伯格平衡定律的假设检验
  + 假设得到1000个样本的有基因型的数据，首先要对基因型数据进行质量检测，哈迪文伯格平衡是质检的重要标准之一。
  + 使用卡方检验进行检测（见练习题x.x）



#### Linkage disequilibrium 连锁不平衡
+ 定义
  + 连锁平衡：如果两个位点处于连锁平衡状态，说明它们在每一代中都是完全独立遗传的。那么AB单倍型的频率P(AB) = P(A) * P(B)
  + 连锁不平衡：指在某一个群体中，不同座位上两个基因同时遗传的频率明显高于预期的随机频率的现象。D(AB) = P(AB) - P(A) * P(B)
</br>
+ 连锁不平衡的原理
  + The mixing time of two original populations was insufficient to produce complete randomization
  + Two mutations are too close to be separated through recombination
  + Some alleles in linkage are advantage in selection
</br>
+ R语言绘制连锁不平衡图
一块一块的就是连锁，颜色越深 连锁越强

#### Heritability 遗传力
+ 定义
  + 用于估计由遗传变异引起的表型性状的变异程度
  + 度量遗传对表型的贡献程度
  + **`Phenotype (P) = Genotype (G) + Environment (E)`**
  + **`Var(P) = Var(G) + Var(E) + 2 Cov(G,E)`**
+ 表型方差的组成


+ 广义遗传力


+ 狭义遗传力

+ 一些容易误解的概念


</br>
+ 估计遗传力的方法
  + 相关性/回归的方法：（1）父母-子代回归 （2）近亲比较 （3）双胞胎研究 （详细解释）
  + 方差分析的方法：（1）ANOVA
  + 线性混合模型：（1）基于家系的研究 （2）基于SNP的研究




---
## 练习题
#### 使用R语言计算等位基因和基因型的频率

#### 使用R语言绘制HWE图

#### 