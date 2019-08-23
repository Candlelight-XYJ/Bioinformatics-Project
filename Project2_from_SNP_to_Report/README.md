[toc]

##  Introduction to data format

CHROM：染色体编号
POS：检测的SNP位点在基因组上的位置
ID： SNP的rsID
REF：野生基因型(如果ref发生了突变，说明这个人这个位点有变异，如果基因型
和ref一致，那么针对参考基因组而言，这个人没有这个snp突变，就没有这种病
的风险)
CHIA-3：基因型 

#### 关于质控
从果壳拿到的数据已经经过了质控，因此我们直接进行下一步的基因型数据填充即可

---

## Step1.trios phasing
trios phasing的目的是区分母系和父系遗传的等位基因，即，一个人的某个基因型是AT，通过phasing的步骤，就可以了解A和T分别遗传自母亲还是父亲
【此步骤对后续祖源分析有用】
trios phasing分为两类：有家系数据和没有家系数据

+ **`有家系数据`**
如果能够获取家庭成员的基因检测数据，那么就可以以父母的数据作为背景，进行phasing，这样准确率是更高的

+ **`无家系数据`**
但是大多数情况无法获取用户的父母数据，因此，只能以大人群队列的数据作为参考
由于目前也没有东亚人群的数据库，暂定 *HapMap* 和 *1000 Genomes*的数据作为reference population

+ **`使用IMPUTES2完成phasing`**
```bash
# imputes2的phase选项可用于完成phasing
# Example
./impute2 \
 -phase \
 -m ./Example/example.chr22.map \
 -g ./Example/example.chr22.study.gens \
 -int 20.4e6 20.5e6 \
 -Ne 20000 \
 -o ./Example/example.chr22.phasing.impute2
```

---

## Step2.缺失基因型推断（missing genotype imputation）
这里仍使用IMPUTES2完成imputation工作

```bash
# Example
./impute2 \
 -m ./Example/example.chr22.map \
 -h ./Example/example.chr22.1kG.haps \
 -l ./Example/example.chr22.1kG.legend \
 -g ./Example/example.chr22.study.gens \
 -strand_g ./Example/example.chr22.study.strand \
 -int 20.4e6 20.5e6 \
 -Ne 20000 \
 -o ./Example/example.chr22.one.phased.impute2
```

IMPUTES可以完成多个种类的imputation工作，具体使用指南可以参考：[IMPUTES2手册](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#examples)

---

## Step3.风险评分(Genetic risk score)
PRS(polygenic risk score)是一个比较大的概念，称作多基因风险评分，它还有个名字叫做 Genetic Risk Score(GRS)。GRS分为多个种类：简单GRS，加权GRS等等，目前基因检测公司所用的风险预测模型大多都是加权GRS(wGRS) ， 因此我们在计算得时候也采用这个模型。
#### wGRS模型定义
+ wGRS=ΣβiSi(βi 为第 i 个 SNPs 的权重， Si 为第 i 个 SNPs)。该算法认为每个风险等位基因对疾病的影响不同，通过给每个风险等位基因赋予一个相应的权重来显示不同SNPs 对疾病的影响程度不同 
+ βi 值来源于已有的GWAS研究中的OR值(离散表型为OR值，例如单眼皮或者双眼皮，逻辑回归)或者β值(连续表型为β，例如身高体重，线性回归)

#### 一个计算的例子
以我们运动基因<肩袖损伤可能性>这个项目为例，它给出的文献是 Genome-wide association study identifies a locus associated with rotator cuff injury
假设现在有一个人，它的这个位点的携带基因型是AA，我们要计算他的风险概率，(其实就可以视作我们已有OR值后，利用逻辑回归的公式反推概率P)
+ **`第一步，设定计算标准：`**
含有风险等位基因纯合子（有两个高风险等位基因）——记2分

杂合子——记1分

没有风险等位基因——记0分


+ **`第二步，查看GWAS文献中给出的统计系数（OR值）`**
**文献给出了风险snp和对应的OR值，风险碱基是A**


+ **`第三步，计算GRS`**
  + 因为这个人携带2个风险碱基(AA)，因此 GRS=ΣβiSi=1.25*2=2.5
  + 逻辑回归预测公式： p=1/(1+e^(-ΣβiSi))， 大写的P代表pheno, 小写的p代表概率

  + 计算风险概率 = 1/(1+2.71828^(-2.5)) = 0.9241417 ≈ 92.4%
  + **注意** GWAS研究中，可能还给了性别，年龄等协变量的OR值，我们在计算GRS时，如果可以获得个体样本的年龄，性别信息，也可以纳入模型一起计算。为简化计算，此处我们没有纳入常量B0。

+ **`总结`**
（1）疾病风险预测的模型有很多，我们此处用的是最基本的logistic回归模型。可以看到SNP如果越多，纳入的项目越多，那么我们在预测一个人的风险的时候，会更加综合和平均。
（2）有些GWAS研究会自己构建的PRS模型，进行预测验证，这个时候如果要追求预测准确，最佳的方法是，我们直接拿到GWAS研究者所构建的预测模型来做预测。但是通常情况下，以运动基因为例，GWAS研究比较少，构建了预测模型并给出源码的 就更少了。
（3） **我们有好多运动基因其实没拿到OR值(因为没有GWAS研究去研究他们和表型的关联性)，因此此时我们预测风险高低的时候，就用最简单的SNP加和(ΣSi，分越高 风险越大)**

#### **`关于PRSice`**


> 参考链接

---

## Step4.祖源分析(Ancestry)
+ 虽然现在有一些公开的祖源分析网站，但是大多是国外的网站，所以对于中国人的祖源分析没有太多的参考意义
+ **`WeGene`** 公司使用的是基于美国加州大学洛杉矶分校Admixture工具改进的算法来完成的祖源相似性比较
+ 我们目前暂时使用R包radmixture来做祖源分析，但是我们使用的是公共参考数据库(上海马普所的一个老师好像建立了一个汉族人群队列的数据库，之后可以探究使用)，因此准确性有待提高

