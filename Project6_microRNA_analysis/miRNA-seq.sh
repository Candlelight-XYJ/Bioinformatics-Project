## Using trinity to assembly transcriptome

#PBS -N TrinityForRnaseq
#PBS -q big 
#PBS -l walltime=12:00:00  
#PBS -j oe
cd /home/zhaoll/liuyn/2-clean_data/
 ~/miniconda3/bin/Trinity --seqType fq --SS_lib_type RF  \
           --left fastp_animinal_S482-01-T02_good_1.fq.gz,fastp_animinal_S482-01-T01_good_1.fq.gz,fastp_animinal_S482-01-T03_good_1.fq.gz \
           --right fastp_animinal_S482-01-T02_good_2.fq.gz,fastp_animinal_S482-01-T01_good_2.fq.gz,fastp_animinal_S482-01-T03_good_2.fq.gz \
           --CPU 2 --max_memory 1G
           --output /home/zhaoll/liuyn/3-ref/trinity_out_dir

###########
## miRNA ##
###########

#!/bin/bash
#PBS -N build_index
#PBS -q blade
#PBS -l nodes=node7:ppn=4
#PBS -l walltime=12:00:00
#PBS -j oe

cd /home/yjxiang/NGS-project/1-miRNA/
~/miniconda3/bin/bowtie-build Monochamus_version1.fa ma


## mapping
srun -p husn -c 4 -J bowtie sh bowtie_aglin.sh &

cd /home/yjxiang/NGS-project/1-miRNA/LKY-Ep/
~/miniconda3/bin/bowtie -f -a -p 4 -m 20 -v 0 /home/yjxiang/NGS-project/1-miRNA/ma clean.fa ma_mi_res


####################
## 1-鉴定已知miRNA ##
####################

# blastn 与GenBank, Rfam数据库比对,去除冗余的非miRNA(如 pirna, rrna等)

## 1. 提交miRNA到Blast官网,检索 genbank数据库

~/miniconda3/bin/update_blastdb.pl nt > log
#/home/miniconda2/bin/update_blastdb.pl nt > log 

## 2. Rfam数据库

# step1, download infernal
wget eddylab.org/infernal/infernal-1.1.2.tar.gz
tar xf infernal-1.1.2.tar.gz
cd infernal-1.1.2
./configure
make
make install
# 测试软件运行正常
make check
# 二进制软件在src文件夹中，可以直接使用


# step2, download 
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin

# step3, cmpress to index the Rfam.cm file
/home/yjxiang/tools/infernal-1.1.2/src/cmpress Rfam.cm

# step4, 估计基因组大小
/home/yjxiang/tools/infernal-1.1.2/easel/miniapps/esl-seqstat /home/yjxiang/NGS-project/1-miRNA/Monochamus_version1.fa

# step5, use cmscan 注释基因组中的RNA
cd /home/yjxiang/NGS-project/1-miRNA/LKY-Ep/
/home/yjxiang/tools/infernal-1.1.2/src/cmscan -Z 169.173472 \
--cpu 4 \
--cut_ga --rfam --nohmmonly --tblout /home/yjxiang/NGS-project/1-miRNA/LKY_Ep_ma-transRfam.tblout \
--fmt 2 \
--clanin /home/yjxiang/tools/Rfam.clanin /home/yjxiang/tools/Rfam.cm /home/yjxiang/NGS-project/1-miRNA/LKY-Ep/clean.fa > LKY_Ep.cmscan

srun -p husn -c 4 -J cmscan sh cmscan.sh &

## 以下检索步骤已被淘汰
# step2, download rfam scan perl script
#wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/tools/rfam_scan-1.0.4.pl
# step3, donwload rfam database
#wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
# step4, use rfam 
#rfam_scan.pl -blastdb /opt/biosoft/rfam/Rfam.fasta /opt/biosoft/rfam/Rfam.cm genome.fasta -o rfam.gff3


## 3. blast miRbase
wget ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz
wget ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz

# 筛选出六足昆虫类的所有mirna

#cut -d ' ' -f 1 listofHexapoda.txt | uniq > 3char.Hexapoda.list

#cd /home/yjxiang/NGS-project/1-miRNA/ref/
#for fn in `cat listofHexapoda.txt`;   #listofHexapoda.txt is a filename
#do
#echo ${fn}
#grep ${fn} /home/database/miRBase/hairpin.fa > ${fn}.loop.fa 
#done

## use perl script to select sequences from all sequences
#~/miniconda3/bin/perl /home/yjxiang/NGS-project/1-miRNA/selectFa.pl loop.perlInput /home/database/miRBase/hairpin.fa > hexapoda.loop.seq



# 1.从hairpin.fasta文件构建数据库
srun -p husn -c 1 -J buildblastDB sh buildblastDB.sh &

cd /home/yjxiang/NGS-project/1-miRNA/
~/miniconda3/bin/makeblastdb -in /home/yjxiang/NGS-project/1-miRNA/ref/hairpin.fa -parse_seqids \
-title "hairpin" \
-out "hairpin" \
-dbtype nucl

# 2.检索鉴定
# -task blastn-short -word_size 7 -evalue 1
srun -p husn -c 4 -J runBlast sh runBlast.sh &
cd /home/yjxiang/NGS-project/1-miRNA/LKY-Ep/
# cd /home/yjxiang/NGS-project/1-miRNA/NJ-Ep/
# cd /home/yjxiang/NGS-project/1-miRNA/LKY-Mg/
# cd /home/yjxiang/NGS-project/1-miRNA/NJ-Mg/
~/miniconda3/bin/blastn -db /home/yjxiang/NGS-project/1-miRNA/hairpin \
-task blastn-short \
-ungapped \
-word_size 15 \
-query clean.fa \
-num_threads 4 \
-evalue 1 \
-outfmt 7 \
-out LKY_Ep_blastRes
# -out NJ_Ep_blastRes
# -out LKY_Mg_blastRes
# -out NJ_Mg_blastRes


## 法二，使用bowtie

cd /home/yjxiang/NGS-project/1-miRNA/LKY-Ep/
~/miniconda3/bin/bowtie -f -a -p 4 -m 20 -v 0 /home/yjxiang/NGS-project/1-miRNA/hairpin.fa clean.fa LKY_Ep_bowtieRes

# 从结果中筛选出六足昆虫的haripin比对序列
# sed 's/原字符串/新字符串/g' 文件
sed 's/_stem-loop/z/g' hexapodaLoopSeq_forBlastDB.fa > hexapodaLoopSeq_xres.fa

# 从结果中筛选出六足昆虫的mature比对序列
# 1. 根据六足昆虫的名单，取出mature里对应物种的序列名称
cd /home/yjxiang/NGS-project/1-miRNA/ref/
for fn in `cat listofHexapoda.txt`;   #listofHexapoda.txt is a filename
do
echo ${fn}
grep ${fn} /home/database/miRBase/mature.fa > ${fn}.mature.fa 
done
# 2. R脚本根据>名称 提取序列
setwd("E:/学习资料存放处/13-Project/IOZ-miRNA/")
options(stringsAsFactors = F)
genelist <- read.table("hexapoda.mature.txt", sep = "\t")
allfasta <- read.table("mature.fa",sep ="\n")
allres <- c()
for( i in 1:nrow(genelist)){
  print(i)
  # i=1
  # i=2
  # i=105
  tmprowname <- grep(genelist[i,1], allfasta$V1)
  kk = tmprowname
  resStr = paste0(allfasta[kk,],"\n")
  if(length(kk)!=0){
  for(j in kk:(kk+50)){   ## kk+50 表示 > 之后的序列预计在50行之内
    if(substr(allfasta[j+1,1],1,1) != ">"){
      print(j)
      resStr <- paste0(resStr,allfasta[j+1,])       
    }else{
      break
    }
  }
    }
  allres <- c(allres, resStr)
  
}  

# 3. 处理列名格式，使其保持50个字符内 (blast构建数据库，必需序列名称在50字符以内)
rm(list=ls())
hexapodaSeq <- read.table("hexapoda.mature.seq",header=F,sep="\n")
library(tidyverse)
split_res=separate(hexapodaSeq,col=V1, into=c("mir","MI","spe1","spe2","other","zz"),sep=" ")

split_res$MI[is.na(split_res$MI)] <- ""
split_res$spe1[is.na(split_res$spe1)] <- ""
split_res$spe2[is.na(split_res$spe2)] <- ""
res <- paste0(split_res$mir,split_res$spe1,"_",split_res$spe2)
#res <- paste0(split_res$mir,"_",split_res$MI,"_",split_res$spe)
replaceRes <- gsub("_","",res)
write.table(replaceRes, "hexapodaLoopSeq_forBlastDB.fa",quote = F,row.names = F)

# 4. 构建六足昆虫的blast数据库

## loop
cd /home/yjxiang/NGS-project/1-miRNA/
~/miniconda3/bin/makeblastdb -in /home/yjxiang/NGS-project/1-miRNA/ref/hexapodaLoopSeq_forBlastDB.fa -parse_seqids \
-title "hexapodaLoop" \
-out "hexapodaLoop" \
-dbtype nucl

## mature
cd /home/yjxiang/NGS-project/1-miRNA/
~/miniconda3/bin/makeblastdb -in /home/yjxiang/NGS-project/1-miRNA/ref/hexapodaMatureSeq_forBlastDB.fa -parse_seqids \
-title "hexapodaMature" \
-out "hexapodaMature" \
-dbtype nucl

# 5. 搜索鉴定已知micro RNA 
# -task blastn-short -word_size 7 -evalue 1
srun -p husn -c 4 -J runBlast sh runBlast.sh &
cd /home/yjxiang/NGS-project/1-miRNA/LKY-Ep/
# cd /home/yjxiang/NGS-project/1-miRNA/NJ-Ep/
# cd /home/yjxiang/NGS-project/1-miRNA/LKY-Mg/
# cd /home/yjxiang/NGS-project/1-miRNA/NJ-Mg/
~/miniconda3/bin/blastn -db /home/yjxiang/NGS-project/1-miRNA/hairpin \
-task blastn-short \
-ungapped \
-word_size 10 \
-query clean.fa \
-num_threads 4 \
-evalue 1 \
-outfmt 7 \
-out LKY_Ep_blastRes
# -out NJ_Ep_blastRes
# -out LKY_Mg_blastRes
# -out NJ_Mg_blastRes



## 4.miRdeep 预测新mircroRNA & 定量

# 4.1 数据预处理
srun -p husn -c 4 -J runMapper sh mirDeepmapper.sh &
cd /home/yjxiang/NGS-project/1-miRNA/LKY-Ep/
~/miniconda3/bin/mapper.pl LKY_Ep.fa -c -j -l 15 -m -p /home/yjxiang/NGS-project/1-miRNA/ref/ma \
-s LKY_Ep_collapsed.fa \
-t LKY_Ep_vs_ma.arf \
-v -o 4 \

cd /home/yjxiang/NGS-project/1-miRNA/NJ-Ep/
~/miniconda3/bin/mapper.pl NJ_Ep.fa -c -j -l 15 -m -p /home/yjxiang/NGS-project/1-miRNA/ref/ma \
-s NJ_Ep_collapsed.fa \
-t NJ_Ep_vs_ma.arf \
-v -o 4 \

cd /home/yjxiang/NGS-project/1-miRNA/LKY-Mg/
~/miniconda3/bin/mapper.pl LKY_Mg.fa -c -j -l 15 -m -p /home/yjxiang/NGS-project/1-miRNA/ref/ma \
-s LKY_Mg_collapsed.fa \
-t LKY_Mg_vs_ma.arf \
-v -o 4 \

cd /home/yjxiang/NGS-project/1-miRNA/NJ-Mg/
~/miniconda3/bin/mapper.pl NJ_Mg.fa -c -j -l 15 -m -p /home/yjxiang/NGS-project/1-miRNA/ref/ma \
-s NJ_Mg_collapsed.fa \
-t NJ_Mg_vs_ma.arf \
-v -o 4 \



# 自己安装mirDeep2
git clone https://github.com/rajewsky-lab/mirdeep2.git



# 4.2 novel mirna
srun -p husn -c 1 -J runNovel sh mirDeepNovel.sh &
cd /home/yjxiang/NGS-project/1-miRNA/LKY-Ep/
~/home/yjxiang/tools/mirdeep2/miRDeep2.pl LKY_Ep_collapsed.fa \
/home/yjxiang/NGS-project/1-miRNA/ref/Monochamus_version1.fa \
LKY_Ep_vs_ma.arf \
none \
/home/yjxiang/NGS-project/1-miRNA/ref/hexapodaMatureSeq_forBlastDB.fa \
/home/yjxiang/NGS-project/1-miRNA/ref/hexapodaLoopSeq_forBlastDB.fa 2>LKY_Ep.novel.report.log

cd /home/yjxiang/NGS-project/1-miRNA/LKY-Ep/
~/miniconda3/bin/perl /home/yjxiang/tools/mirdeep2/bin/miRDeep2.pl LKY_Ep_collapsed.fa /home/yjxiang/NGS-project/1-miRNA/ref/Monochamus_version1.fa LKY_Ep_vs_ma.arf none none none

#conda remove mirdeep2

# 4.3 定量
srun -p husn -c 1 -J runQuantify sh mirDeepQuantify.sh &

cd /home/yjxiang/NGS-project/1-miRNA/LKY-Ep/
/home/yjxiang/tools/mirdeep2/bin/quantifier.pl -p /home/yjxiang/NGS-project/1-miRNA/ref/hexapodaLoopSeq_forBlastDB.fa \
-m /home/yjxiang/NGS-project/1-miRNA/ref/hexapodaMatureSeq_forBlastDB.fa -r LKY_Ep_collapsed.fa \
-y now \

cd /home/yjxiang/NGS-project/1-miRNA/NJ-Ep/
/home/yjxiang/tools/mirdeep2/bin/quantifier.pl -p /home/yjxiang/NGS-project/1-miRNA/ref/hexapodaLoopSeq_forBlastDB.fa \
-m /home/yjxiang/NGS-project/1-miRNA/ref/hexapodaMatureSeq_forBlastDB.fa -r NJ_Ep_collapsed.fa \
-y now \

cd /home/yjxiang/NGS-project/1-miRNA/LKY-Mg/
/home/yjxiang/tools/mirdeep2/bin/quantifier.pl -p /home/yjxiang/NGS-project/1-miRNA/ref/hexapodaLoopSeq_forBlastDB.fa \
-m /home/yjxiang/NGS-project/1-miRNA/ref/hexapodaMatureSeq_forBlastDB.fa -r LKY_Mg_collapsed.fa \
-y now \

cd /home/yjxiang/NGS-project/1-miRNA/NJ-Mg/
/home/yjxiang/tools/mirdeep2/bin/quantifier.pl -p /home/yjxiang/NGS-project/1-miRNA/ref/hexapodaLoopSeq_forBlastDB.fa \
-m /home/yjxiang/NGS-project/1-miRNA/ref/hexapodaMatureSeq_forBlastDB.fa -r NJ_Mg_collapsed.fa \
-y now \


# 5. mirna 靶基因预测
#软件预测都选取结合能量值低于-20 kcal/mol 且种子序列区不含 GU 对的基因，RNAhybrid 预测时设置helix 和loop 长度都在 10 以内

# 靶基因预测
srun -p husn -c 1 -J runTarget sh target.sh > targetGenelog &
## rnahybrid
cd /home/yjxiang/NGS-project/1-miRNA/targetGenes/
/home/yjxiang/tools/RNAhybrid-2.1.2/src/RNAhybrid -b 1 -e -20 -f 10 -v 2 -u 6 -s 3utr_fly \
-t /home/yjxiang/NGS-project/1-miRNA/ref/Monochamus_version1.fa -q /home/yjxiang/NGS-project/1-miRNA/diffmiRNA/LKY_Ep_Vs_LKY_Mg.seq > LKY_Ep_Vs_LKY_Mg_targetGenes

# LKY_Ep_Vs_NJ_Ep.seq
# LKY_Ep_Vs_NJ_Mg.seq
# NJ_Ep_Vs_LKY_Mg.seq
# NJ_Mg_Vs_LKY_Mg.seq
# NJ_Mg_Vs_NJ_Ep.seq

## miranda
~/miniconda3/bin/miranda /home/yjxiang/NGS-project/1-miRNA/diffmiRNA/LKY_Ep_Vs_NJ_Ep.seq /home/yjxiang/NGS-project/1-miRNA/ref/Monochamus_version1.fa \
-en -20 \
-out LKY_Ep_Vs_NJ_Ep_miranda_targetGenes

~/miniconda3/bin/miranda /home/yjxiang/NGS-project/1-miRNA/diffmiRNA/LKY_Ep_Vs_NJ_Mg.seq /home/yjxiang/NGS-project/1-miRNA/ref/Monochamus_version1.fa \
-en -20 \
-out LKY_Ep_Vs_NJ_Mg_miranda_targetGenes

~/miniconda3/bin/miranda /home/yjxiang/NGS-project/1-miRNA/diffmiRNA/NJ_Ep_Vs_LKY_Mg.seq /home/yjxiang/NGS-project/1-miRNA/ref/Monochamus_version1.fa \
-en -20 \
-out NJ_Ep_Vs_LKY_Mg_miranda_targetGenes

~/miniconda3/bin/miranda /home/yjxiang/NGS-project/1-miRNA/diffmiRNA/NJ_Mg_Vs_LKY_Mg.seq /home/yjxiang/NGS-project/1-miRNA/ref/Monochamus_version1.fa \
-en -20 \
-out NJ_Mg_Vs_LKY_Mg_miranda_targetGenes

~/miniconda3/bin/miranda /home/yjxiang/NGS-project/1-miRNA/diffmiRNA/NJ_Mg_Vs_NJ_Ep.seq /home/yjxiang/NGS-project/1-miRNA/ref/Monochamus_version1.fa \
-en -20 \
-out NJ_Mg_Vs_NJ_Ep_miranda_targetGenes
