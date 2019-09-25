#公共conda地址
/home/miniconda2/bin/
# my
~/miniconda3/bin/conda

#######################
## 公共数据库进行搜索 ##  
#######################
srun -p husn -c 1 -J blastdb sh blast.sh  &

## 下载公共数据库

~/miniconda3/bin/update_blastdb.pl nt > log
#/home/miniconda2/bin/update_blastdb.pl nt > log 

## 搜索
/home/miniconda2/bin/blastn -db nt -query nt.fsa -out results.out


########################
## 自构建数据库进行搜索 ##  
########################

srun -p husn -c 1 -J createdb sh createDB.sh  &

## 从fasta文件构建数据库
cd /home/yjxiang/skills/1-local_blast/
~/miniconda3/bin/makeblastdb -in miRNA_novel.fa -parse_seqids \
-title "Ma_mirna" \
-dbtype 'mirna'

## 搜索
srun -p husn -c 4 -J runBlast sh runBlast.sh &
cd /home/yjxiang/skills/1-local_blast/
~/miniconda3/bin/blastn -db miRNA_novel \
-query xianchong.fa \
-num_threads 4 \
-evalue 1e-6 \  ## 10^-6
-outfmt 7 \
-out queryResult.txt

################################
## microRNA短序列需要注意的参数 ##  
################################

# -task blastn-short -word_size 7 -evalue 1
srun -p husn -c 4 -J runBlast sh runBlast.sh &
cd /home/yjxiang/skills/1-local_blast/
~/miniconda3/bin/blastn -db tianniu \
-task blastn-short \
-word_size 20 \
-query xianchong \
-num_threads 4 \
-evalue 1 \ 
-outfmt 7 \
-out queryResult.txt

