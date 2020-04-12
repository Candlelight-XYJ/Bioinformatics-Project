## step1 download genome data
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz

## step2 create samtools index
~/miniconda3/bin/samtools faidx E.coli_K12_MG1655.fa

# 有了索引之后，可以快速提取fasta某特定区域的序列
# 例如提取整个ｆａｓｔａ中序列ID含有NC_000913.3的ｆａｓｔａ序列
samtools faidx E.coli_K12_MG1655.fa NC_000913.3　>　testt

## step3 download sequencing data
# 首先下载安装sratoolkit
# 转sra为fastq文件
/home/yjxiang/tools/sratoolkit.2.9.6-1-centos_linux64/bin/fastq-dump --split-files SRR9723002
# 使用bgzip 压缩.fq为.gz文件
~/miniconda3/bin/bgzip -f SRR9723002_1.fastq
~/miniconda3/bin/bgzip -f SRR9723002_2.fastq

## step4 quality control 略

## step5 mapping
# build index for genome
/home/miniconda2/bin/bwa index E.coli_K12_MG1655.fa

# mapping
cd /home/yjxiang/NGS-project/5-GATK_practice/input/
time /home/miniconda2/bin/bwa mem -t 4 -R '@RG\tID:foo\tPL:illumina\tSM:E.coli_K12' \
E.coli_K12_MG1655.fa SRR9723002_1.fastq.gz SRR9723002_2.fastq.gz | ~/miniconda3/bin/samtools view -Sb - > ../output/E_coli_K12.bam && echo "** bwa mapping done **"

# samtools sort bam file
cd /home/yjxiang/NGS-project/5-GATK_practice/output/
time /home/miniconda2/bin/samtools sort -@ 4 -m 4G \
-O bam -o E_coli_K12.sorted.bam E_coli_K12.bam && echo "** BAM sort done"

# flag PCR replicates
cd /home/yjxiang/NGS-project/5-GATK_practice/output/
time /home/software/install/GATK/GATK-4.0.10.1/gatk-4.0.10.1/gatk MarkDuplicates \
-I E_coli_K12.sorted.bam -O E_coli_K12.sorted.markdup.bam -M E_coli_K12.sorted.markdup_metrics.txt && echo "** markdup done **"

# index bam 
time /home/miniconda2/bin/samtools index /home/yjxiang/NGS-project/5-GATK_practice/output/E_coli_K12.sorted.markdup.bam && echo "** index done **"


## step6 Variant calling
# create dictionary
cd /home/yjxiang/NGS-project/5-GATK_practice/input/
/home/software/install/GATK/GATK-4.0.10.1/gatk-4.0.10.1/gatk CreateSequenceDictionary \
-R E.coli_K12_MG1655.fa \
-O E.coli_K12_MG1655.dict \
&& echo "** dict done **"

# 生成中间文件gvcf
time /home/software/install/GATK/GATK-4.0.10.1/gatk-4.0.10.1/gatk HaplotypeCaller \
 -R /home/yjxiang/NGS-project/5-GATK_practice/input/E.coli_K12_MG1655.fa \
 --emit-ref-confidence GVCF \
 -I /home/yjxiang/NGS-project/5-GATK_practice/output/E_coli_K12.sorted.markdup.bam \
 -O /home/yjxiang/NGS-project/5-GATK_practice/output/E_coli_K12.g.vcf && echo "** gvcf done **"

# 通过gvcf检测变异
time /home/software/install/GATK/GATK-4.0.10.1/gatk-4.0.10.1/gatk GenotypeGVCFs \
 -R /home/yjxiang/NGS-project/5-GATK_practice/input/E.coli_K12_MG1655.fa \
 -V /home/yjxiang/NGS-project/5-GATK_practice/output/E_coli_K12.g.vcf \
 -O /home/yjxiang/NGS-project/5-GATK_practice/output/E_coli_K12.vcf && echo "** vcf done **"


