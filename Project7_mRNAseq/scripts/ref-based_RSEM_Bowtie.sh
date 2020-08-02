#### Step1 Building Index
/home/xiangyj/tool/RSEM-1.2.25/rsem-prepare-reference -p 4 \
--gtf /home/xiangyj/ref/XXX.chr.gtf \
--bowtie2 \
--bowtie2-path ~/miniconda3/bin \
/home/xiangyj/ref/refGenome \
/home/xiangyj/index/refIndex

#### Step2 Mapping and quant
cd /home/ningjing/C5_cleandata/china/
for fn in C5SX{1,2,3};
do
/home/xiangyj/tool/RSEM-1.2.25/rsem-calculate-expression --paired-end -p 4 \
    --bowtie2 \
    --bowtie2-path ~/miniconda3/bin \
    --append-names --no-bam-output \
    ${fn}_1.clean.fq.gz \
    ${fn}_2.clean.fq.gz \
    /home/xiangyj/index/refIndex \
    /home/rsem/${fn}_rsem_out
done


#### Step3 extract feature matrix
cd /home/rsem/
/home/xiangyj/tool/RSEM-1.2.25/rsem-generate-data-matrix *_rsem_out > featureMatrix.txt