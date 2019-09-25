#http://zzz.bwh.harvard.edu/plink/pimputation.shtml

#### quality control 

## 0.准备文件
#1.遗传图谱 genetic map文件, download : ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz
2.将果壳的 genotype 数据转换为ped格式
./vcftools --vcf input_data.vcf --plink --chr 1 --out output_in_plink
# Family ID, Individual ID, Paternal ID, Maternal ID, Sex (1=male; 2=female; other=unknown), Phenotype


## 1.转换坐标
#果壳给的数据已经是hg19坐标系的了, 所以我们不用转换

## 2.将Genotype文件转换为IMPUTE格式
#使用Gtool进行转换
gtool -P --ped example/example.ped --map example/example.map --og example/out.gen --os example/out.sample

#### phasing
# imputes2的phase选项可用于完成phasing
# Example
./impute2 \
 -phase \
 -m ./Example/example.chr22.map \
 -g ./Example/example.chr22.study.gens \
 -int 20.4e6 20.5e6 \
 -Ne 20000 \
 -o ./Example/example.chr22.phasing.impute2


#### imputation
