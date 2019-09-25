## 1.转换坐标
# 果壳给的数据已经是hg19坐标系的了, 所以我们不用转换

## 2.将.ped和.map文件转换为IMPUTE格式
# 使用Gtool进行转换
/root/tools/gtool -P --ped data/test.ped --map data/test.map --og data/out.gen --os data/out.sample

## 3.使用IMPUTE2进行phasing
# imputes2的phase选项可用于完成phasing
# Example
./impute2 \
 -phase \
 -m data/test.map \
 -g data/test.gen \
 -int 0 50000 \
 -Ne 20000 \
 -o data/test.phasing.impute2

## 4.使用IMPUTE2进行Imputation
# 下载用于imputation的参考数据集
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz

# Imputation
./impute2 \
 -use_prephased_g \
 -m data/test.map \
 -h ./Example/example.chr22.1kG.haps \
 -l ./Example/example.chr22.1kG.legend \
 -known_haps_g data/test.phasing.impute2_haps \
 -strand_g ./Example/example.chr22.study.strand \
 -int 0 50000 \
 -Ne 20000 \
 -o data/imputated_result \
 




