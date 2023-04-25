
# the first pipeline 
# Step1: get the chrM bam file 
samtools view -b possorted_bam.bam chrM  > possorted_bam.chrM.bam
samtools index possorted_bam.chrM.bam

# Step2: picard remove the replicates
picard MarkDuplicates CREATE_INDEX=true ASSUME_SORTED=true \
VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true \
INPUT=possorted_bam.chrM.bam \
OUTPUT=possorted_bam.chrM.rmdup.bam  \
METRICS_FILE=possorted_bam.chrM.rmdup.metrics

# Step3: split bam by barcode 
# split_bam_byBC.sh
mkdir split_bam_byBC
nohup samtools view possorted_bam.chrM.rmdup.bam| cut -f 12 | tr "\t" "\n"  | grep  "^CB:Z:"  | cut -d ':' -f 3 | sort | uniq | while read S; do samtools view -h possorted_bam.chrM.rmdup.bam |  awk -v tag="CB:Z:$S" '($0 ~ /^@/ || index($0,tag)>0)' | samtools view -@8 -b > ./split_bam_byBC/${S}.bam ; done &

# Step4: for every barcode
mkdir split_bam_byBC_snv

ls ./split_bam_byBC/*.bam| cut -d "/" -f 3|cut -d "." -f 1 | while read id ; 
do 
samtools mpileup -l /md01/nieyg/ref/10X/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/star/chrMNameLength.txt \
-q 30 -Q 30 -f /md01/nieyg/ref/10X/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa  \
-x ./split_bam_byBC/$id.bam > ./split_bam_byBC_snv/$id.mpileup
VarScanr pileup2snp ./split_bam_byBC_snv/$id.mpileup --min-var-freq 0.01  --min-reads2 2 >  ./split_bam_byBC_snv/$id.mpileup.snv
/md01/jinxu/bin/pileup_inf_rj.pl ./split_bam_byBC_snv/$id.mpileup > ./split_bam_byBC_snv/$id.mpileup.counts
done

#
mkdir split_bam_gatk

ls ./split_bam_byBC/*.bam| cut -d "/" -f 3|cut -d "." -f 1 | while read id ; 
do 
samtools index ./split_bam_byBC/$id.bam 
gatk HaplotypeCaller \
-R /md01/nieyg/ref/10X/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
-L chrM -I ./split_bam_byBC/$id.bam  -O ./split_bam_gatk/$id.gvcf --emit-ref-confidence GVCF

gatk GenotypeGVCFs -R /md01/nieyg/ref/10X/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
-V ./split_bam_gatk/$id.gvcf -O ./split_bam_gatk/$id.vcf

gatk SelectVariants -V ./split_bam_gatk/$id.vcf -O ./split_bam_gatk/$id.snp.vcf --select-type-to-include SNP

done



gatk HaplotypeCaller \
-R /md01/nieyg/ref/10X/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
-L chrM \
-I AAACATCGAAACATCG-1.bam \
-O AAACATCGAAACATCG-1.bam.HC.gvcf \
--emit-ref-confidence GVCF ## 生成中间文件gvcf

##通过gvcf检测变异, -V 添加上步得到的gvcf
gatk GenotypeGVCFs -R /md01/nieyg/ref/10X/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
-V AAACATCGAAACATCG-1.bam.HC.gvcf \
-O AAACATCGAAACATCG-1.bam.HC.vcf


## 提取SNP
gatk SelectVariants -V AAACATCGAAACATCG-1.bam.HC.vcf -O AAACATCGAAACATCG-1.bam.HC.snp.vcf --select-type-to-include SNP


# combin gvcf script 
echo 'gatk CombineGVCFs \' > combine_gatk_gvcf.sh
echo '-R /md01/nieyg/ref/10X/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \'  >> combine_gatk_gvcf.sh
ls ./split_bam_gatk/*.gvcf|while read id;
do 
  echo '-V '${id}' \' >> combine_gatk_gvcf.sh;
done;
echo '-O gatk_combined.vcf' >> combine_gatk_gvcf.sh;


#gatk CombineGVCFs -V ./split_bam_gatk/GCCACATACCGAAGTA-1.gvcf ./split_bam_gatk/GCCACATACCGTGAGA-1.gvcf  -O gatk_combined.gvcf -R  /md01/nieyg/ref/10X/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa




## 提取INDEL
gatk SelectVariants -V AAACATCGAAACATCG-1.bam.HC.vcf -O  AAACATCGAAACATCG-1.bam.HC.indel.vcf --select-type-to-include INDEL

gatk VariantFiltration -O test.snp.fil.vcf.temp -V test.snp.vcf --filter-expression 'QUAL < 30.0 || QD < 2.0 || FS > 60.0 ||  SOR > 4.0' \
　　　　--filter-name lowQualFilter --cluster-window-size 10  --cluster-size 3 --missing-values-evaluate-as-failing

## 参数
-O 输出filt.vcf文件
-V 输入vcf文件
--filter-expression 过滤条件, VCF INFO 信息
--cluster-window-size 以10个碱基为一个窗口
--cluster-size 10个碱基为窗口，若存在3以上个则过滤
--filter-name 被过滤掉的SNP不会删除，而是给一个标签， 比如 Filter
--missing-values-evaluate-as-failing 当筛选标准比较多的时候，可能有一些位点没有筛选条件当中的一条或几条，例如下面的这个表达式；QUAL < 30.0 || QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 并不一定所有位点都有这些信息，这种情况下GATK运行的时候会报很多WARNING信息，用这个参数可以把这些缺少某些FLAG的位点也给标记成没有通过筛选的。
## QualByDepth(QD): 变异位点可信度除以未过滤的非参考read数
## FisherStrand (FS): Fisher精确检验评估当前变异是strand bias的可能性，这个值在0-60间
# RMSMappingQuality (MQ): 所有样本中比对质量的平方根
# MappingQualityRankSumTest (MQRankSum): 根据REF和ALT的read的比对质量来评估可信度
# ReadPosRankSumTest (ReadPosRankSum) : 通过变异在read的位置来评估变异可信度，通常在read的两端的错误率比较高
# StrandOddsRatio (SOR) : 综合评估strand bias的可能性



## 根据FILTER那列信息进行筛选
grep PASS test.snp.fil.vcf.temp >  test.snp.fil.vcf





