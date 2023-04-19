
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



# Cromwell
#cp /md01/nieyg/software/gatk-master/scripts/mitochondria_m2_wdl/ExampleInputsMitochondriaPipeline.json test.json
#
#cromwell run /md01/nieyg/software/gatk-master/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl \
# --inputs test.json --labels udocker.config
#
#  default = "Local"
#  providers {
#    Local {
#      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
#      config {
#        runtime-attributes = """
#        String? docker
#        String? docker_user
#        """
#        submit-docker = """docker run --rm ${ "--user " + docker_user } -v ${cwd}:${docker_cwd} -i ${docker} /bin/bash ${docker_cwd}/execution/script"""
#        .
#        .
#        .
#
#
## Run the Mitochondrial short variant discovery (SNVs + Indels) step by step 
## 1. Subset to keep only the reads mapped to the mitochondrial genome
##    SubsetBamToChrM：
#gatk PrintReads \
#      -R /data/R02/nieyg/software/gatk-master/scripts/mitochondria_m2_wdl/ref/hg38/Homo_sapiens_assembly38.fasta \
#      -L chrM \
#      --read-filter MateOnSameContigOrNoMappedMateReadFilter \
#      --read-filter MateUnmappedAndUnmappedReadFilter \
#      -I /md01/nieyg/project/spatial_ATAC_mmhs/data/rawdata/Processed_data/human/Hippo/SRR18745417/CellRanger/SRR18745417/outs/possorted_bam.bam \
#      --read-index /md01/nieyg/project/spatial_ATAC_mmhs/data/rawdata/Processed_data/human/Hippo/SRR18745417/CellRanger/SRR18745417/outs/possorted_bam.bam.bai \
#      -O ./possorted_bam_chrM.bam
#
## 2. Revert the ChrM mapped reads from an aligned BAM to an unaligned BAM file
##    RevertSam: Removes alignment information while retaining recalibrated base qualities and original alignment tags
#picard  RevertSam \
#    INPUT=./possorted_bam_chrM.bam \
#    OUTPUT_BY_READGROUP=false \
#    OUTPUT=possorted_bam_chrM_revert.bam \
#    VALIDATION_STRINGENCY=LENIENT \
#    ATTRIBUTE_TO_CLEAR=FT \
#    ATTRIBUTE_TO_CLEAR=CO \
#    SORT_ORDER=queryname \
#    RESTORE_ORIGINAL_QUALITIES=false
#
## 3. Align the unmapped BAM file with the reference aligned BAM and shifted reference aligned BAM
##    AlignAndCall.AlignAndCall
## 1) AlignToMt:
#
#    # set the bash variable needed for the command-line
#    bash_ref_fasta=/data/R02/nieyg/software/gatk-master/scripts/mitochondria_m2_wdl/ref/hg38/Homo_sapiens_assembly38.chrM.fasta
#    picard   SamToFastq \
#      INPUT=possorted_bam_chrM_revert.bam \
#      FASTQ=/dev/stdout \
#      INTERLEAVE=true \
#      NON_PF=true | \
#    bwa mem -K 100000000 -p -v 3 -t 2 -Y /data/R02/nieyg/software/gatk-master/scripts/mitochondria_m2_wdl/ref/hg38/Homo_sapiens_assembly38.chrM.fasta /dev/stdin - 2> >(tee possorted_bam_chrM_aligntoMT.bwa.stderr.log >&2) | \
#    picard  MergeBamAlignment \
#      VALIDATION_STRINGENCY=SILENT \
#      EXPECTED_ORIENTATIONS=FR \
#      ATTRIBUTES_TO_RETAIN=X0 \
#      ATTRIBUTES_TO_REMOVE=NM \
#      ATTRIBUTES_TO_REMOVE=MD \
#      ALIGNED_BAM=/dev/stdin \
#      UNMAPPED_BAM=possorted_bam_chrM_revert.bam \
#      OUTPUT=mba.bam \
#      REFERENCE_SEQUENCE=/data/R02/nieyg/software/gatk-master/scripts/mitochondria_m2_wdl/ref/hg38/Homo_sapiens_assembly38.chrM.fasta \
#      PAIRED_RUN=true \
#      SORT_ORDER="unsorted" \
#      IS_BISULFITE_SEQUENCE=false \
#      ALIGNED_READS_ONLY=false \
#      CLIP_ADAPTERS=false \
#      MAX_RECORDS_IN_RAM=2000000 \
#      ADD_MATE_CIGAR=true \
#      MAX_INSERTIONS_OR_DELETIONS=-1 \
#      PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
#      PROGRAM_RECORD_ID="bwamem" \
#      PROGRAM_GROUP_VERSION="0.7.17-r1188" \
#      PROGRAM_GROUP_COMMAND_LINE="bwa mem -K 100000000 -p -v 3 -t 2 -Y $bash_ref_fasta" \
#      PROGRAM_GROUP_NAME="bwamem" \
#      UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
#      ALIGNER_PROPER_PAIR_FLAGS=true \
#      UNMAP_CONTAMINANT_READS=true \
#      ADD_PG_TAG_TO_READS=false
#
#    picard  MarkDuplicates \
#      INPUT=mba.bam \
#      OUTPUT=possorted_bam_chrM_mba_md.bam \
#      METRICS_FILE=possorted_bam_chrM_mba_md_duplicate_metrics.txt \
#      VALIDATION_STRINGENCY=SILENT \
#      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
#      ASSUME_SORT_ORDER="queryname" \
#      CLEAR_DT="false" \
#      ADD_PG_TAG_TO_READS=false
#
#    picard SortSam \
#      INPUT=possorted_bam_chrM_mba_md.bam \
#      OUTPUT=possorted_bam_chrM_mba_md.realigned.bam \
#      SORT_ORDER="coordinate" \
#      CREATE_INDEX=true \
#      MAX_RECORDS_IN_RAM=300000
#
## 2) AlignToShiftedMt:
#
#    bash_ref_fasta=/data/R02/nieyg/software/gatk-master/scripts/mitochondria_m2_wdl/ref/hg38/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta
#    picard   SamToFastq \
#      INPUT=possorted_bam_chrM_revert.bam \
#      FASTQ=/dev/stdout \
#      INTERLEAVE=true \
#      NON_PF=true | \
#    bwa mem -K 100000000 -p -v 3 -t 2 -Y /data/R02/nieyg/software/gatk-master/scripts/mitochondria_m2_wdl/ref/hg38/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta /dev/stdin - 2> >(tee ~{possorted_bam_chrM}.shifted_by_8000_bases.bwa.stderr.log >&2) | \
#    picard   MergeBamAlignment \
#      VALIDATION_STRINGENCY=SILENT \
#      EXPECTED_ORIENTATIONS=FR \
#      ATTRIBUTES_TO_RETAIN=X0 \
#      ATTRIBUTES_TO_REMOVE=NM \
#      ATTRIBUTES_TO_REMOVE=MD \
#      ALIGNED_BAM=/dev/stdin \
#      UNMAPPED_BAM= possorted_bam_chrM_revert.bam \
#      OUTPUT=possorted_bam_chrM_mba.shifted_by_8000_bases.bam \
#      REFERENCE_SEQUENCE=/data/R02/nieyg/software/gatk-master/scripts/mitochondria_m2_wdl/ref/hg38/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta \
#      PAIRED_RUN=true \
#      SORT_ORDER="unsorted" \
#      IS_BISULFITE_SEQUENCE=false \
#      ALIGNED_READS_ONLY=false \
#      CLIP_ADAPTERS=false \
#      MAX_RECORDS_IN_RAM=2000000 \
#      ADD_MATE_CIGAR=true \
#      MAX_INSERTIONS_OR_DELETIONS=-1 \
#      PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
#      PROGRAM_RECORD_ID="bwamem" \
#      PROGRAM_GROUP_VERSION="0.7.17-r1188" \
#      PROGRAM_GROUP_COMMAND_LINE="bwa mem -K 100000000 -p -v 3 -t 2 -Y /data/R02/nieyg/software/gatk-master/scripts/mitochondria_m2_wdl/ref/hg38/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta" \
#      PROGRAM_GROUP_NAME="bwamem" \
#      UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
#      ALIGNER_PROPER_PAIR_FLAGS=true \
#      UNMAP_CONTAMINANT_READS=true \
#      ADD_PG_TAG_TO_READS=false
#
#    picard  MarkDuplicates \
#      INPUT=possorted_bam_chrM_mba.shifted_by_8000_bases.bam \
#      OUTPUT=possorted_bam_chrM_mba_md.shifted_by_8000_bases.bam \
#      METRICS_FILE=possorted_bam_chrM_mba_md_duplicate_metrics.shifted_by_8000_bases.txt \
#      VALIDATION_STRINGENCY=SILENT \
#      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
#      ASSUME_SORT_ORDER="queryname" \
#      CLEAR_DT="false" \
#      ADD_PG_TAG_TO_READS=false
#
#    picard SortSam \
#      INPUT=possorted_bam_chrM_mba_md.shifted_by_8000_bases.bam \
#      OUTPUT=possorted_bam_chrM_mba_md.shifted_by_8000_bases.realigned.bam \
#      SORT_ORDER="coordinate" \
#      CREATE_INDEX=true \
#      MAX_RECORDS_IN_RAM=300000
#
## 3) CollectWgsMetrics:
#    picard CollectWgsMetrics \
#      INPUT= possorted_bam_chrM_mba_md.realigned.bam \
#      VALIDATION_STRINGENCY=SILENT \
#      REFERENCE_SEQUENCE=/data/R02/nieyg/software/gatk-master/scripts/mitochondria_m2_wdl/ref/hg38/Homo_sapiens_assembly38.chrM.fasta \
#      OUTPUT=metrics.txt \
#      USE_FAST_ALGORITHM=true \
#      INCLUDE_BQ_HISTOGRAM=true \
#      THEORETICAL_SENSITIVITY_OUTPUT=theoretical_sensitivity.txt
#
#    R --vanilla <<CODE
#      df = read.table("metrics.txt",skip=6,header=TRUE,stringsAsFactors=FALSE,sep='\t',nrows=1)
#      write.table(floor(df[,"MEAN_COVERAGE"]), "mean_coverage.txt", quote=F, col.names=F, row.names=F)
#      write.table(df[,"MEDIAN_COVERAGE"], "median_coverage.txt", quote=F, col.names=F, row.names=F)
#    CODE
#
#
#CallMt
#
#
#CallShiftedMt
#
#LiftoverAndCombineVcfs
#
#MergeStats
#
#InitialFilter
#
#SplitMultiAllelicsAndRemoveNonPassSites
#
#GetContamination
#
#FilterContamination
#CoverageAtEveryBase
#
#
#
#
#
#
#

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





