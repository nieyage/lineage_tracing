# After get the alignment bam file

# Step1: get the bam file 
samtools index atac_possorted_bam.bam
# Step2: picard remove the replicates
picard MarkDuplicates CREATE_INDEX=true ASSUME_SORTED=true \
VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true \
INPUT=atac_possorted_bam.bam \
OUTPUT=atac_possorted_bam_rmdup.bam  \
METRICS_FILE=atac_possorted_bam.rmdup.metrics

# Step3: make re-align target regions

picard AddOrReplaceReadGroups \
      I=atac_possorted_bam_rmdup.bam \
      O=atac_possorted_bam_rmdup.AddOrReplaceReadGroups.bam \
      RGID=4 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=20
samtools index atac_possorted_bam_rmdup.AddOrReplaceReadGroups.bam

java -Xmx8g -jar /md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar -R /md01/nieyg/ref/hard-mask/mm10_hard_masked/fasta/genome.fa  -T RealignerTargetCreator  -nt 32 -I atac_possorted_bam_rmdup.AddOrReplaceReadGroups.bam  -o atac_possorted_bam_rmdup.target.intervals

java -Xmx100g -jar /md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar -R /md01/nieyg/ref/hard-mask/mm10_hard_masked/fasta/genome.fa  -T IndelRealigner -filterNoBases -maxReads 10000000  -I atac_possorted_bam_rmdup.AddOrReplaceReadGroups.bam -targetIntervals atac_possorted_bam_rmdup.target.intervals -o atac_possorted_bam_rmdup.realign.bam

# Step4: mpileup get file for snv calling

samtools mpileup  -l /md01/nieyg/ref/hard-mask/genome_modify/chrM.len -f /md01/nieyg/ref/hard-mask/mm10_hard_masked/fasta/genome.fa -q 10 -Q 10 atac_possorted_bam_rmdup.realign.bam >  atac_possorted_bam_rmdup.mpileup

# make count table for each allele 
/md01/nieyg/software/ATAC_mito_sc-master/src/pileup_inf_rj.pl atac_possorted_bam_rmdup.mpileup > atac_possorted_bam_rmdup.count

# Step5: varscan call somatic mutation in chrM
varscan  pileup2snp atac_possorted_bam_rmdup.mpileup --min-var-freq 0.000001  --min-reads2 2 > atac_possorted_bam_rmdup_all.snv

varscan  mpileup2snp atac_possorted_bam_rmdup.mpileup --variants SNP --p-value 0.99 --output-vcf 1 --min-coverage 3 --min-var-freq 0.0005  --min-reads2 2 > masked_atac_possorted_bam_rmdup_filtered.vcf


# filter somatic mutation  --step6
# STEP 1: filter mtSNV
# $17>=1 && $18>=1 ( make sure both positive and negative strand were called 
# Minor  Allele frequency more than 0.01 
# remove G->T and C->A, "N",based on previous empirical evidence. 
# sequence depth >= 20 for each cell, the sequence depth criteria should be adjust according to the data. 

awk '$17>=1 && $18>=1' atac_possorted_bam_rmdup_all.snv | sed '1,1d' |awk '$6/($6+$5)>=0.10' |awk '!($3=="G" && $6=="T")' |awk '!($3=="C" && $6=="A")' |awk '$3!="N"'   | awk ' $6+$5>=20' >>atac_possorted_bam_rmdup.mtDNA.snv.filter

awk '{print $1"\t"$2"\t"$3"\t"$19}' AR3-0.5-7-2.mtDNA.snv.filter|sort -k2n -u >AR3-0.5-7-2.mtDNA.snv.filter.cut.uniq
