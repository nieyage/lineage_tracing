# After get the alignment bam file

# Step1: get the bam file 
samtools index P7-15lab_L2.sorted.bam 
# Step2: picard remove the replicates
picard MarkDuplicates CREATE_INDEX=true ASSUME_SORTED=true \
VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true \
INPUT=P7-15lab_L2.sorted.bam \
OUTPUT=P7-15lab_L2_sorted_rmdup.bam  \
METRICS_FILE=P7-15lab_L2_sorted.rmdup.metrics

# Step3: make re-align target regions
cd ../7_mutation
picard AddOrReplaceReadGroups \
      I=../2_bam/P7-15lab_L2_sorted_rmdup.bam \
      O=P7-15lab_L2_sorted_rmdup.AddOrReplaceReadGroups.bam \
      RGID=4 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=20
samtools index P7-15lab_L2_sorted_rmdup.AddOrReplaceReadGroups.bam

java -Xmx8g -jar /md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar -R /md01/nieyg/ref/bowtie2/mm10/mm10.fa  -T RealignerTargetCreator  -nt 32 -I P7-15lab_L2_sorted_rmdup.AddOrReplaceReadGroups.bam  -o P7-15lab_L2_sorted_rmdup.bam.target.intervals

java -Xmx100g -jar /md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar -R /md01/nieyg/ref/hard-mask/genome_modify/old_genome.fa  -T IndelRealigner -filterNoBases -maxReads 10000000  -I P7-15lab_L2_sorted_rmdup.AddOrReplaceReadGroups.bam -targetIntervals P7-15lab_L2_sorted_rmdup.bam.target.intervals -o P7-15lab_L2_sorted_rmdup.realign.bam

# Step4: mpileup get file for snv calling

samtools mpileup  -l /md01/nieyg/ref/hard-mask/genome_modify/chrM.len -f /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa -q 10 -Q 10 ../2_bam/P7-15lab_L2_sorted_rmdup.realign.bam >  P7-15lab_L2.mpileup

samtools mpileup  -l /md01/nieyg/ref/hard-mask/genome_modify/chrM.len -f /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa -q 10 -Q 10 ../2_bam/P7-7Va_L2_sorted_rmdup.realign.bam >  P7-7Va_L2.mpileup


samtools mpileup  -l /md01/nieyg/ref/hard-mask/genome_modify/chrM.len -f /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa -q 10 -Q 10 ../2_bam/AR3-1-5_sorted_rmdup.realign.bam >  AR3-1-5.mpileup

samtools mpileup  -l /md01/nieyg/ref/hard-mask/genome_modify/chrM.len -f /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa -q 10 -Q 10 ../2_bam/AR3-0.5-7-2_sorted_rmdup.realign.bam > AR3-0.5-7-2.mpileup


# make count table for each allele 
/md01/nieyg/software/ATAC_mito_sc-master/src/pileup_inf_rj.pl P7-15lab_L2.mpileup > P7-15lab_L2.count
/md01/nieyg/software/ATAC_mito_sc-master/src/pileup_inf_rj.pl   P7-7Va_L2.mpileup   > P7-7Va_L2.count
/md01/nieyg/software/ATAC_mito_sc-master/src/pileup_inf_rj.pl     AR3-1-5.mpileup     > AR3-1-5.count
/md01/nieyg/software/ATAC_mito_sc-master/src/pileup_inf_rj.pl AR3-0.5-7-2.mpileup > AR3-0.5-7-2.count

# Step5: varscan call somatic mutation in chrM
varscan  pileup2snp P7-15lab_L2.mpileup --min-var-freq 0.000001  --min-reads2 2 > P7-15lab_L2_all.snv
varscan  pileup2snp   P7-7Va_L2.mpileup --min-var-freq 0.000001  --min-reads2 2 >   P7-7Va_L2_all.snv
varscan  pileup2snp     AR3-1-5.mpileup --min-var-freq 0.000001  --min-reads2 2 >     AR3-1-5_all.snv
varscan  pileup2snp AR3-0.5-7-2.mpileup --min-var-freq 0.000001  --min-reads2 2 > AR3-0.5-7-2_all.snv

varscan  pileup2snp P7-15lab_L2.mpileup --min-coverage 3 --min-var-freq 0.0005  --min-reads2 2 > P7-15lab_L2_filtered.snv
varscan  pileup2snp   P7-7Va_L2.mpileup --min-coverage 3 --min-var-freq 0.0005  --min-reads2 2 >   P7-7Va_L2_filtered.snv
varscan  pileup2snp     AR3-1-5.mpileup --min-coverage 3 --min-var-freq 0.0005  --min-reads2 2 >     AR3-1-5_filtered.snv
varscan  pileup2snp AR3-0.5-7-2.mpileup --min-coverage 3 --min-var-freq 0.0005  --min-reads2 2 > AR3-0.5-7-2_filtered.snv

# filter somatic mutation  --step6
# STEP 1: filter mtSNV
# $17>=1 && $18>=1 ( make sure both positive and negative strand were called 
# Minor  Allele frequency more than 0.01 
# remove G->T and C->A, "N",based on previous empirical evidence. 
# sequence depth >= 20 for each cell, the sequence depth criteria should be adjust according to the data. 

awk '$17>=1 && $18>=1' AR3-0.5-7-2_filtered.snv | sed '1,1d' |awk '$6/($6+$5)>=0.10' |awk '!($3=="G" && $6=="T")' |awk '!($3=="C" && $6=="A")' |awk '$3!="N"'   | awk ' $6+$5>=20' >>AR3-0.5-7-2.mtDNA.snv.filter

awk '{print $1"\t"$2"\t"$3"\t"$19}' AR3-0.5-7-2.mtDNA.snv.filter|sort -k2n -u >AR3-0.5-7-2.mtDNA.snv.filter.cut.uniq

awk '$17>=1 && $18>=1' AR3-1-5_filtered.snv | sed '1,1d' |awk '$6/($6+$5)>=0.10' |awk '!($3=="G" && $6=="T")' |awk '!($3=="C" && $6=="A")' |awk '$3!="N"'   | awk ' $6+$5>=20' >>AR3-1-5.mtDNA.snv.filter
awk '{print $1"\t"$2"\t"$3"\t"$19}' AR3-1-5.mtDNA.snv.filter|sort -k2n -u >AR3-1-5.mtDNA.snv.filter.cut.uniq

awk '$17>=1 && $18>=1' P7-15lab_L2_filtered.snv | sed '1,1d' |awk '$6/($6+$5)>=0.10' |awk '!($3=="G" && $6=="T")' |awk '!($3=="C" && $6=="A")' |awk '$3!="N"'   | awk ' $6+$5>=20' >>P7-15lab_L2.mtDNA.snv.filter
awk '{print $1"\t"$2"\t"$3"\t"$19}' P7-15lab_L2.mtDNA.snv.filter|sort -k2n -u >P7-15lab_L2.mtDNA.snv.filter.cut.uniq

awk '$17>=1 && $18>=1' P7-7Va_L2_filtered.snv | sed '1,1d' |awk '$6/($6+$5)>=0.10' |awk '!($3=="G" && $6=="T")' |awk '!($3=="C" && $6=="A")' |awk '$3!="N"'   | awk ' $6+$5>=20' >>P7-7Va_L2.mtDNA.snv.filter
awk '{print $1"\t"$2"\t"$3"\t"$19}' P7-7Va_L2.mtDNA.snv.filter|sort -k2n -u >P7-7Va_L2.mtDNA.snv.filter.cut.uniq


varscan  mpileup2snp P7-15lab_L2.mpileup --variants SNP --p-value 0.99 --output-vcf 1 --min-coverage 3 --min-var-freq 0.0005  --min-reads2 2 > P7-15lab_L2_filtered.vcf
varscan  mpileup2snp   P7-7Va_L2.mpileup --variants SNP --p-value 0.99 --output-vcf 1 --min-coverage 3 --min-var-freq 0.0005  --min-reads2 2 >   P7-7Va_L2_filtered.vcf
varscan  mpileup2snp     AR3-1-5.mpileup --variants SNP --p-value 0.99 --output-vcf 1 --min-coverage 3 --min-var-freq 0.0005  --min-reads2 2 >     AR3-1-5_filtered.vcf
varscan  mpileup2snp AR3-0.5-7-2.mpileup --variants SNP --p-value 0.99 --output-vcf 1 --min-coverage 3 --min-var-freq 0.0005  --min-reads2 2 > AR3-0.5-7-2_filtered.vcf

