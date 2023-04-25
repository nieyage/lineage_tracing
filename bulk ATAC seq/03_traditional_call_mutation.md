# After get the alignment bam file

# Step1: get the bam file 
samtools index AR3-1-5.sorted.bam 
# Step2: picard remove the replicates
picard MarkDuplicates CREATE_INDEX=true ASSUME_SORTED=true \
VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true \
INPUT=AR3-1-5.sorted.bam \
OUTPUT=AR3-1-5_sorted_rmdup.bam  \
METRICS_FILE=AR3-1-5_sorted.rmdup.metrics

# Step3: make re-align target regions
cd ../7_mutation
picard AddOrReplaceReadGroups \
      I=../2_bam/AR3-1-5_sorted_rmdup.bam \
      O=AR3-1-5_sorted_rmdup.AddOrReplaceReadGroups.bam \
      RGID=4 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=20
samtools index AR3-1-5_sorted_rmdup.AddOrReplaceReadGroups.bam

java -Xmx8g -jar /md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar -R /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa  -T RealignerTargetCreator  -nt 32 -I AR3-1-5_sorted_rmdup.AddOrReplaceReadGroups.bam  -o AR3-1-5_sorted_rmdup.bam.target.intervals

java -Xmx100g -jar /md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar -R /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa -T IndelRealigner -filterNoBases -maxReads 10000000  -I AR3-1-5_sorted_rmdup.AddOrReplaceReadGroups.bam -targetIntervals AR3-1-5_sorted_rmdup.bam.target.intervals -o AR3-1-5_sorted_rmdup.realign.bam

# Step4: mpileup get file for snv calling

samtools mpileup  -l /md01/nieyg/ref/hard-mask/genome_modify/chrM.len -f /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa -q 10 -Q 10 AR3-1-5_sorted_rmdup.realign.bam >  AR3-1-5.mpileup

# make count table for each allele 
/md01/nieyg/software/ATAC_mito_sc-master/src/pileup_inf_rj.pl AR3-1-5.mpileup > AR3-1-5.count

# Step5: varscan call somatic mutation in chrM
varscan  pileup2snp AR3-1-5.mpileup --min-var-freq 0.000001  --min-reads2 2 > AR3-1-5.snv
