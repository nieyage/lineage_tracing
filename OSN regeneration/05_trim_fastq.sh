# trim fastq files
## PBS configure
#PBS -N trim_fastq
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=6
#PBS -l mem=10G

cd /md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/rawdata/scATAC

cat OSN_scATAC.input|while read id;

do echo $id
arr=($id)
sample=${arr[0]}
input1=${arr[1]}
input2=${arr[2]}

java -jar /md01/nieyg/ori/biosoft/package/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
 -phred33 $input1 $input2 ./trim_data/$input1 ./trim_data/$sample_unpaired_R1.fq.gz ./trim_data/$input2 ./trim_data/$sample_unpaired_R3.fq.gz LEADING:35 TRAILING:35 CROP:100 MINLEN:35 AVGQUAL:30 -threads 6
done

#!/bin/bash

# List of sample names to process
#sample_names=("S1" "S2" "S3" "S4" )
sample_names=("S5" "S6" "S7" "S8" "S9")

for sample_name in "${sample_names[@]}"; do

    # Decompress R1 and R2 FASTQ files
    gunzip -c "/data/R01/yangjw28/data/OSN/jointpolgF1met5d/scATAC/test/A-52-scATAC_${sample_name}_L001_I2_001.fastq.gz" > "A-52-scATAC_${sample_name}_L001_I2_001.fastq"
    gunzip -c "A-52-scATAC_${sample_name}_L001_R1_001.fastq.gz" > "A-52-scATAC_${sample_name}_L001_R1_001.fastq"

    # Use fastq_pair to pair the reads
    fastq_pair "A-52-scATAC_${sample_name}_L001_I2_001.fastq" "A-52-scATAC_${sample_name}_L001_R1_001.fastq"

    # Remove temporary files
    rm "A-52-scATAC_${sample_name}_L001_I2_001.fastq"
    rm "A-52-scATAC_${sample_name}_L001_R1_001.fastq"

    # Rename the paired file
    mv "A-52-scATAC_${sample_name}_L001_I2_001.fastq.paired.fq" "A-52-scATAC_${sample_name}_L001_I2_001.fastq"

    # Compress the paired file
    gzip "A-52-scATAC_${sample_name}_L001_I2_001.fastq"
done


# cellranger for trimed data

# Step1: bwa realign the chrM bam file 
nohup samtools view -b atac_possorted_bam.bam chrM  > OSN_100bp_possorted_chrM.bam &
nohup samtools index OSN_100bp_possorted_chrM.bam &

nohup samtools collate -Oun128 OSN_100bp_possorted_chrM.bam | samtools fastq -OT RG,CB,CR,CY,TR,TQ - \
  | bwa mem -pt8 -CH <(samtools view -H OSN_100bp_possorted_chrM.bam|grep ^@RG) /md01/nieyg/ref/hard-mask/bwa_index/assembly_hardmask.fa - \
  | samtools sort -@4 -m4g -o OSN_100bp_chrM_bwa.bam - &

conda activate r4-base

# in R 
library(Rsamtools)
library(ggplot2)
bam_file <- "OSN_100bp_chrM_bwa.bam"  # 替换为你的 BAM 文件路径
bam <- scanBam(bam_file)
mapping_qualities <- as.data.frame(as(bam$qual, "integer"))
pdf("./OSN_100bp_chrM_bwa_Mapping Quality Distribution.pdf",width=6,height=4)
ggplot(mapping_qualities, aes(x = V1)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(title = "Mapping Quality Distribution", x = "Mapping Quality", y = "Frequency") +
  theme_minimal()
dev.off()


###### not rmdup ######

# Step3: Filter chrM reads, q >60,PE for SNP calling  
# remove reads mapped with different chrom

samtools view -b -q  60 -f 0x2 OSN_100bp_chrM_bwa.bam -o OSN_100bp_chrM_bwa_chrMQ60PE.bam chrM
samtools index OSN_100bp_chrM_bwa_chrMQ60PE.bam
samtools view -H OSN_100bp_chrM_bwa_chrMQ60PE.bam > header.sam
samtools view  -@ 12 OSN_100bp_chrM_bwa_chrMQ60PE.bam |awk '{if ($7 == "=") print $0}' | grep -v "MM:i:1" | cat header.sam - | samtools view -Sb -@ 12  - > OSN_100bp_chrM_bwa_chrMQ60PE.unique.bam
samtools index OSN_100bp_chrM_bwa_chrMQ60PE.unique.bam
samtools sort -@ 12 OSN_100bp_chrM_bwa_chrMQ60PE.unique.bam -o OSN_100bp_chrM_bwa_chrMQ60PE_filtered.bam
samtools index OSN_100bp_chrM_bwa_chrMQ60PE_filtered.bam

# Step4: make re-align target regions
java -Xmx8g -jar /md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar \
-R /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa  \
-T RealignerTargetCreator  -nt 32 \
-I OSN_100bp_chrM_bwa_chrMQ60PE_filtered.bam \
-o OSN_100bp_chrM_bwa_chrMQ60PE_filtered.bam.target.intervals

nohup java -Xmx100g -jar /md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar \
-R /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa  \
-T IndelRealigner -filterNoBases -maxReads 10000000  \
-I OSN_100bp_chrM_bwa_chrMQ60PE_filtered.bam \
-targetIntervals OSN_100bp_chrM_bwa_chrMQ60PE_filtered.bam.target.intervals \
-o OSN_100bp_chrM_bwa_chrMQ60PE_filtered.realign.bam &

samtools view -h OSN_100bp_chrM_bwa_chrMQ60PE_filtered.realign.bam | awk '{if ($0 !~ /NM:i:[2-9]/ && $0 !~ /NM:i:[0-9][0-9]/) print $0}' | samtools view -b -o OSN_100bp_chrM_bwa_chrMQ60PE_filtered.realign_filtered_unique_mismatches.bam -


# mgatk ppl
#PBS -N OSN_not_rmdup
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=20
#PBS -l mem=64G
source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked_add_trimed/outs/mgatk
source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate r4-base
source /md01/nieyg/software/venv3/bin/activate
mgatk tenx -i /md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked_add_trimed/outs/OSN_100bp_chrM_bwa_chrMQ60PE_filtered.realign_filtered_unique_mismatches.bam \
-o  /md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked_add_trimed/outs/mgatk \
-n OSN_100bp_chrM_bwa -g /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask_chrM.fasta \
-c 10 -bt CB -b /md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked/mgatk_for_filtered_cells/barcode.tsv



