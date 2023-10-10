# test the Mitosort 
# clone MitoSort repository 
git clone https://github.com/tangzhj/MitoSort.git

# create a conda environment and install all the required python packages
conda env create -f /md01/nieyg/software/MitoSort/MitoSort_env.yaml
conda activate MitoSort

picard CreateSequenceDictionary R=genome.fa O=genome.dict

## PBS configure 
#PBS -N mitosort_C4
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=16
#PBS -l mem=32G

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/mito_sort_out
#source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
#conda activate MitoSort

## realign MT reads 
python /md01/nieyg/software/MitoSort/MitoSort_pipeline.py mt-realign \
-b /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/outs/possorted_bam.bam \
-f /md01/nieyg/ref/hard-mask/mm10_hard_masked/fasta/genome.fa \
--gatk_path /md01/nieyg/software/MitoSort/GenomeAnalysisTK_3.5-0.jar \
-o /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/mito_sort_out

## generate SNP matrices
python /md01/nieyg/software/MitoSort/MitoSort_pipeline.py generate-snp-matrix \
-b /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/mito_sort_out/MitoSort/BAM/possorted_chrM_realign.bam \
-f /md01/nieyg/ref/hard-mask/mm10_hard_masked/fasta/genome.fa \
-c /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/outs/singlecell.csv \
-m /md01/nieyg/ref/hard-mask/genome_modify/chrM.len \
--varscan_path /md01/nieyg/software/MitoSort/VarScan.v2.3.7.jar \
-o /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/mito_sort_out

# test pipeline 
python /md01/nieyg/software/MitoSort/MitoSort_pipeline.py mt-realign \
-b /md01/nieyg/software/MitoSort/test_data/MitoSort_test_data/test_DOGMAseq_atac_possorted_chrM.bam \
-f /md01/nieyg/ref/10X/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
--gatk_path /md01/nieyg/software/MitoSort/GenomeAnalysisTK_3.5-0.jar \
-o .

python /md01/nieyg/software/MitoSort/MitoSort_pipeline.py generate-snp-matrix \
-b /md01/nieyg/software/MitoSort/test_data/MitoSort_test_data/MitoSort/BAM/possorted_chrM_realign.bam \
-f /md01/nieyg/ref/10X/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
-c /md01/nieyg/software/MitoSort/test_data/MitoSort_test_data/test_DOGMAseq_barcode.txt \
-m /md01/nieyg/software/MitoSort/data/hg38_chrM.bed \
--varscan_path /md01/nieyg/software/MitoSort/VarScan.v2.3.7.jar \
-o /md01/nieyg/software/MitoSort/test_data/MitoSort_test_data/

samtools view -s 0.2 -b /md01/nieyg/software/MitoSort/test_data/MitoSort_test_data/MitoSort/BAM/possorted_chrM_realign.bam > possorted_chrM_realign_0.2.bam
samtools mpileup -l /md01/nieyg/software/MitoSort/data/hg38_chrM.bed -q 30 -Q 30 \
-f /md01/nieyg/ref/10X/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
-x possorted_chrM_realign_0.2.bam -o possorted_chrM_realign_0.2.mpileup

python /md01/nieyg/software/MitoSort/MitoSort_pipeline.py generate-snp-matrix \
-b /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/mgatk/mito_sort_out/qc_bam/mgatk_merged.bam \
-f /md01/nieyg/ref/10X/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
-c /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/outs/singlecell.csv \
-m /md01/nieyg/ref/hard-mask/genome_modify/chrM.len \
--varscan_path /md01/nieyg/software/MitoSort/VarScan.v2.3.7.jar \
-o /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/mito_sort_out

python /md01/nieyg/software/MitoSort/MitoSort_pipeline.py generate-snp-matrix \
-b /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/mito_sort_out/MitoSort/BAM/possorted_chrM_realign.bam \
-f /md01/nieyg/ref/hard-mask/mm10_hard_masked/fasta/genome.fa \
-c /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/mgatk/barcode.tsv \
-m /md01/nieyg/ref/hard-mask/genome_modify/chrM.len \
--varscan_path /md01/nieyg/software/MitoSort/VarScan.v2.3.7.jar \
-o /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/mito_sort_out



## PBS configure 
#PBS -N mitosort_C4
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=16
#PBS -l mem=32G

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/mgatk_realign2/qc_bam


python /md01/nieyg/software/MitoSort/MitoSort_pipeline.py generate-snp-matrix \
-b /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/mgatk_realign2/qc_bam/possorted_chrM_realign.bam \
-f /md01/nieyg/ref/hard-mask/mm10_hard_masked/fasta/genome.fa \
-c /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/mgatk/barcode.tsv \
-m /md01/nieyg/ref/hard-mask/genome_modify/chrM.len \
--varscan_path /md01/nieyg/software/MitoSort/VarScan.v2.3.7.jar \
-o /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/mgatk_realign2/qc_bam


## PBS configure 
#PBS -N mitosort_C5
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=10
#PBS -l mem=32G

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/mgatk_realign/qc_bam

python /md01/nieyg/software/MitoSort/MitoSort_pipeline.py generate-snp-matrix \
-b /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/mgatk_realign/qc_bam/possorted_chrM_realign.bam \
-f /md01/nieyg/ref/hard-mask/mm10_hard_masked/fasta/genome.fa \
-c /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/mgatk/barcode.tsv \
-m /md01/nieyg/ref/hard-mask/genome_modify/chrM.len \
--varscan_path /md01/nieyg/software/MitoSort/VarScan.v2.3.7.jar \
-o /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/mgatk_realign/qc_bam

