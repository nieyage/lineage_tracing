
conda activate r4-base
source /md01/nieyg/software/venv3/bin/activate
samtools faidx /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa chrM > assembly_hardmask_chrM.fasta


## PBS configure 
#PBS -N mutaion
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=32
#PBS -l mem=32G

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/mgatk/test
source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate r4-base
source /md01/nieyg/software/venv3/bin/activate
mgatk tenx -i /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/outs/atac_possorted_bam.bam \
-o /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/mgatk/ \
-n plogF1met5d_masked -g /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask_chrM.fasta \
-c 32 -bt CB -b /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/mgatk/barcodes.tsv 

#######################
######## test #########
#######################

nohup mgatk tenx -i /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/outs/atac_possorted_bam.bam \
-o /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/mgatk/test \
-n plogF1met5d_masked -g /md01/nieyg/ref/hard-mask/mm10_hard_masked/fasta/genome.fa \
-c 32 -bt CB -b /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/mgatk/test.tsv &

nohup mgatk check -i /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/outs/atac_possorted_bam.bam \
-o /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/mgatk/test4 \
-n plogF1met5d_masked -g /md01/nieyg/ref/hard-mask/mm10_hard_masked/fasta/genome.fa \
-c 32 -bt CB -b /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/mgatk/test.tsv & 





/md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/outs/atac_possorted_bam.bam --cell-barcodes test.tsv --cores 10 --out-bam test.bam
samtools sort -@ 12 -O bam -o test_sorted.bam test.bam
samtools index test_sorted.bam 

mgatk check -i test_sorted.bam \
-o /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/mgatk/test4 \
-n plogF1met5d_masked -g /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask_chrM.fasta -c 32 -bt CB \
-b /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/mgatk/test.tsv


Given the .fa or .fasta file that was part of your reference genome for whichever alignment tool, 
one can extract only the mitochondrial genome using the following command for chrM 
(update your mitochondrial genome name accordingly-- other examples may be MT).
 For example:


mgatk tenx -i test_barcode.bam -n bc1 -o test -bt CB -b test_barcodes.txt -c 32





## PBS configure 
#PBS -N mitosort_OSN
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=16
#PBS -l mem=32G

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked_add/outs/mitosort

## realign MT reads 
python /md01/nieyg/software/MitoSort/MitoSort_pipeline.py mt-realign \
-b /md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked_add/outs/atac_possorted_bam.bam \
-f /md01/nieyg/ref/hard-mask/mm10_hard_masked/fasta/genome.fa \
--gatk_path /md01/nieyg/software/MitoSort/GenomeAnalysisTK_3.5-0.jar \
-o /md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked_add/outs/mitosort

## generate SNP matrices
python /md01/nieyg/software/MitoSort/MitoSort_pipeline.py generate-snp-matrix \
-b /md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked_add/outs/mitosort/MitoSort/BAM/possorted_chrM_realign.bam \
-f /md01/nieyg/ref/hard-mask/mm10_hard_masked/fasta/genome.fa \
-c /md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked_add/outs/per_barcode_metrics.csv \
-m /md01/nieyg/ref/hard-mask/genome_modify/chrM.len \
--varscan_path /md01/nieyg/software/MitoSort/VarScan.v2.3.7.jar \
-o /md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked_add/outs/mitosort


write.table(colnames(OSN_last),"/md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked_add/outs/mgatk/barcode.tsv",row.names=F,col.names=F)

sed -i 's/"//g' barcode.tsv

## PBS configure 
#PBS -N mgatk_OSN
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=32
#PBS -l mem=32G

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked_add/outs/mgatk

source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate r4-base
source /md01/nieyg/software/venv3/bin/activate
mgatk tenx -i /md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked_add/outs/atac_possorted_bam.bam \
-o /md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked_add/outs/mgatk \
-n OSN_regeneration -g /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask_chrM.fasta \
-c 32 -bt CB -b /md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked_add/outs/mgatk/barcode.tsv











