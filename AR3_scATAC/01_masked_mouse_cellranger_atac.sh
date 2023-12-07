#PBS -N masked_mouse
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=32
#PBS -l mem=30G

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/scATAC/

cellranger-atac count --id=AR3_C4_scATAC \
                        --reference=/md01/nieyg/ref/hard-mask/mm10_hard_masked \
                        --fastqs=/md01/nieyg/project/lineage_tracing/scATAC/rawdata/AR3_C4 \
                        --sample=AR3_C4 \
                        --localcores=32 \
                        --localmem=64

#PBS -N masked_mouse
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=32
#PBS -l mem=30G
source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/scATAC/
cellranger-atac count --id=AR3_C5_scATAC \
                        --reference=/md01/nieyg/ref/hard-mask/mm10_hard_masked \
                        --fastqs=/md01/nieyg/project/lineage_tracing/scATAC/rawdata/AR3_C5 \
                        --sample=AR3_C5 \
                        --localcores=32 \
                        --localmem=64


# mgatk 
## PBS configure 
#PBS -N mgatk_C4
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=16
#PBS -l mem=32G

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/scATAC/AR3_C4_scATAC/mgatk
source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate r4-base
source /md01/nieyg/software/venv3/bin/activate
mgatk tenx -i /md01/nieyg/project/lineage_tracing/scATAC/AR3_C4_scATAC/outs/possorted_bam.bam \
-o /md01/nieyg/project/lineage_tracing/scATAC/AR3_C4_scATAC/mgatk \
-n AR3_C4_scATAC -g /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask_chrM.fasta \
-c 16 -bt CB -b /md01/nieyg/project/lineage_tracing/scATAC/AR3_C4_scATAC/outs/filtered_peak_bc_matrix/barcodes.tsv 

## PBS configure 
#PBS -N mgatk_C5
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=16
#PBS -l mem=32G

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/scATAC/AR3_C5_scATAC/mgatk
source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate r4-base
source /md01/nieyg/software/venv3/bin/activate
mgatk tenx -i /md01/nieyg/project/lineage_tracing/scATAC/AR3_C5_scATAC/outs/possorted_bam.bam \
-o /md01/nieyg/project/lineage_tracing/scATAC/AR3_C5_scATAC/mgatk \
-n AR3_C5_scATAC -g /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask_chrM.fasta \
-c 16 -bt CB -b /md01/nieyg/project/lineage_tracing/scATAC/AR3_C5_scATAC/outs/filtered_peak_bc_matrix/barcodes.tsv 


# 2023-06-28
# cell ranger
# AR3_C4 add to 500G
#PBS -N masked_mouse
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=32
#PBS -l mem=30G
source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/scATAC/
cellranger-atac count --id=AR3_C5_scATAC_add500G \
                        --reference=/md01/nieyg/ref/hard-mask/mm10_hard_masked \
                        --fastqs=/data/R03/zhangwx/rawdata/mouse/Yage-mouse/AR3/AR3_C4/ \
                        --sample=AR3_C4 \
                        --localcores=32 \
                        --localmem=64

# ref:https://github.com/caleblareau/mgatk/wiki/Circular-mtDNA
# Handling the circular mtDNA
# 哺乳动物mtDNA是圆形的，因此许多来自ATAC-seq和其他检测的读数会在连接点附近出现对齐问题（在d环中；人类16569<->1的参考位置）。
/md01/nieyg/software/gmap-2019-09-12/util/gmap_build -d mm10 /md01/nieyg/ref/hard-mask/mm10_hard_masked/fasta/genome.fa -c chrM -D .

gsnap --nthreads 4 --gunzip -D $GSNAP_DIR -d hg38 $fastq1 $fastq2 -A sam | \
     samtools view -bS - | samtools sort -@ 4 - -o "${sample}.st.bam"
     

less singlecell.csv |cut -d "," -f 1 > barcode.tsv


sed -i 's/"//g' barcode.tsv
## PBS configure 
#PBS -N mgatk_C4 for all barcode 
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=32
#PBS -l mem=32G

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/mgatk

source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate r4-base
source /md01/nieyg/software/venv3/bin/activate
mgatk tenx -i /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/outs/possorted_bam.bam \
-o /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/mgatk \
-n AR3_C4_scATAC -g /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask_chrM.fasta \
-c 32 -bt CB -b /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/mgatk/barcode.tsv


## PBS configure 
#PBS -N mgatk_C5
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=10
#PBS -l mem=32G

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/mgatk

source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate r4-base
source /md01/nieyg/software/venv3/bin/activate
mgatk tenx -i /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/outs/possorted_bam.bam \
-o /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/mgatk \
-n AR3_C5_scATAC -g /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask_chrM.fasta \
-c 10 -bt CB -b /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/mgatk/barcode.tsv




# mgatk realign done!

# Step1: merge the all qc_bam files and sort the last bam file 
nohup samtools merge merged.qc.bam barcodes*qc.bam &
nohup samtools sort -@8 merged.qc.bam -o merged.qc.sorted.bam & 









