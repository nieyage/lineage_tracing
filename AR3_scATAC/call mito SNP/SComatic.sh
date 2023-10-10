conda create -n SComatic -c bioconda python=3.7 r-base=4.2 samtools datamash bedtools
conda activate SComatic
pip install -r requirements.txt
Rscript r_requirements_install.R
gunzip PoNs/PoN.scRNAseq.hg38.tsv.gz
gunzip PoNs/PoN.scATACseq.hg38.tsv.gz 

# Step 0: get the metadata
data<- as.data.frame(combined@meta.data)
data$barcode<- rownames(data)
data$barcode<- gsub("AR3_C4_","",data$barcode)
data$barcode<- gsub("AR3_C5_","",data$barcode)

AR3_C4_data<- data[data$orig.ident=="AR3_C4",]
AR3_C5_data<- data[data$orig.ident=="AR3_C5",]

write.table(AR3_C4_data[,c("barcode","Annotation")],"/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/SComatic/01_SplitBamCellTypes/AR3_C4_metadata.tsv",row.names=F,col.names=F,sep="\t")
write.table(AR3_C5_data[,c("barcode","Annotation")],"/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/SComatic/01_SplitBamCellTypes/AR3_C5_metadata.tsv",row.names=F,col.names=F,sep="\t")

sed -i 's/"//g' AR3_C4_metadata.tsv
sed -i 's/"//g' AR3_C5_metadata.tsv

# ref: https://github.com/cortes-ciriano-lab/SComatic/blob/main/docs/faqs.md#1-are-the-scomatic-parameters-for-scatac-seq-data-the-same-as-for-scrna-seq-data
# 1. Are the SComatic parameters for scATAC-seq data the same as for scRNA-seq data?
# No, they are different. ATAC-seq data is a DNA-based approach. Therefore, there are differences in the way of processing the bam files. These are the parameters that you should change in comparison to the scRNA-seq workflow:
# 
# Step 1: remove the --max_nM and --max_NH parameters, set the mapping quality filter to --min_MQ 30 .
# Step 2: set the mapping quality filter to --min_mq 30
# Step 4.2: Remove the --editing parameter, and --pon altered to point at the scATACseq PoN provided in this GitHub repo (or custom one)


# Step 1: Splitting alignment file into cell-type-specific bams
## PBS configure 
#PBS -N cell-type-specific-bams
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=32G

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/SComatic/01_SplitBamCellTypes
source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate SComatic

SCOMATIC=/md01/nieyg/software/SComatic/
output_dir=/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/SComatic/01_SplitBamCellTypes
mkdir -p $output_dir

python $SCOMATIC/scripts/SplitBam/SplitBamCellTypes.py \
--bam /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/outs/possorted_bam.bam \
--meta AR3_C4_metadata.tsv --id AR3_C4_scATAC  --min_MQ 30 \
--outdir $output_dir

## PBS configure 
#PBS -N cell-type-specific-bams
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=32G

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/SComatic/01_SplitBamCellTypes
source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate SComatic

SCOMATIC=/md01/nieyg/software/SComatic/
output_dir=/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/SComatic/01_SplitBamCellTypes
mkdir -p $output_dir

python $SCOMATIC/scripts/SplitBam/SplitBamCellTypes.py \
--bam /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/outs/possorted_bam.bam \
--meta AR3_C5_metadata.tsv --id AR3_C5_scATAC --max_nM 5 --min_MQ 30 \
--outdir $output_dir

# Step 2: Collecting base count information

## PBS configure 
#PBS -N 02_BaseCellCounts
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=10
#PBS -l mem=32G
source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/SComatic/
source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate SComatic
REF=/md01/nieyg/ref/hard-mask/mm10_hard_masked/fasta/genome.fa
SCOMATIC=/md01/nieyg/software/SComatic/
output_dir=/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/SComatic/
output_dir1=$output_dir/01_SplitBamCellTypes
output_dir2=$output_dir/02_BaseCellCounts
mkdir -p $output_dir2
for bam in $(ls -d $output_dir1/*bam);do
  # Cell type
  cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')
  # Temp folder
  temp=$output_dir2/temp_${cell_type}
  mkdir -p $temp
  # Command line to submit to cluster
  python $SCOMATIC/scripts/BaseCellCounter/BaseCellCounter.py --bam $bam \
    --ref $REF \
    --chrom chrM \
    --out_folder $output_dir2 \
    --min_mq 30 \
    --id AR3_C5_scATAC \
    --tmp_dir $temp \
    --nprocs 10

  rm -rf $temp
done


# Step2: test 
REF=/md01/nieyg/ref/hard-mask/mm10_hard_masked/fasta/genome.fa
SCOMATIC=/md01/nieyg/software/SComatic/
output_dir=/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/SComatic/
output_dir1=$output_dir/01_SplitBamCellTypes
output_dir2=$output_dir/02_BaseCellCounts
  # Cell type
  cell_type="FB"
  # Temp folder
  temp=$output_dir2/temp_test_${cell_type}

  # Command line to submit to cluster
  python $SCOMATIC/scripts/BaseCellCounter/BaseCellCounter.py --bam ./01_SplitBamCellTypes/AR3_C4_scATAC.FB.bam \
    --ref $REF \
    --chrom all --min_dp 1 \
    --out_folder $output_dir2 \
    --min_bq 30 \
    --id AR3_C4_scATAC \
    --tmp_dir $temp \
    --nprocs 20


# Step 3: Merging base count matrices
## PBS configure 
#PBS -N 03_BaseCellCountsMerged
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=10
#PBS -l mem=32G
source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/SComatic/01_SplitBamCellTypes
source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate SComatic
sample=AR3_C5
SCOMATIC=/md01/nieyg/software/SComatic/
output_dir=/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/SComatic/
output_dir1=$output_dir/01_SplitBamCellTypes
output_dir2=$output_dir/02_BaseCellCounts
output_dir3=$output_dir/03_BaseCellCountsMerged
mkdir -p $output_dir3

python $SCOMATIC/scripts/MergeCounts/MergeBaseCellCounts.py --tsv_folder ${output_dir2} \
  --outfile ${output_dir3}/${sample}.BaseCellCounts.AllCellTypes.tsv

# Step 4: Detection of somatic mutations

## PBS configure 
#PBS -N 04_VariantCalling
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=10
#PBS -l mem=32G
source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/SComatic/01_SplitBamCellTypes
source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate SComatic
sample=AR3_C4
SCOMATIC=/md01/nieyg/software/SComatic/
output_dir=/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/SComatic/
output_dir1=$output_dir/01_SplitBamCellTypes
output_dir2=$output_dir/02_BaseCellCounts
output_dir3=$output_dir/03_BaseCellCountsMerged
mkdir -p $output_dir3
# Step 4.1
output_dir4=$output_dir/04_VariantCalling
mkdir -p $output_dir4
REF=/md01/nieyg/ref/hard-mask/mm10_hard_masked/fasta/genome.fa

python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step1.py \
          --infile ${output_dir3}/${sample}.BaseCellCounts.AllCellTypes.tsv \
          --outfile ${output_dir4}/${sample} \
          --ref $REF


REF=/md01/nieyg/ref/hard-mask/mm10_hard_masked/fasta/genome.fa
SCOMATIC=/md01/nieyg/software/SComatic/
output_dir=/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/SComatic/
output_dir1=$output_dir/01_SplitBamCellTypes
output_dir2=$output_dir/02_BaseCellCounts

  # Command line to submit to cluster
  python $SCOMATIC/scripts/BaseCellCounter/BaseCellCounter.py --bam $output_dir1/AR3_C4_scATAC.B.bam \
    --ref $REF \
    --chrom chrM \
    --out_folder $output_dir2 \
    --min_bq 30 \
    --id AR3_C4_scATAC \
    --nprocs 10





