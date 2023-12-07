# 4.0 preprocessing: generate mitochondrial genotype matrixes

# 4.0.1 After alignment, new bam files consisting of MtDNA records, which were extracted from the alignment result with SAMtools, were obtained.
#samtools view -bS {STAR_OUT_SAM} > {OUT_BAM}
#samtools sort {INPUT_BAM} -o {OUT_SORTED_BAM}
#samtools index {OUT_SORTED_BAM}

nohup samtools view -@ 10 -b /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/outs/possorted_bam.bam  chrM > possorted_bam_chrM.bam &
nohup samtools view -@ 10 -b /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/outs/possorted_bam.bam  chrM > possorted_bam_chrM.bam &
ln -s /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/LINEAGE/possorted_bam_chrM.bam AR3_C5_scATAC_possorted_bam_chrM.bam
ln -s /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/LINEAGE/possorted_bam_chrM.bam AR3_C4_scATAC_possorted_bam_chrM.bam

# ref: https://github.com/songjiajia2018/ppl
# Generate the SNP matrix 
# 1. For 10x data:
samtools view -h -@ 10 -o AR3_C4_chrM.sam AR3_C4_scATAC_possorted_bam_chrM.bam
samtools view -h -@ 10 -o AR3_C5_chrM.sam AR3_C5_scATAC_possorted_bam_chrM.bam

## PBS configure 
#PBS -N LINEAGE
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=5
#PBS -l mem=32G
source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/LINEAGE
source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate r4-base

python /md01/nieyg/software/ppl/ppl2_run.py \
-t 1 --qalign 0 --qbase 0 \
--reference /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask_chrM.fasta \
--split-sam --input-filelist -p -m -r \
--outdir /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/LINEAGE \
--name AR3_ppl10x_hardmask_chrM \
--input 10x_split_chrM_bam_file.list 

Rscript /md01/nieyg/software/ppl/mtMatrix.R processed.MAE_mito.rds
Rscript /md01/nieyg/software/ppl/Script/cells.R processed.MAE_mito.rds annotation 10
Rscript /md01/nieyg/software/ppl/extract_af.R processed.MAE_mito.rds



# split bam file by barcode (+samplename)
nohup samtools view -@ 10 -bS AR3_C4_chrM.sam > AR3_C4_scATAC_possorted_bam_chrM.bam &
nohup samtools view -@ 10 -bS AR3_C5_chrM.sam > AR3_C5_scATAC_possorted_bam_chrM.bam &

scsplit -i AR3_C4_scATAC_possorted_bam_chrM.bam  -o AR3_C4_chrM_bam
scsplit -i AR3_C5_scATAC_possorted_bam_chrM.bam  -o AR3_C5_chrM_bam


# Perform a pileup for all bam samples
# 指定包含拆分后BAM文件的文件夹路径
input_directory="AR3_C4_chrM_bam"  # 替换为实际的文件夹路径
prefix="AR3_C4_"
# 循环遍历文件夹中的所有BAM文件
for bam_file in ${input_directory}/*.bam; do
    # 从BAM文件名中提取barcode名称，假设文件名包含barcode信息
    barcode=$(basename "$bam_file" .bam)  # 提取文件名并去掉扩展名
    barcode=${barcode##*_}  # 使用下划线分隔的barcode情况，提取最后的部分
    new_barcode="${prefix}${barcode}"
    # 在这里运行您的命令，使用"$bam_file"引用当前BAM文件，"$barcode"引用当前barcode名称
    python /md01/nieyg/software/ppl/Script/01_pileup_counts.py $bam_file processed_data/$new_barcode 16571 0 $new_barcode 0
done

sh 02_merge_pileup_counts.sh processed_data BM0820_example
Rscript 03_makeRDS.R processed_data mito_fastas/hg19.fasta

# ref :https://github.com/jon-xu/scSplit

git clone https://github.com/jon-xu/scSplit

freebayes -f /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask_chrM.fasta -iXu -C 2  AR3_C4_scATAC_possorted_bam_chrM.bam > AR3_C4_scATAC_possorted_bam_chrM.vcf

scSplit count -v filtered.vcf -i filtered_rmdup_sorted.bam -b barcodes.tsv -r ref_filtered.csv -a alt_filtered.csv -o .













# get the input file form mgatk?

# 4.1 mode1: run parallel iterative optimization
dir=/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/LINEAGE/AR3_ppl10x_hardmask_chrM_notsplitsam

data<- readRDS("/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/LINEAGE/AR3_ppl10x_hardmask_chrM_notsplitsam/processed.MAE_mito.rds")


library(LINEAGE)

data("TF1_clones")
data=TF1_clones$data
label=TF1_clones$rlabel    #reference clone labels
result=lineage(data = data, repeats = 30, thread = 20)

# The inferred clone labels are embedded in the returned result
# as result$label. We also provide lineage_tree function to trace # the lineage tree of the returned result.
hc=lineage_tree(result)
str(hc, max.level = 1)

plots0=traceplot(result, label)  #plots with reference clone labels
plots=traceplot(result, result$label)  #plots with inferred clone labels
# 2d visualization cluster results with clonal labels
# colored in reference clone labels 
print(plots0$d2d)
# or inferred labels
print(plots$d2d)
# Heatmap results with markers information across cells and color bar is
# colored in reference clone labels
print(plots$heatmap)
# The recommended result is embedded in the returned result as result$best.
best=list(result=result$best, plots=plots)





data=TF1_clones$data
result=lineage(data = data, repeats = 30, thread = 20)
#    data: A dataframe containing mitochondrial variants frequency
#          matrix, where a column represents a single cell and a row
#          represents variants frequency of a specific mitochondrial
#          genotype, and two other columns "altAllele" and "refAllele".
#          "altAllele" column represents the mutant allele; "refAllele"
#          column represents the reference allele. Required.
#
# repeats: Number of iterations. Default:30.
#
#  thread: Number of threads to use. NULL for non-parallel iterative
#          optimization. Default:10.


# 4.3 trace the lineage tree

# The inferred clone labels are embedded in the returned result
# as result$label. We also provide lineage_tree function to trace # the lineage tree of the returned result.
hc=lineage_tree(result)
str(hc, max.level = 1)

result <- hc
# 4.4 visualization of the result
plots=traceplot(result, result$labels)  #plots with inferred clone labels

# 4.2 mode2: run non-parallel iterative optimization
data("TF1_clones")
data=TF1_clones$data
result=lineage(data = data, repeats = 30, thread = NULL)


# 2d visualization cluster results with clonal labels
# colored in reference clone labels 
print(plots0$d2d)
# or inferred labels
print(plots$d2d)
# Heatmap results with markers information across cells and color bar is
# colored in reference clone labels
print(plots$heatmap)

# The recommended result is embedded in the returned result as result$best.
best=list(result=result$best, plots=plots)

mutation_matrix<- GetAssayData(combined,assay="alleles")

