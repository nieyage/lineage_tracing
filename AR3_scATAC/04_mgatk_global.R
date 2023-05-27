library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
combined<- readRDS("./03_all_celltype/AR3_integrated_all_celltype_annotated.rds")
AR3_C4_last<- readRDS("./03_all_celltype/AR3_C4_scATAC.rds")
AR3_C5_last<- readRDS("./03_all_celltype/AR3_C5_scATAC.rds")

write.table(colnames(AR3_C4_last),"~/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC/mgatk_for_filtered_cells/barcode.tsv",row.names=F,col.names=F)
write.table(colnames(AR3_C5_last),"~/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC/mgatk_for_filtered_cells/barcode.tsv",row.names=F,col.names=F)

# load mgatk output
AR3_C4_mito.data <- ReadMGATK(dir = "~/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC/mgatk_for_filtered_cells/final/")
AR3_C5_mito.data <- ReadMGATK(dir = "~/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC/mgatk_for_filtered_cells/final/")

# create an assay
AR3_C4_mito_counts<- AR3_C4_mito.data$counts
colnames(AR3_C4_mito_counts)<- paste("AR3_C4_",colnames(AR3_C4_mito_counts),sep="")
AR3_C5_mito_counts<- AR3_C5_mito.data$counts
colnames(AR3_C5_mito_counts)<- paste("AR3_C5_",colnames(AR3_C5_mito_counts),sep="")
all_counts<- cbind(AR3_C4_mito_counts,AR3_C5_mito_counts)

AR3_mito <- CreateAssayObject(counts = all_counts)
# Subset to cell present in the scATAC-seq assat
AR3_mito <- subset(AR3_mito, cells = colnames(combined))

# add assay and metadata to the seurat object
combined[["mito"]] <- AR3_mito

rownames(AR3_C4_mito.data$depth) <- paste("AR3_C4_",rownames(AR3_C4_mito.data$depth),sep="")
rownames(AR3_C5_mito.data$depth) <- paste("AR3_C5_",rownames(AR3_C5_mito.data$depth),sep="")
depth_all<- rbind(AR3_C4_mito.data$depth,AR3_C5_mito.data$depth)

combined <- AddMetaData(combined, metadata = depth_all, col.name = "mtDNA_depth")
pdf("./03_all_celltype/all_cell_type_mtDNA_depth.pdf",width=16,height=4)
VlnPlot(crc, "mtDNA_depth", pt.size = 0.1) + scale_y_log10()
dev.off()

# filter cells based on mitochondrial depth

crc <- subset(combined, mtDNA_depth >= 10)
crc
variable.sites <- IdentifyVariants(crc, assay = "mito", refallele = AR3_C4_mito.data$refallele)
pdf("./03_all_celltype/04_mgatk/VariantPlot_global.pdf")
VariantPlot(variants = variable.sites)
dev.off()

# Establish a filtered data frame of variants based on this processing
high.conf <- subset(
	variable.sites, 
	subset = n_cells_conf_detected >= 5 &
    strand_correlation >= 0.65 &
    vmr > 0.01
)

high.conf[,c(1,2,5)]
crc <- AlleleFreq(
  object = crc,
  variants = high.conf$variant,
  assay = "mito"
)
crc[["alleles"]]
DefaultAssay(crc) <- "alleles"

# all celltype mutation 


# major celltype specific mutation 
DoHeatmap(crc, features = rownames(crc), slot = "data", disp.max = 1) + scale_fill_viridis_c()


# detail celltype specific mutation 


pdf("./03_all_celltype/04_mgatk/Variant_Featureplot.pdf")
alleles.view <- c("12889G>A", "16147C>T", "9728C>T", "9804G>A")
FeaturePlot(
  object = crc,
  features = alleles.view,
  order = TRUE,
  cols = c("grey", "darkred"),
  ncol = 4
) & NoLegend()
dev.off()

#

DefaultAssay(tf1) <- "alleles"
tf1 <- FindClonotypes(tf1)

saveRDS(crc,"./03_all_celltype/04_mgatk/all_cell_type_mgatk.rds")

#
#crc$log10_mtDNA_depth <- log10(crc$mtDNA_depth)
#FeaturePlot(crc, features = c("pct_reads_in_DNase"))
#FeaturePlot(crc, features = c("log10_mtDNA_depth"))
#source("variant_calling.R") 
#
## mgatk_se is the Summarized Experiment .rds file
## That is automatically produced from running
## The mgatk CLI python package 
#
#mut_se <- call_mutations_mgatk(mgatk_se)



