library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
combined<- readRDS("./03_all_celltype/03_recall_peak/AR3_integrated_all_celltype_annotated_recall_peak.rds")
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
Idents(combined)<- combined$detail_anno
pdf("./03_all_celltype/04_mgatk/all_cell_type_mtDNA_depth.pdf",width=16,height=8)
VlnPlot(combined, c("mtDNA_depth","nCount_mito","nFeature_mito"), pt.size = 0,ncol=1) + scale_y_log10()
dev.off()
saveRDS(crc,"./03_all_celltype/04_mgatk/all_cell_type_mgatk.rds")
