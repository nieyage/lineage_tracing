library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
combined<- readRDS("./03_all_celltype/AR3_integrated_all_celltype_annotated.rds")
# detailed QC to make sure the celltype is truly exist 

# feature plot 
pdf("./03_all_celltype/all_cell_type_QC.pdf",width=16,height=16)
FeaturePlot(
  object = combined,features = c("nCount_ATAC","nFeature_ATAC","total","duplicate","unmapped","lowmapq",
  	"mitochondrial","passed_filters","is__cell_barcode","TSS_fragments","peak_region_fragments",
  	"nucleosome_signal","TSS.enrichment","pct_reads_in_peaks","nCount_RNA","nFeature_RNA"),
  order = FALSE,
  cols = c("grey", "darkred"),
  ncol = 4
) & NoLegend()
dev.off()

# vlnplot 
Idents(combined)<- combined$detail_anno
pdf("./03_all_celltype/all_cell_type_QC_vlnplot.pdf",width=16,height=10)
VlnPlot(object = combined,features = c("nCount_ATAC","nFeature_ATAC","total","duplicate"), ncol = 1, pt.size = 0)
VlnPlot(object = combined,features = c("unmapped","lowmapq","mitochondrial","passed_filters"), ncol = 1, pt.size = 0)
VlnPlot(object = combined,features = c("is__cell_barcode","TSS_fragments","peak_region_fragments","nucleosome_signal"), ncol = 1, pt.size = 0)
VlnPlot(object = combined,features = c("TSS.enrichment","pct_reads_in_peaks","nCount_RNA","nFeature_RNA"), ncol = 1, pt.size = 0)

dev.off()

