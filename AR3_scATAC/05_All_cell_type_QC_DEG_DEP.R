library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
combined<- readRDS("./03_all_celltype/03_recall_peak/AR3_integrated_all_celltype_annotated_recall_peak.rds")
# detailed QC to make sure the celltype is truly exist 
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
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
VlnPlot(object = combined,features = c("is_cell_barcode","TSS_fragments","peak_region_fragments","nucleosome_signal"), ncol = 1, pt.size = 0)
VlnPlot(object = combined,features = c("TSS.enrichment","pct_reads_in_peaks","nCount_RNA","nFeature_RNA"), ncol = 1, pt.size = 0)
dev.off()

# Major celltype 
# Find DEG between clusters(by Gene Score)
DefaultAssay(combined)<-"RNA"
combined <- ScaleData(combined,features=rownames(combined))
markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)
table(markers$cluster)
write.csv(markers,"./03_all_celltype/05_DEG_and_DEP/MajorCelltype_FindAllMarkers_gene.csv")
# verification
signif_markers <- markers[markers$p_val_adj<0.05,] 
top5 <- signif_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC);
top5_Avg <- AverageExpression(combined,features=top5$gene,assays = "RNA")
library(pheatmap)
count=t(scale(t(top5_Avg$RNA),scale = T,center = T))
pdf("./03_all_celltype/05_DEG_and_DEP/MajorCelltype_FindAllMarkers_heatmap_top5_avg.pdf",width=8,height=10)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "firebrick3"))(100),show_rownames=T,show_colnames=T)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "#E41A1C"))(100),show_rownames=T,show_colnames=T)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),show_rownames=T,show_colnames=T)
dev.off();
pdf("./03_all_celltype/05_DEG_and_DEP/MajorCelltype_FindAllMarkers_heatmap_top5_avg.pdf",width=15,height=15)
DoHeatmap(object = combined,features=top5$gene,label=T, group.colors =myUmapcolors,
  disp.min = -1,disp.max = 1,size = 2,group.by = "Annotation") + scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"))
DoHeatmap(object = combined,features=markers$gene,label=T, group.colors =myUmapcolors,
  size = 2,group.by = "Annotation") + scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"))+NoLegend()
dev.off();

# GO and KEGG for Major cell type DEG 
library(org.Mm.eg.db)
library(biomaRt)
library(dplyr)
library(goseq)
library(DOSE)
library(stringr)
library(clusterProfiler)
library(GenomicRanges)
library(AnnotationDbi)

# GO:BP for marker genes in Major celltype 
all_ego<- data.frame()
pdf("./03_all_celltype/05_DEG_and_DEP/Majorcelltype_GO_BP.pdf")
for (i in levels(combined)){
  gene<- markers[markers$cluster==i,]$gene;
  id<- mapIds(x = org.Mm.eg.db,
               keys = gene,
               keytype = "SYMBOL",
               column = "ENTREZID");
  ego<-enrichGO(
    gene= id,
    OrgDb= org.Mm.eg.db,
    ont  = "BP",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE,
    pool = FALSE)
  data<-as.data.frame(ego)
  if(!nrow(data)==0){
      p<-barplot(ego, main=i,showCategory=20);
      print(p);
    write.csv(ego,paste0(i,"_GO_BP.csv"))}
  ego<-ego[ego$pvalue<0.05,]
  ego$celltype=i
  all_ego<-rbind(all_ego,ego)
}

dev.off()

write.csv(all_ego,"./03_all_celltype/05_DEG_and_DEP/Majorcelltype_GO_BP.csv")

top5_term <- all_ego %>% group_by(celltype) %>% top_n(n = 5, wt = -log10(p.adjust));
library(ggplot2)
pdf("./03_all_celltype/05_DEG_and_DEP/Majorcelltype_GO_BP.pdf",width=14,height=20)
p <- ggplot(top5_term,aes(y=Count,x=Description,fill=-log10(p.adjust))) + 
      geom_bar(stat="identity",position = "dodge") +
      facet_grid(celltype~.,scales = "free",space = "free") + 
      coord_flip() + 
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            strip.text.y = element_text(size = 14),
            legend.position="right",
            legend.title = element_text(size=18),
            legend.text = element_text(size=14),
            axis.text.x = element_text(size=14),
            axis.text.y = element_text(size=18),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
p
dev.off()

# Find differentially accessible peaks between clusters
DefaultAssay(combined)<- "peaks";
da_peaks <- FindAllMarkers(  object = combined, test.use = 'LR',logfc.threshold = 0.1, min.pct = 0.05, latent.vars = 'peak_region_fragments')
head(da_peaks)
write.csv(da_peaks,"./03_all_celltype/05_DEG_and_DEP/MajorCelltype_FindAllMarkers_peaks.csv")

peak_Avg <- AverageExpression(combined,features=da_peaks$gene,assays = "peaks")
count=t(scale(t(peak_Avg$peaks),scale = T,center = F))
pdf("./03_all_celltype/05_DEG_and_DEP/DEP_MajorCelltype_FindAllMarkers_heatmap_avg.pdf",width=15,height=15)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "firebrick3"))(100),show_rownames=F,show_colnames=T)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "#E41A1C"))(100),show_rownames=F,show_colnames=T)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),show_rownames=F,show_colnames=T)
dev.off();
pdf("./03_all_celltype/05_DEG_and_DEP/DEP_MajorCelltype_FindAllMarkers_heatmap.pdf",width=15,height=15)
DoHeatmap(object = combined,features=da_peaks$gene,label=T, group.colors =myUmapcolors,size = 2,group.by = "Annotation") + scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"))+NoLegend()
dev.off();

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
# add motif information
mouse <- AddMotifs(
  object = combined,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)
combined <- RunChromVAR(
  object = combined,
  genome = BSgenome.Mmusculus.UCSC.mm10
)
DefaultAssay(combined) <- 'chromvar'
differential.activity <- FindAllMarkers(object = combined, only.pos = TRUE,  mean.fxn = rowMeans,  fc.name = "avg_diff")
write.csv(differential.activity,"./03_all_celltype/05_DEG_and_DEP/MajorCelltype_FindAllMarkers_TF_chromvar.csv")


#closest_genes_da_peaks <- ClosestFeature(pbmc, regions = rownames(da_peaks))
#
## test enrichment
#enriched.motifs <- FindMotifs(
#  object = mouse_brain,
#  features = top.da.peak
#)
#MotifPlot(
#  object = mouse_brain,
#  motifs = head(rownames(enriched.motifs))
#)