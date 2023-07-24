library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)

combined<- readRDS("./03_all_celltype/04_mgatk/all_cell_type_mgatk_alleles.rds")
Idents(combined)<- combined$detail_anno
DefaultAssay(combined) <- "alleles"
combined <- FindVariableFeatures(object = combined)

combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- RunUMAP(object = combined, reduction.key = "UMAP_alleles_",dims = 1:30)
combined <- FindClonotypes(combined,resolution = 3)

saveRDS(combined,"./03_all_celltype/04_mgatk/all_cell_type_mgatk_alleles_FindClonetype.rds")

# plot mito Umap 
pdf("./03_all_celltype/04_mgatk/alleles-cluster-Umap.pdf",width=18,height=6)
p0<- DimPlot(combined,group.by = 'Annotation', label = TRUE,repel = TRUE) + NoLegend()
p1<- DimPlot(combined,group.by = 'detail_anno', label = TRUE,repel = TRUE) + NoLegend()
p2<- DimPlot(combined, label = TRUE,repel = TRUE) + NoLegend()
p0+p1+p2
dev.off()
pdf("./03_all_celltype/04_mgatk/alleles-cluster-VF-heatmap.pdf",width=12,height=6)
DoHeatmap(combined, features = VariableFeatures(combined), slot = "data", disp.max = 0.1) +
  scale_fill_viridis_c()
dev.off()
# find top3 markers(mutations) 
markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)
table(markers$cluster)
write.csv(markers,"./03_all_celltype/04_mgatk/MajorClonotypes_FindAllMarkers_mutations.csv")
# verification
signif_markers <- markers[markers$p_val_adj<0.05,] 
top5 <- signif_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC);

pdf("./03_all_celltype/04_mgatk/MajorClonotypes_FindAllMarkers_mutations_heatmap_top5_avg.pdf",width=15,height=15)
DoHeatmap(object = combined,features=top5$gene,label=T, group.colors =myUmapcolors,
  disp.min = -1,disp.max = 1,size = 2,group.by = "Annotation") + scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"))
DoHeatmap(object = combined,features=markers$gene,label=T, group.colors =myUmapcolors,
  size = 2,group.by = "Annotation") + scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"))+NoLegend()
dev.off();
