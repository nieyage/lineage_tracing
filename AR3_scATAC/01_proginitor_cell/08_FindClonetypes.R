library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
library(future)
plan("multicore", workers = 12)
plan()
options(future.globals.maxSize = 70 * 1024 ^ 3) # for 50 Gb RAM

combined<- readRDS("./03_all_celltype/04_mgatk/all_cell_type_mgatk_alleles.rds")
DefaultAssay(combined) <- "alleles"

AR3_C4_mito.data <- ReadMGATK(dir = "~/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC/mgatk_for_filtered_cells/final/")
crc <- combined
variable.sites <- IdentifyVariants(crc, assay = "mito", refallele = AR3_C4_mito.data$refallele)
write.csv(variable.sites,"./03_all_celltype/variable_sites_info.csv")
variable.sites<- read.csv("./03_all_celltype/variable_sites_info.csv")
summary(variable.sites$n_cells_conf_detected)
pdf("./04_proginitor_cells/alleles_n_cells_conf_detected_density.pdf",width=10,height=6)
plot(density(variable.sites$n_cells_conf_detected),xlim=c(0,100))
dev.off()


high.conf <- subset(
    variable.sites,
    subset = n_cells_conf_detected >= 50 &
    strand_correlation >= 0.65 &
    vmr > 0.01
)
  crc <- AlleleFreq(
    object = combined,
    variants = high.conf$variant,
    assay = "mito"
  )
  DefaultAssay(crc) <- "alleles"
  crc<- FindVariableFeatures(crc)

crc <- FindClonotypes(crc,features = VariableFeatures(crc)[1:200],resolution = 1,assay="alleles",metric = "cosine",algorithm = 3)

markers <- FindAllMarkers(crc, only.pos = TRUE, min.pct = 0, logfc.threshold = 0.1)

library(dittoSeq)
crc$Clonotypes<- Idents(crc)

pdf("./04_proginitor_cells/3_top200_heatmap_addcelltype.pdf",width=20,height=10)
DoHeatmap(crc, features = VariableFeatures(crc)[1:200], slot = "data", disp.max = 0.1) +scale_fill_viridis_c()
dittoHeatmap(crc, VariableFeatures(crc)[1:200], annot.by = c("Clonotypes", "Annotation"))
dev.off()

# for all 
crc <- FindClonotypes(combined,resolution = 1,assay="alleles",metric = "cosine",algorithm = 3)
crc$Clonotypes<- Idents(crc)
pdf("./04_proginitor_cells/1_all_mutation_heatmap.pdf",width=20,height=10)
DoHeatmap(crc, features = VariableFeatures(crc)[1:200], slot = "data", disp.max = 0.1) +scale_fill_viridis_c()
DoHeatmap(crc, features = VariableFeatures(crc)[1:200], slot = "data") +scale_fill_viridis_c()
dev.off()


# at least 125 cells
combined <- FindClonotypes(combined,resolution = 1,assay="alleles",metric = "cosine",k = 20,algorithm = 3)
# top 100 features 

high.conf <- variable.sites[order(variable.sites$n_cells_conf_detected,decreasing=T),]

# top 10 variant plot 
DefaultAssay(combined) <- "alleles"
alleles.view <- rownames(high.conf[1:9,])
top9_vf<-  VariableFeatures(combined)[1:9]

pdf("./04_proginitor_cells/top9_variant_umap_plot.pdf",width=10,height=10)
FeaturePlot(
  object = combined,
  features = alleles.view,
  order = TRUE,
  cols = c("grey", "darkred"),
  ncol = 3
) & NoLegend()
FeaturePlot(
  object = combined,
  features = top9_vf,
  order = TRUE,
  cols = c("grey", "darkred"),
  ncol = 3
) & NoLegend()
dev.off()


# top 100 variant plot 
combined <- FindClonotypes(combined,features = VariableFeatures(combined)[1:1000],resolution = 1,assay="alleles",metric = "cosine",k = 20,algorithm = 3)
table(Idents(combined))


pdf("./04_proginitor_cells/top9_variant_umap_plot.pdf",width=10,height=10)
DoHeatmap(combined, features = VariableFeatures(combined)[1:100], slot = "data", disp.max = 0.1) +
  scale_fill_viridis_c()
dev.off()

saveRDS(combined,"./04_proginitor_cells/all_cell_type_mgatk_alleles_FindClonetype.rds")



# call mutation by mitosort 

data<- read.csv("/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/mgatk_realign2/qc_bam/possorted_chrM_realign.snv",sep="\t")
variant <- paste0(as.character(data$Position), data$Ref, ">", data$VarAllele)
variant<- variant[which(variant%in%rownames(combined))]
combined <- FindClonotypes(combined,features = variant,resolution = 1,assay="alleles",algorithm = 3)
table(Idents(combined))
pdf("./04_proginitor_cells/4_mutation_by_mitosort_heatmap.pdf",width=10,height=10)
DoHeatmap(combined, features = VariableFeatures(combined)[1:100], slot = "data", disp.max = 0.1) +scale_fill_viridis_c()
DoHeatmap(combined, features = variant, slot = "data", disp.max = 0.1) +scale_fill_viridis_c()
DoHeatmap(combined, features = VariableFeatures(combined)[1:100], slot = "data") +scale_fill_viridis_c()
DoHeatmap(combined, features = variant, slot = "data") +scale_fill_viridis_c()
dev.off()









# cluster like ATAC
high.conf <- subset(
    variable.sites,
    subset = n_cells_conf_detected >= 50 &
    strand_correlation >= 0.65 &
    vmr > 0.01
)
  crc <- AlleleFreq(
    object = combined,
    variants = high.conf$variant,
    assay = "mito"
  )
DefaultAssay(crc) <- "alleles"
crc<- FindVariableFeatures(crc)
combined<- crc
DefaultAssay(combined) <- "alleles"
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 90)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, reduction = "lsi",reduction.key = "UMAPalleles_", dims = 1:10)
combined <- FindClonotypes(combined,resolution = 1)


# plot mito Umap 
pdf("./03_all_celltype/04_mgatk/2_alleles-cluster-Umap.pdf",width=18,height=6)
p0<- DimPlot(combined,group.by = 'Annotation', label = TRUE,repel = TRUE) + NoLegend()
p1<- DimPlot(combined,group.by = 'detail_anno', label = TRUE,repel = TRUE) + NoLegend()
p2<- DimPlot(combined, label = TRUE,repel = TRUE) + NoLegend()
p0+p1+p2
dev.off()
pdf("./03_all_celltype/04_mgatk/2_alleles-cluster-VF-heatmap.pdf",width=12,height=6)
DoHeatmap(combined, features = VariableFeatures(combined)[1:100], slot = "data", disp.max = 0.1) +scale_fill_viridis_c()
DoHeatmap(combined, features = VariableFeatures(combined)[1:100], slot = "data") +scale_fill_viridis_c()
DoHeatmap(combined, features = VariableFeatures(combined), slot = "data", disp.max = 0.1) +scale_fill_viridis_c()
dev.off()
# find top3 markers(mutations) 
markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)
table(markers$cluster)
write.csv(markers,"./03_all_celltype/04_mgatk/MajorClonotypes_FindAllMarkers_mutations.csv")
# verification
signif_markers <- markers[markers$p_val_adj<0.05,] 
top5 <- signif_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC);


pdf("./03_all_celltype/04_mgatk/MajorClonotypes_FindAllMarkers_mutations_heatmap_top5_avg.pdf",width=15,height=15)
DoHeatmap(object = combined,features=top5$gene,slot = "data",label=T, size = 2) 
DoHeatmap(object = combined,features=signif_markers$gene,slot = "data",label=T, size = 2) + scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"))+NoLegend()
dev.off();



# for LINEAGE pipeline 
# get the matrix (all sites and all ATCG) # the last two columu is altAllele and refAllele
AR3_C4_mito.data <- ReadMGATK(dir = "~/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC/mgatk_for_filtered_cells/final/")

mutation_matrix<- GetAssayData(combined,assay="mito")

A_fwd <- mutation_matrix[grep("^A.*fwd$",rownames(mutation_matrix),value=TRUE),]
A_rev <- mutation_matrix[grep("^A.*rev$",rownames(mutation_matrix),value=TRUE),]
T_fwd <- mutation_matrix[grep("^T.*fwd$",rownames(mutation_matrix),value=TRUE),]
T_rev <- mutation_matrix[grep("^T.*rev$",rownames(mutation_matrix),value=TRUE),]
C_fwd <- mutation_matrix[grep("^C.*fwd$",rownames(mutation_matrix),value=TRUE),]
C_rev <- mutation_matrix[grep("^C.*rev$",rownames(mutation_matrix),value=TRUE),]
G_fwd <- mutation_matrix[grep("^G.*fwd$",rownames(mutation_matrix),value=TRUE),]
G_rev <- mutation_matrix[grep("^G.*rev$",rownames(mutation_matrix),value=TRUE),]

rownames(A_fwd)<- gsub("-fwd","",rownames(A_fwd))
rownames(T_fwd)<- gsub("-fwd","",rownames(T_fwd))
rownames(C_fwd)<- gsub("-fwd","",rownames(C_fwd))
rownames(G_fwd)<- gsub("-fwd","",rownames(G_fwd))
rownames(A_rev)<- gsub("-rev","",rownames(A_rev))
rownames(T_rev)<- gsub("-rev","",rownames(T_rev))
rownames(C_rev)<- gsub("-rev","",rownames(C_rev))
rownames(G_rev)<- gsub("-rev","",rownames(G_rev))

A_mut <- A_fwd+A_rev
T_mut <- T_fwd+T_rev
C_mut <- C_fwd+C_rev
G_mut <- G_fwd+G_rev

A_mut <- as.data.frame(A_mut)
T_mut <- as.data.frame(T_mut)
C_mut <- as.data.frame(C_mut)
G_mut <- as.data.frame(G_mut)
A_mut$altAllele<- "A"
T_mut$altAllele<- "T"
C_mut$altAllele<- "C"
G_mut$altAllele<- "G"
A_mut$refAllele<- AR3_C4_mito.data$refallele$ref
T_mut$refAllele<- AR3_C4_mito.data$refallele$ref
C_mut$refAllele<- AR3_C4_mito.data$refallele$ref
G_mut$refAllele<- AR3_C4_mito.data$refallele$ref

all_mut<- rbind(A_mut,T_mut,C_mut,G_mut)

label<- data.frame(barcode=rownames(combined),Label=combined$Annotation)

result=lineage(data = all_mut, repeats = 30, thread = 20)
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













