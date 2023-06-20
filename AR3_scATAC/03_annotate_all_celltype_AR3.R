library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)

#################################
# annotate all cell type in AR3 #
#################################

combined<- readRDS("./03_all_celltype/combined_AR3_scATAC.rds")
combined$orig.ident<- gsub("_................-1","",colnames(combined))

DefaultAssay(combined) <- 'ATAC'
pdf("./03_all_celltype/DepthCor.pdf")
DepthCor(combined)
dev.off()

combined <- RunUMAP(object = combined, reduction = 'lsi', dims = 2:30)
combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:30)
combined <- FindClusters(object = combined,resolution = 2, verbose = FALSE, algorithm = 3)

pdf("./03_all_celltype/Unannotated_allcelltype_UMAP.pdf",width=6,height=5)
DimPlot(combined, label = T, repel = TRUE,  reduction = "umap",group.by = "seurat_clusters")
DimPlot(combined, label = F, repel = TRUE,  reduction = "umap",group.by = "orig.ident")
dev.off()
# compute gene activities
gene.activities <- GeneActivity(combined)

# add the gene activity matrix to the Seurat object as a new assay
combined[['RNA']] <- CreateAssayObject(counts = gene.activities)
combined <- NormalizeData(
  object = combined,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined$nCount_RNA)
)
library(scCustomize)

markerGenes  <- c(
"Tnnt2","Myl2","Actc1","Nppa","Myh6","Atp2a2","Acta1", # cardiomyocyte
"Ccl3","Clec4a1","Rgs1","Fcgr1","Cd14","Csf1r","Cd163","Cd68","Itgam","Lgals3","Mrc1", # Macrophages
"Cd3e","Cd3d","Cd8a","Cd8b1","Igfbp4","Lat","Itk","Cd3g", # T cell
"Cd22","Cd79a","Cd79b","Mzb1","Ly6d","Ms4a1", # B cell
"Cd74","Cd83","Cd86","Flt3","Cd209a","Ccl17", # DC
"Plac8", # monocyte
"S100a9", "S100a8",#Granulocy
"Nkg7","Gzma", #NK
"Chl1","Kcna6", # Glial
"Myct1","Cdh5","Ednrb","Aqp7","Emcn","Madcam1","Epas1","Flt1","Tie1","Fabp4","Esam","Kdr","Tek", # endothelial
"Dpt","Col3a1","Mmp2","Col1a2","Fstl1","Gsn","Thy1","Pdgfra", # fibroblast
"Abcc9","Rgs4","Ano1","Acta2","Myh11", # smoothmuscle
"Msln","Krt19","Plet1","Prr15","Slc9a3r1","Lrrn4","Slc39a8","Krt8","Anxa8",# epicardial
"Gfpt2","Hk2","Ddr2","Slc1a7","Adamtsl3","Layn","Cd248","Mcam","Fgfr1" # Pericyte
)

label<- c(rep("CM",7),rep("MP",11),rep("T",8),rep("B",6),rep("DC",6),rep("monocyte",1),rep("Granulocy",2),rep("NK",2),rep("Glial",2),rep("EC",13),
	rep("FB",8),rep("SMC",5),rep("Epi",9),rep("Pericyte",9))
DefaultAssay(combined)<-"RNA"
pdf("./03_all_celltype/AR3_cluster-annotation-all_celltype.pdf",width=20,height=8)
p<-DotPlot(combined, features = markerGenes,dot.scale = 3)+theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p
p&scale_x_discrete(labels=label)
dev.off()
pdf("./03_all_celltype/AR3_cluster-annotation-all_celltype_cluster.pdf",width=20,height=8)
Clustered_DotPlot(seurat_object = combined, flip=T,features = markerGenes,k=11)
dev.off()

# make the tree 

# make the trans dist tree 
object <- combined
embeddings <- Embeddings(object = object, reduction = "lsi")[,1:50]
data.dims <- lapply(X = levels(x = object), FUN = function(x) {
    cells <- WhichCells(object = object, idents = x)
    if (length(x = cells) == 1) {
        cells <- c(cells, cells)
    }
    temp <- colMeans(x = embeddings[cells, ])
})
data.dims <- do.call(what = "cbind", args = data.dims)
colnames(x = data.dims) <- levels(x = object)
library(lsa)
cosine_dist <- as.dist(1-cosine(data.dims))
data.tree <- ape::as.phylo(x = hclust(d = cosine_dist))
library(ggtree);

pdf("./03_all_celltype/AR3-combined-tree-cosine.pdf",width=6,height=6)
ggtree(data.tree,layout = "circular") + geom_tiplab()+ geom_treescale()
dev.off()

library(Seurat)
# Load the pre-processed scRNA-seq data for PBMCs
scRNA <- readRDS("/md01/Chenzh275/ori/Data/scATAC/cardiacmyocyte_regeneration/integrative_analysis/cmc_sct.rds")
DefaultAssay(scRNA) <- "RNA" # Create dummy new assay to demo switching default assays
scRNA<- NormalizeData(scRNA)
scRNA<- ScaleData(scRNA)
# first annotation 
pdf("./03_all_celltype/integrated_scRNA-cluster-annotation-all_celltype.pdf",width=20,height=8)
p<-DotPlot(scRNA, features = markerGenes,dot.scale = 3)+theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p
p&scale_x_discrete(labels=label)
dev.off()

CR_features<- c("Myh6","Cdh5","Pdgfra","Msln","Myh11","Kcnj8",
	"Adgre1","Itgal","Naaa","S100a9","Cd3g","Ms4a1","Plp1")
pdf("./03_all_celltype/integrated_scRNA-cluster-annotation-CR.pdf",width=20,height=8)
p<-DotPlot(scRNA, features = CR_features,dot.scale = 3)+theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p
dev.off()

pdf("./03_all_celltype/integrated_scRNA-cluster-Umap.pdf",width=12,height=6)
p1<- DimPlot(scRNA,group.by = 'cell_types', label = TRUE,repel = TRUE) + NoLegend()
p2<- DimPlot(scRNA, label = TRUE,repel = TRUE) + NoLegend()
p1+p2
dev.off()


# modify the annotation in scRNA 
# 11-- pericyte
# 19-- Epi
# 24-- Glial 
cell_types<- scRNA$cell_types
seurat_clusters<- Idents(scRNA)
data<- data.frame(cell_types,seurat_clusters)
for (i in 1:nrow(data)){
	if(data[i,2]=="11"){data[i,1]="Pericyte"}
	if(data[i,2]=="24"){data[i,1]="Glial"}
	if(data[i,2]=="19"){data[i,1]="Epicardial"}
}
scRNA$re_annotation<- data$cell_types
Idents(scRNA)<- scRNA$re_annotation
DefaultAssay(scRNA) <- "RNA" # Create dummy new assay to demo switching default assays
pdf("./03_all_celltype/integrated_scRNA-reannotation-Umap.pdf",width=12,height=6)
p1<- DimPlot(scRNA,group.by = 're_annotation', label = TRUE,repel = TRUE) + NoLegend()
p2<- DimPlot(scRNA, label = TRUE,repel = TRUE) + NoLegend()
p1+p2
dev.off()

scRNA <- FindVariableFeatures(
  object = scRNA,
  nfeatures = 5000
)

transfer.anchors <- FindTransferAnchors(
  reference = scRNA,
  query = combined,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = scRNA$re_annotation,
  weight.reduction = combined[['lsi']],
  dims = 2:30
)

combined <- AddMetaData(object = combined, metadata = predicted.labels)
plot1 <- DimPlot(object = scRNA, group.by = 're_annotation', label = TRUE,repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
plot2 <- DimPlot(  object = combined, group.by = 'predicted.id',  label = TRUE,  repel = TRUE)  + ggtitle('scATAC-seq')
plot3 <- DimPlot(  object = combined, group.by = 'seurat_clusters',  label = TRUE,  repel = TRUE)+ NoLegend()  + ggtitle('scATAC-seq')

pdf("./03_all_celltype/integrative_analysis_scRNA.pdf",width=17,height=5)
plot1 + plot2 + plot3
dev.off()

Idents(combined)<-combined$seurat_clusters
#####further annotation########
combined <- RenameIdents(
  object = combined,
  '0' = 'MP_DC',
  '1' = 'EC',
  '2' = 'MP_DC',
  '3' = 'EC',
  '4' = 'CM',
  '5' = 'MP_DC',
  '6' = 'FB',
  '7' = 'FB',
  '8' = 'CM',
  '9' = 'FB',
  '10' = 'FB',
  '11' = 'FB',
  '12' = 'CM',
  '13' = 'FB',
  '14' = 'Pericyte',
  '15' = 'EC',
  '16' = 'FB',
  '17' = 'FB',
  '18' = 'Epi',
  '19' = 'MP_DC',
  '20' = 'T',
  '21' = 'Gra',
  '22' = 'Epi',
  '23' = 'MP_DC',
  '24' = 'FB',
  '25' = 'FB',
  '26' = 'SMC',
  '27' = 'FB',
  '28' = 'SMC',
  '29' = 'EC',
  '30' = 'B',
  '31' = 'EC',
  '32' = 'Glial'
  )
combined@meta.data$Annotation<-Idents(combined)
table(combined$Annotation,combined$orig.ident)

combined$Annotation<-factor(combined$Annotation,levels=c("CM",'EC','FB','Epi','SMC','Pericyte',"MP_DC","T","B","Gra","Glial"))
Idents(combined)<-combined$Annotation;
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
combined$detail_anno<- paste(combined$seurat_clusters,combined$Annotation,sep=":")

pdf("./03_all_celltype/Annotated_allcelltype_UMAP.pdf",width=6,height=5)
DimPlot(combined, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap",group.by = "Annotation")
dev.off()
pdf("./03_all_celltype/Annotated_allcelltype_UMAP_detail.pdf",width=7,height=5)
DimPlot(combined, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap",group.by = "detail_anno")
dev.off()
object <- combined
Idents(object)<- object$detail_anno
embeddings <- Embeddings(object = object, reduction = "lsi")[,1:50]
data.dims <- lapply(X = levels(x = object), FUN = function(x) {
    cells <- WhichCells(object = object, idents = x)
    if (length(x = cells) == 1) {
        cells <- c(cells, cells)
    }
    temp <- colMeans(x = embeddings[cells, ])
})
data.dims <- do.call(what = "cbind", args = data.dims)
colnames(x = data.dims) <- levels(x = object)
library(lsa)
cosine_dist <- as.dist(1-cosine(data.dims))
data.tree <- ape::as.phylo(x = hclust(d = cosine_dist))
library(ggtree);

pdf("./03_all_celltype/AR3-combined-tree-cosine_detail_anno.pdf",width=6,height=6)
ggtree(data.tree,layout = "circular") + geom_tiplab()+ geom_treescale()
dev.off()


# Save rds have annotation information 
DefaultAssay(combined) <- "RNA"
saveRDS(combined,"./03_all_celltype/AR3_integrated_all_celltype_annotated.rds")

# Marker gene for different cell types (color on UMAP ) 
combined<- readRDS("./03_all_celltype/AR3_integrated_all_celltype_annotated.rds")
pdf('./03_all_celltype/All_Marker_gene_FeaturePlot_uamp.pdf', width=5.5, height=5)
for (i in 1:length(markerGenes)){
  p<-FeaturePlot(combined,order=F, reduction = 'umap',max.cutoff = 10, features = markerGenes[i], ncol = 1)
  print(p)
}
dev.off()

# call peak 

DefaultAssay(combined)<-"ATAC"
peak<-CallPeaks(
       combined,
       group.by = "detail_anno",
       macs2.path = "/public/home/nieyg/biosoft/conda/bin/macs2",
       broad = FALSE,
       format = "BED",
       fragment.tempdir = tempdir(),
       effective.genome.size = 2.7e+09,
       outdir="./03_all_celltype/recall_peak/",
       combine.peaks=TRUE
)
macs2_counts <- FeatureMatrix(
     fragments = Fragments(combined),
     features = peak,
     cells = colnames(combined)
     )     
#macs2_counts<-macs2_counts[-which(rownames(macs2_counts)=="GroupUN243-311766-311965"),]
# create a new assay using the MACS2 peak set and add it to the Seurat object
combined[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  #fragments = fragpath[[i]],
  fragments = Fragments(combined),
  annotation = Annotation(combined)
)

##Track for Marker genes promoters
Idents(combined)<-combined$detail_anno
library(BSgenome.Mmusculus.UCSC.mm10)
DefaultAssay(combined) <- "peaks"
# first compute the GC content for each peak
combined <- RegionStats(combined, genome = BSgenome.Mmusculus.UCSC.mm10)
#Annotation(combined)$tx_id <- Annotation(combined)$gene_name 

# link peaks to genes
combined <- LinkPeaks(
  object = combined,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = markerGenes
)
######Visulize track and RNA exp######
idents.plot <- Idents(combined)

pdf("./03_all_celltype/Marker_gene-peaktrack-RNAexp.pdf",height=16,width=8)
for(i in markerGenes){
  print(i)
  p1 <- CoveragePlot(
  object = combined,
  region = i,
  features = i,
  expression.assay = "RNA",
  idents = idents.plot,
  extend.upstream = 500,
  annotation=TRUE,
  extend.downstream = 500
)
print(p1)}
dev.off()

saveRDS(combined,"./03_all_celltype/03_recall_peak/AR3_integrated_all_celltype_annotated_recall_peak.rds")

