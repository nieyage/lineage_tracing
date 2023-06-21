library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
set.seed(1234)

OSN<- readRDS("./03_All_celltype/OSN_ATAC_all_celltype.rds")
# compute gene activities
gene.activities <- GeneActivity(OSN)
# add the gene activity matrix to the Seurat object as a new assay
OSN[['RNA']] <- CreateAssayObject(counts = gene.activities)
OSN <- NormalizeData(
  object = OSN,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(OSN$nCount_RNA)
)

obj.integrated <- readRDS("/md01/yangjw28/projects/OSN_injury_download/obj_10_35_1.5.rds20222.0222")
# #Assigning cell type identity to clusters定义/合并群-----------------------------
new.cluster.ids <- c("HBC", "immatureOSNs","HBC","HBC","Sustentacular","Sustentacular","HBC","HBC","Sustentacular","actGBC","HBC","actGBC","HBC","immatureOSNs", "respiratoryHBC","respiratoryHBC","matureOSNs","INP","GBC","actHBC","actGBC","HBC","GBC","MicrovillarCell", "BrushCell","HBC","respiratoryHBC","BrushCell","respiratoryHBC")
names(new.cluster.ids) <- levels(obj.integrated)
obj.integrated <- RenameIdents(obj.integrated, new.cluster.ids)
lineage_df <- data.frame(new.cluster.ids=as.vector(Idents(obj.integrated)),row.names=names(Idents(obj.integrated )))
obj.integrated <- AddMetaData(obj.integrated,lineage_df)
obj.integrated$Annotation<- Idents(obj.integrated)
scRNA<- obj.integrated
pdf("./04_scRNA_published/integrated_scRNA-cluster-Umap.pdf",width=12,height=6)
p1<- DimPlot(scRNA,group.by = 'Annotation', label = TRUE,repel = TRUE) + NoLegend()
p2<- DimPlot(scRNA,group.by = 'seurat_clusters',  label = TRUE,repel = TRUE) + NoLegend()
p1+p2
dev.off()
# first annotation 
DefaultAssay(scRNA) <- 'RNA'
features <- c(
              "Omp","Gng13",#mOSN
              "Emx1",#atypical OSN
              "Nqo1","Ncam2",
              "Gap43","Gng8",#iOSN
              "Sox2","Ermn","Cyp2g1","Cbr2",#SUS
              "Ephx1",# SUS V 
              "Sult1c1", # SUS D
              "Neurog1","Neurod1",#INP
              "Ascl1","Kit" ,#GBC
              "Mki67","Top2a",#cycling
              "Trp63","Krt5","Krt14",#HBC
              "Cxcl14","Meg3", #olfactory HBCs were distinguished from respiratory HBCs
              "Reg3g", #respiratory
              "Ascl3",#olfactory microvillar cells
              "Krt18",
              "Cftr",# MV lonocyte
              "Trpm5" #MV brush cell
            
)
label<- c(
  "mOSN","mOSN","atypical_OSN",rep("iOSN",4),rep("SUS",4),"SUS_V","SUS_D","INP","INP",
  rep("GBC",2),rep("cycling",2),rep("HBC",3),rep("olfactory_HBC",2),"respiratory","olfactory_MC","MVloncyte","MVloncyte",
  "MV_brush"
  )
pdf("./04_scRNA_published/integrated_scRNA-cluster-annotation-all_celltype.pdf",width=16,height=8)
p<-DotPlot(scRNA, features = features,dot.scale = 3)+theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p
p&scale_x_discrete(labels=label)
dev.off()

allen_rna <- FindVariableFeatures(
  object = scRNA,
  nfeatures = 5000
)
DefaultAssay(OSN)<-"RNA"
transfer.anchors <- FindTransferAnchors(
  reference = allen_rna,
  query = OSN,
  reduction = 'cca',
  dims = 1:50
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = allen_rna$Annotation,
  weight.reduction = OSN[['lsi']],
  dims = 2:30
)

OSN <- AddMetaData(object = OSN, metadata = predicted.labels)
plot1 <- DimPlot(allen_rna, group.by = 'Annotation', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
plot2 <- DimPlot(OSN, group.by = 'predicted.id', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
pdf("./04_scRNA_published/integrated_scRNA-annotated_cluster.pdf",width=10,height=5)
plot1 + plot2
dev.off()

###annotation 
DefaultAssay(OSN) <- 'RNA'
pdf("./03_All_celltype/cluster_dotplot_OSN_all.pdf")
p<-DotPlot(OSN, features = features,dot.scale = 3)+theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p
p&scale_x_discrete(labels=label)
dev.off()
object <- OSN
Idents(object)<- object$predicted.id
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
pdf("./03_All_celltype/OSN-tree-cosine_detail_anno.pdf",width=6,height=6)
ggtree(data.tree,layout = "circular") + geom_tiplab()+ geom_treescale()
dev.off()

# add annotation info in OSN 
Idents(OSN)<-OSN$predicted.id
#####further annotation########
OSN <- RenameIdents(
  object = OSN,
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
OSN@meta.data$Annotation<-Idents(OSN)
table(OSN$Annotation,OSN$orig.ident)

Idents(OSN)<-OSN$Annotation;
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

pdf("./03_all_celltype/Annotated_allcelltype_UMAP.pdf",width=6,height=5)
DimPlot(OSN, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap",group.by = "Annotation")
dev.off()

saveRDS(scRNA,"./04_scRNA_published/scRNA_published.rds")


# Marker gene for different cell types (color on UMAP ) 
pdf('./03_All_celltype/All_Marker_gene_FeaturePlot_uamp.pdf', width=5.5, height=5)
for (i in 1:length(features)){
  p<-FeaturePlot(OSN,order=F, reduction = 'umap',max.cutoff = 10, features = features[i], ncol = 1)
  print(p)
}
dev.off()

# call peak 
DefaultAssay(OSN)<-"ATAC"
peak<-CallPeaks(
       OSN,
       group.by = "Annotation",
       macs2.path = "/public/home/nieyg/biosoft/conda/bin/macs2",
       broad = FALSE,
       format = "BED",
       fragment.tempdir = tempdir(),
       effective.genome.size = 2.7e+09,
       outdir="./03_All_celltype/01_recall_peak/",
       combine.peaks=TRUE
)
macs2_counts <- FeatureMatrix(
     fragments = Fragments(OSN),
     features = peak,
     cells = colnames(OSN)
     )     
#macs2_counts<-macs2_counts[-which(rownames(macs2_counts)=="GroupUN243-311766-311965"),]
# create a new assay using the MACS2 peak set and add it to the Seurat object
OSN[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  #fragments = fragpath[[i]],
  fragments = Fragments(OSN),
  annotation = Annotation(OSN)
)

##Track for Marker genes promoters
Idents(OSN)<-OSN$Annotation
library(BSgenome.Mmusculus.UCSC.mm10)
DefaultAssay(OSN) <- "peaks"
# first compute the GC content for each peak
OSN <- RegionStats(OSN, genome = BSgenome.Mmusculus.UCSC.mm10)
# link peaks to genes
OSN <- LinkPeaks(
  object = OSN,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = features
)
######Visulize track and RNA exp######
idents.plot <- Idents(OSN)
pdf("./03_All_celltype/Marker_gene-peaktrack-RNAexp.pdf",height=16,width=8)
for(i in features){
  print(i)
  p1 <- CoveragePlot(
  object = OSN,
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
saveRDS(OSN,"./03_All_celltype/OSN_all_celltype_annotated_recall_peak.rds")






