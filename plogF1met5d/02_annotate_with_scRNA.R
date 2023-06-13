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
OSN[['gene.activities']] <- CreateAssayObject(counts = gene.activities)

OSN <- NormalizeData(
  object = OSN,
  assay = 'gene.activities',
  normalization.method = 'LogNormalize',
  scale.factor = median(OSN$nCount_RNA)
)

## Load the pre-processed scRNA-seq data
raw_counts <- read.csv("/md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/03_publish_data/GSE157068_countMatrix_k5.csv",check.names=FALSE)


allen_rna <- FindVariableFeatures(
  object = allen_rna,
  nfeatures = 5000
)

transfer.anchors <- FindTransferAnchors(
  reference = allen_rna,
  query = OSN,
  reduction = 'cca',
  dims = 1:40
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = allen_rna$subclass,
  weight.reduction = OSN[['lsi']],
  dims = 2:30
)

OSN <- AddMetaData(object = OSN, metadata = predicted.labels)
plot1 <- DimPlot(allen_rna, group.by = 'subclass', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
plot2 <- DimPlot(OSN, group.by = 'predicted.id', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
plot1 + plot2

###annottaion 
DefaultAssay(OSN) <- 'RNA'
pdf("./03_All_celltype/cluster_dotplot_OSN_all.pdf")
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
              "Ascl3",,#olfactory microvillar cells
              "Krt18",
              "Cftr",# MV lonocyte
              "Trpm5" #MV brush cell
            
)
DotPlot(OSN, features = features,dot.scale = 3) + RotatedAxis() + theme(axis.text.x=element_text(size=8))
dev.off()
saveRDS(OSN, file = "./03_All_celltype/OSN_ATAC_GeneActivity_all_celltype.rds")


