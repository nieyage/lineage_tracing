library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
set.seed(1234)

counts <- Read10X_h5("/md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked/outs/filtered_feature_bc_matrix.h5")

OSN_assay <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = '/md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked/outs/atac_fragments.tsv.gz',
  min.cells= 10, min.features=200
  )
metadata <- read.csv(
  file = "/md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked/outs/per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)
OSN <- CreateSeuratObject(
  counts = OSN_assay,
  assay = 'ATAC',
  project = 'OSN_regeneration',
  meta.data= metadata
  )

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(OSN) <- annotations
OSN <- NucleosomeSignal(object = OSN)
OSN <- TSSEnrichment(OSN, fast = FALSE)

pdf("./01_qc/QC_NS_TSS.pdf",width=8,height=4)
OSN$NS <- ifelse(OSN$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = OSN, group.by = 'NS',region = "chr1-1-195471971")
OSN$high.tss <- ifelse(OSN$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(OSN, group.by = 'high.tss') + NoLegend()
dev.off()

OSN$pct_reads_in_peaks <- OSN$atac_peak_region_fragments / OSN$atac_fragments * 100
pdf("./01_qc/QC_NS_TSS_density.pdf",width=8,height=4)
plot(density(OSN@meta.data$nucleosome_signal),xlim=c(0,4))
plot(density(OSN@meta.data$TSS.enrichment),xlim=c(0,20))
plot(density(OSN@meta.data$nCount_ATAC),xlim=c(0,200000))
plot(density(OSN@meta.data$nFeature_ATAC),xlim=c(0,50000))
plot(density(OSN@meta.data$pct_reads_in_peaks),xlim=c(0,100))
plot(density(OSN@meta.data$atac_peak_region_fragments),xlim=c(0,20000))
dev.off()

#OSN$blacklist_ratio <- OSN$blacklist_region_fragments / OSN$atac_peak_region_fragments
pdf("./01_qc/QC_Vln_plot.pdf",width=8,height=6)
VlnPlot(
  object = OSN,
  features = c('pct_reads_in_peaks', 'peak_region_fragments','TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 4
)
dev.off()

# filter cells 
OSN <- subset(
  x = OSN,
  subset = atac_peak_region_fragments > 1000 &
    atac_peak_region_fragments < 100000 &
    pct_reads_in_peaks > 40 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    nCount_ATAC > 1000 &
    nCount_ATAC < 100000  
)

# ArchR 
library(ArchR)
addArchRThreads(threads = 1) 
addArchRGenome("mm10")
inputFiles <- "/md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked/outs/atac_fragments.tsv.gz"
names(inputFiles)<-"OSN"

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)

projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "02_OSN_ArchR",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
getAvailableMatrices(projHeme1)
head(projHeme1$cellNames)
df <- getCellColData(projHeme1, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(
    x = df[,1], 
    y = df[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = projHeme1, addDOC = FALSE)

p1 <- plotFragmentSizes(ArchRProj = projHeme1)
p2 <- plotTSSEnrichment(ArchRProj = projHeme1)
plotPDF(p1,p2, name = "./02_OSN_ArchR/QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = projHeme1, outputDirectory = "Save-ProjHeme1-for-filterdoublets", load = FALSE)

projHeme2 <- filterDoublets(projHeme1)
filterDoublets<- gsub("OSN#","",projHeme2$cellNames)
last<-intersect(colnames(OSN),filterDoublets)
OSN_last<- subset(OSN,cells=last)

OSN <- RunTFIDF(OSN_last)
OSN <- FindTopFeatures(OSN, min.cutoff = 'q0')
OSN <- RunSVD(object = OSN)
pdf("./03_All_celltype/DepthCor.pdf")
DepthCor(OSN)
dev.off()
OSN <- RunUMAP(
  object = OSN,
  reduction = 'lsi',
  dims = 2:30
)
OSN <- FindNeighbors(
  object = OSN,
  reduction = 'lsi',
  dims = 2:30
)
OSN <- FindClusters(
  object = OSN,
  algorithm = 3,
  resolution = 1.2,
  verbose = FALSE
)
table(OSN$seurat.clusters)
pdf("./03_All_celltype/OSN_ATAC_Umap.pdf",width=6,height=6)
DimPlot(object = OSN, label = TRUE) + NoLegend()
dev.off()
saveRDS(OSN, file = "./03_All_celltype/OSN_ATAC_all_celltype.rds")

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

## Load the pre-processed scRNA-seq data
#allen_rna <- readRDS("../vignette_data/allen_OSN.rds")
#allen_rna <- FindVariableFeatures(
#  object = allen_rna,
#  nfeatures = 5000
#)
#
#transfer.anchors <- FindTransferAnchors(
#  reference = allen_rna,
#  query = OSN,
#  reduction = 'cca',
#  dims = 1:40
#)
#
#predicted.labels <- TransferData(
#  anchorset = transfer.anchors,
#  refdata = allen_rna$subclass,
#  weight.reduction = OSN[['lsi']],
#  dims = 2:30
#)
#
#OSN <- AddMetaData(object = OSN, metadata = predicted.labels)
#plot1 <- DimPlot(allen_rna, group.by = 'subclass', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
#plot2 <- DimPlot(OSN, group.by = 'predicted.id', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
#plot1 + plot2

###annottaion 
DefaultAssay(OSN) <- 'RNA'
pdf("./03_All_celltype/cluster_dotplot_OSN_all.pdf")
features <- c(
              "Omp","Gng13",#成熟OSN
              "Nqo1","Ncam2",
              "Gap43","Gng8",#不成熟OSN
              "Sox2","Ermn","Cyp2g1",#支持细胞
              "Neurog1","Neurod1",#INP
              "Ascl1","Kit" ,#GBC
              "Mki67","Top2a",#cycling
              "Trp63","Krt5","Krt14",#HBC
              "Cxcl14","Meg3", #olfactory HBCs were distinguished from respiratory HBCs
              "Reg3g", #respiratory
              "Ascl3","Cftr",#olfactory microvillar cells
              "Krt18","Trpm5" #brush cell
            
)
DotPlot(OSN, features = features,dot.scale = 3) + RotatedAxis() + theme(axis.text.x=element_text(size=8))
dev.off()
saveRDS(OSN, file = "./03_All_celltype/OSN_ATAC_GeneActivity_all_celltype.rds")


