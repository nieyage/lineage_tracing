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
plot(density(OSN@meta.data$peak_region_fragments),xlim=c(0,100000))
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
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 40 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    nCount_ATAC > 1000 &
    nCount_ATAC < 100000  
)

# doubletFinder 
  dataset.name <-"OSN"
  sratDecontx <- OSN
  # Compute expected doublet rate
  cellN=nrow(sratDecontx@meta.data)
  expDoubletRate = (cellN*0.0008 + 0.0527)*0.01
  normalizationMethod='SCTransform'
  sweep.res.list_scData <- paramSweep_v3(sratDecontx, 
                                         PCs = 1:50, 
                                         sct = normalizationMethod == 'SCTransform', 
                                         num.cores = 4) #num.cores = 4
  sweep.stats_scData <- summarizeSweep(sweep.res.list_scData, GT = FALSE)
  bcmvn_scData <- find.pK(sweep.stats_scData)
  bcmvn_scData$pK <- as.numeric(as.character(bcmvn_scData$pK))
  pK1=bcmvn_scData$pK[bcmvn_scData$BCmetric==max(bcmvn_scData$BCmetric)]
  print(head(pK1))
  # PLOT: pK selection
  p1=ggplot(data=bcmvn_scData, 
            aes(x=pK, y=BCmetric, group=2)) +
    geom_line(color="blue")+
    geom_point()+
    geom_vline(xintercept=pK1, linetype="dashed", color = "red")+
    labs(title="pK Selection",x="pK", y = "BCmvn")+
    theme_classic()
  filename <- paste0("./01_doubletfinder/2a_doubletfinder_pkselection_","OSN",".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p1)
  # More doublet finder
  pK1=as.numeric(as.character( pK1 ))
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations <- sratDecontx@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)   
  nExp_poi <- round(expDoubletRate*nrow(sratDecontx@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  sratDecontx <- doubletFinder_v3( sratDecontx, PCs = sratDecontx@commands$RunUMAP.SCT.pca$dims,
                                   pN = 0.25, pK = pK1, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
  sratDecontx@meta.data$DoubletFinder =  sratDecontx@meta.data[,grep('DF.classifications', colnames( sratDecontx@meta.data))]
  # PLOT: Doublet Finder graphs
  p2 <- FeatureScatter(sratDecontx, feature1 = "nCount_ATAC", feature2 = "percent.mt", group.by = 'DoubletFinder')
  p3 <- FeatureScatter(sratDecontx, feature1 = "nCount_ATAC", feature2 = "nFeature_RNA", group.by = 'DoubletFinder')
  g = arrangeGrob(p2,p3, ncol = 2)
  filename <- paste0("./01_doubletfinder/2b_doubletfinder_ScatterPlots_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,g)
  rm(filename, g)
  # PLOT: Violin Plots
  p4 <- VlnPlot(sratDecontx, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), group.by = 'DoubletFinder', pt.size = 0)
  filename <- paste0("./01_doubletfinder/2c_doubletfinder_DoubletFinder-ViolinPlots_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p4)
  objList2[[i]] <- sratDecontx

# create the geneactivate matrix 

# ArchR 
library(ArchR)
addArchRThreads(threads = 16) 
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
projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "03_OSN_ArchR",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
getAvailableMatrices(projHeme1)
head(projHeme1$cellNames)
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

p
plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = projHeme1, addDOC = FALSE)

p1 <- plotFragmentSizes(ArchRProj = projHeme1)
p1
p2 <- plotTSSEnrichment(ArchRProj = projHeme1)
p2
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = projHeme1, outputDirectory = "Save-ProjHeme1", load = FALSE)
projHeme2 <- filterDoublets(projHeme1)




doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)



OSN
OSN <- RunTFIDF(OSN)
OSN <- FindTopFeatures(OSN, min.cutoff = 'q0')
OSN <- RunSVD(object = OSN)

DepthCor(OSN)
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

DimPlot(object = OSN, label = TRUE) + NoLegend()
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
DefaultAssay(OSN) <- 'RNA'
FeaturePlot(
  object = OSN,
  features = c('Sst','Pvalb',"Gad2","Neurod6","Rorb","Syt6"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
# Load the pre-processed scRNA-seq data
allen_rna <- readRDS("../vignette_data/allen_OSN.rds")
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




