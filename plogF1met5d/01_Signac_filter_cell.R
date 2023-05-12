library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
set.seed(1234)

counts <- Read10X_h5("/md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/outs/filtered_feature_bc_matrix.h5")

OSN_assay <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = '/md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/outs/atac_fragments.tsv.gz',
  min.cells= 10, min.features=200
  )
metadata <- read.csv(
  file = "/data/R02/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/outs/per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)
OSN <- CreateSeuratObject(
  counts = OSN_assay,
  assay = 'ATAC',
  project = 'plogF1met5d_masked',
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

pdf("./QC_NS_TSS_density.pdf",width=8,height=4)
plot(density(OSN@meta.data$nucleosome_signal),xlim=c(0,0.6))
plot(density(OSN@meta.data$TSS.enrichment),xlim=c(0,20))
plot(density(OSN@meta.data$nCount_ATAC),xlim=c(0,10000))
plot(density(OSN@meta.data$nFeature_ATAC),xlim=c(0,10000))
dev.off()

pdf("./QC_NS_TSS.pdf",width=8,height=4)
OSN$NS <- ifelse(OSN$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = OSN, group.by = 'NS',region = "chr1-1-195471971")
OSN$high.tss <- ifelse(OSN$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(OSN, group.by = 'high.tss') + NoLegend()
dev.off()

OSN$pct_reads_in_peaks <- OSN$atac_peak_region_fragments / OSN$atac_fragments * 100
OSN$blacklist_ratio <- OSN$blacklist_region_fragments / OSN$atac_peak_region_fragments
pdf("./QC_Vln_plot.pdf",width=8,height=6)
VlnPlot(
  object = OSN,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

# filter cells 

OSN <- subset(
  x = OSN,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.025 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
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




