
# 1. load in signac 
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)

# combine peaksets 
# read in peak sets
peaks.AR3_C4_scATAC <- read.table(
  file = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.AR3_C5_scATAC <- read.table(
  file = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
# convert to genomic ranges
gr.AR3_C4_scATAC <- makeGRangesFromDataFrame(peaks.AR3_C4_scATAC)
gr.AR3_C5_scATAC <- makeGRangesFromDataFrame(peaks.AR3_C5_scATAC)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.AR3_C4_scATAC, gr.AR3_C5_scATAC))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# load metadata
md.AR3_C4_scATAC <- read.table(
  file = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

md.AR3_C5_scATAC <- read.table(
  file = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

# perform an initial filtering of low count cells
md.AR3_C4_scATAC <- md.AR3_C4_scATAC[md.AR3_C4_scATAC$passed_filters > 500, ]
md.AR3_C5_scATAC <- md.AR3_C5_scATAC[md.AR3_C5_scATAC$passed_filters > 500, ]

# create fragment objects
frags.AR3_C4_scATAC <- CreateFragmentObject(
  path = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC/outs/fragments.tsv.gz",
  cells = rownames(md.AR3_C4_scATAC)
)
frags.AR3_C5_scATAC <- CreateFragmentObject(
  path = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC/outs/fragments.tsv.gz",
  cells = rownames(md.AR3_C5_scATAC)
)

AR3_C4_scATAC.counts <- FeatureMatrix(
  fragments = frags.AR3_C4_scATAC,
  features = combined.peaks,
  cells = rownames(md.AR3_C4_scATAC)
)
AR3_C5_scATAC.counts <- FeatureMatrix(
  fragments = frags.AR3_C5_scATAC,
  features = combined.peaks,
  cells = rownames(md.AR3_C5_scATAC)
)

AR3_C4_assay <- CreateChromatinAssay(AR3_C4_scATAC.counts, fragments = frags.AR3_C4_scATAC)
AR3_C4 <- CreateSeuratObject(AR3_C4_assay, assay = "ATAC", meta.data=md.AR3_C4_scATAC)

AR3_C5_assay <- CreateChromatinAssay(AR3_C5_scATAC.counts, fragments = frags.AR3_C5_scATAC)
AR3_C5 <- CreateSeuratObject(AR3_C5_assay, assay = "ATAC", meta.data=md.AR3_C5_scATAC)


# 2. calculate the qc feature 

# 3. doubletfinder(ArchR or DoubletFinder?)

# 4. filter cells (low quality and doublets)

# 5. intergrated object and analyze 

# 6. annotate all cell type in AR3 

