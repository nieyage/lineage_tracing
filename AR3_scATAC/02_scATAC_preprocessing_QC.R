#####################
# 1. load in signac #
#####################
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
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

###############################
# 2. calculate the qc feature #
###############################
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
objList <- c(AR3_C4,AR3_C5)

# calculate the score of NS and TSS
for (i in seq_len(length(objList))) {
    Annotation(objList[[i]]) <- annotations # add the gene information to the object
    objList[[i]] <- NucleosomeSignal(objList[[i]])
    objList[[i]] <- TSSEnrichment(objList[[i]],fast=FALSE)
    }

# plot TSS and fragrament distribution plot 
  pdf("./01_qc/TSS_distribution_C4.pdf")
  objList[[1]]$high.tss<-ifelse(objList[[1]]$TSS.enrichment > 2, 'High', 'Low')
  TSS<-TSSPlot(objList[[1]], group.by = 'high.tss') + NoLegend()+ labs(title = "NE")
  objList[[1]]$NS <- ifelse(objList[[1]]$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  Frag<-FragmentHistogram(object = objList[[1]], group.by = 'NS',region = "chr1-1-195471971")+ labs(title = "NE")
  print(TSS);
  print(Frag);
  dev.off();
  pdf("./01_qc/TSS_distribution_C5.pdf")
  objList[[2]]$high.tss<-ifelse(objList[[2]]$TSS.enrichment > 2, 'High', 'Low')
  TSS<-TSSPlot(objList[[2]], group.by = 'high.tss') + NoLegend()+ labs(title = "Nurse")
  objList[[2]]$NS <- ifelse(objList[[2]]$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  Frag<-FragmentHistogram(object = objList[[2]], group.by = 'NS',region = "chr1-1-195471971")+ labs(title = "Nurse")
  print(TSS);
  print(Frag);
  dev.off();

objList[[1]]$pct_reads_in_peaks <- objList[[1]]$peak_region_fragments / objList[[1]]$passed_filters * 100
objList[[1]]$blacklist_ratio <- objList[[1]]$blacklist_region_fragments / objList[[1]]$peak_region_fragments
pdf("./01_qc/QC_NS_TSS_density_C4.pdf",width=8,height=4)
plot(density(objList[[1]]@meta.data$nucleosome_signal),xlim=c(0,4))
plot(density(objList[[1]]@meta.data$TSS.enrichment),xlim=c(0,20))
plot(density(objList[[1]]@meta.data$nCount_ATAC),xlim=c(0,200000))
plot(density(objList[[1]]@meta.data$nFeature_ATAC),xlim=c(0,50000))
plot(density(objList[[1]]@meta.data$pct_reads_in_peaks),xlim=c(0,100))
plot(density(objList[[1]]@meta.data$blacklist_ratio),xlim=c(0,20000))
dev.off()

objList[[2]]$pct_reads_in_peaks <- objList[[2]]$peak_region_fragments / objList[[2]]$passed_filters * 100
objList[[2]]$blacklist_ratio <- objList[[2]]$blacklist_region_fragments / objList[[2]]$peak_region_fragments
pdf("./01_qc/QC_NS_TSS_density_C5.pdf",width=8,height=4)
plot(density(objList[[2]]@meta.data$nucleosome_signal),xlim=c(0,4))
plot(density(objList[[2]]@meta.data$TSS.enrichment),xlim=c(0,20))
plot(density(objList[[2]]@meta.data$nCount_ATAC),xlim=c(0,200000))
plot(density(objList[[2]]@meta.data$nFeature_ATAC),xlim=c(0,50000))
plot(density(objList[[2]]@meta.data$pct_reads_in_peaks),xlim=c(0,100))
plot(density(objList[[2]]@meta.data$blacklist_ratio),xlim=c(0,20000))
dev.off()

sample<-c("AR3_C4","AR3_C5")
for (i in seq_len(length(objList))) {
  # plot QC plot 
  pdf(file = paste("./01_qc/QC_before_Vlnplot_",sample[i],".pdf", sep = ""),width=10,height=6)
  qc<-VlnPlot(object = objList[[i]],
            features = c('pct_reads_in_peaks', 'peak_region_fragments','TSS.enrichment', 'nucleosome_signal'),
            ncol = 3,
            pt.size = 0.01
          )
  print(qc)
  dev.off();
   }

#############################################
# 3. doubletfinder(ArchR) #
#############################################
library(ArchR)
addArchRThreads(threads = 1) 
addArchRGenome("mm10")
inputFiles <- c("/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC/outs/fragments.tsv.gz",
  "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC/outs/fragments.tsv.gz")
names(inputFiles)<- c("AR3_C4","AR3_C5")
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
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)
projHeme2 <- filterDoublets(projHeme1)

filterDoublets_C4<- projHeme2$cellNames[grep("AR3_C4",projHeme2$cellNames)]
filterDoublets_C4<- gsub("AR3_C4#","",filterDoublets_C4)
filterDoublets_C5<- projHeme2$cellNames[grep("AR3_C5",projHeme2$cellNames)]
filterDoublets_C5<- gsub("AR3_C5#","",filterDoublets_C5)



##############################################
# 4. filter cells (low quality and doublets) #
##############################################

# filter out low quality cells
# To remove doublets,select different cutoff#####
# NE 
  objList2<-c()
  objList2[[1]]<-subset(x=objList[[1]],
   subset = nCount_ATAC < 100000 &nCount_ATAC > 1000 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    peak_region_fragments > 1000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 25 
    )
  print(objList2[[1]])
# Nurse 
  objList2[[2]]<-subset(x=objList[[2]],
   subset = nCount_ATAC < 100000 &nCount_ATAC > 1000 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    peak_region_fragments > 1000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 25 
    )
  print(objList2[[2]])

# Subset identified singlets
last_C4<-intersect(colnames(objList2[[1]]),filterDoublets_C4)
AR3_C4_last<- subset(objList2[[1]],cells=last_C4)
last_C5<-intersect(colnames(objList2[[2]]),filterDoublets_C5)
AR3_C5_last<- subset(objList2[[2]],cells=last_C5)
saveRDS(AR3_C4_last,"./03_all_celltype/AR3_C4_scATAC.rds")
saveRDS(AR3_C5_last,"./03_all_celltype/AR3_C5_scATAC.rds")

#####################################
# 5. intergrated object and analyze #
#####################################
combined <- merge(
  x = AR3_C4_last,
  y = AR3_C5_last,
  add.cell.ids = c("AR3_C4", "AR3_C5")
)
combined[["ATAC"]]
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
saveRDS(combined,"./03_all_celltype/combined_AR3_scATAC.rds")

####################################
# 6. annotate all cell type in AR3 #
####################################

