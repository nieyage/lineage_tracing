library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
OSN<- readRDS("./03_All_celltype/OSN_all_celltype_annotated_recall_peak.rds")

write.table(colnames(OSN),"~/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked/mgatk_for_filtered_cells/barcode.tsv",row.names=F,col.names=F)

# load mgatk output
mito.data <- ReadMGATK(dir = "/md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked/mgatk_for_filtered_cells/final/")
# create an assay
mito_counts<- mito.data$counts
mito <- CreateAssayObject(counts = mito_counts)
# Subset to cell present in the scATAC-seq assat
mito <- subset(mito, cells = colnames(OSN))
# add assay and metadata to the seurat object
OSN[["mito"]] <- mito

OSN <- AddMetaData(OSN, metadata = mito.data$depth, col.name = "mtDNA_depth")
Idents(OSN)<- OSN$Annotation;

pdf("./03_All_celltype/02_mgatk/all_cell_type_mtDNA_depth.pdf",width=16,height=8)
VlnPlot(OSN, c("mtDNA_depth","nCount_mito","nFeature_mito"), pt.size = 0,ncol=1) + scale_y_log10()
dev.off()
saveRDS(OSN,"./03_All_celltype/02_mgatk/all_cell_type_mgatk.rds")

OSN<- readRDS("./03_All_celltype/02_mgatk/all_cell_type_mgatk.rds")
# 1. mtDNA coverage among celltype 
# load mgatk output
cov.data <- read.table(gzfile("/md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked/mgatk_for_filtered_cells/final/OSN.coverage.txt.gz"),sep=",")

# add prefix in barcode 
cov.data$celltype<- OSN$Annotation[match(cov.data$V2,rownames(OSN@meta.data))]
# Gental smooth of just 1 bp for plot aesthetics
library(roll)
smooth <- 5
cov.data$V3<- roll_mean(cov.data$V3, smooth)

mdf <- cov.data[,-2]
colnames(mdf)<- c("pos","coverage","celltype")
mdf[which(is.na(mdf$coverage)),]$coverage=0

df<- aggregate(mdf$coverage,by=list(mdf$pos,mdf$celltype),mean)
colnames(df)<- c("pos","celltype","coverage")

# Mean coverage X
aggregate(mdf$coverage,by=list(mdf$celltype),mean)
           Group.1         x
1              HBC  7.435169
2           actGBC  5.048911
3              GBC  6.581020
4     immatureOSNs  7.833693
5    Sustentacular 16.292197
6   respiratoryHBC 11.432890
7        BrushCell 10.179746
8       matureOSNs  7.395644
9           actHBC  9.620186
10 MicrovillarCell 11.431367
11             INP 10.242323

# the coverage plot 
library(ComplexHeatmap)
library(GenomeInfoDb)
library(BuenColors) 
# Theme to remove any of the other plot riff raff
xxtheme <-   theme(
  axis.line = element_blank(),
  axis.ticks.y = element_blank(),   
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.margin = unit(c(0, 0, 0, 0), "cm"),
  plot.margin = unit(c(-0.35, -0.35, -0.35, -0.35), "cm")) 
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
# Visualize the rolled means

P1 <- ggplot(df, aes(x = pos, y = coverage, color = celltype)) + 
  geom_line() +  expand_limits(y = c(-5, 4)) +
  pretty_plot(fontsize = 8)  + scale_color_manual(values = myUmapcolors[1:11]) +
  coord_polar(direction = 1) + labs(x = "", y = "Coverage") + scale_y_log10() #+
  #theme(legend.position = "none")
cowplot::ggsave2(P1, file = "./03_All_celltype/02_mgatk/All_celltype_rollMean_coverage.pdf", width = 4, height = 4)

# IdentifyVariants
# filter cells based on mitochondrial depth
mito.data <- ReadMGATK(dir = "/md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked/mgatk_for_filtered_cells/final/")
crc <- OSN
variable.sites <- IdentifyVariants(crc, assay = "mito", refallele = mito.data$refallele)
pdf("./03_All_celltype/02_mgatk/VariantPlot_global.pdf")
VariantPlot(variants = variable.sites)
dev.off()


# Establish a filtered data frame of variants based on this processing
high.conf <- subset(
  variable.sites, 
  subset = n_cells_conf_detected >= 2 &
    strand_correlation >= 0.65 &
    vmr > 0.01
)
high.conf<- high.conf[order(high.conf$mean,decreasing=T),]
high.conf[,c(1,2,5)]
# global VAF distribution of mutation 
# the mutation freq barplot 
OSN <- AlleleFreq(
    object = OSN,
    variants = high.conf$variant,
    assay = "mito"
)

DefaultAssay(OSN)<- "alleles"
data<- as.numeric(GetAssayData(OSN))
data<- data[data!=0]
data<- data.frame(variance=data)
data$variance<- data$variance*100;
for(i in 1:nrow(data)){
  if(data$variance[i] < 0.1){data$VAF[i]<- "0.0-0.1%"}
  if(data$variance[i] >0.1  && data$variance[i] < 0.5){data$VAF[i]<- "0.1-0.5%"}
  if(data$variance[i] >0.5  && data$variance[i] < 1){data$VAF[i]<- "0.5-1%"}
  if(data$variance[i] >1    && data$variance[i] < 90){data$VAF[i]<- "1-90%"}
  if(data$variance[i] >90   && data$variance[i] < 100){data$VAF[i]<- "90-100%"}
}

pdf("./03_All_celltype/02_mgatk/All_celltype_mutation_freq.pdf",width=5,height=3)
ggplot(data=data,aes(x=VAF))+
geom_histogram(stat="count",fill="#795885",color="#795885",alpha=0.8)+ 
ylab("mtDNA mutations")+theme_bw()+xlab("% VAF")+
ggtitle("All celltype")
ggplot(data=data,aes(x=variance))+
geom_histogram(fill="#795885",color="#795885",alpha=0.8)+ 
ylab("mtDNA mutations")+theme_bw()+xlab("% VAF")+
ggtitle("All celltype")
dev.off()

# VAF distribution of mutation in difference celltype 
# the mutation freq barplot 
Idents(OSN)<- OSN$Annotation
pdf("./03_All_celltype/02_mgatk/celltype_mutation_freq.pdf",width=5,height=3)
for (j in levels(OSN)){
  print(j);
  obj<- subset(OSN,idents=j)
  variable.sites <- IdentifyVariants(obj, assay = "mito", refallele = mito.data$refallele)
  high.conf <- subset(
  variable.sites, 
  subset = n_cells_conf_detected >= 5 &
    strand_correlation >= 0.65 &
    vmr > 0.01
  )
  p1 <- VariantPlot(variants = variable.sites,concordance.threshold=0.65)
  print(p1)
  if(nrow(high.conf)>1){
  crc <- AlleleFreq(
    object = obj,
    variants = high.conf$variant,
    assay = "mito"
  )
  DefaultAssay(crc) <- "alleles"
  data<- GetAssayData(crc)
  d_varition1=as.vector(data)
  d_varition1=d_varition1[which(d_varition1!=0)]
  data<- data.frame(variance=d_varition1)
  data$variance<- data$variance;
  p=ggplot(data,aes(x=variance))+
  geom_histogram(
                 binwidth = 0.1,
                 fill="#69b3a2",##69b3a2
                 color="#e9ecef",##e9ecef
                 alpha=0.9,
                 breaks=seq(0,1,0.1))+ 
  theme_bw()+
  labs(x="Frequence",y="Count",title=j)+
  theme(plot.title = element_text(size = 18, vjust = 0.5, hjust = 0.5, face = "bold"),
        axis.text.x = element_text(size=18,colour="black", face = "bold"),          
        axis.text.y = element_text(size=18,colour="black", face = "bold"),
        axis.title.x = element_text(size=14, face = "bold"), 
        axis.title.y = element_text(size=14, face = "bold"),
        panel.background = element_blank(),
        line = element_line(size=1),
        axis.line = element_line(size =1.0,colour = "black"),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank())+scale_x_continuous(limits = c(0,1),
                     breaks = seq(0,1,0.2))
    print(p)
}}
dev.off()

Idents(OSN)<- OSN$Annotation
OSN$mito_reads_rate<- (OSN$atac_mitochondrial_reads/OSN$atac_raw_reads)*100
library(scCustomize)

pdf("./03_All_celltype/02_mgatk/all_cell_type_mtDNA_depth.pdf",width=20,height=10)
VlnPlot(OSN, c("mito_reads_rate","mtDNA_depth","nCount_mito","nFeature_mito"), pt.size = 0,ncol=1) + scale_y_log10()
Stacked_VlnPlot(seurat_object = OSN, features = c("mito_reads_rate","mtDNA_depth","nCount_mito","nFeature_mito","nCount_alleles","nFeature_alleles"), x_lab_rotate = TRUE,
    colors_use = myUmapcolors)
dev.off()
saveRDS(OSN,"./03_All_celltype/02_mgatk/all_cell_type_mgatk.rds")

DefaultAssay(OSN) <- "mito"
OSN <- FindClonotypes(OSN)



OSN2 <- AlleleFreq(OSN, variants = high.conf$variant, assay = "mito")
OSN[["alleles"]]
DefaultAssay(OSN) <- "alleles"

DoHeatmap(OSN, features = rownames(high.conf), disp.max = 0.1) +scale_fill_viridis_c()




fastqs,sample,library_type
/md01/pengdj3/WT_mNeu/workdir/GEX_fastq,mNeu12-RNA,Gene Expression
/md01/pengdj3/WT_mNeu/workdir/ATAC_fastq,mNeu12-ATAC,Chromatin Accessibility






