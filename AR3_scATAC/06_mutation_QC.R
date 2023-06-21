library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
combined<- readRDS("./03_all_celltype/04_mgatk/all_cell_type_mgatk.rds")
# 1. mtDNA coverage among celltype 

# load mgatk output
AR3_C4_cov.data <- read.table(gzfile("~/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC/mgatk_for_filtered_cells/final/AR3_C4_scATAC.coverage.txt.gz"),sep=",")
AR3_C5_cov.data <- read.table(gzfile("~/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC/mgatk_for_filtered_cells/final/AR3_C5_scATAC.coverage.txt.gz"),sep=",")

# add prefix in barcode 
AR3_C4_cov.data$V2<- paste("AR3_C4_",AR3_C4_cov.data$V2,sep="")
AR3_C5_cov.data$V2<- paste("AR3_C5_",AR3_C5_cov.data$V2,sep="")
cov.data<- rbind(AR3_C4_cov.data,AR3_C5_cov.data)
cov.data$celltype<- combined$Annotation[match(cov.data$V2,rownames(combined@meta.data))]
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
       CM 65.924376
       EC 35.275374
       FB 34.965598
      Epi 35.550709
      SMC 26.057016
 Pericyte 33.800104
    MP_DC 15.968654
        T 19.030238
        B 23.596686
      Gra  9.830597
    Glial 21.240215

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
P1 <- ggplot(df, aes(x = pos, y = log2(coverage+1), color = celltype)) + 
  geom_line() +  expand_limits(y = c(-5, 4)) +
  pretty_plot(fontsize = 8)  + scale_color_manual(values = myUmapcolors[1:11]) +
  coord_polar(direction = 1) + labs(x = "", y = "log2 Coverage") + scale_y_log10() +
  theme(legend.position = "none") + xxtheme 
cowplot::ggsave2(P1, file = "./03_all_celltype/04_mgatk/rollMean_coverage.pdf", width = 4, height = 4)

P1 <- ggplot(df, aes(x = pos, y = log2(coverage+1), color = celltype)) + 
  geom_line() +  expand_limits(y = c(-5, 4)) +
  pretty_plot(fontsize = 8)  + scale_color_manual(values = myUmapcolors[1:11]) +
  coord_polar(direction = 1) + labs(x = "", y = "log2 Coverage") + scale_y_log10() +
  theme(legend.position = "none")
cowplot::ggsave2(P1, file = "./03_all_celltype/04_mgatk/rollMean_coverage2.pdf", width = 4, height = 4)

# IdentifyVariants
# filter cells based on mitochondrial depth
#crc <- subset(combined, mtDNA_depth >= 10)
AR3_C4_mito.data <- ReadMGATK(dir = "~/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC/mgatk_for_filtered_cells/final/")
crc <- combined
variable.sites <- IdentifyVariants(crc, assay = "mito", refallele = AR3_C4_mito.data$refallele)
pdf("./03_all_celltype/04_mgatk/VariantPlot_global.pdf")
VariantPlot(variants = variable.sites)
dev.off()

# all celltype mutation 
# No homoplasmic mutation?---> sequence depth not enough!

# Establish a filtered data frame of variants based on this processing
high.conf <- subset(
	variable.sites, 
	subset = n_cells_conf_detected >= 5 &
    strand_correlation >= 0.65 &
    vmr > 0.01
)
high.conf<- high.conf[order(high.conf$mean,decreasing=T),]
high.conf[,c(1,2,5)]
# global VAF distribution of mutation 
# the mutation freq barplot 
DefaultAssay(combined)<- "alleles"
data<- as.numeric(GetAssayData(combined))
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

pdf("./03_all_celltype/04_mgatk/All_celltype_mutation_freq.pdf",width=5,height=3)
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
Idents(combined)<- combined$Annotation
pdf("./03_all_celltype/04_mgatk/celltype_mutation_freq.pdf",width=5,height=3)
for (j in levels(combined)){
	print(j);
	obj<- subset(combined,idents=j)
	variable.sites <- IdentifyVariants(obj, assay = "mito", refallele = AR3_C4_mito.data$refallele)
  high.conf <- subset(
  variable.sites, 
  subset = n_cells_conf_detected >= 5 &
    strand_correlation >= 0.5 &
    vmr > 0.01
  )
  p1 <- VariantPlot(variants = variable.sites,concordance.threshold=0.5)
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


Idents(combined)<- combined$detail_anno
combined$mito_reads_rate<- (combined$mitochondrial/combined$total)*100
library(scCustomize)

pdf("./03_all_celltype/04_mgatk/all_cell_type_mtDNA_depth.pdf",width=20,height=10)
VlnPlot(combined, c("mito_reads_rate","mtDNA_depth","nCount_alleles","nFeature_alleles"), pt.size = 0,ncol=1) + scale_y_log10()
Stacked_VlnPlot(seurat_object = combined, features = c("mito_reads_rate","mtDNA_depth","nCount_mito","nFeature_mito","nCount_alleles","nFeature_alleles"), x_lab_rotate = TRUE,
    colors_use = myUmapcolors)
dev.off()




crc <- AlleleFreq(
  object = crc,
  variants = high.conf$variant,
  assay = "mito"
)
crc[["alleles"]]
DefaultAssay(crc) <- "alleles"


pdf("CM_feature_all_celltype.pdf",height=4,width=16)
DefaultAssay(obj) <- "alleles"
alleles.view <- c("6785G>A", "9581C>T", "6886C>T", "12892C>T")
FeaturePlot(
  object = obj,
  features = alleles.view,
  order = TRUE,
  cols = c("grey", "darkred"),
  ncol = 4
) 
FeaturePlot(
  object = combined,
  features = alleles.view,
  order = TRUE,
  cols = c("grey", "darkred"),
  ncol = 4
) 
dev.off()







CM<- subset(combined,idents="CM")
  variable.sites <- IdentifyVariants(CM, assay = "mito", refallele = AR3_C4_mito.data$refallele)
  high.conf <- subset(
  variable.sites, 
  subset = n_cells_conf_detected >= 5 &
    strand_correlation >= 0.5 &
    vmr > 0.01
  )
  high.conf<- high.conf[order(high.conf$mean,decreasing=T),]
high.conf[1:4,c(1,2,5)]

  crc <- AlleleFreq(
    object = CM,
    variants = high.conf$variant,
    assay = "mito"
  )
  DefaultAssay(crc) <- "alleles"

  data<- GetAssayData(crc)
  d_varition1=as.vector(data)
  d_varition1=d_varition1[which(d_varition1!=0)]
  data<- data.frame(variance=d_varition1)
  data$variance<- data$variance
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



