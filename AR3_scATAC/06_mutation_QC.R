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
  coord_polar(direction = 1) + labs(x = "", y = "Coverage") + scale_y_log10() +
  theme(legend.position = "none") + xxtheme 


  cowplot::ggsave2(P1, file = "./03_all_celltype/04_mgatk/rollMean_coverage.pdf", width = 4, height = 4)
P1 <- ggplot(df, aes(x = pos, y = coverage, color = celltype)) + 
  geom_line() +  expand_limits(y = c(-5, 4)) +
  pretty_plot(fontsize = 8)  + scale_color_manual(values = myUmapcolors[1:11]) +
  coord_polar(direction = 1) + labs(x = "", y = "Coverage") + scale_y_log10() +
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

high.conf[,c(1,2,5)]
# global VAF distribution of mutation 
# the mutation freq barplot 
data<- variable.sites
data$variance<- data$variance*100;
#for(i in 1:nrow(data)){
#	if(data$variance[i] < 0.1){data$VAF[i]<- "0.0-0.1%"}
#	if(data$variance[i] >0.1  && data$variance[i] < 0.5){data$VAF[i]<- "0.1-0.5%"}
#	if(data$variance[i] >0.5  && data$variance[i] < 1){data$VAF[i]<- "0.5-1%"}
#	if(data$variance[i] >1    && data$variance[i] < 90){data$VAF[i]<- "1-90%"}
#	if(data$variance[i] >90   && data$variance[i] < 100){data$VAF[i]<- "90-100%"}
#}

pdf("./03_all_celltype/04_mgatk/All_celltype_mutation_freq.pdf",width=5,height=3)
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
	data<- variable.sites
    data$variance<- data$variance*100;
#for(i in 1:nrow(data)){
#	if(data$variance[i] < 0.1){data$VAF[i]<- "0.0-0.1%"}
#	if(data$variance[i] >0.1  && data$variance[i] < 0.5){data$VAF[i]<- "0.1-0.5%"}
#	if(data$variance[i] >0.5  && data$variance[i] < 1){data$VAF[i]<- "0.5-1%"}
#	if(data$variance[i] >1    && data$variance[i] < 90){data$VAF[i]<- "1-90%"}
#	if(data$variance[i] >90   && data$variance[i] < 100){data$VAF[i]<- "90-100%"}
#}
    p<- ggplot(data=data,aes(x=variance))+
    geom_histogram(fill="#795885",color="#795885",alpha=0.8)+ 
    ylab("mtDNA mutations")+theme_bw()+xlab("% VAF")+ggtitle(j)
    print(p)
}
dev.off()

Idents(combined)<- combined$detail_anno
combined$mito_reads_rate<- (combined$mitochondrial/combined$total)*100
library(scCustomize)

pdf("./03_all_celltype/04_mgatk/all_cell_type_mtDNA_depth.pdf",width=16,height=10)
#VlnPlot(combined, c("mito_reads_rate","mtDNA_depth","nCount_mito","nFeature_mito","nCount_alleles","nFeature_alleles"), pt.size = 0,ncol=1) + scale_y_log10()
Stacked_VlnPlot(seurat_object = combined, features = c("mito_reads_rate","mtDNA_depth","nCount_mito","nFeature_mito","nCount_alleles","nFeature_alleles"), x_lab_rotate = TRUE,
    colors_use = myUmapcolors)
dev.off()

# visulization signature of mutation 
library(MutationalPatterns)
library(BSgenome)
head(available.genomes())
ref_genome <- "BSgenome.Mmusculus.UCSC.mm10"
library(ref_genome, character.only = TRUE)


# Annotate with called variants
called_variants <- rowData(combined)$variant
ref_all_long$called <- ref_all_long$variant %in% called_variants

# Compute changes in expected/observed
total <- dim(ref_all_long)[1]
total_called <- sum(ref_all_long$called)
prop_df <- ref_all_long %>% group_by(three_plot, group_change, strand) %>%
  summarize(observed_prop_called = sum(called)/total_called, expected_prop = n()/total, n = n()) %>%
  mutate(fc_called = observed_prop_called/expected_prop)

prop_df$change_plot <- paste0(prop_df$group_change, "_", prop_df$three_plot)

# Visualize
p1 <- ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_called)) +
  geom_bar(stat = "identity", position = "dodge") + pretty_plot(fontsize = 8) + L_border() + 
  #theme(axis.title.x=element_blank(), axis.text.x =element_blank())+
  scale_fill_manual(values= c("firebrick", "dodgerblue3")) +
  theme(legend.position = "top",axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.5)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_hline(yintercept = 1, linetype =2, color = "black") +
  labs(x = "Change in nucleotide", y = "Substitution Rate (Expected / Observed)")
cowplot::ggsave2(p1, file = "./all_mito_signature.pdf", width = 4, height = 2.4)




mut_mat <- read.table("your_matrix.txt", sep = "\t", header = T)
makeGRangesFromDataFrame
GRangesList

grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome,group="circular", type = c("snv"))
time <- c(rep("AR3", 2), rep("AR7", 2))
#snv_grl <- get_mut_type(grl, type = "snv")
muts <- mutations_from_vcf(grl[[1]])
head(muts, 12)
library("gridExtra")
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
head(mut_mat)
pdf("./All_Sample_96_profile.pdf",width=12,height=6)
plot_96_profile(mut_mat[, 1:4])
dev.off()


#Strand bias profile
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
genes_hg19 <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
strand <- mut_strand(grl[[1]], genes_hg19)
head(strand, 10)
mut_mat_s <- mut_matrix_stranded(grl, ref_genome, genes_hg19)
mut_mat_s[1:5, 1:5]

pdf("./All_Sample_192_profile.pdf",width=12,height=6)
plot_192_profile(mut_mat_s[, 1:4])
dev.off()

# crc <- AlleleFreq(
#   object = crc,
#   variants = high.conf$variant,
#   assay = "mito"
# )
# crc[["alleles"]]
# DefaultAssay(crc) <- "alleles"
# crc <- FindClonotypes(crc)
# pdf("./03_all_celltype/04_mgatk/VariantPlot_global_heatmap.pdf",width=10,height=10)
# DoHeatmap(crc, features = VariableFeatures(crc), slot = "data", disp.max = 0.1) +  scale_fill_viridis_c()
# dev.off()

# pdf("./03_all_celltype/04_mgatk/Variant_Featureplot.pdf")
# alleles.view <- c("12889G>A", "16147C>T", "9728C>T", "9804G>A")
# FeaturePlot(
#   object = crc,
#   features = alleles.view,
#   order = TRUE,
#   cols = c("grey", "darkred"),
#   ncol = 4
# ) & NoLegend()
# dev.off()

#

#
#crc$log10_mtDNA_depth <- log10(crc$mtDNA_depth)
#FeaturePlot(crc, features = c("pct_reads_in_DNase"))
#FeaturePlot(crc, features = c("log10_mtDNA_depth"))
#source("variant_calling.R") 
#
## mgatk_se is the Summarized Experiment .rds file
## That is automatically produced from running
## The mgatk CLI python package 
#
#mut_se <- call_mutations_mgatk(mgatk_se)



