library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)

combined<- readRDS("./03_all_celltype/04_mgatk/all_cell_type_mgatk_alleles.rds")

# # of alleles and # of cells
# identify by MAESTER

high.conf <- subset(
    variable.sites,
    subset = n_cells_conf_detected >= 5 &
    strand_correlation >= 0.65 &
    vmr > 0.01
)

Allele_matrix<- as.matrix(combined@assays$alleles@counts)

# plot the Mean coverage per base, Number of mtDNA UMIs, Total number of UMIs
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

data <- combined@meta.data
library(ggpubr)
library(cowplot)
pdf("./07_identify_mutation/Mean_coverage_per_base_celltype_mtDNA_total_UMI.pdf",width=12,height=4)
p1 <- ggboxplot(df, x="celltype", y="coverage", bxp.errorbar = T,width=0.6, notch = F)+
theme(legend.position="none",axis.text.x = element_text(angle = 45, hjust = 1))+
ylab("Mean coverage per base")
p2 <- ggboxplot(data, x="Annotation", y="mitochondrial", bxp.errorbar = T,width=0.6, notch = F)+
theme(legend.position="none",axis.text.x = element_text(angle = 45, hjust = 1))+
ylab("Number of mtDNA reads")
p3 <- ggboxplot(data, x="Annotation", y="nCount_ATAC", bxp.errorbar = T,width=0.6, notch = F)+
theme(legend.position="none",axis.text.x = element_text(angle = 45, hjust = 1))+
ylab("Total number of ATAC reads")
p4 <- ggboxplot(data, x="Annotation", y="nFeature_alleles", bxp.errorbar = T,width=0.6, notch = F)+
theme(legend.position="none",axis.text.x = element_text(angle = 45, hjust = 1))+
ylab("# of alleles")
p1+p2+p3
p4
dev.off()

# THE VAF distribution
DefaultAssay(combined)<- "alleles"
  data<- GetAssayData(combined)
  d_varition1=as.vector(data)
  d_varition1=d_varition1[which(d_varition1!=0)]
  data<- data.frame(variance=d_varition1)
  data$variance<- data$variance;
pdf("./03_all_celltype/04_mgatk/THE VAF distribution.pdf",width=8,height=8)
  p=ggplot(data,aes(x=variance))+
  geom_histogram(
                 binwidth = 0.01,
                 fill="#69b3a2",##69b3a2
                 color="#e9ecef",##e9ecef
                 alpha=0.9,
                 breaks=seq(0,1,0.01))+ 
  theme_bw()+
  labs(x="Frequence",y="Count",title="Allcelltype")+
  theme(plot.title = element_text(size = 18, vjust = 0.5, hjust = 0.5, face = "bold"),
        axis.text.x = element_text(size=18,angle = 90,colour="black", face = "bold"),          
        axis.text.y = element_text(size=18,colour="black", face = "bold"),
        axis.title.x = element_text(size=14, face = "bold"), 
        axis.title.y = element_text(size=14, face = "bold"),
        panel.background = element_blank(),
        line = element_line(size=1),
        axis.line = element_line(size =1.0,colour = "black"),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank())+scale_x_continuous(limits = c(0,1),
                     breaks = seq(0,1,0.05))
    print(p)

   p=ggplot(data,aes(x=variance))+
  geom_histogram(
                 binwidth = 0.01,
                 fill="#69b3a2",##69b3a2
                 color="#e9ecef",##e9ecef
                 alpha=0.9,
                 breaks=seq(0,0.2,0.01))+ 
  theme_bw()+
  labs(x="Frequence",y="Count",title="Allcelltype")+
  theme(plot.title = element_text(size = 18, vjust = 0.5, hjust = 0.5, face = "bold"),
        axis.text.x = element_text(size=18,angle = 90,colour="black", face = "bold"),          
        axis.text.y = element_text(size=18,colour="black", face = "bold"),
        axis.title.x = element_text(size=14, face = "bold"), 
        axis.title.y = element_text(size=14, face = "bold"),
        panel.background = element_blank(),
        line = element_line(size=1),
        axis.line = element_line(size =1.0,colour = "black"),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank())+scale_x_continuous(limits = c(0,0.2),
                     breaks = seq(0,1,0.01))
    print(p)

dev.off()

# Calculate the number of cells that exceed VAF thresholds 0, 1, 5, 10, 20 (3 minutes)
Allele_matrix<- as.matrix(combined@assays$alleles@counts)*100
high.conf <- high.conf %>%
    mutate(n0 = apply(Allele_matrix, 1, function(x) sum(x == 0))) %>%
    mutate(n1 = apply(Allele_matrix, 1, function(x) sum(x > 1))) %>%
    mutate(n5 = apply(Allele_matrix, 1, function(x) sum(x > 5))) %>%
    mutate(n10 = apply(Allele_matrix, 1, function(x) sum(x > 10))) %>%
    mutate(n20 = apply(Allele_matrix, 1, function(x) sum(x > 20)))

pdf("./07_identify_mutation/n0_rate of mutation.pdf",width=5,height=3)
ggplot(high.conf, aes(x = n0/ncol(Allele_matrix))) +
  geom_histogram(stat="bin",fill="grey",color="grey",alpha=0.8)+ 
  labs( x = "n0", y = "# of mutations")+theme_bw()
dev.off()


high.conf_filter <- high.conf %>% filter(mean_coverage > 5,n0 > 0.9*ncol(Allele_matrix))


# all mutation: 3988;
# all cell: 25006
# Selection of informative variants

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Test different variant selection thresholds #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Specify variant selection thresholds to test. voi = variant of interest
conditions.tib <- tibble(min_clone_size = rep(1:50, 4),
                         min_vaf = rep(c(1, 5, 10, 20), each = 50),
                         vois = NA,
                         n_vois = NA,
                         cells = NA,
                         n_cells = NA,
                         transitions = NA)
vois.ls <- vector(mode = "list", length = nrow(conditions.tib))
cells.ls <- vector(mode = "list", length = nrow(conditions.tib))

# Quality threshold, including selection of variants that are absent from 90% of cells
Allele_matrix<- as.matrix(combined@assays$alleles@counts*100)[high.conf_filter$variant,]
library(stringr)
# Fill in conditions.tib
for (x in 1:nrow(conditions.tib)) {
    print(x)
    min_clone_size <- conditions.tib$min_clone_size[x]
    min_vaf <- conditions.tib$min_vaf[x]

    # Define variants of which the number of cells exceeding min_vaf is higher than min_clone_size
    count_greater_than_10 <- apply(Allele_matrix, 1, function(row) sum(row >= min_vaf))
    voi.ch<- names(which(count_greater_than_10 > min_clone_size))

    ## Exclude bad variants (as well as 1583_A>G based on downstream analysis)
    #voi.ch <- setdiff(voi.ch, c("1583_A>G", blocklist.var))

    # Which cells are positive for at least one of the variants?
    # af_subset.dm <- af.dm[voi.ch,]
    af_subset.dm<- Allele_matrix[voi.ch,]
    positive_cells <- colnames(af_subset.dm)[apply(af_subset.dm, 2, function(col) any(col >= 1))]
    
    # Add information to summary table
    conditions.tib[x,"n_vois"] <- length(voi.ch)
    conditions.tib[x,"n_cells"] <- length(positive_cells)
    # Transitions vs. transversions
    conditions.tib[x,"transitions"] <- mean( str_count(voi.ch, "G>A|A>G|C>T|T>C") )
    # Save variants and cells for these parameters
    vois.ls[[x]] <- voi.ch
    cells.ls[[x]] <- positive_cells
}
conditions.tib$vois <- vois.ls
conditions.tib$cells <- cells.ls
data<- conditions.tib %>%
    mutate(min_vaf = str_c(min_vaf, "%")) %>%
    mutate(min_vaf = factor(min_vaf, levels = c("1%", "5%", "10%", "20%"))) 

# Visualize
pdf("./07_identify_mutation/All_mutation_using_two_selection_criteria_Scatterplot.pdf", width = 5, height = 5)
data%>%
    ggplot(aes(x = min_clone_size, y = n_cells, color = min_vaf)) + #, size = n_vois
    geom_hline(yintercept = ncol(Allele_matrix)) +
    geom_point() +
    coord_cartesian(ylim = c(0, ncol(Allele_matrix))) +
    ylab("Number of cells with informative variant") +
    xlab("Minimum clone size") +
    theme_bw() +
    guides(color = guide_legend(override.aes = list(size = 3),
                                title = "Minimum VAF")) +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black"))
dev.off()

pdf("./07_identify_mutation//n_mutation_using_two_selection_criteria_Scatterplot.pdf", width = 5, height = 5)
data%>%
    ggplot(aes(x = min_clone_size, y = n_vois, color = min_vaf)) + #, size = n_vois
    geom_hline(yintercept = nrow(Allele_matrix)) +
    geom_point() +
    coord_cartesian(ylim = c(0, nrow(Allele_matrix))) +
    ylab("Number of informative variant") +
    xlab("Minimum clone size") +
    theme_bw() +
    guides(color = guide_legend(override.aes = list(size = 3),
                                title = "Minimum VAF")) +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black"))
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Deeper analysis of four variant selection thresholds #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Here we look into more detail into four parameter combinations:
# Select variants present in at least 5,10,20,50 cells with a VAF of >10% or >5%
conditions_subset.tib <- conditions.tib %>% filter(min_clone_size %in% c(5,10,20,50), min_vaf %in% c("5","10"))

# Pick one of the threshold settings
a <- 1
a <- 2
a <- 3
a <- 4 

# load mgatk output
AR3_C4_cov.data <- read.table(gzfile("~/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC/mgatk_for_filtered_cells/final/AR3_C4_scATAC.coverage.txt.gz"),sep=",")
AR3_C5_cov.data <- read.table(gzfile("~/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC/mgatk_for_filtered_cells/final/AR3_C5_scATAC.coverage.txt.gz"),sep=",")

# add prefix in barcode 
AR3_C4_cov.data$V2<- paste("AR3_C4_",AR3_C4_cov.data$V2,sep="")
AR3_C5_cov.data$V2<- paste("AR3_C5_",AR3_C5_cov.data$V2,sep="")
cov.data<- rbind(AR3_C4_cov.data,AR3_C5_cov.data)
cov.data$celltype<- combined$Annotation[match(cov.data$V2,rownames(combined@meta.data))]
library(reshape2)

cov.data_matrix <- dcast(cov.data, V1 ~ V2, value.var = "V3")
metadata.tib <- as_tibble(combined@meta.data, rownames = "cell")

myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")


library(dittoSeq)
library(ComplexHeatmap)
library(circlize)
# this is the one we ended up using for the paper
# Or run all four
for (a in 1:8) {
print(a)
voi.ch <- conditions_subset.tib$vois[[a]]

# List cell IDs that are positive for each of the vois --------------------------------------------
positive_cells.ls <- list()
for (v in voi.ch) {
    # Determine cells with an appreciable VAF
    current_cells.ch <- colnames(Allele_matrix)[Allele_matrix[v,]>1]
    # Save cell IDs for positive cells
    positive_cells.ls[[v]] <- current_cells.ch
}

# Make a tibble of cells marked by each voi
positive_cells.tib <- as_tibble(bind_rows(lapply(positive_cells.ls, function(x) data.frame(cell = x)), .id = "variant")[,2:1]) %>%
    mutate(variant = factor(variant, levels = voi.ch))
cov.mat <- cov.data_matrix[,-1] # location X cell, value is cov;

positive_cells.tib$cov <- apply(positive_cells.tib, 1, function(x) { cov.mat[as.numeric(gsub("[^0-9]+", " ", x[2])),x[1]] } )

# Add cell type columns.
#positive_cells.tib$cell_types <- apply(positive_cells.tib, 1, function(x) { cov.mat[as.numeric(cutf(x[2], d = "_")),x[1]] } )
positive_cells.tib <- positive_cells.tib %>%
    left_join(dplyr::select(as_tibble(combined@meta.data, rownames = "cell"), cell, Annotation,detail_anno))

# Unlike the full analysis of the clonal hematopoiesis sample (4.3_Variants_Of_Interest.R), I am not grouping
# similar variants together into groups (clones), but treating each variant separately.
# Heatmap -----------------------------------------------------------------------------------------

# Prepare matrix of variants of interest in cells that are positive for at least one
af_voi.mat <- Allele_matrix[voi.ch,]
af_subset.mat <- af_voi.mat[,apply(af_voi.mat, 2, function(x) sum(x > 1) > 0)]

# Customize column order. This is different from the strategy for K562 subclones.
plot_order.mat <- af_subset.mat
for (x in rev(voi.ch)) { plot_order.mat <- plot_order.mat[,order(-plot_order.mat[x,])] }

# Set upper VAF limit for visualization clarity
# if (a %in% c(1, 2)) { plot_order.mat[plot_order.mat > 20] <- 20 }

# Add annotation bar
anno.tib <- tibble(cell = colnames(plot_order.mat)) %>% left_join(metadata.tib, by = "cell")  %>% dplyr::select(Annotation,detail_anno)
Annotation<- levels(combined$Annotation);
col_ls<- myUmapcolors[1:11]
names(col_ls)<- Annotation;
combined$detail_anno<- as.factor(combined$detail_anno)
detail_anno<- levels(combined$detail_anno);
col_ls2<- myUmapcolors[1:length(detail_anno)]
names(col_ls2)<- detail_anno

ha <- HeatmapAnnotation(df = data.frame(anno.tib), col = list(Annotation = col_ls,detail_anno=col_ls2))

# Plot
hm <- Heatmap(plot_order.mat,
              col = colorRamp2(seq(0, round(max(plot_order.mat)), length.out = 9),
                               c("#FCFCFC","#FFEDB0","#FFDF5F","#FEC510","#FA8E24","#F14C2B","#DA2828","#BE2222","#A31D1D")),
              show_row_names = ifelse(nrow(plot_order.mat) < 50, T, F),
              show_column_names = F,
              cluster_columns = F,
              cluster_rows = F,
              row_names_gp = gpar(fontsize = 10),
              name = "AF",
              heatmap_legend_param = list(border = "#000000", grid_height = unit(10, "mm")),
              top_annotation = ha,
              border = T,
              width = unit(200, "mm"),
              height = unit(100, "mm"),
              use_raster = T,
              raster_quality = 5);
hm_50 <- Heatmap(plot_order.mat,
              col = colorRamp2(seq(0, 50, length.out = 9),
                               c("#FCFCFC","#FFEDB0","#FFDF5F","#FEC510","#FA8E24","#F14C2B","#DA2828","#BE2222","#A31D1D")),
              show_row_names = ifelse(nrow(plot_order.mat) < 50, T, F),
              show_column_names = F,
              cluster_columns = F,
              cluster_rows = F,
              row_names_gp = gpar(fontsize = 10),
              name = "AF",
              heatmap_legend_param = list(border = "#000000", grid_height = unit(10, "mm")),
              top_annotation = ha,
              border = T,
              width = unit(200, "mm"),
              height = unit(100, "mm"),
              use_raster = T,
              raster_quality = 5)
hm_20 <- Heatmap(plot_order.mat,
              col = colorRamp2(seq(0, 20, length.out = 9),
                               c("#FCFCFC","#FFEDB0","#FFDF5F","#FEC510","#FA8E24","#F14C2B","#DA2828","#BE2222","#A31D1D")),
              show_row_names = ifelse(nrow(plot_order.mat) < 50, T, F),
              show_column_names = F,
              cluster_columns = F,
              cluster_rows = F,
              row_names_gp = gpar(fontsize = 10),
              name = "AF",
              heatmap_legend_param = list(border = "#000000", grid_height = unit(10, "mm")),
              top_annotation = ha,
              border = T,
              width = unit(200, "mm"),
              height = unit(100, "mm"),
              use_raster = T,
              raster_quality = 5)
hm_10 <- Heatmap(plot_order.mat,
              col = colorRamp2(seq(0, 10, length.out = 9),
                               c("#FCFCFC","#FFEDB0","#FFDF5F","#FEC510","#FA8E24","#F14C2B","#DA2828","#BE2222","#A31D1D")),
              show_row_names = ifelse(nrow(plot_order.mat) < 50, T, F),
              show_column_names = F,
              cluster_columns = F,
              cluster_rows = F,
              row_names_gp = gpar(fontsize = 10),
              name = "AF",
              heatmap_legend_param = list(border = "#000000", grid_height = unit(10, "mm")),
              top_annotation = ha,
              border = T,
              width = unit(200, "mm"),
              height = unit(100, "mm"),
              use_raster = T,
              raster_quality = 5)
pdf(str_c("./07_identify_mutation/variant_selection_thresholds", a ,"_Heatmap.pdf"), width = 12, height = 6)
print(hm)
print(hm_50)
print(hm_20)
print(hm_10)
dev.off()

}

# Save for Supplemental Figure 11b
write.csv(dplyr::select(conditions_subset.tib[1:8], min_clone_size, min_vaf, n_vois, n_cells, transitions, ARI), file = "./07_identify_mutation/variant_selection_thresholds_results.csv")

# Select variants present in at least 20 cells with a VAF of >=10%
# Select variants present in at least 50 cells with a VAF of >=10%

# do find clone type by function FindClonotypes()
DefaultAssay(combined) <- "alleles"
voi.ch <- conditions_subset.tib$vois[[8]]
#crc<- subset(combined,cells=unique(unlist(positive_cells.ls)))
crc <- FindClonotypes(combined,features = voi.ch,resolution = 2,assay="alleles",metric = "cosine",algorithm = 3)
table(Idents(crc))

crc$Clonotypes<- Idents(crc)
pdf(str_c("./07_identify_mutation/variant_selection_thresholds_Findclonetype", a ,"_Heatmap.pdf"), width = 12, height = 6)
DoHeatmap(crc, features = voi.ch, slot = "data", disp.max = 0.1) +scale_fill_viridis_c()
DoHeatmap(crc, features = voi.ch, slot = "data", disp.max = 0.2) +scale_fill_viridis_c()
DoHeatmap(crc, features = voi.ch, slot = "data", disp.max = 0.5) +scale_fill_viridis_c()
DoHeatmap(crc, features = voi.ch, slot = "data", disp.max = 1) +scale_fill_viridis_c()
dev.off();




















