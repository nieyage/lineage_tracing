library(Signac)
library(Seurat)
library(ComplexHeatmap)
library(GenomeInfoDb)
library(ggplot2)
library(BuenColors) 
library(roll)
"%ni%" <- Negate("%in%")
P7_bulk_ATAC<- readRDS("./P7_bulk_ATAC.rds")
mut_se_signac <- readRDS("./P7_bulk_ATAC.signac.rds")

source("/data/R02/nieyg/project/lineage_tracing/bulk_ATAC/7_mgatk_2/variant_calling.R") 

# mgatk_se is the Summarized Experiment .rds file
# That is automatically produced from running
# The mgatk CLI python package 

mut_se <- call_mutations_mgatk(P7_bulk_ATAC)

misc_df <- data.frame(rowData(mut_se))
pdf("./P7_bulk_ATAC_VMR_strand.pdf",width=6,height=4)
ggplot(misc_df,aes(x = strand_correlation, y = log10(vmr))) +
  geom_point() +
  labs(color = "HQ", x = "Pearson correlation (strand)", y = "log VMR") +
  geom_vline(xintercept = 0.65, color = "black", linetype = 2) +
  geom_hline(yintercept = -2, color = "black", linetype = 2) +
  pretty_plot(fontsize = 12) + L_border() +
  theme(legend.position = "bottom")
dev.off()

# mtDNA QC 


# rollMean_coverage 
# Import the after
mcols(P7_bulk_ATAC)$refAllele <- toupper(as.character(mcols(P7_bulk_ATAC)$refAllele))
SE <- P7_bulk_ATAC

# Process the after
top2k <- tail(head(sort(colData(SE)$depth, decreasing = TRUE),2000),1)
cov_new <- rowMeans((assays(SE)[["coverage"]]))
mean(cov_new)

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

# Gental smooth of just 1 bp for plot aesthetics

smooth <- 1
df<- data.frame(
  pos = roll_mean(1:nrow(P7_bulk_ATAC@assays@data$coverage), smooth),
  P7_7Va = roll_mean(P7_bulk_ATAC@assays@data$coverage[,1], smooth),  
  P7_15lab = roll_mean(P7_bulk_ATAC@assays@data$coverage[,2], smooth)
)

mdf <- reshape2::melt(df, id.var = "pos")

# Visualize the rolled means
P1 <- ggplot(mdf, aes(x = pos, y = value, color = variable)) + 
  geom_line() +  expand_limits(y = c(-5, 4)) +
  pretty_plot(fontsize = 8)  + scale_color_manual(values = c("firebrick", "dodgerblue4")) +
  coord_polar(direction = 1) + labs(x = "", y = "log2 Coverage") + scale_y_log10() +
  theme(legend.position = "none") + xxtheme 

cowplot::ggsave2(P1, file = "./rollMean_coverage.pdf", width = 4, height = 4)

# A custom script was used for calling of bulk mtDNA mutations with base coverage >20 , variant coverage >10 ,and strand concordance >0.5 based on available documentation.  
mean_base_coverage<- data.frame(position=1:nrow(P7_bulk_ATAC@assays@data$coverage), mean_base_coverage=cov_new)
data_last<- rowData(mut_se)[,c(1,3,12)]
data_last$mean_base_coverage <- mean_base_coverage[match(data_last$position,mean_base_coverage$position),]$mean_base_coverage

misc_df <- data.frame(rowData(mut_se))
pdf("./P7_bulk_ATAC_VMR_strand.pdf",width=6,height=4)
# Make the standard variant calling plot
ggplot(misc_df, aes(x = mean_coverage, y = log10(mean), color = log10(mean) > -2 & mean_coverage > 20000)) +
  geom_point() + scale_color_manual(values = c("black", "firebrick")) +
  labs(color = "HQ", x = "mean_coverage", y = "log mutation frequency") +
  #pretty_plot(fontsize = 8) + L_border() +
  theme(legend.position = "bottom") + 
  geom_vline(xintercept = 20000, linetype = 2) +
  geom_hline(yintercept = -2, linetype = 2) + theme(legend.position = "none")
dev.off()