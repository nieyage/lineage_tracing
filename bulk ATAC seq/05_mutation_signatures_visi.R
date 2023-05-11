library(Signac)
library(Seurat)
library(ComplexHeatmap)
library(GenomeInfoDb)
library(ggplot2)
library(BuenColors) 
library(roll)

# % mtDNA reads
sample<- c("AR3-0.5-7-2","AR3-1-5","P7-15lab_L2","P7-7Va_L2")
mtDNA_rate<- c(57.50,63.10,55.22,43.47)

data<- data.frame(sample,mtDNA_rate)

pdf("./mtDNA_reads.pdf",width=4,height=2)
p<-ggplot(data = data, aes_string(x = "sample", y = "mtDNA_rate", 
        fill = "sample")) +  xlab(" ") + ylab("% mtDNA reads") + 
        scale_fill_manual(values = c("#896493" ,"#896493" ,"#896493","#896493")) + 
        geom_bar( stat = "identity",alpha=0.6, width = 0.6) +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 25,vjust = 0.5,hjust = 0.5));
p
#add gene number in plot 
p+geom_text(aes(label = mtDNA_rate), size = 3, hjust = 0.5, vjust = 3) 
dev.off();

# caculate average depth 
awk 'BEGIN{sum=0;num=0}''{sum+=$4;num++}''END{print sum/num}'  count  
sample<- c("AR3-0.5-7-2","AR3-1-5","P7-15lab_L2","P7-7Va_L2")
mt_aver_depth<- c(5963.74,6460.89,7181.27,7347.2)
data<- data.frame(sample,mt_aver_depth)
pdf("./mt_aver_depth.pdf",width=4,height=3)
p<-ggplot(data = data, aes_string(x = "sample", y = "mt_aver_depth", 
        fill = "sample")) +  xlab(" ") + ylab("mtDNA_aver_depth") + 
        scale_fill_manual(values = c("#75B7B0" ,"#75B7B0" ,"#75B7B0","#75B7B0")) + 
        geom_bar( stat = "identity",alpha=0.6, width = 0.6) +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 25,vjust = 0.5,hjust = 0.5));
p
#add gene number in plot 
p+geom_text(aes(label = mt_aver_depth), size = 3, hjust = 0.5, vjust = 3) 
dev.off();

# the coverage plot 
data1<- read.table("./7_mutaion/AR3-0.5-7-2.count")
data2<- read.table("./7_mutaion/AR3-1-5.count")
data3<- read.table("./7_mutaion/P7-15lab_L2.count")
data4<- read.table("./7_mutaion/P7-7Va_L2.count")

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


smooth <- 5
df<- data.frame(
  pos = roll_mean(data1$V2, smooth),
  AR3_0.5_7_2 = roll_mean(data1$V4, smooth),  
  AR3_1_5 = roll_mean(data2$V4, smooth),
  P7_15lab_L2=roll_mean(data3$V4, smooth),
  P7_7Va_L2=roll_mean(data4$V4, smooth)
)

mdf <- reshape2::melt(df, id.var = "pos")

# Visualize the rolled means
P1 <- ggplot(mdf, aes(x = pos, y = value, color = variable)) + 
  geom_line() +  expand_limits(y = c(-5, 4)) +
  pretty_plot(fontsize = 8)  + scale_color_manual(values = c("#795885", "#7B96AF","#CAC7A6","#FAEF5E")) +
  coord_polar(direction = 1) + labs(x = "", y = "log2 Coverage") + scale_y_log10() +
  theme(legend.position = "none") + xxtheme 
  cowplot::ggsave2(P1, file = "./rollMean_coverage.pdf", width = 4, height = 4)

P1 <- ggplot(mdf, aes(x = pos, y = value, color = variable)) + 
  geom_line() +  expand_limits(y = c(-5, 4)) +
  pretty_plot(fontsize = 8)  + scale_color_manual(values = c("#795885", "#7B96AF","#CAC7A6","#FAEF5E")) +
  coord_polar(direction = 1) + labs(x = "", y = "log2 Coverage") + scale_y_log10() +
  theme(legend.position = "none")
cowplot::ggsave2(P1, file = "./rollMean_coverage2.pdf", width = 4, height = 4)
# the mutation freq barplot 

# the coverage plot 
data1<- read.table("./7_mutaion/AR3-0.5-7-2.mtDNA.snv.filter")
data2<- read.table("./7_mutaion/AR3-1-5.mtDNA.snv.filter")
data3<- read.table("./7_mutaion/P7-15lab_L2.mtDNA.snv.filter")
data4<- read.table("./7_mutaion/P7-7Va_L2.mtDNA.snv.filter")
data1$sample<- "AR3-0.5-7-2"
data2$sample<- "AR3-1-5"
data3$sample<- "P7-15lab"
data4$sample<- "P7-7Va"
#data<- rbind(data1,data2,data3,data4)
data1$V7<- as.numeric(sub("%","",data1$V7))
data2$V7<- as.numeric(sub("%","",data2$V7))
data3$V7<- as.numeric(sub("%","",data3$V7))
data4$V7<- as.numeric(sub("%","",data4$V7))
pdf("mutation_freq.pdf",width=5,height=3)
ggplot(data=data1,aes(x=V7))+geom_histogram(fill="#795885",color="#795885",alpha=0.8)+ ylab("mtDNA mutations")+theme_bw()+xlab("% VAF")+ggtitle("AR3-0.5-7-2")
ggplot(data=data2,aes(x=V7))+geom_histogram(fill="#7B96AF",color="#7B96AF",alpha=0.8)+ ylab("mtDNA mutations")+theme_bw()+xlab("% VAF")+ggtitle("AR3-1-5")
ggplot(data=data3,aes(x=V7))+geom_histogram(fill="#CAC7A6",color="#CAC7A6",alpha=0.8)+ ylab("mtDNA mutations")+theme_bw()+xlab("% VAF")+ggtitle("P7-15lab")
ggplot(data=data4,aes(x=V7))+geom_histogram(fill="#FAEF5E",color="#FAEF5E",alpha=0.8)+ ylab("mtDNA mutations")+theme_bw()+xlab("% VAF")+ggtitle("P7-7Va")
dev.off()

library(data.table)
library(dplyr)
library(BuenColors)

# Simple reverse complement function
reverse_complement <- function(s){
  chartr("ATGC","TACG",s)
}

# Process 3 digit signature based on letters

 count<- fread("./7_mutaion/P7-7Va_L2.count")
ref_all<- count[,2:3]
colnames(ref_all) <- c("pos", "ref")
ref_all$ref <- toupper(ref_all$ref)
l <- as.character(ref_all$ref)

# Gs happen to be at the first and last position
ref_all$three <- paste0(c("G", l[-length(l)]), l, c(l[-1], "G"))

# Remove Ns
ref_all <- ref_all[!grepl("N", ref_all$three),]

# Make every possible mutation
ref_all_long <- rbind(ref_all,ref_all, ref_all,ref_all)
ref_all_long$alt <- rep(c("A", "C", "G", "T"), each = dim(ref_all)[1])
ref_all_long <- ref_all_long[ref_all_long$ref != ref_all_long$alt,]

# add some meta data
ref_all_long$variant <- paste0(as.character(ref_all_long$pos), ref_all_long$ref, ">", ref_all_long$alt)
ref_all_long$change <- paste0(ref_all_long$ref, ref_all_long$alt)
ref_all_long$change_rc <- reverse_complement(paste0(ref_all_long$ref, ref_all_long$alt))

# A/G rich strand is "heavy" -- https://en.wikipedia.org/wiki/Heavy_strand
table(ref_all$ref) # so the reference strand is light (more C/T)
# Change to C/T as ref allele
ref_all_long$strand <- ifelse(ref_all_long$ref %in% c("C","T"), "L", "H")


ref_all_long$rc3 <- reverse_complement(ref_all_long$three)
ref_all_long$three_plot <- ifelse(ref_all_long$strand == "L", ref_all_long$three, ref_all_long$rc3)
ref_all_long$group_change <- ifelse(ref_all_long$strand == "L", ref_all_long$change, ref_all_long$change_rc)

# Annotate with called variants
#called_variants <- rowData(readRDS("./filtered_mitoSE_CD34-500.rds"))$variant
called_variants<- fread("./7_mutaion/AR3-0.5-7-2.mtDNA.snv.filter.cut.uniq") 
colnames(called_variants)<-c("chrom","Position","Ref","VarAllele")
called_variants$change<- paste(called_variants$Ref,called_variants$VarAllele,sep=">")
called_variants<- paste(called_variants$Position,called_variants$change,sep="")
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
  theme(axis.title.x=element_blank(), axis.text.x =element_blank())+
  scale_fill_manual(values= c("firebrick", "dodgerblue3")) +
  theme(legend.position = "top") +
  scale_y_continuous(expand = c(0,0)) +
  geom_hline(yintercept = 1, linetype =2, color = "black") +
  labs(x = "Change in nucleotide", y = "Substitution Rate (Expected / Observed)")+ggtitle("AR3-1-5.test")
cowplot::ggsave2(p1, file = "./AR3-1-5.test_all_mito_signature.pdf", width = 4, height = 2)



BiocManager::install("MutationalPatterns")
library(MutationalPatterns)
library(BSgenome)
head(available.genomes())
ref_genome <- "BSgenome.Mmusculus.UCSC.mm10"
library(ref_genome, character.only = TRUE)
vcf_files <- list.files(path="./7_mutaion/",pattern = ".vcf", full.names = TRUE)
sample_names <- c(  "AR3-0.5-7-2", "AR3-1-5", "P7-15lab_L2", "P7-7Va_L2")
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



