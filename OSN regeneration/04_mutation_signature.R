library(data.table)
library(dplyr)
library(BuenColors)
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
OSN<- readRDS("./06_mgatk/all_cell_type_mgatk_add_alleles.rds")

# Simple reverse complement function
reverse_complement <- function(s){
  chartr("ATGC","TACG",s)
}
# Process 3 digit signature based on letters
ref_all <- fread("/md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked_add/outs/mgatk/final/chrM_refAllele.txt")
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
ref_all_long$strand <- ifelse(ref_all_long$ref %in% c("C","T"), "L", "H")

# Change to C/T as ref allele
ref_all_long$rc3 <- reverse_complement(ref_all_long$three)
ref_all_long$three_plot <- ifelse(ref_all_long$strand == "L", ref_all_long$three, ref_all_long$rc3)
ref_all_long$group_change <- ifelse(ref_all_long$strand == "L", ref_all_long$change, ref_all_long$change_rc)

# Annotate with called variants

mito.data <- ReadMGATK(dir = "/md01/nieyg/project/lineage_tracing/OSN_regeneration/00_data/plogF1met5d_masked_add/outs/mgatk/final/")
crc <- OSN
variable.sites <- IdentifyVariants(crc, assay = "mito", refallele = mito.data$refallele)
high.conf <- subset(
  variable.sites, 
  subset = n_cells_conf_detected >= 5 &
    strand_correlation >= 0.65 &
    vmr > 0.01
)

called_variants <- rownames(high.conf)
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
  scale_fill_manual(values= c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom",
    axis.text.x =element_text(angle=90,hjust=1,size=2)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_hline(yintercept = 1, linetype =2, color = "black") +
  labs(x = "Change in nucleotide", y = "Substitution Rate (Expected / Observed)")
cowplot::ggsave2(p1, file = "./06_mgatk/all_mito_signature.pdf", width = 4, height = 2.4)

