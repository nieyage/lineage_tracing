library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
library(future)
plan("multicore", workers = 12)
options(future.globals.maxSize = 70 * 1024 ^ 3) # for 50 Gb RAM
combined<- readRDS("./03_all_celltype/04_mgatk/all_cell_type_mgatk_alleles.rds")
Idents(combined)<- combined$Annotation
CM<- subset(combined,idents="CM")
summary(CM$passed_filters)

library(ggplot2)
library(dplyr)
library(hrbrthemes)
data<- data.frame(CM_ATAC_reads = CM$passed_filters)
pdf("./06_CM_ploidy/CM_DNA_content_distribution.pdf",width=15,height=6)
data %>%
  filter( CM_ATAC_reads< 40000) %>%
  ggplot( aes(x=CM_ATAC_reads)) +
    geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
    ggtitle("CM DNA content distribution") +
    theme_bw()
dev.off()



