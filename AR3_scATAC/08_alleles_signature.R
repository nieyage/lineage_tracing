library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)

combined<- readRDS("./03_all_celltype/04_mgatk/all_cell_type_mgatk_alleles.rds")

# # of alleles and # of cells

high.conf <- subset(
    variable.sites,
    subset = n_cells_conf_detected >= 5 &
    strand_correlation >= 0.65 &
    vmr > 0.01
)

Allele_matrix<- as.matrix(combined@assays$alleles@counts)

Allele_matrix[Allele_matrix != 0] <- 1
data<- as.data.frame(table(colSums(Allele_matrix)))
colnames(data)<- c("alleles_number","cell_number")
data$alleles_number<- as.integer(data$alleles_number)
pdf("./03_all_celltype/04_mgatk/of alleles and of cells_for_alleles_matrix.pdf",width=8,height=3)
ggplot(data, aes(x = alleles_number, y = cell_number)) +
  geom_bar(stat = "identity", fill = "grey", color = "grey") +
  labs( x = "# of alleles", y = "# of cells")+theme_bw()+
  scale_x_continuous(breaks = seq(0, 600, by = 20))
dev.off()

# nFeature_alleles
data2 <- as.data.frame(table(combined@meta.data$nFeature_alleles))
colnames(data2)<- c("alleles_number","cell_number")
data2$alleles_number<- as.integer(data2$alleles_number)
pdf("./03_all_celltype/04_mgatk/of alleles and of cells_fornFeature_alleles.pdf",width=8,height=3)
ggplot(data2, aes(x = alleles_number, y = cell_number)) +
  geom_bar(stat = "identity", fill = "grey", color = "grey") +
  labs( x = "# of alleles", y = "# of cells")+theme_bw()+
  scale_x_continuous(breaks = seq(0, 600, by = 20))
dev.off()

# mitochondrial vs cell number
data2 <- as.data.frame(combined@meta.data$mitochondrial)
colnames(data2)<- c("mitochondrial")

pdf("./03_all_celltype/04_mgatk/mitochondrial of cells.pdf",width=8,height=3)
ggplot(data2, aes(x = mitochondrial)) +
  geom_histogram(stat="bin",fill="grey",color="grey",alpha=0.8)+ 
  labs( x = "mitochondrial", y = "# of cells")+theme_bw()
dev.off()
# mito_reads_rate vs cell number
combined$mito_reads_rate<- (combined$mitochondrial/combined$total)*100

data2 <- as.data.frame(combined@meta.data$mito_reads_rate)
colnames(data2)<- c("mito_reads_rate")

pdf("./03_all_celltype/04_mgatk/mito_reads_rate of cells.pdf",width=8,height=3)
ggplot(data2, aes(x = mito_reads_rate)) +
  geom_histogram(stat="bin",fill="grey",color="grey",alpha=0.8)+ 
  labs( x = "mito_reads_rate", y = "# of cells")+theme_bw()
dev.off()

# the alleles_counts high and low in each cell type 

data<- combined@meta.data
data$alleles_count_type<- ifelse(data$nCount_alleles>=3,"high","low")
cell_number <- as.data.frame(table(data$Annotation,data$alleles_count_type))#计算各组样本不同细胞群比例

Cellratio <- prop.table(table(data$Annotation,data$alleles_count_type), margin =1)#计算各组样本不同细胞群比例
Cellratio<- as.data.frame(Cellratio)
pdf("./03_all_celltype/04_mgatk/Cellratio of alleles_counts.pdf",width=8,height=3)
p1<- ggplot(Cellratio) + 
  geom_bar(aes(x =Var1, y= Freq, fill = Var2),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Cell type',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
p2<- ggplot(cell_number) + 
  geom_bar(aes(x =Var1, y= Freq, fill = Var2),stat = "identity",position="stack",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Cell type',y = 'cell number')+theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
p1+p2
dev.off()

# mutation feature VMR and # of alleles
variable.sites<- read.csv("./03_all_celltype/variable_sites_info.csv")
high.conf <- subset(
    variable.sites,
    subset = n_cells_conf_detected >= 5 &
    strand_correlation >= 0.65 &
    vmr > 0.01
)

pdf("./03_all_celltype/04_mgatk/VMR_alleles_number.pdf",width=8,height=3)
ggplot(high.conf, aes(x = vmr)) +
  geom_histogram(stat="bin",fill="grey",color="grey",alpha=0.8)+ 
  labs( x = "vmr", y = "# of alleles")+theme_bw()
dev.off()
# coverage vs VMR 

pdf("./03_all_celltype/04_mgatk/coverage vs VMR.pdf",width=8,height=3)
ggplot(high.conf, aes(x = mean_coverage, y = vmr)) +
  geom_bar(stat = "identity", fill = "grey", color = "grey") +
  labs( x = "mean_coverage", y = "vmr")+theme_bw()+
  scale_x_continuous(breaks = seq(0, 600, by = 20))
dev.off()




















