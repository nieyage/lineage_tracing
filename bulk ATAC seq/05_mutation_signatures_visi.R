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


smooth <- 1
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
# filter somatic mutation  --step6
# STEP 1: filter mtSNV
# $17>=1 && $18>=1 ( make sure both positive and negative strand were called 
# Minor  Allele frequency more than 0.01 
# remove G->T and C->A, "N",based on previous empirical evidence. 
# sequence depth >= 20 for each cell, the sequence depth criteria should be adjust according to the data. 

awk '$17>=1 && $18>=1' $file | sed '1,1d' |awk '$6/($6+$5)>=0.10' |awk '!($3=="G" && $6=="T")' |awk '!($3=="C" && $6=="A")' |awk '$3!="N"'  AR3-0.5-7-2.snv | awk ' $6+$5>=20' >>AR3-0.5-7-2.mtDNA.snv.filter

awk '{print $1"\t"$2"\t"$3"\t"$19}' AR3-0.5-7-2.mtDNA.snv.filter|sort -k2n -u >AR3-0.5-7-2.mtDNA.snv.filter.cut.uniq

awk '$17>=1 && $18>=1' $file | sed '1,1d' |awk '$6/($6+$5)>=0.10' |awk '!($3=="G" && $6=="T")' |awk '!($3=="C" && $6=="A")' |awk '$3!="N"'  AR3-1-5.snv | awk ' $6+$5>=20' >>AR3-1-5.mtDNA.snv.filter
awk '{print $1"\t"$2"\t"$3"\t"$19}' AR3-1-5.mtDNA.snv.filter|sort -k2n -u >AR3-1-5.mtDNA.snv.filter.cut.uniq

awk '$17>=1 && $18>=1' $file | sed '1,1d' |awk '$6/($6+$5)>=0.10' |awk '!($3=="G" && $6=="T")' |awk '!($3=="C" && $6=="A")' |awk '$3!="N"'  P7-15lab_L2.snv | awk ' $6+$5>=20' >>P7-15lab_L2.mtDNA.snv.filter
awk '{print $1"\t"$2"\t"$3"\t"$19}' P7-15lab_L2.mtDNA.snv.filter|sort -k2n -u >P7-15lab_L2.mtDNA.snv.filter.cut.uniq

awk '$17>=1 && $18>=1' $file | sed '1,1d' |awk '$6/($6+$5)>=0.10' |awk '!($3=="G" && $6=="T")' |awk '!($3=="C" && $6=="A")' |awk '$3!="N"'  P7-7Va_L2.snv | awk ' $6+$5>=20' >>P7-7Va_L2.mtDNA.snv.filter
awk '{print $1"\t"$2"\t"$3"\t"$19}' P7-7Va_L2.mtDNA.snv.filter|sort -k2n -u >P7-7Va_L2.mtDNA.snv.filter.cut.uniq

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

