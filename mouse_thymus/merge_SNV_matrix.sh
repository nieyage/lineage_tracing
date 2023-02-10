cp -r /data/R04/chenbzh5/bottleneck2/mouse_thymus/mk_realign/SNV_filter_* .

#For site_infor.txt:
sed -i 's/chrM://g' ./*/site_infor.txt
sed -i 's/:/_/g' ./*/site_infor.txt
sed -i 's/ //g' ./*/site_infor.txt

#add barcode and total_rate
# add at first row
sed -i '1 i barcode' ./*/site_infor.txt
# add at last row
sed -i '$ a total_rate' ./*/site_infor.txt


# For fre_matrix.txt: 
cat ./*/fre_matrix.txt | awk '{print NF}'|uniq
cat ./*/site_infor.txt | wc -l 

# remove the first line 
sed -i '1d' ./*/fre_matrix.txt
# row > col

awk '{for(i=1;i<=NF;i=i+1){a[NR,i]=$i}}END{for(j=1;j<=NF;j++){str=a[1,j];for(i=2;i<=NR;i++){str=str " " a[i,j]}print str}}' fre_matrix.txt > trans_fre_matrix.txt
paste site_infor.txt  trans_fre_matrix.txt >  last_matrix.txt

#R 
library(data.table)
library(stringr)
library(tidyverse) 
df_empty <- data.frame()

library(plyr)

allDat <- lapply(list.files('./',pattern='SNV_filter'),function(f){
  print(f);
  tmp=read.table(file.path(f,"last_matrix.txt"), header = T);
  sample=gsub("SNV_filter_","",f);
  colnames(tmp)<- paste(sample,colnames(tmp),sep="-")
  colnames(tmp)[1]<-"barcode" 
  tmp_cell_mutation_rate<-as.numeric(tmp[nrow(tmp),2:ncol(tmp)]);
  as.data.frame(tmp_cell_mutation_rate) %>%
  #filter( price<300 ) %>%
  ggplot( aes(x=tmp_cell_mutation_rate)) +
    geom_histogram( binwidth=0.01,fill="#69b3a2", color="#e9ecef", alpha=0.9) +
    ggtitle("tmp_cell_mutation_rate") +
    theme(
      plot.title = element_text(size=15)
    )
  #df_empty=merge(df_empty, tmp,by = "barcode")  
  return(tmp)
})

raw.data <- do.call(rbind.fill, allDat)
merge<- raw.data[!duplicated(raw.data$barcode),]
merge<-merge[-which(merge$barcode=="total_rate"),]
rownames(merge)<-merge$barcode
merge<-merge[,-1]

SNV<-rownames(merge)
#trans NA to 0 
merge<-as.data.frame(lapply(merge,as.numeric));
merge[is.na(merge)] <- 0;
rownames(merge)<-SNV
write.csv(merge,"all_cell_SNV.csv");

merge<-read.csv("all_cell_SNV.csv",header=T,row.names=1)

everysite_mutation_rate<-rowSums(merge)
everycell_mutation_rate<-colSums(merge)



# Libraries
library(tidyverse)
library(hrbrthemes)
# plot
pdf("everycell_mutation_rate.pdf")
p <- as.data.frame(everycell_mutation_rate) %>%
  #filter( price<300 ) %>%
  ggplot( aes(x=everycell_mutation_rate)) +
    geom_histogram( binwidth=0.01,fill="#69b3a2", color="#e9ecef", alpha=0.9) +
    ggtitle("everycell_mutation_rate") +
    theme(
      plot.title = element_text(size=15)
    )
p
dev.off()
pdf("everysite_mutation_rate.pdf")
p <- as.data.frame(everysite_mutation_rate) %>%
  #filter( price<300 ) %>%
  ggplot( aes(x=everysite_mutation_rate)) +
    geom_histogram( binwidth=0.01,fill="#69b3a2", color="#e9ecef", alpha=0.9) +
    ggtitle("everysite_mutation_rate") +
    theme(
      plot.title = element_text(size=15)
    )
p
dev.off()

names(everycell_mutation_rate)<-colnames(merge)

data<-merge[,names(everycell_mutation_rate[everycell_mutation_rate > 0])]

names(everysite_mutation_rate)<- rownames(merge)

data<-data[names(everysite_mutation_rate[everysite_mutation_rate > 0.1]),]

#for heatmap 
mycolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF" )


library(pheatmap)
celltype<-sapply(strsplit(colnames(data),split="\\."),"[[",1)
annotation_col = data.frame(
  celltype
)
rownames(annotation_col) = factor(colnames(data))
ann_colors= list(
  celltype = c("Baso"=mycolors[1], "C2"  =mycolors[2], "C6"  =mycolors[3], "DC1" =mycolors[4], "DC2" =mycolors[5],       
"DC3" =mycolors[6], "NCL1"=mycolors[7], "NCL2"=mycolors[8], "T0"  =mycolors[9], "T1"  =mycolors[10],       
"T11" =mycolors[11], "T3"  =mycolors[12], "T4"  =mycolors[13], "TEC1"=mycolors[14], "TEC2"=mycolors[15],       
"Tprogenitor"=mycolors[16]))




#bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
#random select 500 cells to show 
cell<-sample(colnames(data),500)

data_heatmap<-data[,grep("Tprogenitor",colnames(data))]
data_heatmap<-data_heatmap[which(rowSums(data_heatmap)>0),]
bk <- c(seq(0,0.1,by=0.001))
pdf("SNV-heatmap.pdf",width=20,height=20)
pheatmap(data,cluster_cols = T,cluster_rows = T,
         color =colorRampPalette(colors = c("white","red"))(100),
         #legend_breaks=seq(-2,2,1),
         breaks=bk,
         #cellwidth = 10, cellheight = 10,
         annotation_col = annotation_col, 
         #annotation_row = annotation_row,
         #cutree_rows=2,
         annotation_colors = ann_colors,
         show_rownames=F,
         show_colnames=F)

  trans_dist <- 1-cosine(as.matrix(data))
  #col<-my47colors[1:length(unique(barcode_label_pheatmap$OR))]
  #names(col)<-unique(barcode_label_pheatmap$OR)
  #ann_colors= list(
  #  OR = col)
  trans_dist[is.na(trans_dist)]<-0
  pheatmap(trans_dist,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         color = colorRampPalette(c("#4B0C1A", "#C25B3F","#F1ECEC", "#3082BD","#1C214F"))(100),
           annotation_col = annotation_col,
           annotation_row = annotation_col,
           annotation_colors = ann_colors,
           annotation_legend = TRUE,
           show_rownames=F,
           show_colnames=F
  )
dev.off() 
color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
#make tree 

library(lsa)
dist <- as.dist(t(data_heatmap))
data.tree <- ape::as.phylo(x = hclust(d = dist))
library(ggtree);

pdf("cell-tree.pdf",width=20,height=20)
ggtree(data.tree) + geom_tiplab()+ geom_treescale()
dev.off()

pdf("SNV-heatmap-celltype.pdf",width=20,height=20)
for (i in unique(annotation_col$celltype)){
  print(i)
  data_heatmap<-data[,grep(paste(i,".",sep=""),colnames(data))]
  data_heatmap<-data_heatmap[which(rowSums(data_heatmap)>0),]
  bk <- c(seq(0,0.1,by=0.001))
  p <- pheatmap(data_heatmap,cluster_cols = T,cluster_rows = T,
            color =colorRampPalette(colors = c("white","red"))(100),
            #legend_breaks=seq(-2,2,1),
            breaks=bk,
            #cellwidth = 10, cellheight = 10,
            annotation_col = annotation_col, 
            #annotation_row = annotation_row,
            #cutree_rows=2,
            annotation_colors = ann_colors,
            show_rownames=F,
            show_colnames=F)
  print(p);
  #cosine dist 
  trans_dist <- 1-cosine(as.matrix(data_heatmap))
  trans_dist[is.na(trans_dist)]<-0
  p2<-pheatmap(trans_dist,
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           color = colorRampPalette(c("#4B0C1A", "#C25B3F","#F1ECEC", "#3082BD","#1C214F"))(100),
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           annotation_legend = TRUE,
           show_rownames=F,
           show_colnames=F
  )
  print(p2)
}
dev.off()

