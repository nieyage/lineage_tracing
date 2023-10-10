install.packages("devtools")
devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")

library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)

plan("multicore", workers = 1)
options(future.globals.maxSize = 30 * 1024 ^ 3) # for 50 Gb RAM

combined<- readRDS("/Users/fraya/Documents/project/heart regeneration lineage tracing/01_AR3_scATAC/AR3_integrated_all_celltype_annotated.rds")

# Building trajectories with Monocle 3
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(patchwork)

# We can convert the Seurat object to a CellDataSet object using the as.cell_data_set() function from SeuratWrappers and build the trajectories using Monocle 3. 
DefaultAssay(combined)<-"ATAC"
#combined.cds <- importCDS(combined)
combined.cds <- as.cell_data_set(combined)
combined.cds <- cluster_cells(cds = combined.cds, reduction_method = "UMAP")
combined.cds <- learn_graph(combined.cds, use_partition = TRUE)
combined.cds <- order_cells(combined.cds)

plot_cells(cds = combined.cds,label_leaves=FALSE,color_cells_by = "pseudotime",show_trajectory_graph = TRUE)
plot_cells(cds = combined.cds,color_cells_by = "Annotation",show_trajectory_graph = TRUE)
plot_cells(cds = combined.cds,color_cells_by = "detail_anno",show_trajectory_graph = TRUE)
plot_cells(cds = combined.cds,color_cells_by = "passed_filters",show_trajectory_graph = TRUE)

combined <- AddMetaData(
  object = combined,
  metadata = combined.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "all_celltype_monocle"
)



DefaultAssay(combined)<-"ATAC"
Idents(combined)<- combined$Annotation
FB<- subset(combined,idents=c("FB","Epi","SMC","Pericyte"))
#combined.cds <- importCDS(combined)
combined.cds <- as.cell_data_set(FB)
combined.cds <- cluster_cells(cds = combined.cds, reduction_method = "UMAP")
combined.cds <- learn_graph(combined.cds, use_partition = TRUE)
combined.cds <- order_cells(combined.cds)

plot_cells(cds = combined.cds,label_leaves=FALSE,color_cells_by = "pseudotime",show_trajectory_graph = TRUE)
plot_cells(cds = combined.cds,color_cells_by = "Annotation",show_trajectory_graph = TRUE)
plot_cells(cds = combined.cds,color_cells_by = "detail_anno",show_trajectory_graph = TRUE)
plot_cells(cds = combined.cds,color_cells_by = "passed_filters",show_trajectory_graph = TRUE)
