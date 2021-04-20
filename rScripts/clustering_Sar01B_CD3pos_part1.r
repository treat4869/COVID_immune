setwd("/projectbig/jupyternotebook/AG_Panzer/COVID/")

library(Seurat)
library(dplyr)
#library(tidyverse)
library(here)
#library(readxl)
#library(future)
#library(Matrix)
library(Signac)
library(ggplot2)
#library(sctransform)
library("xlsx")
library(future)
#plan("multiprocess", workers = 20)
options(future.globals.maxSize = 20 * 1000 * 1024^2)#20GB
options(repr.matrix.max.rows=600, repr.matrix.max.cols=200)
RhpcBLASctl::blas_set_num_threads(1)

samples_ID="Sar01B_CD3pos"

seu<-readRDS(paste0("Prefiltering_objects/",samples_ID,"_filtered.rds"))
seu

seu <- NormalizeData(object = seu,
                     normalization.method = "LogNormalize",
                     scale.factor = 10000)

seu<-FindVariableFeatures(seu, 
        selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seu), 10)
top10

#DefaultAssay(object = seu) <- "RNA"
seu <- ScaleData(seu,features = rownames(seu))

seu <- RunPCA(object = seu, features = VariableFeatures(object = seu),
              verbose = F)

DimPlot(object = seu, pt.size = 0.1,reduction = 'pca')

# Examine and visualize PCA results 
print(seu[["pca"]], dims = 1:5, nfeatures = 5)



#seu <- RunTSNE(object = seu, dims = 1:30)
seu <- RunUMAP(object = seu, dims = 1:30)

#DimPlot(object = seu,
#        reduction = 'tsne',label = F, 
#        pt.size = 0.1)+ theme(aspect.ratio=1)
DimPlot(object = seu,
        reduction = 'umap',label = F, 
        pt.size = 0.1)+ theme(aspect.ratio=1)

seu <- FindNeighbors(object = seu, dims = 1:30)

seu <- FindClusters(object = seu, resolution = 0.1)
table(Idents(seu))
#DimPlot(object = seu, reduction = 'tsne',label = TRUE, 
#        pt.size = 0.1,label.size = 6, repel = T) + theme(aspect.ratio=1)
DimPlot(object = seu, reduction = 'umap',label = TRUE, 
        pt.size = 0.1,label.size = 6, repel = T) + theme(aspect.ratio=1)

seu <- FindClusters(object = seu, resolution = 0.2)
table(Idents(seu))
#DimPlot(object = seu, reduction = 'tsne',label = TRUE, 
#        pt.size = 0.1,label.size = 6, repel = T) + theme(aspect.ratio=1)
DimPlot(object = seu, reduction = 'umap',label = TRUE, 
        pt.size = 0.1,label.size = 6, repel = T) + theme(aspect.ratio=1)

seu <- FindClusters(object = seu, resolution = 0.3)
table(Idents(seu))
#DimPlot(object = seu, reduction = 'tsne',label = TRUE, 
#        pt.size = 0.1,label.size = 6, repel = T) + theme(aspect.ratio=1)
DimPlot(object = seu, reduction = 'umap',label = TRUE, 
        pt.size = 0.1,label.size = 6, repel = T) + theme(aspect.ratio=1)

seu <- FindClusters(object = seu, resolution = 0.4)
table(Idents(seu))
#DimPlot(object = seu, reduction = 'tsne',label = TRUE, 
#        pt.size = 0.1,label.size = 6, repel = T) + theme(aspect.ratio=1)
DimPlot(object = seu, reduction = 'umap',label = TRUE, 
        pt.size = 0.1,label.size = 6, repel = T) + theme(aspect.ratio=1)

seu <- FindClusters(object = seu, resolution = 0.6)
table(Idents(seu))
#DimPlot(object = seu, reduction = 'tsne',label = TRUE, 
#        pt.size = 0.1,label.size = 6, repel = T) + theme(aspect.ratio=1)
DimPlot(object = seu, reduction = 'umap',label = TRUE, 
        pt.size = 0.1,label.size = 6, repel = T) + theme(aspect.ratio=1)



######################

seu <- FindClusters(object = seu, resolution = 0.2)
table(Idents(seu))
#DimPlot(object = seu, reduction = 'tsne',label = TRUE, 
#        pt.size = 0.1,label.size = 6, repel = T) + theme(aspect.ratio=1)
DimPlot(object = seu, reduction = 'umap',label = TRUE, 
        pt.size = 0.1,label.size = 6, repel = T) + theme(aspect.ratio=1)

Idents(seu) <- seu@meta.data$seurat_clusters
ident_use = seu@meta.data$seurat_clusters

# nCells in total 
paste0(ncol(GetAssayData(object = seu, slot = "scale.data"))," cells in total")
# show percentage of nCells from each cluster 
cluster_nCell<-as.data.frame(table(ident_use))
colnames(cluster_nCell)<- c("cluster","nCells")
cluster_nCell$percent <-  round((cluster_nCell$nCells / sum(cluster_nCell$nCells)*100),2)
cluster_nCell

marker_qc<-c("nFeature_RNA","nCount_RNA","frac.mito","frac.ribo","CD3_count")
for (n in marker_qc){
    print(VlnPlot(object = seu, features = n,
                  #group.by = "sample",
                  x.lab.rot=T, 
                  #size.x.use = 5,
                  pt.size = 0.01
                       )+NoLegend()
    )
        print(VlnPlot(object = seu, features = n,
                  #group.by = "sample",
                  x.lab.rot=T, 
                  #size.x.use = 5,
                  pt.size = 0
                       )+NoLegend()
    )
    }

# find markers for every cluster compared to all remaining cells, report only the positive ones
seu.markers <- FindAllMarkers(object = seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
#logfc.threshold=0.25 (default) instead of old version thresh.use=0.25
head(seu.markers)
dim(seu.markers)

write.xlsx(seu.markers, file = paste0("Prefiltering_objects/",samples_ID,"_cluster_marker.xlsx"), 
           sheetName = "cluster_marker",col.names = TRUE, row.names = TRUE, append = FALSE)



top3 <- seu.markers %>% group_by(cluster) %>% top_n(3, avg_logFC)
top5 <- seu.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
top10 <- seu.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

top10

plt1<-DoHeatmap(object = seu, features = top5$gene,
         size = 3, angle=0)+NoLegend()
plt1

plt1<-DotPlot(seu, features = unique(top5$gene),
              dot.scale = 4
              #scale.by = "size"
             ) + coord_flip()+
theme(#strip.background = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 10,angle = 45,hjust = 1),
    axis.text.y = element_text(size = 10),
     legend.position = "right",
     #legend.spacing = unit(0, "mm"),
     legend.direction = "vertical",
        legend.text = element_text(size=5),
        legend.key.width = unit(2, "mm"),
        legend.key.height = unit(2, "mm"),
        legend.box.spacing = unit(1, "mm"),
        legend.margin = margin(2),
        legend.title = element_text(size = 7,angle = 90)
)
plt1

# Search for known marker genes in seu
leukos <- c("PTPRC") 
Tcells <- c("CD3G","CD3D","CD3E")
CD4<-c("CD4")
CD8 <- c("CD8A", "PRF1", "GZMB")
Memory<-c("CCR7","SELL","KLF2","CD69","S1PR1")
NK <- c("KLRC1")
Th1 <- c("TBX21", "IFNG", "LTA") 
Th2 <- c("GATA3", "IL4", "IL5", "IL13") 
Th17 <- c("RORC","IL17A","IL17F", "IL23R","CCR6") 
Tregs <- c("FOXP3", "IL2RA", "CTLA4")
Tr1<-c("IL10","ITGA2","LAG3","HAVCR2")#,"Ahr","Irf4","Prdm1","Maf")
Tgd<- c("TRDC","TCRG")
Prolif<-c("STMN1","MKI67")
#Tgd_Scart<-c("SCART2")
#T_APC<-c("C1QA","C1QB")
#others<-c("ZCCHC12","KLRC3","VCAM1")

known_markers<-list(
leukos,
Tcells,
CD4,
CD8,
Memory,
NK,
Th1,
Th2,
Th17,
Tregs,
Tr1,
Tgd,
Prolif
#Tgd_Scart,
#T_APC,
#others   
)
known_markers

marker_gene_list<-known_markers
length(unlist(marker_gene_list))
marker_gene_list_expressed <- intersect(unlist(marker_gene_list), rownames(GetAssayData(seu)))
length(marker_gene_list_expressed)
setdiff(unlist(marker_gene_list),marker_gene_list_expressed)



genes_to_plot<- marker_gene_list_expressed
length(genes_to_plot)

plt1<-DotPlot(seu, features = genes_to_plot,
              dot.scale = 4
              #scale.by = "size"
             ) + coord_flip()+
theme(#strip.background = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 10,angle = 45,hjust = 1),
    axis.text.y = element_text(size = 10),
     legend.position = "right",
     #legend.spacing = unit(0, "mm"),
     legend.direction = "vertical",
        legend.text = element_text(size=5),
        legend.key.width = unit(2, "mm"),
        legend.key.height = unit(2, "mm"),
        legend.box.spacing = unit(1, "mm"),
        legend.margin = margin(2),
        legend.title = element_text(size = 7,angle = 90)
)
plt1




####################tbc!!!!!!!!!!!!!!!!!!!!!!!!!!!

#saveRDS(seu, paste0("Prefiltering_objects/",samples_ID,"_clustering_part1.rds"))


