# 2- scRNAseq

## 2-1 One-Sample Analysis

```R
sample2='Cattle_lung'
pbmc.data_2=readRDS("Cattle.rds")
data_2 <- pbmc.data_2@assays$RNA@counts

pbmc.data_2 <- CreateSeuratObject(counts = data_2, project = sample2) %>% 
  NormalizeData(verbose = FALSE) %>% 
  ScaleData(verbose = FALSE) %>% 
  FindVariableFeatures(verbose = FALSE)

resolution=0.5
pbmc.data_2 <- RunPCA(pbmc.data_2, verbose = FALSE, npcs = 30)
pbmc.data_2 <- FindNeighbors(pbmc.data_2,dims = 1:20)
pbmc.data_2 <- FindClusters(pbmc.data_2, resolution = resolution)
pbmc.data_2 <- RunUMAP(pbmc.data_2, dims = 1:20, check_duplicates = FALSE)

DimPlot(pbmc.data_2, reduction = "umap",label=T)
DimPlot(pbmc.data_2, reduction = "umap")

logFC_filter=1
pbmc.markers_2 <- FindAllMarkers(object = pbmc.data_2)
sig_marker_2=subset(pbmc.markers_2,pbmc.markers_2$p_val_adj<0.05)
filter_marker_2=subset(sig_marker_2,sig_marker_2$avg_log2FC > logFC_filter)
filter_marker_1=subset(sig_marker_2,sig_marker_2$avg_log2FC < -logFC_filter)
filter_marker=rbind(filter_marker_2,filter_marker_1)

saveRDS(pbmc.data_2, file = paste(sample2,resolution,"_CCA.rds",sep=''))
```

## 2-2 Multi-sample CCA integrated analysis

```R
pbmc.data_1=readRDS("cattle.rds")
pbmc.data_2=readRDS("yak.rds")
pbmc.data_1@project.name="Cattle"
pbmc.data_2@project.name="Yak"

ob.list <- list()
ob.list
ob.list[['Cattle']] <- pbmc.data_1
ob.list[['Yak']] <- pbmc.data_2

pancreas.anchors <- FindIntegrationAnchors(object.list = ob.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
pancreas.integrated@meta.data$orig.ident <- as.factor(pancreas.integrated@meta.data$orig.ident)
```

```R
resolution=0.1

pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
pancreas.integrated <- FindNeighbors(pancreas.integrated, reduction = "pca", dims = 1:20)
pancreas.integrated <- FindClusters(pancreas.integrated, resolution = resolution)

DimPlot(pancreas.integrated, reduction = "umap", group.by = "orig.ident",pt.size = 0.000000001)
DimPlot(pancreas.integrated, reduction = "umap", group.by = "seurat_clusters")

pbmc.markers_1 <- FindAllMarkers(object =pancreas.integrated)
sig_marker_2=subset(pbmc.markers_1,pbmc.markers_1$p_val_adj<0.05)
filter_marker_2=subset(sig_marker_2,sig_marker_2$avg_log2FC > 1)
filter_marker_1=subset(sig_marker_2,sig_marker_2$avg_log2FC < -1)
filter_marker=rbind(filter_marker_2,filter_marker_1)
write.table(filter_marker,paste0(resolution,"_CCA_marker.csv"))

saveRDS(pancreas.integrated, file = paste("CCA",resolution,".rds",sep=''))
```

## 2-3 Correlation analysis

```R
pbmc.data_1=readRDS("cattle.rds")
pbmc.data_2=readRDS(".Yak.rds")

AverageExp_pbmc.data_2<-as.data.frame(AverageExpression(pbmc.data_2,features=feature))
AverageExp_pbmc.data_1<-as.data.frame(AverageExpression(pbmc.data_1,features=feature))

coorda<-corr.test(AverageExp_pbmc.data_2,AverageExp_pbmc.data_1,method="pearson")
```