library(Seurat)
data <-Read10X(data.dir="./filtered_feature_bc_matrix")
hm <- CreateSeuratObject(counts =data, project = "hm")
rm(data)
library(SC3)
library(scater)
library(SingleCellExperiment)
sce <- SingleCellExperiment(
    assays = list(
        counts = as.matrix(hm[["RNA"]]@counts),
        logcounts = log2(as.matrix(hm[["RNA"]]@counts) + 1)
    )
)
rm(hm)
rowData(sce)$feature_symbol =rownames(sce)
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
sce <- runPCA(sce)
sce <- sc3(sce, ks = 11, biology = TRUE)
saveRDS(sce, file = "SC3.rds")
sce <- runUMAP(sce)
plotUMAP(
  sce,
  colour_by = "sc3_13_clusters"
)