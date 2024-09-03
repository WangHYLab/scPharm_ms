# Title    : data_process.R
# Author   : wanglab
# Time     : 2024.7

# construct null distribution
library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)
library(future)

set.seed(0)
## normal lung cell
Donor.list = list.files("../null", pattern = "raw_gene_bc_matrices_h5.h5")
data.list = lapply(Donor.list, function(x){
  path = paste("../null/", x, sep = "")
  project = substr(x, start = 12, stop = 19)
  mat = Read10X_h5(path)
  mat = CreateSeuratObject(counts = mat,
                           project = project,
                           min.cells = 3,
                           min.features = 200)
  # add prefix to cell barcode
  mat = RenameCells(mat, add.cell.id = project)
  return(mat)
})
names(data.list) = paste("Donor0", 1:8, sep = "")

NL.list = list.dirs("../null", full.names = F, recursive = F)
data.list2 = lapply(NL.list, function(x){
  path = paste("../null/", x, sep = "")
  mat = Read10X(path)
  mat = CreateSeuratObject(counts = mat,
                           project = x,
                           min.cells = 3,
                           min.features = 200)
  # add prefix to cell barcode
  mat = RenameCells(mat, add.cell.id = x)
  return(mat)
})
names(data.list2) = NL.list

data.list = c(data.list, data.list2)
rm(data.list2)
# merge seurat objects
scRNA = merge(data.list[[1]], data.list[2:length(data.list)])
### vlnplot
theme.set1 = theme(axis.title.x = element_blank())
plot.features = c("nFeature_RNA", "nCount_RNA")
group = "orig.ident"
# plot vlnplot of "nfeature_RNA" and "nCount_RNA" before filter
plots = list()
for (i in seq_along(plot.features)){
  plots[[i]] = VlnPlot(scRNA, group.by = group, pt.size = 0,
                       features = plot.features[i])+
    theme.set1 +
    NoLegend()
}
violin = wrap_plots(plots = plots, nrow = 2)
dir.create("../null/QC")
ggsave("../null/QC/vlnplot_before_qc.pdf", plot = violin, width = 6, height = 8)

### set threshold for filtering
maxGene = 6000
maxUMI = 25000

# remove terrible object
scRNA = subset(scRNA, subset = orig.ident != "Donor_08")

# filter with threshold
scRNA = subset(scRNA, subset = nCount_RNA < maxUMI & nFeature_RNA < maxGene)
plots = list()
for (i in seq_along(plot.features)){
  plots[[i]] = VlnPlot(scRNA, group.by = group, pt.size = 0,
                       features = plot.features[i])+
    theme.set1 +
    NoLegend()
}
violin = wrap_plots(plots = plots, nrow = 2)
ggsave("../null/QC/vlnplot_after_qc.pdf", plot = violin, width = 6, height = 8)
# reduction and cluster
scRNA = NormalizeData(scRNA) %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData()
scRNA = RunPCA(scRNA, verbose = F)
ElbowPlot(scRNA, ndims = 50)
pc.num = 1:30
scRNA = scRNA %>% RunTSNE(dims=pc.num) %>% RunUMAP(dims=pc.num)
scRNA = FindNeighbors(scRNA, dims = pc.num) %>% FindClusters()
p = DimPlot(scRNA, group.by = "orig.ident", raster = F)
ggsave("../null/UAMP_Samples.pdf", p, width = 8, height = 6)
p = DimPlot(scRNA, group.by = "orig.ident", split.by = "orig.ident",
            ncol = 5, raster = F)
ggsave("../null/UMAP_Samples_Split.pdf", p, width = 22, height = 12)
saveRDS(scRNA, file = "../null/scRNA.rds")

# normal breast cell ------------------------------------------------------
Donor.list = list.files("../null/breast/")
data.list = lapply(Donor.list, function(x){
  path = paste("../null/breast/", x, sep = "")
  project = substr(x, start = 1, stop = 13)
  mat = Read10X_h5(path)
  mat = CreateSeuratObject(counts = mat,
                           project = project,
                           min.cells = 3,
                           min.features = 200)
  # add prefix to cell barcode
  mat = RenameCells(mat, add.cell.id = project)
  return(mat)
})
names(data.list) = paste("Donor0", c(1:5, "6_10_CMG", "6_10_HN", 11), sep = "")

# merge object
scRNA = merge(data.list[[1]], data.list[2:length(data.list)])
### vlnplot
theme.set1 = theme(axis.title.x = element_blank())
plot.features = c("nFeature_RNA", "nCount_RNA")
group = "orig.ident"
# vlnplot before filter
plots = list()
for (i in seq_along(plot.features)){
  plots[[i]] = VlnPlot(scRNA, group.by = group, pt.size = 0,
                       features = plot.features[i])+
    theme.set1 +
    NoLegend()
}
violin = wrap_plots(plots = plots, nrow = 3)
dir.create("../null/breast/QC")
ggsave("../null/breast/QC/vlnplot_before_qc.pdf", plot = violin, width = 9, height = 8)
### set threshold
maxGene = 7500
maxUMI = 25000
# filter and ploting
scRNA = subset(scRNA, subset = nCount_RNA < maxUMI & nFeature_RNA < maxGene & nFeature_RNA > 200)
plots = list()
for (i in seq_along(plot.features)){
  plots[[i]] = VlnPlot(scRNA, group.by = group, pt.size = 0,
                       features = plot.features[i])+
    theme.set1 +
    NoLegend()
}
violin = wrap_plots(plots = plots, nrow = 3)
ggsave("../null/breast/QC/vlnplot_after_qc.pdf", plot = violin, width = 9, height = 8)
# reduction and cluster
scRNA = NormalizeData(scRNA) %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData()
scRNA = RunPCA(scRNA, verbose = F)
ElbowPlot(scRNA, ndims = 50)
pc.num = 1:30
scRNA = scRNA %>% RunTSNE(dims=pc.num) %>% RunUMAP(dims=pc.num)
scRNA = FindNeighbors(scRNA, dims = pc.num) %>% FindClusters()
p = DimPlot(scRNA, group.by = "orig.ident", raster = F)
ggsave("../null/breast/UAMP_Samples.pdf", p, width = 8, height = 6)
p = DimPlot(scRNA, group.by = "orig.ident", split.by = "orig.ident",
            ncol = 3, raster = F)
ggsave("../null/breast/UMAP_Samples_Split.pdf", p, width = 22, height = 12)
saveRDS(scRNA, file = "../null/breast/scRNA.rds")

# normal skin cell --------------------------------------------------------
Donor.list = list.files("../null/skin/GSE182861/", pattern = "bc_matrix.h5")
data.list = lapply(Donor.list, function(x){
  path = paste("../null/skin/GSE182861/", x, sep = "")
  project = substr(x, start = 12, stop = 16)
  mat = Read10X_h5(path)
  mat = CreateSeuratObject(counts = mat[[1]],
                           project = project,
                           min.cells = 3,
                           min.features = 200)
  # add prefix to cell barcode
  mat = RenameCells(mat, add.cell.id = project)
  return(mat)
})
names(data.list) = paste("Donor0", 1:2, sep = "")

NL.list = list.dirs("../null/skin/GSE151177/", full.names = F, recursive = F)
data.list2 = lapply(NL.list, function(x){
  path = paste("../null/skin/GSE151177/", x, sep = "")
  mat = Read10X(path)
  mat = CreateSeuratObject(counts = mat,
                           project = x,
                           min.cells = 3,
                           min.features = 200)
  # add prefix to cell barcode
  mat = RenameCells(mat, add.cell.id = x)
  return(mat)
})
names(data.list2) = NL.list

data.list = data.list2
rm(data.list2)
scRNA = merge(data.list[[1]], data.list[2:length(data.list)])
theme.set1 = theme(axis.title.x = element_blank())
plot.features = c("nFeature_RNA", "nCount_RNA")
group = "orig.ident"
plots = list()
for (i in seq_along(plot.features)){
  plots[[i]] = VlnPlot(scRNA, group.by = group, pt.size = 0,
                       features = plot.features[i])+
    theme.set1 +
    NoLegend()
}
violin = wrap_plots(plots = plots, nrow = 2)
dir.create("../null/skin/QC")
ggsave("../null/skin/QC/vlnplot_before_qc.pdf", plot = violin, width = 6, height = 8)

# set threshold
maxGene = 4500
maxUMI = 35000
# filter and ploting
scRNA = subset(scRNA, subset = nCount_RNA < maxUMI & nFeature_RNA < maxGene & nFeature_RNA > 200)
plots = list()
for (i in seq_along(plot.features)){
  plots[[i]] = VlnPlot(scRNA, group.by = group, pt.size = 0,
                       features = plot.features[i])+
    theme.set1 +
    NoLegend()
}
violin = wrap_plots(plots = plots, nrow = 2)
ggsave("../null/skin/QC/vlnplot_after_qc.pdf", plot = violin, width = 6, height = 8)
# reduction and cluster
scRNA = NormalizeData(scRNA) %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData()
scRNA = RunPCA(scRNA, verbose = F)
ElbowPlot(scRNA, ndims = 50)
pc.num = 1:30
scRNA = scRNA %>% RunTSNE(dims=pc.num) %>% RunUMAP(dims=pc.num)
scRNA = FindNeighbors(scRNA, dims = pc.num) %>% FindClusters()
p = DimPlot(scRNA, group.by = "orig.ident", raster = F)
ggsave("../null/skin/UAMP_Samples.pdf", p, width = 8, height = 6)
p = DimPlot(scRNA, group.by = "orig.ident", split.by = "orig.ident",
            ncol = 3, raster = F)
ggsave("../null/skin/UMAP_Samples_Split.pdf", p, width = 22, height = 12)
saveRDS(scRNA, file = "../null/scRNA.rds")


# integrate objects of three healthy tissues
scRNA = readRDS("../null/scRNA.rds") # read merged objects of different tissues from corresponding paths
cellinfo = subset(scRNA@meta.data, select = "orig.ident")
scRNA = CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = cellinfo)
### noramlization with SCTransform
scRNA = SCTransform(scRNA, conserve.memory = TRUE)
### pca
scRNA = RunPCA(scRNA, npcs = 50, verbose = F)
### integrate
scRNA = RunHarmony(scRNA, group.by.vars="orig.ident", assay.use="SCT",
                   max.iter.harmony = 20)
pc.num = 1:30
scRNA = RunTSNE(scRNA, reduction="harmony", dims=pc.num) %>%
  RunUMAP(reduction="harmony",dims=pc.num)
p = DimPlot(scRNA, group.by = "orig.ident")
ggsave("../null/skin/UMAP_Samples_harmony.pdf", p, width = 8, height = 6)
p = DimPlot(scRNA, group.by = "orig.ident", split.by = "orig.ident",
            ncol = 3, raster = F)
ggsave("../null/skin/UMAP_Samples_harmony__Split.pdf", p, width = 22, height = 12)
saveRDS(scRNA, file = "../null/skin/scRNA_SCT_harmony.rds") # save the integrated object

# run scPharm and construct the null distribution

scPharmIdentify.result = scPharmIdentify(object, type = "cellline", cancer = 'LUAD', cores = 8, assay = "RNA") # cancer = 'LUAD' or 'BRCA' or "SKCM
null = list(lung_scPharmIdentify.result, breast_scPharmIdentify.result, skin_scPharmIdentify.result)
null = lapply(null, function(data) {
  data = data@meta.data
  data = data[,seq(6,595,2)]
})
null = lapply(null, function(data) {
  data = data %>% gather(key = "Group", value = "NES")
})
null = bind_rows(null)
null$Group = "null"
null = null[null$NES != 0,]