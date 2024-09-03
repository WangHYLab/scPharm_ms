# Title    : data_process.R
# Author   : wanglab
# Time     : 2024.7

### 评估数据处理及scPharm
library(Seurat)
library(ggplot2)
library(patchwork)
library(CellID)
library(purrr)
library(readxl)
library(clusterProfiler)
library(tidyverse)
library(copykat)
library(tictoc)

Read_data = function(ID){
  ## read single cell data and create seuratobject
  message(paste("...load raw data:", ID, "..."))
  if (ID == "GSE136246"){
    path.list = list.dirs("../verify_data/GSE136246", full.names = F, recursive = F)
    data.list = lapply(path.list, function(x){
      files = list.files(paste("../verify_data/GSE136246", x, sep = "/"))
      if (length(files) > 1){
        data.x.list = lapply(files, function(y){
          counts = read.table(paste("../verify_data/GSE136246", x, y, sep = "/"), header = TRUE, row.names = 1, check.names = F)
          message(paste("...load sample:", y, "..."))
          counts = t(counts)
          counts = CreateSeuratObject(counts = counts,
                                      project = x,
                                      min.cells = 3,
                                      min.features = 200)
          counts = RenameCells(counts, add.cell.id = x)
          return(counts)
        })
        data = merge(data.x.list[[1]], data.x.list[2:length(data.x.list)])
      }else{
        counts = read.table(paste("../verify_data/GSE136246", x, files[1], sep = "/"),
                            header = TRUE, row.names = 1, check.names = F)
        message(paste("...load sample:", files[1], "..."))
        counts = CreateSeuratObject(counts = counts,
                                    project = x,
                                    min.cells = 3,
                                    min.features = 200)
        counts = RenameCells(counts, add.cell.id = x)
        data = counts
      }
      return(data)
    })
    names(data.list) = path.list
  }else if (ID == "GSE171145") {
    path.list = list.dirs(paste("../verify_data", ID, sep = "/"),
                          full.names = F, recursive = F)
    data.list = lapply(path.list, function(x){
      message(paste("...load sample:", x, "..."))
      counts = read.table(paste("../verify_data/GSE171145", x, "counts.tsv.gz", sep = "/"), header = TRUE, row.names = 1, check.names = F)
      data = CreateSeuratObject(counts = counts,
                                project = x,
                                min.cells = 3,
                                min.features = 200)
      data = RenameCells(data, add.cell.id = x)
      return(data)
    })
    names(data.list) = path.list
  }else{
    path.list = list.dirs(paste("../verify_data", ID, sep = "/"),
                          full.names = F, recursive = F)
    data.list = lapply(path.list, function(x){
      message(paste("...load sample:", x, "..."))
      data = Read10X(paste("../verify_data", ID, x, sep = "/"))
      data = CreateSeuratObject(counts = data,
                                project = x,
                                min.cells = 3,
                                min.features = 200)
      data = RenameCells(data, add.cell.id = x)
      return(data)
    })
    names(data.list) = path.list
  }
  return(data.list)
}

vlnplot = function(scRNA_list, ID, qc_status){
  scRNA = merge(scRNA_list[[1]],scRNA_list[2:length(scRNA_list)])
  plot.features = c("nFeature_RNA", "nCount_RNA")
  group = "orig.ident"
  # violin plot
  theme.set1 = theme(axis.title.x = element_blank())
  plots = list()
  for (i in seq_along(plot.features)){
    plots[[i]] = VlnPlot(scRNA, group.by = group, pt.size = 0,
                         features = plot.features[i])+
      theme.set1 +
      NoLegend()
  }
  violin = wrap_plots(plots = plots, nrow = 2)
  if (!dir.exists(paste("../verify_data/Rplot", ID, sep = "/"))){
    dir.create(paste("../verify_data/Rplot", ID, sep = "/"))
  }
  if (qc_status == "before"){
    ggsave(paste("../verify_data/Rplot", ID, "vlnplot_before_qc.pdf",sep = "/"),
           plot = violin, width = 6, height = 8)
  }else{
    ggsave(paste("../verify_data/Rplot", ID, "vlnplot_after_qc.pdf",sep = "/"),
           plot = violin, width = 6, height = 8)
  }
}

qc = function(scRNA_list, ID, maxGene, maxUMI, minGene = 200){
  out.list = lapply(seq(1,length(scRNA_list),1), function(i){
    scRNA = scRNA_list[[i]]
    out = subset(scRNA, subset = nCount_RNA < maxUMI[i] & nFeature_RNA < maxGene[i] & nFeature_RNA > minGene)
    })
  names(out.list) = names(scRNA_list)
  return(out.list)
}
id = 'GSE161529_HER2'
seurat_list = Read_data(id)
vlnplot(seurat_list, id, qc_status = "before")
#  some threshold for qc
{
  # GSE161529_HER2
  # maxGene = c(7000,3000,4000,5000,4000,5500),
  # maxUMI = c(40000,10000,10000,20000,10000,25000)
  # GSE161529_ER
  # maxGene = c(3000,1500,2500,4000,3500,3000,2500,2500,2500,2000,
  #             4000,3500,4000,2500),
  # maxUMI = c(7500,5000,7500,15000,15000,10000,10000,5000,7500,5000,
  #            7500,12000,12000,7500))
  # GSE171145
  # maxGene = c(2500,5500,6500,1500,6000,5000,5000,6000,2500),
  # maxUMI = c(10000,20000,35000,5000,25000,20000,20000,20000,10000)
}
seurat_list = qc(scRNA_list = seurat_list, id,
                 maxGene = c(7000,3000,4000,5000,4000,5500),
                 maxUMI = c(40000,10000,10000,20000,10000,25000))
vlnplot(seurat_list, id, qc_status = "after")

# GSE134839
{
  D0 <- read.table(
    "verify_data/D0.txt",
    row.names = 1,
    header = T)
  D1 <- read.table(
    "verify_data/D1.txt",
    row.names = 1,
    header = T)
  D2 <- read.table(
    "verify_data/D2.txt",
    row.names = 1,
    header = T)
  D4 <- read.table(
    "verify_data/D4.txt",
    row.names = 1,
    header = T)
  D9 <- read.table(
    "verify_data/D9.txt",
    row.names = 1,
    header = T)
  D11 <- read.table(
    "verify_data/D11.txt",
    row.names = 1,
    header = T)


  Seurat.d0 <- CreateSeuratObject(
    counts = D0,
    project = "D0",
    min.cells = 3,
    min.features = 200)
  Seurat.d1 <- CreateSeuratObject(
    counts = D1,
    project = "D1",
    min.cells = 3,
    min.features = 200)
  Seurat.d2 <- CreateSeuratObject(
    counts = D2,
    project = "D2",
    min.cells = 3,
    min.features = 200)

  Seurat.d4 <- CreateSeuratObject(
    counts = D4,
    project = "D4",
    min.cells = 3,
    min.features = 200)
  Seurat.d9 <- CreateSeuratObject(
    counts = D9,
    project = "D9",
    min.cells = 3,
    min.features = 200)

  Seurat.d11 <- CreateSeuratObject(
    counts = D11,
    project = "D11",
    min.cells = 3,
    min.features = 200)

  data<-merge(Seurat.d0,y=c(Seurat.d1,Seurat.d2,Seurat.d4,Seurat.d9,Seurat.d11))
  data[['percent.mt']] = PercentageFeatureSet(data, pattern = "^MT-")
  # violin plot
  VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  # filter
  data = subset(data, subset = nFeature_RNA > 800 & nFeature_RNA < 7500)
  data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 3000)
  # all.genes <- rownames(data)
  data <- ScaleData(data,features=VariableFeatures(data))
  data <- RunPCA(data, verbose = FALSE)
  ElbowPlot(data,ndims = 50)
  data <- RunUMAP(data, features = VariableFeatures(data))
  data <- FindNeighbors(data, dims = 1:10, verbose = FALSE)
  data <- FindClusters(data, verbose = FALSE)
  data@active.ident<-factor(data$orig.ident,levels = c("D0","D1","D2","D4","D9","D11"))
  data$orig.ident<-factor(data$orig.ident,levels = c("D0","D1","D2","D4","D9","D11"))
}

### single cell data of cancer cell lines
# data precess ------------------------------------------------------------
data = readRDS("../verify_data/SCP542/CCLE_scRNAseq_CPM.RDS")
proc.data = list()
for (type in unique(CCLE_metadata$cancer_type)) {
  if (sum(CCLE_metadata[CCLE_metadata$cancer_type == type,3]) >= 1000) {
    type.cellline = rownames(CCLE_metadata[CCLE_metadata$cancer_type == type,])
    type.data = data[type.cellline]
    type.object.list = map2(type.data, names(type.data), function(object,name) {
      object = suppressWarnings(CreateSeuratObject(counts = object,
                                                   project = name,
                                                   min.cells = 3,
                                                   min.features = 200))
      object = RenameCells(object, add.cell.id = name)
      return(object)
    })
    type.object = merge(type.object.list[[1]], type.object.list[2:length(type.object.list)])
    proc.data[[type]] = type.object
    message(paste0(type, "-done"))
  }else {
    warning("cell not enough")
  }
}
saveRDS(proc.data, file = "../verify_data/SCP542/cancer_type_sc.rds")
cancer_sc = readRDS("../verify_data/SCP542/cancer_type_sc.rds")
# remove some cancer type have not sufficient cell lines and cells
cancer_sc = cancer_sc[c(3:5, 7:14)]
cancer_tcgaid = c("UCEC", "LIHC", "LGG", "BLCA",
                  "HNSC", "OV", "COREAD", "STAD",
                  "KIRC", "ESCA", "PAAD")
names(cancer_tcgaid) = names(cancer_sc)

# Run scPharm
scPharmIdentify.result = scPharmIdentify(object, type = "tissue", cancer = "", cores = 8, assay = "RNA")
scPharmDr.result = scPharmDr(scPharmIdentify.result)
scPharmCombo.reuslt = scPharmCombo(scPharmIdentify.result, scPharmDr.result)
scPharmDse.result = scPharmDse(scPharmIdentify.result)