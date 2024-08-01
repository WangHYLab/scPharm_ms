# Title    : Supplementary Figure 6
# Author   : Wanglab
# Time     : 2024.7

library(tidyr)
library(ggplot2)
library(ggalluvial)
library(ggthemes)
library(Cairo)
library(Seurat)

sample = c("AH0319","MH0001","MH0025","MH0029-7C","MH0029-9C","MH0032","MH0040","MH0042","MH0043-T","MH0064-T",
           "MH0151","MH0163","MH0173-T","PM0360"  )
names(sample) = sample


# SuppleFigure 6b ---------------------------------------------------------------
SuppleFigure6b <- function (sample) {
  data.1 = lapply(sample, function(sam) {
    # Load the single-cell pharmacology result object from the specified path
    rank = readRDS(paste0("./scPharm/result/", sam, "_scPharm_object_nmcs_50_nfs_200.rds"))
    return(rank)
  })
  # Loop through each element in the data.1 list to generate a PDF chart for each
  for (sam in names(data.1)) {
  meta.data = data.1[[sam]]@meta.data
  meta.data = meta.data[meta.data$cell.label == "tumor", c("scPharm_label_1200_Fulvestrant", "scPharm_label_1816_Fulvestrant")]
  meta.data = meta.data[order(meta.data$scPharm_label_1200_Fulvestrant),]
  colnames(meta.data) = c("Fulvestrant_low","Fulvestrant_high")
  meta.data = meta.data %>% gather(key = "dosage", value = "cluster")
  meta.data$cell = c(rep(1:(nrow(meta.data)/2),2))
  CairoPDF(paste("./Figure/supfig3/alluvial_", sam,"_1112.pdf", sep = ""), width = 2.2, height = 4)
  alluvial = ggplot(meta.data,
                    aes(x = dosage, stratum = cluster, alluvium = cell, y =cell,
                        fill = cluster, label = NULL)) +
    geom_col(width = 0.3,color=NA) +
    scale_x_discrete(limits = c("Fulvestrant_low","Fulvestrant_high"), expand = c(.1, .1)) +
    geom_flow(alpha = 0.6, linewidth = 0) +
    geom_stratum(alpha = 1, width = .3, color = "white", linewidth = 0) +
    scale_fill_manual(values = c("#B5B5B5", "#FF0000", "#0000FF")) +
    # geom_text(stat = "stratum", size = 8) +
    theme_map()+
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 13),
          axis.text.x = element_text(hjust = 0.1, size = 12, angle = -60)) +
    ggtitle(sam)
  print(alluvial)
  dev.off()
  }
}
# statistic test for SuppleFigure 6b
c.switch = function(label) {
  # transform the label to 0,1,2
  if (label == "other") {
    return(0)
  }else if (label == "resistant") {
    return(1)
  }else {
    return(2)
  }
}
STforSF6b <- function (sample) {
  data.list = lapply(sample, function(sam) {
    # Load the single-cell pharmacology result object from the specified path
    rank = readRDS(paste0("./scPharm/result/", sam, "_scPharm_object_nmcs_50_nfs_200.rds"))
    return(rank)
  })
  signif = c()
  for (sam in names(data.list)) {
    meta.data = data.list[[sam]]@meta.data
    meta.data = meta.data[meta.data$cell.label == "tumor",c("scPharm_label_1200_Fulvestrant", "scPharm_label_1816_Fulvestrant")]
    meta.data = meta.data[order(meta.data$scPharm_label_1200_Fulvestrant),]
    meta.data = meta.data[(meta.data$scPharm_label_1200_Fulvestrant != "sensitive" | meta.data$scPharm_label_1816_Fulvestrant != "sensitive"),]
    sample = sample(1:nrow(meta.data), nrow(meta.data))
    meta.data$H0 = meta.data$scPharm_label_1200_Fulvestrant[sample]
    meta.data = meta.data[meta.data$scPharm_label_1200_Fulvestrant == 'other',]
    test.data = data.frame(matrix(0, 1853, 2))
    colnames(test.data) = c("H0","H1")
    for (i in 1:nrow(meta.data)) {
      test.data[i, "H1"] = c.switch(meta.data$scPharm_label_1816_Fulvestrant[i])
      test.data[i, "H0"] = c.switch(meta.data$H0[i])
    }
    pval = wilcox.test(test.data$H1, test.data$H0, alternative = "greater", paired = T)$p.value
    signif = c(signif, pval)
  }
  return(signif)
}

SuppleFigure6b(sample)