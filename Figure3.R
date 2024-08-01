# Title    : Figure 3
# Author   : Wanglab
# Time     : 2024.7
library(Cairo)
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(GSEAbase)
library(ggpubr)
library(RcolorBrewer)
library(harmony)


sample = c("AH0308", "MH0031", "MH0069", "MH0161", "MH0176", "PM0337")
names(sample) = sample

# figure 3a ---------------------------------------------------------------

Figure3a <- function (sample) {
  fig3.data.1 = lapply(sample, function(sam) {
    # Load the single-cell pharmacology result object from the specified path
    rank = readRDS(paste0("./scPharm/result/", sam, "_scPharm_object_nmcs_50_nfs_200.rds"))
    return(rank)
  })
  # Loop through each element in the fig3.data.1 list to generate a PDF chart for each
  for (i in 1:6) {
    meta.data = fig3.data.1[[i]]@meta.data
    meta.data = meta.data[,c("cell.label", "scPharm_nes_1032_Afatinib")]
    # Filter cells where the label is "tumor"
    meta.data = meta.data[meta.data$cell.label == "tumor",]
    colnames(meta.data) = c("Group", "NES")
    # Create a PDF file using the CairoPDF function, specifying filename, width, and height
    CairoPDF(paste0("./Figure/fig3a/fig3a_afatinib_", names(fig3.data.1)[i],"_20231117.pdf"), width = 2, height = 1.4)
    # Set graphic parameters including outer margins, inner margins, and x-axis tick label positioning
    par(oma = c(1.2, 3, 0.2, 1), mar = c(1, 0, 1.4, 0), xpd = TRUE)
    plot(density(null$NES), col = "#999999", main = "",
        lwd = 1.5, lty = 2, cex.axis = 1.2, font.main = 1, tcl = -0.25, ylab = "", xlab = "",
        bty = "l", las = 1, ylim = c(0, 1.2))
    lines(density(meta.data$NES), col = "#FF0000", lwd=1.5)
    dev.off()
  }
}
Figure3a(sample)


# figure 3b ---------------------------------------------------------------

Figure3b <- function (sample) {
  fig3.data.1 = lapply(sample, function(sam) {
    # Load the single-cell pharmacology result object from the specified path
    rank = readRDS(paste0("./scPharm/result/", sam, "_scPharm_object_nmcs_50_nfs_200.rds"))
    return(rank)
  })
  # Loop through each element in the fig3.data.1 list to generate a PDF chart for each
  for (i in 1:6) {
    meta.data = fig3.data.1[[i]]@meta.data
    meta.data = meta.data[,c("cell.label", "scPharm_nes_1032_Afatinib")]
    # Filter cells where the label is "adjacent"
    meta.data = meta.data[meta.data$cell.label == "adjacent",]
    colnames(meta.data) = c("Group", "NES")
    # Create a PDF file using the CairoPDF function, specifying filename, width, and height
    CairoPDF(paste0("./Figure/fig3b/fig3b_afatinib_", names(fig3.data.1)[i],"_20231117.pdf"), width = 2, height = 1.4)
    # Set graphic parameters including outer margins, inner margins, and x-axis tick label positioning
    par(oma = c(1.2, 3, 0.2, 1), mar = c(1, 0, 1.4, 0), xpd = TRUE)
    plot(density(null$NES), col = "#999999", main = "",
        lwd = 1.5, lty = 2, cex.axis = 1.2, font.main = 1, tcl = -0.25, ylab = "", xlab = "",
        bty = "l", las = 1, ylim = c(0, 1.2))
    lines(density(meta.data$NES), col = "#FF0000", lwd=1.5)
    dev.off()
  }
}
Figure3b(sample)


# figure 3c ---------------------------------------------------------------

Figure3c <- function () {
  # Load the pharmacology result sample MH0176
  MH0176 = readRDS("./scPharm/result/MH0176_scPharm_object_nmcs_50_nfs_200.rds")
  MH0176 = subset(MH0176, subset = cell.label == "tumor")
  umap = data.frame(MH0176@reductions$umap@cell.embeddings)
  umap$label = MH0176$scPharm_label_1032_Afatinib
  CairoPDF("./Figure/fig3c_dimplot.pdf", width = 4.5, height = 3.6)
  DimPlot(MH0176, group.by = "scPharm_label_1032_Afatinib",pt.size = 0.1,
          cols = c("grey", "red", "blue"))+
    theme(
      aspect.ratio = 1,
      axis.line = element_line(arrow = arrow(type = "closed")))+
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL)
  dev.off()
  # Extract the label information about afatinib from the MH0176 object
  fig3c_pie = table(MH0176$scPharm_label_1032_Afatinib)
  fig3c_pie = as.data.frame(fig3c_pie)
  colnames(fig3c_pie) = c("group","value")
  CairoPDF("./Figure/fig3c_pie.pdf", width = 4, height = 3.6)
  ggplot(fig3c_pie, aes(x="", y = value, fill = group))+#数据
    geom_bar(width = 1, stat = "identity",color="white")+#绘制柱状图
    coord_polar('y')+#变为极坐标
    theme_void()+#主题
    scale_fill_manual(values = c("grey", "red", "blue"))+#自定义颜色
    geom_text(aes(y = sum(value)-cumsum(value)+value/2,
                  label = value), size=6, color = "black")+#标签
    theme(legend.text = element_text(size = 14),
          legend.title =  element_blank(),
          legend.position = "right")
  dev.off()
}
Figure3c()

# figure 3d ---------------------------------------------------------------

Figure3d <- function (sample) {
  # Load the drug prioritization of HER2 positive BRCA samples from the specified path
  fig3.data.2 = lapply(sample, function(sam) {
    rank = read.csv(paste0("./scPharm/result/", sam, "_scPharm_drug_rank.csv"))
    return(rank)
  })
  fig3.data.4 = data.frame(matrix(0, nrow = nrow(fig3.data.2[[1]]), ncol = length(fig3.data.2)))
  colnames(fig3.data.4) = names(fig3.data.2)
  for (sam in names(fig3.data.2)) {
    data = fig3.data.2[[sam]]
    fig3.data.4[data[data$DRUG_ID == 1032,]$Rank, sam] = 1
    fig3.data.4[data[data$DRUG_ID == 1549,]$Rank, sam] = 2
    fig3.data.4[data[data$DRUG_ID == 1558,]$Rank, sam] = 3
  }
  fig3d = Heatmap(t(fig3.data.4),
                  # row_split = 1:4,
                  # rect_gp = gpar(col = 'white', lwd = 0.1),
                  # row_split = group,
                  # row_gap = unit(1, "mm"),
                  col = colorRamp2(c(0,1,2,3), c("grey", "#E41A1C", "#377EB8", "#4DAF4A")),
                  # border = "black",
                  column_title = "Rank of all drugs(Dr)",
                  column_title_gp = gpar(fontsize = 14),
                  cluster_rows = F,
                  cluster_columns = F,
                  show_row_names = T,
                  show_column_names = F,
                  row_names_gp = gpar(fontsize = 12),
                  row_names_side = "left",
                  row_title = "",
                  show_heatmap_legend = T,
                  name = "Drug",
                  heatmap_legend_param = list(
                    color_bar = "discrete",
                    at = 0:3,
                    labels = c("other","Afatinib", "Sapitinib", "Lapatinib"),
                    legend_gp = gpar(fill = c("grey","#E41A1C", "#377EB8", "#4DAF4A")),
                    labels_gp = gpar(col = "#000000", fontsize = 12),
                    title_gp = gpar(fontsize = 12)
                  )
                  # width = unit(15, "cm"), height = unit(1.5, "cm")
  )
  CairoPDF("./Figure/fig3d_20231211.pdf", width = 10, height = 1.6)
  print(fig3d)
  dev.off()
}
Figure3d(sample)

# figure 3e ---------------------------------------------------------------

Figure3e_h <- function (sample) {
  # Load the single-cell pharmacology result object from the specified path
  fig3.data5 = lapply(sample, function(sam) {
    rank = readRDS(paste0("./scPharm/result/", sam, "_scPharm_object_nmcs_50_nfs_200.rds"))
    return(rank)
  })
  # Load the KEGG pathway data
  kegg.gmt = getGmt("../kegg_pathway.gmt")
  # Select tumor cell from object
  tumor.list = lapply(fig3.data5, function(sample) {
    data = subset(sample, subset = cell.label == "tumor")
    return(data)
  })
  # integrate the data
  merge.data = merge(tumor.list[[1]], tumor.list[2:6])
  cellinfo = subset(merge.data@meta.data, select = c("orig.ident",
                                                     "scPharm_label_1032_Afatinib",
                                                     "scPharm_label_1549_Sapitinib",
                                                     "scPharm_label_1558_Lapatinib"))
  merge.data = CreateSeuratObject(merge.data@assays$RNA@counts, meta.data = cellinfo)
  # SCT normalization
  merge.data = SCTransform(merge.data, conserve.memory = TRUE)
  # pca
  merge.data = RunPCA(merge.data, npcs = 50, verbose = F)
  # integration
  merge.data = RunHarmony(merge.data, group.by.vars="orig.ident", assay.use="SCT",
                    max.iter.harmony = 20)
  # ElbowPlot(data)
  merge.data = RunTSNE(merge.data, reduction="harmony", dims=1:15) %>%
    RunUMAP(reduction="harmony",dims=1:15)
  merge.data = merge.data %>% FindNeighbors() %>%  FindClusters(resolution = 0.3)
  # summary signaling pathways to plot
  pathway = c("ErbB signaling pathway",
              "MAPK signaling pathway",
              "Ras signaling pathway",
              "mTOR signaling pathway")
  # plot
  for (p in pathway) {
    # calculate pathway activity
    expr = GetAssayData(merge.data, slot = "data")
    expr = expr[which(rownames(expr) %in% kegg.gmt[[p]]@geneIds),]
    expr = as.data.frame(t(expr))
    expr$mean_expr = rowMeans(expr)
    col.name = str_replace_all(p, " ", "_")
    merge.data = AddMetaData(merge.data, expr[,'mean_expr'], col.name = "target_pathway_activity")
    CairoPDF(paste("./Figure/fig3/", col.name, "_20231117.pdf", sep = ""), width = 4.6, height = 3)
    viol_data = merge.data@meta.data[,c("orig.ident", "target_pathway_activity")]
    viol_data$orig.ident = factor(viol_data$orig.ident, levels = c("AH0308","MH0031","MH0069","PM0337","MH0176","MH0161"))
    p <- ggplot(viol_data, aes(x = orig.ident,
                               y = as.numeric(target_pathway_activity),
                               fill = orig.ident))+
      geom_violin(linewidth=0.4)+
      geom_boxplot(width=.2, col="black", fill="white", outlier.color = "#999999", outlier.size = 0.8,)+
      scale_fill_manual(values = c("#D4C2AD","#BFC7D9","#66796B",
                                   "#CDECEF","#9CBCB7","#CDD5C6"))+
      stat_compare_means(aes(group = orig.ident), label = "p.signif", method.args = list(alternative = "greater"),
                         label.y = max(viol_data[,2])+0.01,
                         size=4.5, ref.group = "MH0161")+
      theme_classic()+
      theme(plot.title = element_text(hjust = 0.5, size = 13),
            legend.text = element_text(size = 12),
            legend.title =  element_text(size = 12),
            axis.text.y = element_text(size = 12,colour = "#000000"),
            axis.text.x = element_text(size = 12,colour = "#000000", angle = 45, hjust = 1),
            text = element_text(size = 12),
            axis.title.y = element_text(size = 13, angle = 90),
            plot.margin = unit(c(0.2,0.3,0.2,0.3),'cm'),
            legend.position = "none"
      )+
      ylim(0, max(viol_data[,2])+0.02)+
      xlab(NULL)+
      ylab("Pathway activity")+
      ggtitle(p)
    print(p)
    dev.off()
  }
}
Figure3e_h(sample)