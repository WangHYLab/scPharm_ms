# Title    : Supplementary Figure 4
# Author   : Wanglab
# Time     : 2024.7

library(Cairo)
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(ggpubr)
library(RcolorBrewer)

sample = c("AH0319","MH0001","MH0025","MH0029-7C","MH0029-9C","MH0032","MH0040","MH0042","MH0043-T","MH0064-T",
           "MH0151","MH0163","MH0173-T","PM0360"  )
names(sample) = sample

# SuppleFigure 4a ---------------------------------------------------------------

SuppleFigure4a <- function (sample) {
  data.1 = lapply(sample, function(sam) {
    # Load the single-cell pharmacology result object from the specified path
    rank = readRDS(paste0("./scPharm/result/", sam, "_scPharm_object_nmcs_50_nfs_200.rds"))
    return(rank)
  })
  # Loop through each element in the data.1 list to generate a PDF chart for each
  for (i in 1:14) {
    meta.data = data.1[[i]]@meta.data
    meta.data = meta.data[,c("cell.label", "scPharm_nes_1925_GDC0810")]
    # Filter cells where the label is "tumor"
    meta.data = meta.data[meta.data$cell.label == "tumor",]
    colnames(meta.data) = c("Group", "NES")
    # Create a PDF file using the CairoPDF function, specifying filename, width, and height
    CairoPDF(paste0("./Figure/supplefig4a_GDC0810_", names(data.1)[i],"_T_20231117.pdf"), width = 2, height = 1.4)
    # Set graphic parameters including outer margins, inner margins, and x-axis tick label positioning
    par(oma = c(1.2, 3, 0.2, 1), mar = c(1, 0, 1.4, 0), xpd = TRUE)
    plot(density(null$NES), col = "#999999", main = "",
        lwd = 1.5, lty = 2, cex.axis = 1.2, font.main = 1, tcl = -0.25, ylab = "", xlab = "",
        bty = "l", las = 1, ylim = c(0, 1.2))
    lines(density(meta.data$NES), col = "#FF0000", lwd=1.5)
    dev.off()
  }
}
SuppleFigure4a(sample)


# SuppleFigure 4b ---------------------------------------------------------------

SuppleFigure4b <- function (sample) {
  data.1 = lapply(sample, function(sam) {
    # Load the single-cell pharmacology result object from the specified path
    rank = readRDS(paste0("./scPharm/result/", sam, "_scPharm_object_nmcs_50_nfs_200.rds"))
    return(rank)
  })
  # Loop through each element in the data.1 list to generate a PDF chart for each
  for (i in 1:14) {
    meta.data = data.1[[i]]@meta.data
    meta.data = meta.data[,c("cell.label", "scPharm_nes_1925_GDC0810")]
    # Filter cells where the label is "adjacent"
    meta.data = meta.data[meta.data$cell.label == "adjacent",]
    colnames(meta.data) = c("Group", "NES")
    # Create a PDF file using the CairoPDF function, specifying filename, width, and height
    CairoPDF(paste0("./Figure/supplefig4b_GDC0810_", names(data.1)[i],"_A_20231117.pdf"), width = 2, height = 1.4)
    # Set graphic parameters including outer margins, inner margins, and x-axis tick label positioning
    par(oma = c(1.2, 3, 0.2, 1), mar = c(1, 0, 1.4, 0), xpd = TRUE)
    plot(density(null$NES), col = "#999999", main = "",
        lwd = 1.5, lty = 2, cex.axis = 1.2, font.main = 1, tcl = -0.25, ylab = "", xlab = "",
        bty = "l", las = 1, ylim = c(0, 1.2))
    lines(density(meta.data$NES), col = "#FF0000", lwd=1.5)
    dev.off()
  }
}
SuppleFigure4b(sample)

# SuppleFigure 4c ---------------------------------------------------------------

SuppleFigure4c <- function (sample) {
  # Load the drug prioritization of ER positive BRCA samples from the specified path
  data.2 = lapply(sample, function(sam) {
    rank = read.csv(paste("./scPharm/result/", sam, "_scPharm_drug_rank.csv", sep = ""))
    return(rank)
  })
  data.3 = data.frame(matrix(0, nrow = nrow(data.2[[1]]), ncol = length(data.2)))
  colnames(data.3) = names(data.2)
  for (sam in names(data.2)) {
    data = data.2[[sam]]
    data.3[data[data$DRUG_ID == 1925,]$Rank, sam] = 1
    data.3[data[data$DRUG_ID == 1200,]$Rank, sam] = 2
    data.3[data[data$DRUG_ID == 1816,]$Rank, sam] = 3
  }
  supfig2b = Heatmap(t(data.3),
                  # row_split = 1:14,
                  # row_gap = unit(1, "mm"),
                  # rect_gp = gpar(col = 'white', lwd = 0.1),
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
                  show_heatmap_legend = T,
                  name = "Drug",
                  heatmap_legend_param = list(
                    color_bar = "discrete",
                    at = 0:3,
                    labels = c("other", "GDC0810", "Fulvestrant_low", "Fulvestrant_high"),
                    legend_gp = gpar(fill = c("grey", "#E41A1C", "#377EB8", "#4DAF4A")),
                    labels_gp = gpar(col = "#000000", fontsize = 12),
                    title_gp = gpar(fontsize = 12)
                  )
                  # width = unit(14, "cm"), height = unit(1.5, "cm")
  )
  CairoPDF("./Figure/supfig4c.pdf", width = 9.8, height = 3.6)
  print(supfig2b)
  dev.off()
}
SuupleFigure4c(sample)


# SuppleFigure 4d ---------------------------------------------------------------

SuppleFigure4d <- function (sample) {
  # Load the drug prioritization of ER positive BRCA samples from the specified path
  data.2 = lapply(sample, function(sam) {
    rank = read.csv(paste("./scPharm/result/", sam, "_scPharm_drug_rank.csv", sep = ""))
    return(rank)
  })
  for (sam in names(data.2)) {
    data = data.2[[sam]][1:30,]
    data$Rank = 31-data$Rank
    data$group = "other"
    for (id in c(1925, 1200, 1816)) {
      if (id %in% data$DRUG_ID) {
        if (id == 1200) {
          data[data$DRUG_ID == id,]$group = paste(data[data$DRUG_ID == id,]$DRUG_NAME, "lowconc", sep = "_")
        }else if (id == 1816) {
          data[data$DRUG_ID == id,]$group = paste(data[data$DRUG_ID == id,]$DRUG_NAME, "highconc", sep = "_")
        }else {
          data[data$DRUG_ID == id,]$group = data[data$DRUG_ID == id,]$DRUG_NAME
        }
      }else {
        next
      }
    }
    col_value = c("GDC0810" = "#E41A1C", "Fulvestrant_lowconc"="#377EB8", "Fulvestrant_highconc"="#4DAF4A", "other"="grey")
    p = ggplot(data, aes(Rank, Score, fill = group))+
      geom_col(width = 0.8)+
      coord_flip()+
      scale_x_continuous(name = "Drug", breaks = seq(30,1,-1), labels = data$DRUG_NAME)+
      scale_fill_manual(values = col_value) +
      theme_bw() +
      theme(legend.position = "none",
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.title = element_text(size = 10,colour = "#000000"),
            axis.text = element_text(size = 10,colour = "#000000"),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_blank())+
      ggtitle(sam)
    # assign(paste("fig3c.", sam, sep = ""), p)
    CairoPDF(paste0("./Figure/supfig4d_", sam, "_1113.pdf"), width = 2.2, height = 5.8)
    print(p)
    dev.off()
  }
}
SuppleFigure4d(sample)