# Title    : Figure 3
# Author   : Wanglab
# Time     : 2024.7

library(Seurat)
library(ggplot2)
library(ggpubr)
library(Cairo)
library(dplyr)
library(ggsci)

# Load the single-cell object processed by scPharm
gse = readRDS("../scPharm/result/GSE134839_scPharm_nmcs_50_nfs_200.rds")
# Load a previously generated null distribution for subsequent using
null = readRDS("null.rds")
meta.data = gse@meta.data

# figure 4b, 4d, 4e ---------------------------------------------------------------

Figure4bde = function(gse, null, meta.data){
  # Figure4b, 4d, 4e
  # 4b
  for (time in unique(meta.data$orig.ident)[2:6]) {
    nes = meta.data[meta.data$orig.ident == time,]
  CairoPDF(paste0("./Figure/GSE134839/", time, "_density.pdf"), width = 3, height = 1.6)
  par(oma = c(1.2, 3, 0.2, 1), mar = c(1, 0, 1.4, 0), xpd = TRUE)
  plot(density(null$NES), col = "#999999", main = "",
       lwd = 1.5, lty = 2, cex.axis = 1.2, font.main = 1, tcl = -0.25, ylab = "", xlab = "",
       bty = "l", las = 1, ylim = c(0, 1), xlim = c(-4,4))
  lines(density(nes$scPharm_nes_1168_Erlotinib), col = "#FF0000", lwd=1.5)
  dev.off()
  }
  # 4d and 4e
  {
    # calculate the ratio of each cell type
    cellratio = prop.table(table(Idents(gse), gse$scPharm_label_1168_Erlotinib), margin = 1)
    cellratio = as.data.frame(cellratio)
    colnames(cellratio) = c('orig.ident', 'celltype', 'Freq')
    cellratio = cellratio[cellratio$orig.ident != "D0",]
    CairoPDF("./Figure/GSE134839/ratio_bar.pdf", width = 4, height = 2.4)
    p6 <- ggplot(cellratio)+
      geom_bar(aes(x = orig.ident, y = Freq, fill =celltype), stat = 'identity', width = 0.6, linewidth = 0.5)+
      theme_classic()+
      labs(x = "Time", y = "Ratio", fill = "")+
      scale_fill_manual(values = c('#D3D3D3', '#FF0000', "#0000FF"))+
      theme(axis.text = element_text(size = 12, colour = "#000000"),
            axis.title = element_text(size = 14),
            axis.text.y = element_text(hjust = -4),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14),
            # plot.margin = unit(c(0.1,0.1,0.1,0.1),'cm'),
      )
    print(p6)
    dev.off()
    CairoPDF("./Figure/GSE134839/resistant_line.pdf", width = 4 ,height = 2.4)
    p7 <- ggplot(cellratio[cellratio$celltype == "resistant",])+
      geom_line(aes(x=orig.ident,y=Freq, group = celltype, color = celltype)) +
      labs(x = "Time", y = "Ratio", fill = "")+
      scale_color_manual(values = c('#FF0000',"#0000FF"))+
      theme_classic()+
      theme(axis.text = element_text(size = 12, colour = "#000000"),
            axis.title = element_text(size = 14),
            axis.text.y = element_text(hjust = -4),
            legend.text = element_text(size = 12),
            # plot.margin = unit(c(0.1,0.1,0.1,0.1),'cm'),
            legend.title = element_blank())
    print(p7)
    dev.off()
  }
}
Figure4bde(gse, null, meta.data)


# figure 4c -----------------------------------------------------------------------

Figure4c = function(gse){
  for (time in c("D1","D2","D4","D9","D11")) {
    CairoPDF(paste0("./Figure/GSE134839/", time, "_dimplot.pdf"), width = 3.4, height = 2)
    DimPlot(gse %>% subset(idents = time), group.by = "scPharm_label_1168_Erlotinib",
            pt.size = 0.05,
            cols = c('#D3D3D3', '#FF0000', "#D3D3D3"))
    dev.off()
  }
}
Figure4c(gse)

# figure 4f -----------------------------------------------------------------------

Figure4f = function(gse){
  # label resistant cells and no resistant cells in D9 and D11
  gse$deg_group = "NOT"
  gse@meta.data[gse@meta.data$orig.ident == "D9" | gse@meta.data$orig.ident == "D11",]$deg_group = "D9D11_non_resistant"
  gse@meta.data[gse@meta.data$time_label == "D9_resistant" | gse@meta.data$time_label == "D11_resistant",]$deg_group = "D9D11_resistant"
  gse$time_label = paste(gse$orig.ident, gse$scPharm_label_1168_Erlotinib, sep = "_")
  gse_cp = gse
  CairoPDF("./Figure/GSE134839/dotplot_D9D11.pdf", width =4, height = 3.6)
  DotPlot(gse_cp, group.by = "deg_group" ,
          features = c("EGFR","FGFR1","MAP2K1","KRAS","NRAS","BRAF","RAF1"),
          idents = c("D9","D11")) +
    ylab("Group") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.title = element_text(size = 8),
          plot.margin = unit(c(0.1,0.4,0.1,0.1),'cm'),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    coord_flip()
  dev.off()
}
Figure4f(gse)