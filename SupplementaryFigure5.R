# Title    : Supplementary Figure 5
# Author   : Wanglab
# Time     : 2024.7

library(Cairo)
library(Seurat)

sample = c("AH0319","MH0001","MH0025","MH0029-7C","MH0029-9C","MH0032","MH0040","MH0042","MH0043-T","MH0064-T",
           "MH0151","MH0163","MH0173-T","PM0360"  )
names(sample) = sample

# SuppleFigure 5a ---------------------------------------------------------------

SuppleFigure5a <- function (sample) {
  data.1 = lapply(sample, function(sam) {
    # Load the single-cell pharmacology result object from the specified path
    rank = readRDS(paste0("./scPharm/result/", sam, "_scPharm_object_nmcs_50_nfs_200.rds"))
    return(rank)
  })
  # Loop through each element in the data.1 list to generate a PDF chart for each
  for (i in 1:14) {
    meta.data = data.1[[i]]@meta.data
    meta.data = meta.data[,c("cell.label", "scPharm_nes_1816_Fulvestrant")]
    # Filter cells where the label is "tumor"
    meta.data = meta.data[meta.data$cell.label == "tumor",]
    colnames(meta.data) = c("Group", "NES")
    # Create a PDF file using the CairoPDF function, specifying filename, width, and height
    CairoPDF(paste0("./Figure/supplefig5a_Fulvestrant_high_", names(data.1)[i],"_T_20231117.pdf"), width = 2, height = 1.4)
    # Set graphic parameters including outer margins, inner margins, and x-axis tick label positioning
    par(oma = c(1.2, 3, 0.2, 1), mar = c(1, 0, 1.4, 0), xpd = TRUE)
    plot(density(null$NES), col = "#999999", main = "",
        lwd = 1.5, lty = 2, cex.axis = 1.2, font.main = 1, tcl = -0.25, ylab = "", xlab = "",
        bty = "l", las = 1, ylim = c(0, 1.2))
    lines(density(meta.data$NES), col = "#FF0000", lwd=1.5)
    dev.off()
  }
}
SuppleFigure5a(sample)


# SuppleFigure 5b ---------------------------------------------------------------

SuppleFigure5b <- function (sample) {
  data.1 = lapply(sample, function(sam) {
    # Load the single-cell pharmacology result object from the specified path
    rank = readRDS(paste0("./scPharm/result/", sam, "_scPharm_object_nmcs_50_nfs_200.rds"))
    return(rank)
  })
  # Loop through each element in the ata.1 list to generate a PDF chart for each
  for (i in 1:14) {
    meta.data = data.1[[i]]@meta.data
    meta.data = meta.data[,c("cell.label", "scPharm_nes_1816_Fulvestrant")]
    # Filter cells where the label is "adjacent"
    meta.data = meta.data[meta.data$cell.label == "adjacent",]
    colnames(meta.data) = c("Group", "NES")
    # Create a PDF file using the CairoPDF function, specifying filename, width, and height
    CairoPDF(paste0("./Figure/supplefig5b_Fulvestant_low_", names(data.1)[i],"_A_20231117.pdf"), width = 2, height = 1.4)
    # Set graphic parameters including outer margins, inner margins, and x-axis tick label positioning
    par(oma = c(1.2, 3, 0.2, 1), mar = c(1, 0, 1.4, 0), xpd = TRUE)
    plot(density(null$NES), col = "#999999", main = "",
        lwd = 1.5, lty = 2, cex.axis = 1.2, font.main = 1, tcl = -0.25, ylab = "", xlab = "",
        bty = "l", las = 1, ylim = c(0, 1.2))
    lines(density(meta.data$NES), col = "#FF0000", lwd=1.5)
    dev.off()
  }
}
SuppleFigure5b(sample)