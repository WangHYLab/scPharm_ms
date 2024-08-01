# Title    : Supplementary Figure 2
# Author   : Wanglab
# Time     : 2024.7

library(Cairo)

sample = c("AH0308", "MH0031", "MH0069", "MH0161", "MH0176", "PM0337")
names(sample) = sample

Function = function(sample,  drug_id, drug_name, cell_label){
  fig3.data.1 = lapply(sample, function(sam) {
    # Load the single-cell pharmacology result object from the specified path
    rank = readRDS(paste0("./scPharm/result/", sam, "_scPharm_object_nmcs_50_nfs_200.rds"))
    return(rank)
  })
  # Loop through each element in the fig3.data.1 list to generate a PDF chart for each
  for (i in 1:6) {
    meta.data = fig3.data.1[[i]]@meta.data
    temp = paste("scPharm_nes", drug_id, drug_name, sep = "_")
    meta.data = meta.data[,c("cell.label", temp)]
    # Filter cells where the label == cell_label
    meta.data = meta.data[meta.data$cell.label == cell_label,]
    colnames(meta.data) = c("Group", "NES")
    # Create a PDF file using the CairoPDF function, specifying filename, width, and height
    CairoPDF(paste0("./Figure/SuppleFig2", drug_name, names(fig3.data.1)[i],"20231117.pdf", sep="_"), width = 2, height = 1.4)
    # Set graphic parameters including outer margins, inner margins, and x-axis tick label positioning
    par(oma = c(1.2, 3, 0.2, 1), mar = c(1, 0, 1.4, 0), xpd = TRUE)
    plot(density(null$NES), col = "#999999", main = "",
        lwd = 1.5, lty = 2, cex.axis = 1.2, font.main = 1, tcl = -0.25, ylab = "", xlab = "",
        bty = "l", las = 1, ylim = c(0, 1.2))
    lines(density(meta.data$NES), col = "#FF0000", lwd=1.5)
    dev.off()
  }
}
# SuppleFigure2a,b-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Function(sample, 1558, "Lapatinib", "tumor")
Function(sample, 1558, "Lapatinib", "adjacent")

# SuppleFigure2c,d-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Function(sample, 1549, "Sapitinib", "tumor")
Function(sample, 1549, "Sapitinib", "adjacent")