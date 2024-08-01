# Title    : Supplementary Figure1
# Author   : Wanglab
# Time     : 2024.7
library(Cairo)
library(dplyr)

# variables: luad, brca, skcm are the NES statics of cells from healthy lung, breast and skin tissue, respectively
luad = luad[,seq(6,595,2)]
luad = luad %>% gather(key = "Group", value = "NES")
luad = luad[luad$NES != 0,]

brca = brca[,seq(6,595,2)]
brca = brca %>% gather(key = "Group", value = "NES")
brca = brca[brca$NES != 0,]

skcm = skcm[,seq(6,595,2)]
skcm = skcm %>% gather(key = "Group", value = "NES")
skcm = skcm[skcm$NES != 0,]

# variables: null is the NES statics of cells from pooled tissues
# SuppleFigure1a---------------------------------------------------------------------
CairoPDF(paste0("./Figure/supfig_null/null_lung.pdf"), width = 2.6, height = 1.8)
par(oma = c(1.2, 3, 0.2, 1), mar = c(1, 0, 1.4, 0), xpd = TRUE)
plot(density(luad$NES), col = "#999999", main = "",
     lwd = 1.5, lty = 1, cex.axis = 1.2, font.main = 1, tcl = -0.25, ylab = "", xlab = "",
     bty = "l", las = 1, ylim = c(0, 1))
dev.off()

# SuppleFigure1b---------------------------------------------------------------------
CairoPDF(paste0("./Figure/supfig_null/null_breast.pdf"), width = 2.6, height = 1.8)
par(oma = c(1.2, 3, 0.2, 1), mar = c(1, 0, 1.4, 0), xpd = TRUE)
plot(density(brca$NES), col = "#999999", main = "",
     lwd = 1.5, lty = 1, cex.axis = 1.2, font.main = 1, tcl = -0.25, ylab = "", xlab = "",
     bty = "l", las = 1, ylim = c(0, 1))
dev.off()

# SuppleFigure1c---------------------------------------------------------------------
CairoPDF(paste0("./Figure/supfig_null/null_skcm.pdf"), width = 2.6, height = 1.8)
par(oma = c(1.2, 3, 0.2, 1), mar = c(1, 0, 1.4, 0), xpd = TRUE)
plot(density(skcm$NES), col = "#999999", main = "",
     lwd = 1.5, lty = 1, cex.axis = 1.2, font.main = 1, tcl = -0.25, ylab = "", xlab = "",
     bty = "l", las = 1, ylim = c(0, 1))
dev.off()

# SuppleFigure1d---------------------------------------------------------------------
CairoPDF(paste0("./Figure/supfig_null/null.pdf"), width = 2.6, height = 1.8)
par(oma = c(1.2, 3, 0.2, 1), mar = c(1, 0, 1.4, 0), xpd = TRUE)
plot(density(null$NES), col = "#999999", main = "",
     lwd = 1.5, lty = 1, cex.axis = 1.2, font.main = 1, tcl = -0.25, ylab = "", xlab = "",
     bty = "l", las = 1, ylim = c(0, 1))
dev.off()