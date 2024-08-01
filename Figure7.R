# Title    : Figure 7
# Author   : Wanglab
# Time     : 2024.7

library(tidyverse)
library(readxl)
library(Cairo)

# split cells from all healthy samples by donor
# variables: luad, brca, skcm are the NES statics of cells from healthy lung, breast and skin tissue, respectively
healthy_lung = split(luad, luad$orig.ident)
healthy_breast = split(brca, brca$orig.ident)
healthy_skin = split(skcm, skcm$orig.ident)

# Extract DRUG_ID and DRUG_NAME of drugs from GDSC2
GDSC <- read_xlsx("../GDSC/GDSC2_fitted_dose_response_24Jul22.xlsx")
GDSC <- GDSC[GDSC$TCGA_DESC == 'LUAD',]
drug_id <- GDSC[,8:11]
drug_id = drug_id[!duplicated(drug_id$DRUG_ID),]
drug_id = drug_id %>%
  filter(!if_all(.fns = is.na))
drug_id = drug_id[order(drug_id$DRUG_ID),]

calculate_side_effect = function(healthy){
  # calculate the side effect score in pooled tissue data
  h = healthy[,seq(6,595,2)]
  copy_h = h
  for (i in 1:295) {
    copy_h[,i] = 'other'
    copy_h[rownames(h[h[,i] < -1.751302,]), i] = "sensitive"
    copy_h[rownames(h[h[,i] > 1.518551,]), i] = "resistant"
  }
  score = data.frame(DRUG_ID = sapply(strsplit(colnames(copy_h), split = "_"), function(x){x[2]}),
                     DRUG_NAME = drug_id[sapply(strsplit(colnames(copy_h), split = "_"), function(x){x[2]}),"DRUG_NAME"],
                     Dse = colMeans(copy_h == "sensitive"))
  return(score)
}


Figure7_heatmap = function(healthy, tissue, healthy_label){
  # calculate the side effect score
  healthy_tissue_label = lapply(healthy, function(h) {
    h = h[,seq(6,595,2)]
    copy_h = h
    # print(dim(copy_h))
    for (i in 1:295) {
      copy_h[,i] = 'other'
      copy_h[rownames(h[h[,i] < -1.751302,]), i] = "sensitive"
      copy_h[rownames(h[h[,i] > 1.518551,]), i] = "resistant"
    }
    score = data.frame(DRUG_ID = sapply(strsplit(colnames(copy_h), split = "_"), function(x){x[2]}),
                       DRUG_NAME = drug_id[sapply(strsplit(colnames(copy_h), split = "_"), function(x){x[2]}),"DRUG_NAME"],
                       Dse = colMeans(copy_h == "sensitive"))
                       # score$Dse = score$Sensi_ratio*(1 - score$Resis_ratio)
                       # score = score[order(-score[,3]),]
                       # score$Rank = seq(1,nrow(score),1)
    return(score)
  })
  # transform the data from list to data.frame
  Dse = bind_cols(healthy_tissue_label)
  Dse = Dse[,c(1,2, seq(3,ncol(Dse),3))]
  rownames(Dse) = paste(Dse$DRUG_ID...1, Dse$DRUG_NAME...2, sep = "_")
  Dse = Dse[,-c(1,2)]
  colnames(Dse) = names(healthy_tissue_label)
  Dse= scale(Dse, center = T, scale = TRUE)
  # label top 5 drug in pooled tissue data
  col_anno = data.frame(name = c(paste(healthy_label$LUAD$DRUG_ID[1:5],healthy_label$LUAD$DRUG_NAME[1:5], sep = "_")))
  col_labels = c(rep(" ", 295))
  names(col_labels) = rownames(Dse)
  col_labels[col_anno$name] = col_anno$name
  CairoPDF(paste("./Figure/supfig7/", tissue,".pdf"), width = 10, height = 3.2)
  Heatmap(t(Dse),
          # col = colorRamp2(c(0,1,2,3), c("grey", "#E41A1C", "#377EB8", "#4DAF4A")),
          column_title = "Drugs",
          column_title_gp = gpar(fontsize = 14),
          cluster_rows = F,
          cluster_columns = T,
          column_dend_height = unit(0,"mm"),
          column_dend_reorder = T,
          show_row_names = T,
          row_labels = paste("L", c(0:9, "a","b","c"), sep = ""),
          show_column_names = T,
          column_labels = col_labels,
          column_names_gp = gpar(fontsize = 2),
          row_names_side = "left",
          row_title = "",
          show_heatmap_legend = T,
          name = "Side\neffect",
          # bottom_annotation = bottom_anno
          # col = c("#0000FF","white","#FF0000"),
          # width = unit(12, "cm"), height = unit(5, "cm")
  )
  dev.off()
}


Figure7_density = function (healthy_label, cancer) {
  rank = healthy_label
  for (i in 1:5) {
    if (cancer == "LUAD") {
      healthy_nes = luad
    } else if (cancer == "BRCA") {
      healthy_nes = brca
    } else {
      healthy_nes = skcm
    }
    CairoPDF(paste0("./Figure/healthy_dr/", cancer,"_",rank[i,1],"_",rank[i,2],".pdf"), width = 2, height = 1.4)
    par(oma = c(1.2, 3, 0.2, 1), mar = c(1, 0, 1.4, 0), xpd = TRUE)
    plot(density(null$NES), col = "#999999", main = rank[i,2], cex.main = 0.8,
         lwd = 1.5, lty = 2, cex.axis = 0.8, font.main = 1, tcl = -0.25, ylab = "", xlab = "",
         bty = "l", las = 1, ylim = c(0, 1.2))
    lines(density(healthy_nes[,paste0("id_", rank[i,1],"_nes")]), col = "#FF0000", lwd=1.5)
    dev.off()
  }
}

# Figure7b -----------------------------------------------------------------------
healthy_label = calculate_side_effect(skcm)
Figure7_heatmap(healthy_skin, "skin", healthy_label)

# Figure7c -----------------------------------------------------------------------
Figure7_density(healthy_label, "SKCM")


# Figure7d -----------------------------------------------------------------------
healthy_label = calculate_side_effect(brca)
Figure7_heatmap(healthy_breast, "breast", healthy_label)

# Figure7e -----------------------------------------------------------------------
Figure7_density(healthy_label, "BRCA")


# Figure7f -----------------------------------------------------------------------
healthy_label = calculate_side_effect(luad)
Figure7_heatmap(healthy_lung, "lung", healthy_label)

# Figure7g -----------------------------------------------------------------------
Figure7_density(healthy_label, "LUAD")
