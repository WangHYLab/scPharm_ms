# Title    : Supplementary Figure 7
# Author   : Wanglab
# Time     : 2024.7

library(Seurat)
library(readxl)
library(ggplot2)
library(ggpubr)
library(Cairo)
library(dplyr)
library(ggsci)


cancer_tcgaid = c("UCEC", "LIHC", "LGG", "BLCA",
                  "HNSC", "OV", "COREAD", "STAD",
                  "KIRC", "ESCA", "PAAD", "LUAD",
                  "BRCA", "SKCM")
# Load scPharm result of single cell data simulated by bulk data of cancer cell lines
cl.scPharm = lapply(cancer_tcgaid, function(type) {
  scPharm = readRDS(paste0("../scPharm/result/", type, "_scPharm_object_nmcs_50_nfs_200.rds"))
  return(scPharm)
})
names(cl.scPharm) = cancer_tcgaid

cl.scatter = function(fig2.data.4, fig2.data.5, cancer_type) {
  # fig2.data4: cancer cell line single cell data
  # fig2.data5: gdsc data
  meta.data = fig2.data.4@meta.data
  drug_id <- fig2.data.5[fig2.data.5$TCGA_DESC == cancer_type,] ## cancer type
  drug_id <- drug_id[,c(8,9,10,11)]
  drug_id <- drug_id[!duplicated(drug_id$DRUG_ID),]
  drug_id = drug_id %>%
    filter(!if_all(.fns = is.na))
  scatter.data = data.frame(matrix(0,1,5))
  colnames(scatter.data) = c("DRUG_ID","DRUG_NAME","MEDIAN_SENSI","MEDIAN_RESIS", "P_VAL1")
  message("calculate p value")
  for (i in seq(1, nrow(drug_id))) {
    x = fig2.data.5[(fig2.data.5$DRUG_ID == drug_id$DRUG_ID[i] & fig2.data.5$TCGA_DESC == cancer_type), c("CELL_LINE_NAME", "AUC")]
    x = na.omit(x)
    x = x[order(x[,2]),]
    colnames(x) = c("cellline", "AUC")
    # x$cellline = sapply(x$cellline, function(name) {return(toupper(gsub("-","",name)))})
    x$label = "other"
    x1 = x
    x1[x1$AUC >= median(x1$AUC), 3] = "AUC_high.50%"
    x1[x1$AUC < median(x1$AUC), 3] = "AUC_low.50%"
    # message(paste0(drug_id$DRUG_ID[i], ":x done"))
    if (paste("scPharm_nes", drug_id$DRUG_ID[i], drug_id$DRUG_NAME[i], sep = "_") %in% colnames(meta.data)) {
      meta.data$orig.ident = rownames(meta.data)
      y = meta.data[, c("orig.ident",paste("scPharm_nes", drug_id$DRUG_ID[i], drug_id$DRUG_NAME[i], sep = "_"))]
      colnames(y) = c("cellline", "response")
      # y$cellline = sapply(strsplit(y$orig.ident, split ="_"), function(x) {return(x[[1]][1])})
      y$group = "other"
      y[which(y$cellline %in% x1[x1$label == "AUC_high.50%",]$cellline),"group"] = "AUC_high.50%"
      y[which(y$cellline %in% x1[x1$label == "AUC_low.50%",]$cellline),"group"] = "AUC_low.50%"
      y = y[y$group != "other",]
      # message(paste0(drug_id$DRUG_ID[i], "y done"))
      # print(paste(sum(y$group == "AUC_low.50%"), sum(y$group == "AUC_high.50%"), sep = "---"))
      if (sum(y$group == "AUC_low.50%") < 3 | sum(y$group == "AUC_high.50%") < 3) {
        print("next2")
        next
      }
      p.val1 = wilcox.test(y[y$group=="AUC_low.50%",]$response, y[y$group=="AUC_high.50%",]$response,
                           alternative = "less") # x compare y
      p.val1 = p.val1$p.value
      scatter.data[nrow(scatter.data)+1,] = c(drug_id$DRUG_ID[i],
                                              drug_id$DRUG_NAME[i],
                                              median(y[y$group == "AUC_low.50%",]$response),
                                              median(y[y$group == "AUC_high.50%",]$response),
                                              p.val1)
    }else {
      print("next1")
      next
    }
  }
  if (nrow(scatter.data) > 1) {
    message("ploting scatter plot")
    scatter.data = scatter.data[-1,]
    scatter.data[,3:5] = apply(scatter.data[,3:5], 2, as.numeric)
    scatter.data$signif = "pval >= 0.1"
    scatter.data[scatter.data$P_VAL1 < 0.1, 6] = "pval < 0.1"
    axis_begin<- -4
    axis_end<-4
    total_ticks<-9
    tick_frame<-data.frame(ticks=seq(axis_begin,
                                     axis_end,
                                     length.out = total_ticks),
                           zero=0)%>%
      subset(ticks != 0)
    label_frame<-data.frame(lab=seq(axis_begin,axis_end),
                            zero=0)%>%
      subset(lab!=0)
    p_axis <- ggplot(tick_frame) +
      # draw axis line
      geom_segment(x=0,xend=0,y=-4.5,yend=4.5, linewidth = 0.4)+
      geom_segment(x=-4.5,xend=4.5,y=0,yend=0, linewidth = 0.4)+
      # x ticks
      geom_segment(data=tick_frame,aes(x=zero,xend=zero+0.1,
                                       y=ticks,yend=ticks))+
      # y ticks
      geom_segment(data=tick_frame,aes(x=ticks,xend=ticks,
                                       y=zero,yend=zero+0.1))+
      # labels
      geom_text(data=label_frame,aes(x=zero-0.3,y=lab,label=lab), size = 3, colour = "#4F4F4F")+
      geom_text(data=label_frame,aes(x=lab,y=zero-0.3,label=lab), size = 3, colour = "#4F4F4F")+
      theme(panel.grid = element_blank(), #次网格线
            panel.border = element_blank(), #边框,
            # panel.border = element_rect(fill = NA, size = 0.8, colour = "#000000", linetype = "solid"),
            axis.title = element_text(size = 12),
            axis.text = element_blank(),
            axis.ticks = element_blank())
    p2 <- p_axis + geom_point(data = scatter.data, aes(x = MEDIAN_SENSI, y = MEDIAN_RESIS,
                                                       col = signif), size = 1)+
      geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "dashed", linewidth = 0.4)+
      theme_bw() +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, size = 14),
            panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_text(size = 13),
            axis.title.y = element_text(size = 13, angle = 90),
            legend.text = element_text(size = 12),
            legend.title =  element_blank(),
            legend.position = c(0.76, 0.12),
            legend.background = element_rect(fill = rgb(1,1,1, alpha = 0.001), colour = NA)
      )+
      xlab("NES (AUC_low)")+
      ylab("NES (AUC_high)")+
      scale_color_manual(values = c("#6A5ACD", "#FD8D3C")) +
      guides(color=guide_legend(override.aes = list(size=2))) +
      ggtitle(paste0(cancer_type, "\n(# of drug : ",nrow(scatter.data),")"))
    CairoPDF(paste0("./Figure/figs7/", cancer_type,".pdf"), width = 3.4, height = 3.6)
    print(p2)
    dev.off()
    return(scatter.data)
  }else {
    return(NULL)
  }
}

for (i in 1:14) {
  print(names(cl.scPharm)[i])
  scatter.data = cl.scatter(fig2.data.4 = cl.scPharm[[i]], fig2.data.5 = gdsc.data, cancer_type = names(cl.scPharm)[i])
}