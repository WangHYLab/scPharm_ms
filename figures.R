library(patchwork)
library(readxl)
library(ComplexHeatmap)

GDSC <- read_xlsx("../GDSC/GDSC2_fitted_dose_response_24Jul22.xlsx")
GDSC <- GDSC[GDSC$TCGA_DESC == 'LUAD',]
drug_id <- GDSC[,8:11]
drug_id = drug_id[!duplicated(drug_id$DRUG_ID),]
drug_id = drug_id %>% 
  filter(!if_all(.fns = is.na))
drug_id = drug_id[order(drug_id$DRUG_ID),]




# figure 2a ---------------------------------------------------------------


fig2.data.1 = readRDS("./scPharm/info/bulkdata.rds") # 细胞系表达谱TPM
fig2.data.2 = read_xlsx("./scPharm/info/GDSC2_fitted_dose_response_24Jul22.xlsx") # GDSC药敏数据

fig2.data.2 = fig2.data.2[(fig2.data.2$TCGA_DESC == "LUAD" & fig2.data.2$DRUG_ID == 2169),]
fig2.data.2 = fig2.data.2[which(fig2.data.2$CELL_LINE_NAME %in% colnames(fig2.data.1)),] #筛选有表达谱的细胞系
fig2.data.2 = fig2.data.2[order(fig2.data.2$AUC, decreasing = T),]
fig2.data.1 = fig2.data.1[,fig2.data.2$CELL_LINE_NAME]
fig2.data.1 = fig2.data.1[rowSums(fig2.data.1) > 0,] # 去除不表达的基因

geneList = apply(fig2.data.1, 1, function(expr) {
  return(cor.test(expr, fig2.data.2$AUC, method = "pearson", alternative = 'two.sided')$estimate)
})
geneList = sort(geneList, decreasing = TRUE) #计算相关性
fig2.data.1 = data.frame(t(scale(t(fig2.data.1), scale = T, center = T)), check.names = F) # z-score标化表达谱
fig2.data.1 = fig2.data.1[rowSums(fig2.data.1 > -0.25) > 1,] # 去除只在某个细胞系中高表达的基因

geneList = geneList[which(names(geneList) %in% rownames(fig2.data.1))]
fig2.data.1 = fig2.data.1[names(geneList),]

fig2.data.1 = fig2.data.1[,-3]
# 列注释，AUC
# phenotype = 1:12
# names(phenotype) = rownames(fig2.data.2)
# fig2a.left.col.ha = HeatmapAnnotation(
#   phenotype = phenotype, # pch = 1, 
#                           # pt_size = unit(c(seq(12,6,-1)/4, seq(6,12,1)/4), "mm"),
#                           # pt_gp = gpar(col = c(rep("red",6), rep("blue",6))),
#                           
#   which = "col",
#   show_legend = T,
#   show_annotation_name = F,
#   simple_anno_size = unit(1.5, "mm"),
#   # height = unit(0.1, "mm"),
#   annotation_legend_param = list(title_gp = gpar(fontsize = 10, fontfamily = "sans"),
#                                  at = c(3, 10),
#                                  labels = c("Resistant", "Sensitive"),
#                                  label_gp = gpar(fontsize = 8,fontfamily = "sans")),
#   col = list(
#     phenotype = colorRamp2(
#       c(1, 12), 
#       c("#045A8D", "#ECE7F2")
#     )
#   )
# )

# 行注释，基因排名
fig2a.left.row.ha = HeatmapAnnotation(
  cor = geneList,
  which = "row",
  show_legend = F,
  show_annotation_name = F,
  simple_anno_size = unit(8, "mm"),
  # width = unit(1, "mm"),
  # annotation_legend_param = list(title_gp = gpar(fontsize = 10, fontfamily = "sans")),
  col = list(
    cor = colorRamp2(
      c(min(geneList), quantile(geneList, 0.1), quantile(geneList, 0.2), 0, quantile(geneList, 0.8),quantile(geneList, 0.9), max(geneList)), 
      c("#313695", "#74ADD1", "#E0F3F8", "#F7F7F7", "#FDE0EF", "#D53E4F", "#9E0142")
    )
  )
)
fig2a.left = Heatmap(fig2.data.1,
                     col = colorRamp2(seq(min(fig2.data.1), max(fig2.data.1), length = 3), c("#8470FF", "#EEEEEE", "red")),
                     cluster_rows = F, 
                     cluster_columns = F,
                     show_row_names = F,
                     show_column_names = F,
                     row_title = "Ranked gene list\n(corr)",
                     row_title_gp = gpar(fontsize = 14, fontfamily = "Arial"),
                     row_title_rot = 90,
                     column_title = "Cell line",
                     column_title_gp = gpar(fontsize = 14, fontfamily = "Arial"),
                     column_title_side = "top",
                     show_heatmap_legend = F,
                     name = "Expr",
                     heatmap_legend_param = list(title_gp = gpar(fontsize = 14, fontfamily = "Arial")),
                     # top_annotation = fig2a.left.col.ha,
                     # right_annotation = fig2a.left.row.ha,
                     width = unit(6, "cm"),
                     height = unit(7, "cm"),
                     use_raster = TRUE
)
CairoPNG("./Figure/fig2a_left.png", units = "in",
         width = 3.6, height = 3.2, dpi = 300)
print(fig2a.left)
dev.off()
CairoPDF("./Figure/fig2a_left.pdf",width = 3, height = 3.4)
print(fig2a.left)
dev.off()


# colnames(fig2.data.3) = "gene"
bulk_scPharm = readRDS("./scPharm/result/luad_bulkdata_scPharm_object_nmcs_20_nfs_500.rds")
cell_line_2169_nes = bulk_scPharm@meta.data[,c("scPharm_label_2169_AZD6482", "scPharm_nes_2169_AZD6482")]
pos_cl = rownames(cell_line_2169_nes[cell_line_2169_nes$scPharm_nes_2169_AZD6482 == max(cell_line_2169_nes$scPharm_nes_2169_AZD6482),])
neg_cl = rownames(cell_line_2169_nes[cell_line_2169_nes$scPharm_nes_2169_AZD6482 == min(cell_line_2169_nes$scPharm_nes_2169_AZD6482),])
zero_cl = "TK10"

show_gene = GetCellGeneSet(bulk_scPharm,
                           reduction = "mca",
                           dims = 1:20,
                           n.features = 200) # [c(pos_cl, zero_cl, neg_cl)]



# fig2.data.3[show_gene[[1]], 1] = 1
# fig2.data.3[show_gene[[2]], 2] = 1
# fig2.data.3[show_gene[[3]], 3] = 1
# colnames(fig2.data.3) = c("pos","neg","zero")

for (i in 1:length(show_gene)) {
  fig2.data.3 = data.frame(matrix(0,nrow = length(geneList), ncol = 1))
  rownames(fig2.data.3) = names(geneList)
  fig2.data.3[show_gene[[i]], 1] = 1
  fig2a.right = Heatmap(fig2.data.3,
                        col = colorRamp2(c(0, 1), c("white", "#FC8D59")),
                        cluster_rows = F, 
                        cluster_columns = F,
                        show_row_names = F,
                        show_column_names = F,
                        show_heatmap_legend = F,
                        name = "Expression",
                        border_gp = gpar(code = "black", lwd = 0.2),
                        width = unit(0.8, "cm"),
                        height = unit(12, "cm")
  )
  CairoPNG(paste("./Figure/fig2a_right/fig2a_blackright_", names(show_gene)[i], ".png" ), units = "in",
           width = 0.5, height = 5.6, dpi = 300)
  print(fig2a.right)
  dev.off()
}








# figure 2b ---------------------------------------------------------------


library(ggpubr)
fig2.data.4 = readRDS("./scPharm/result/luad_bulkdata_scPharm_object_nmcs_20_nfs_500.rds")
fig2.data.5 = read_xlsx("./scPharm/info/GDSC2_fitted_dose_response_24Jul22.xlsx")

# 单药AUC-NES图
fig2b = function(fig2.data.4, fig2.data.5) {
  meta.data = fig2.data.4@meta.data
  drug_id <- fig2.data.5[fig2.data.5$TCGA_DESC == "LUAD",]
  drug_id <- drug_id[,c(8,9,10,11)]
  drug_id <- drug_id[!duplicated(drug_id$DRUG_ID),]
  drug_id = drug_id %>% 
    filter(!if_all(.fns = is.na))
  scatter.data = data.frame(matrix(0,1,5))
  colnames(scatter.data) = c("DRUG_ID","DRUG_NAME","MEDIAN_SENSI","MEDIAN_RESIS", "P_VAL1")
  message("calculate p value")
  for (i in seq(1, nrow(drug_id))) {
    x = fig2.data.5[(fig2.data.5$DRUG_ID == drug_id$DRUG_ID[i] & fig2.data.5$TCGA_DESC == "LUAD"), c("CELL_LINE_NAME", "AUC")]
    x = x[order(x[,2]),]
    y = meta.data[x$CELL_LINE_NAME, paste("scPharm_nes", drug_id$DRUG_ID[i], drug_id$DRUG_NAME[i], sep = "_")]
    x$nes = y
    x = na.omit(x)
    x$label = "other"
    
    x1 = x
    x1[x1$AUC >= median(x1$AUC), 4] = "AUC_high.50%"
    x1[x1$AUC < median(x1$AUC), 4] = "AUC_low.50%"
    colnames(x1) = c("cellline", "AUC", "response","label")
    if (sum(x1$label == "AUC_low.50%") < 3 | sum(x1$label == "AUC_high.50%") < 3) {
      next
    }
    
    p.val1 = wilcox.test(x1[x1$label=="AUC_low.50%",]$response, x1[x1$label=="AUC_high.50%",]$response,
                         alternative = "less")
    p.val1 = p.val1$p.value

    if (drug_id$DRUG_ID[i] == 1079) {
      p1 <- ggplot(x1, aes(x = label,
                         y = response,
                         fill = label))+
        geom_violin(linewidth=0.4,col="black")+
        geom_boxplot(width=.1,col="black",fill="white", linewidth = 0.3)+
        scale_fill_manual(values = c("#77787b","#77787b"))+ #c("#FF0000","#0000FF")
        theme_classic()+
        theme(plot.title = element_text(hjust = 0.3,size = 12),
              axis.text.y = element_text(size = 12,colour = "#000000"),
              axis.text.x = element_text(size = 12,colour = "#000000", angle = 30, hjust = 1),
              axis.title.y = element_text(size = 12),
              plot.margin = unit(c(0.2,0.3,0.2,0.3),'cm'),
              legend.position = "none"
        )+
        ylim(min(x1$response)-0.2, max(x1$response)+0.5)+
        xlab(NULL)+
        ylab("NES")+
        ggtitle(paste0("ABL Inhibitor:",drug_id[drug_id$DRUG_ID == 1079, 2])) +
        stat_compare_means(aes(label = paste0("pval = ", ..p.format..)),
                           method = "wilcox.test",
                           method.args = list(alternative = "greater"),
                           label.x = 1, label.y = max(x1$response)+0.3,
                           size=4.5)
      # CairoPNG("./Figure/fig2b_left.png", units = "in", width = 2.2, height = 2.6, dpi = 300)
      # print(p1)
      # dev.off()
      CairoPDF("./Figure/fig2b_left.pdf", width = 2.8, height = 3.6)
      print(p1)
      dev.off()
    }
    scatter.data[nrow(scatter.data)+1,] = c(drug_id$DRUG_ID[i],
                                            drug_id$DRUG_NAME[i],
                                            median(x1[x1$label == "AUC_low.50%",]$response),
                                            median(x1[x1$label == "AUC_high.50%",]$response),
                                            p.val1)
  }
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
  
  p2 <- p_axis + geom_point(data = scatter.data, aes(x = MEDIAN_SENSI, y = MEDIAN_RESIS, col = signif),
                            size = 1)+
    geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "dashed", linewidth = 0.4)+
    theme_bw() +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 12),
          axis.title.y = element_text(size = 12, angle = -90),
          legend.text = element_text(size = 12),
          legend.title =  element_blank(),
          legend.position = c(0.76, 0.12),
          legend.background = element_rect(fill = rgb(1,1,1, alpha = 0.001), colour = NA)
          )+
    # scale_x_continuous(expand = c(0,0),limits = c(-4,4))+
    # scale_y_continuous(expand = c(0,0),limits = c(-4,4))+
    xlab("NES_median\n(AUC_low.50%)")+
    ylab("NES_median\n(AUC_high.50%)")+
    # scale_color_gradient(low="#6A5ACD", high = "#FD8D3C")+
    scale_color_manual(values = c("#6A5ACD", "#FD8D3C")) +
    guides(color=guide_legend(override.aes = list(size=2))) +
    ggtitle("All drugs(n=295)")
    
  # CairoPNG("./Figure/fig2b_right.png", units = "in", width = 3.6, height = 3.6, dpi = 300)
  # print(p2)
  # dev.off()
  CairoPDF("./Figure/fig2b_right.pdf", width = 3.4, height = 3.6)
  print(p2)
  dev.off()
  return(scatter.data)
}

fig2b.scatter = fig2b(fig2.data.4 = fig2.data.4, fig2.data.5 = fig2.data.5)


# figure 2c ---------------------------------------------------------------


# 单药AUC-NES曲线图
fig2c = function(fig2.data.4, fig2.data.5) {
  meta.data = fig2.data.4@meta.data
  drug_id <- fig2.data.5[fig2.data.5$TCGA_DESC == "LUAD",]
  drug_id <- drug_id[,c(8,9,10,11)]
  drug_id <- drug_id[!duplicated(drug_id$DRUG_ID),]
  drug_id = drug_id %>% 
    filter(!if_all(.fns = is.na))
  scatter.data = data.frame(matrix(0,1,5))
  colnames(scatter.data) = c("DRUG_ID","DRUG_NAME","MEDIAN_SENSI","MEDIAN_RESIS", "P_VAL1")
  message("calculate p value")
  for (i in seq(1, nrow(drug_id))) {
    x = fig2.data.5[(fig2.data.5$DRUG_ID == drug_id$DRUG_ID[i] & fig2.data.5$TCGA_DESC == "LUAD"), c("CELL_LINE_NAME", "AUC")]
    x = x[order(x[,2]),]
    y = meta.data[x$CELL_LINE_NAME, paste("scPharm_nes", drug_id$DRUG_ID[i], drug_id$DRUG_NAME[i], sep = "_")]
    x$nes = y
    x = na.omit(x)
    x$label = "other"
    
    x1 = x
    x1[x1$AUC >= quantile(x1$AUC, 0.8), 4] = "AUC_high.20%"
    x1[x1$AUC <= quantile(x1$AUC, 0.2), 4] = "AUC_low.20%"
    colnames(x1) = c("cellline", "AUC", "response","label")
    if (sum(x1$label == "AUC_low.20%") < 3 | sum(x1$label == "AUC_high.20%") < 3) {
      next
    }
    
    p.val1 = wilcox.test(x1[x1$label=="AUC_low.20%",]$response, x1[x1$label=="AUC_high.20%",]$response,
                         alternative = "less")
    p.val1 = p.val1$p.value
    scatter.data[nrow(scatter.data)+1,] = c(drug_id$DRUG_ID[i],
                                            drug_id$DRUG_NAME[i],
                                            median(x1[x1$label == "AUC_low.20%",]$response),
                                            median(x1[x1$label == "AUC_high.20%",]$response),
                                            p.val1)
    if (drug_id$DRUG_ID[i] == 1079) {
      x1 = x1[x1$label != "other",]
      p1 <- ggplot(x1, aes(x = label,
                           y = response,
                           fill = label))+
        geom_violin(linewidth=0.4, col="black")+
        geom_boxplot(width=.1,col="black",fill="white", linewidth = 0.3)+
        scale_fill_manual(values = c("#77787b","#77787b"))+ #c("#FF0000","#0000FF")
        theme_classic()+
        theme(plot.title = element_text(hjust = 0.3,size = 12),
              axis.text.y = element_text(size = 12,colour = "#000000"),
              axis.text.x = element_text(size = 12,colour = "#000000", angle = 30, hjust = 1),
              axis.title.y = element_text(size = 12),
              plot.margin = unit(c(0.2,0.3,0.2,0.3),'cm'),
              legend.position = "none"
        )+
        ylim(min(x1$response)-0.2, max(x1$response)+0.5)+
        xlab(NULL)+
        ylab("NES")+
        ggtitle(paste0("ABL Inhibitor:",drug_id[drug_id$DRUG_ID == 1079, 2])) +
        stat_compare_means(aes(label = paste0("pval = ", ..p.format..)),
                           method = "wilcox.test",
                           method.args = list(alternative = "greater"),
                           label.x = 1, label.y = max(x1$response)+0.3,
                           size=4.5)

      CairoPDF("./Figure/fig2c_left.pdf", width = 2.8, height = 3.6)
      print(p1)
      dev.off()
    }
    
  }
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
  
  p2 <- p_axis + geom_point(data = scatter.data, aes(x = MEDIAN_SENSI, y = MEDIAN_RESIS, col = signif),
                            size = 1)+
    geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "dashed", linewidth = 0.4)+
    theme_bw() +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 12),
          axis.title.y = element_text(size = 12, angle = -90),
          legend.text = element_text(size = 12),
          legend.title =  element_blank(),
          legend.position = c(0.76, 0.12),
          legend.background = element_rect(fill = rgb(1,1,1, alpha = 0.001), colour = NA)
    )+
    xlab("NES_median\n(AUC_low.20%)")+
    ylab("NES_median\n(AUC_high.20%)")+
    scale_color_manual(values = c("#6A5ACD", "#FD8D3C")) +
    guides(color=guide_legend(override.aes = list(size=2))) +
    ggtitle("All drugs(n=295)")
  
  CairoPDF("./Figure/fig2c_right.pdf", width = 3.4, height = 3.6)
  print(p2)
  dev.off()
  return(scatter.data)
}

fig2c.scatter = fig2c(fig2.data.4 = fig2.data.4, fig2.data.5 = fig2.data.5)


# figure 2d ---------------------------------------------------------------

library(plotrix)

number1 = length(unique(union(fig2b.scatter[fig2b.scatter$P_VAL1 < 0.1,]$DRUG_ID,
                              fig2c.scatter[fig2c.scatter$P_VAL1 < 0.1,]$DRUG_ID)
                        )
                 )
number2 = nrow(fig2b.scatter) - number1

fig2.data.6 = data.frame(cluster = c("included", "excluded"),
                         number = c(number1, number2))
fig2.data.6$ratio = round((fig2.data.6$number / sum(fig2.data.6$number))*100, 2)

CairoPNG("./Figure/fig2d.png", units = "in", width = 4, height = 4, dpi = 300)
# CairoPDF("./Figure/fig2d.pdf", width = 4, height = 3.5)

pielabels = paste(fig2.data.6$ratio, "%", sep = "")
pie3D(fig2.data.6$number, 
      radius=1,
      height = 0.1,
      start = 2.05,
      labels = pielabels, 
      theta = pi/6, 
      labelcex=1,
      explode = 0.06, 
      main = "",
      shade = 0.5,
      col = c("#6A5ACD", "#FD8D3C"),
      border = NULL
      )
dev.off()

analysable = unique(union(fig2b.scatter[fig2b.scatter$P_VAL1 < 0.1,]$DRUG_ID,
                          fig2c.scatter[fig2c.scatter$P_VAL1 < 0.1,]$DRUG_ID))
analysable = fig2c.scatter[which(fig2c.scatter$DRUG_ID %in% analysable),c(1,2)]
write.csv(analysable, file = "./scPharm/info/LUAD_analysable_drug.csv", row.names = F)
# Figure3 -----------------------------------------------------------------

fig3.data.1 = lapply(sample, function(sam) {
  rank = readRDS(paste0("./scPharm/result/", sam, "_scPharm1114_object_nmcs_50_nfs_200.rds"))
  return(rank)
})

fig3.data.1 = fig3.data.1[c("AH0308", "MH0031", "MH0069", "PM0337")]

# figure 3a&b ---------------------------------------------------------------


for (i in 1:6) {
  meta.data = fig3.data.1[[i]]@meta.data
  meta.data = meta.data[,c("cell.label", "scPharm_nes_1558_Lapatinib")]
  meta.data = meta.data[meta.data$cell.label == "tumor",]
  colnames(meta.data) = c("Group", "NES")
  # meta.data = rbind(meta.data, null)
  # p <- ggplot(meta.data, aes(NES, fill=Group, color=Group, linetype=Group))+
  #   geom_density(alpha=0, linewidth = 0.5)+
  #   scale_linetype_manual(values = c("dashed","solid"))+ # "dotted",
  #   scale_color_manual(values = c("#999999",'#FF0000'))+ # "#00FFFF",
  #   geom_rug()+
  #   # ggtitle(plot_title)+
  #   theme_classic()+
  #   # labs(tag = tag.list[i])+
  #   xlab(NULL)+
  #   ylab("Density")+
  #   # ylim(0,ylimit)+
  #   # scale_y_continuous(breaks = 0:5, labels = c("0.0","1.0","2.0","3.0","4.0","5.0"))+
  #   theme(plot.title = element_text(hjust = 0.5,size = 12),
  #         plot.margin = unit(c(0.2,0.1,0.2,0.1),'cm'),
  #         # plot.tag = element_text(family = "sans", face = "bold", size = 26),
  #         # axis.text = element_text(colour = "#000000"),
  #         axis.text = element_text(size = 10, colour = "#000000"),
  #         # title = element_text(family = "sans", face = "bold", size = 14),
  #         axis.title = element_blank(),
  #         legend.text = element_text(size = 12),
  #         legend.title = element_text(size = 12),
  #         legend.position = "right"
  #   )+
  #   ggtitle(paste(names(fig3.data.1)[i],"(HER2+)"))
  CairoPDF(paste0("./Figure/fig3a/fig3a_lapatinib_", names(fig3.data.1)[i],"_20231117.pdf"), width = 2, height = 1.4)
  par(oma = c(1.2, 3, 0.2, 1), mar = c(1, 0, 1.4, 0), xpd = TRUE)
  plot(density(null$NES), col = "#999999", main = "",
       lwd = 1.5, lty = 2, cex.axis = 1.2, font.main = 1, tcl = -0.25, ylab = "", xlab = "",
       bty = "l", las = 1, ylim = c(0, 1.2))
  lines(density(meta.data$NES), col = "#FF0000", lwd=1.5)
  dev.off()
}

# figure 3c ---------------------------------------------------------------


library(GSEABase)
library(ggpubr)

kegg.gmt = getGmt("../kegg_pathway.gmt")
fig3b.data = map2(fig3.data.1, names(fig3.data.1), function(sample, name) {
  assay = GetAssayData(sample, slot = "data")
  expr = assay[which(rownames(assay) %in% kegg.gmt[["ErbB signaling pathway"]]@geneIds),]
  expr = as.data.frame(t(expr))
  expr$mean_expr = rowMeans(expr)
  sample = AddMetaData(sample, expr[,'mean_expr'], col.name = "target_pathway_activity")
  viol_data = sample@meta.data[,c("cell.label", "scPharm_label_1549_Sapitinib", "target_pathway_activity")]
  viol_data = viol_data[viol_data$cell.label == "tumor",c("scPharm_label_1549_Sapitinib", "target_pathway_activity")]
  colnames(viol_data) = c("cell.label", "target_pathway_activity")
  # viol_data[viol_data$cell.label != "sensitive", "cell.label"] = "no.sensitive"
  viol_data = viol_data[viol_data$cell.label != "other",]
  viol_data$sample = name
  return(viol_data)
})

for (i in 1:6) {
  CairoPDF(paste0("./Figure/fig3c/fig3c_sr_sapitinib_", names(fig3b.data)[i], ".pdf"), width = 4, height = 2)
  p <- ggplot(fig3b.data[[i]], aes(x = cell.label,
                                   y = as.numeric(target_pathway_activity),
                                   fill = cell.label)) +
    stat_boxplot(mapping=aes(x=cell.label,y=as.numeric(target_pathway_activity)),
                 geom ="errorbar",                             ##添加箱子的bar为最大、小值
                 width=0.15,position=position_dodge(0.4))+
    geom_boxplot(aes(fill = cell.label), 
                 position=position_dodge(0.4),                 ##因为分组比较，需设组间距
                 width=0.4,                                    ##箱子的宽度
                 outlier.color = "#999999",
                 outlier.size = 0.8,
                 linewidth = 0.3)+
    stat_summary(mapping=aes(group=cell.label),                    ##分组计算的变量
                 fun="mean",                                   ##箱线图添加均值
                 geom="point",shape=23,size=2,fill="white",    ##均值图形的设置
                 position=position_dodge(0.4))+                ##因为分组比较，需设置两组间距
    scale_fill_manual(values = c("#D73027","#4575B4"))+
    stat_compare_means(aes(group = cell.label), label = "p.signif", method = "wilcox.test",
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***","**","*","ns")),
                       label.x = 1.5,
                       label.y = max(fig3b.data[[i]]$target_pathway_activity)+0.01,
                       size=5)+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          axis.text.y = element_text(size = 10, colour = "#000000"),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 12, colour = "#000000"),
          text = element_text(colour = "#000000"),
          plot.margin = unit(c(0.2,0.1,0.2,0.1),'cm'),
          legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size = 12)
    )+
    ylim(0, max(fig3b.data[[i]]$target_pathway_activity)+0.05)+
    xlab(NULL)+
    ylab("")+
    ggtitle(names(fig3b.data)[i])
  # assign(paste("fig3b", i, sep = "_"), p)
  print(p)
  dev.off()
}

fig3b <- (fig3b_1 | fig3b_2 | fig3b_3 | fig3b_4) +
  plot_layout(guides = 'collect')
CairoPDF("./Figure/fig3c/fig3b.pdf", width = 14, height = 2)
print(fig3b)
dev.off()

MH0176 = readRDS("./scPharm/result/MH0176_scPharm1114_object_nmcs_50_nfs_200.rds")
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

# figure 3e ---------------------------------------------------------------


## 天梯图
sample = c("AH0308", "MH0031", "MH0069", "MH0161", "MH0176", "PM0337")
names(sample) = sample
fig3.data.2 = lapply(sample, function(sam) {
  rank = read.csv(paste0("./scPharm/result/", sam, "_scPharm1208_drug_rank.csv"))
  return(rank)
})

for (sam in names(fig3.data.2)) {
  data = fig3.data.2[[sam]][1:30,]
  data$Rank = 31-data$Rank
  data$group = "other"
  for (id in c(1032, 1549, 1558)) {
    if (id %in% data$DRUG_ID) {
      data[data$DRUG_ID == id,]$group = data[data$DRUG_ID == id,]$DRUG_NAME
    }else {
      next
    }
  }
  if (length(unique(data$group)) < 4) {
    col_value = c("#E41A1C", "#4DAF4A", "grey")
  }else {
    col_value = c("#E41A1C", "#4DAF4A", "grey", "#377EB8")
  }
  
  
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
  
  CairoPDF(paste0("./Figure/fig3e/fig3e_", sam, ".pdf"), width = 2.2, height = 5.8)
  print(p)
  dev.off()
}

# fig3c = (fig3c.AH0308 | fig3c.MH0031 | fig3c.MH0069 | fig3c.PM0337) +
#   plot_layout(guides = 'auto') # &
#   # theme(legend.position = "bottom")
# CairoPNG("./Figure/fig3c.png", units = "in", width = 10, height = 6, dpi = 300)
# print(fig3c)
# dev.off()



# figure 3d ---------------------------------------------------------------


library(ComplexHeatmap)

fig3.data.4 = data.frame(matrix(0, nrow = nrow(fig3.data.2[[1]]), ncol = length(fig3.data.2)))
colnames(fig3.data.4) = names(fig3.data.2)
for (sam in names(fig3.data.2)) {
  data = fig3.data.2[[sam]]
  fig3.data.4[data[data$DRUG_ID == 1032,]$Rank, sam] = 1
  fig3.data.4[data[data$DRUG_ID == 1549,]$Rank, sam] = 2
  fig3.data.4[data[data$DRUG_ID == 1558,]$Rank, sam] = 3
}

# group = data.frame(sample = colnames(fig3.data.4), group = c("a","b","c","d"))
# group$group = factor(group$group)

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

for (i in 1:14) {
  print(supfig2.data.2[[i]][supfig2.data.2[[i]]$DRUG_NAME == "Docetaxel", c("Rank","DRUG_ID")])
}
# figure 5a ---------------------------------------------------------------

library(igraph)
library(tidyverse)
library(RColorBrewer)

fig5.data.1 = readRDS("../verify_data/GSE161529_HER2/scPharm_targeted_drug_combo1114.rds")
# location = c("central",
#              "nordic","southwestern",
#              "southeastern","southern","sourthern",
#              "western","northwestern","northeastern","northern",
#              )
set.seed(0)
CairoPDF("./Figure/fig5a20231114.pdf",
         width = 9, height = 8)
par(mfrow = c(2,3))
for (i in 1:6) {
  data = fig5.data.1[[i]][[1]]
  data = data[!duplicated(data$DRUG_NAME),]
  name = data.frame(c(data$DRUG_FIRST, data$DRUG_NAME))
  nodes = name %>%
    distinct() # %>%
  # mutate(location=location[1:(nrow(data)+1)])
  colnames(nodes) = c("drug") # , "location")
  
  edges = data[,c(3,1,4)]
  colnames(edges) = c("from", "to", "weight")
  
  net_pc = graph_from_data_frame(
    d = edges, vertices = nodes,
    directed = T)
  ###计算节点的度
  deg<-degree(net_pc,mode="all")
  ###指定节点的颜色、大小
  V(net_pc)$color<- "#999999" #c("#7B68EE",brewer.pal(n=nrow(edges), name = "Paired"))
  
  ###指定边的颜色
  E(net_pc)$width<-E(net_pc)$weight*2
  E(net_pc)$color<-"#999999" #brewer.pal(n=nrow(edges), name = "Paired")
  E(net_pc)$length<-E(net_pc)$weight*2
  
  l <- layout_as_star(net_pc)
  plot(net_pc,
       vertex.size=c(sum(E(net_pc)$width*2), E(net_pc)$width*5),
       layout = l,
       vertex.frame.color=NA,
       vertex.label.cex=1.2,
       vertex.label.dist=1.5,
       vertex.label.color="black",
       vertex.label.font=1,
       vertex.label.degree=c(rep(-pi/2,5), rep(pi/2,5)),
       edge.arrow.size = 0,
       edge.curved=.1)
}
dev.off()

# figure 5b ---------------------------------------------------------------

fig5.data.2 = readRDS("../verify_data/GSE161529_HER2/scPharm_toxicity_index1114.rds")

for (sam in names(fig5.data.2)) {
  data = fig5.data.2[[sam]][1:10,]
  data$Rank = 11-data$Rank
  data$group = "other"
  # for (id in c(1032, 1549, 1558)) {
  #   if (id %in% data$DRUG_ID) {
  #     data[data$DRUG_ID == id,]$group = data[data$DRUG_ID == id,]$DRUG_NAME
  #   }else {
  #     next
  #   }
  # }
  # if (length(unique(data$group)) < 4) {
  #   col_value = c("#E41A1C", "#4DAF4A", "#B3CDE3")
  # }else {
  #   col_value = c("#E41A1C", "#4DAF4A", "#B3CDE3", "#377EB8")
  # }
  
  
  p = ggplot(data, aes(Rank, R2, fill = group))+
    geom_col(width = 0.8)+
    coord_flip()+
    scale_x_continuous(name = "Drug", breaks = seq(10,1,-1), labels = data$DRUG_NAME)+
    scale_fill_manual(values = c("grey")) +
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
  CairoPDF(paste0("./Figure/fig5b/fig5b_", sam, "_1114.pdf"), width = 2.2, height = 1.8)
  print(p)
  dev.off()
}

# figure 5c ---------------------------------------------------------------

# HER2靶向药毒性排名

fig5.data.3 = data.frame(matrix(0, nrow = nrow(fig5.data.2[[1]]), ncol = length(fig5.data.2)))
colnames(fig5.data.3) = names(fig5.data.2)
for (sam in names(fig5.data.2)) {
  data = fig5.data.2[[sam]]
  fig5.data.3[data[data$DRUG_ID == 1032,]$Rank, sam] = 1
  fig5.data.3[data[data$DRUG_ID == 1549,]$Rank, sam] = 2
  fig5.data.3[data[data$DRUG_ID == 1558,]$Rank, sam] = 3
}

fig5c = Heatmap(t(fig5.data.3),
                # row_split = 1:4,
                col = colorRamp2(c(0,1,2,3), c("grey", "#E41A1C", "#377EB8", "#4DAF4A")),
                border = "black",
                column_title = "Side effect rank of all drugs",
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
                  labels = c("other","Afatinib", "Sapitinib", "Lapatinib"),
                  legend_gp = gpar(fill = c("grey","#E41A1C", "#377EB8", "#4DAF4A")), 
                  labels_gp = gpar(col = "#000000", fontsize = 12),
                  title_gp = gpar(fontsize = 12)
                )
                # width = unit(15, "cm"), height = unit(1.5, "cm")
) 

CairoPDF("./Figure/fig5c.pdf", width = 10, height = 1.8)
print(fig5c)
dev.off()



# 20231212 healthy Dr -----------------------------------------------------

healthy = list(LUAD = luad,
               BRCA = brca,
               SKCM = skcm)

healthy_lung = split(luad, luad$orig.ident)
healthy_breast = split(brca, brca$orig.ident)
healthy_skin = split(skcm, skcm$orig.ident)

healthy_skin_label = lapply(healthy_skin, function(h) {
  h = h[,seq(6,595,2)]
  copy_h = h
  print(dim(copy_h))
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

Dse = bind_cols(healthy_lung_label)

Dse = Dse[,c(1,2, seq(3,ncol(Dse),3))]
rownames(Dse) = paste(Dse$DRUG_ID...1, Dse$DRUG_NAME...2, sep = "_")
Dse = Dse[,-c(1,2)]

colnames(Dse) = names(healthy_lung_label)
Dse= scale(Dse, center = T, scale = TRUE)


col_anno = data.frame(name = c(paste(healthy_label$LUAD$DRUG_ID[1:5],healthy_label$LUAD$DRUG_NAME[1:5], sep = "_")))
col_labels = c(rep(" ", 295))
names(col_labels) = rownames(Dse)
col_labels[col_anno$name] = col_anno$name

CairoPDF("./Figure/supfig7/lung.pdf", width = 10, height = 3.2)
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

for (cancer in names(healthy_label)) {
  rank = healthy_label[[cancer]]
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

test = cbind(bind_cols(healthy_skin_label), bind_cols(healthy_breast_label), bind_cols(healthy_lung_label))

test = test[, c(1,2, seq(3,78,3))]
rownames(test) = paste(test$DRUG_ID...1, test$DRUG_NAME...2, sep = "_")
test = test[,-c(1,2)]

colnames(test) = c(paste("skin", 1:5, sep = ""), paste("brea", 1:8, sep = ""), paste("lung", 1:13, sep = ""))

CairoPDF("./Figure/supfig7/all.pdf", width = 10, height = 8)
Heatmap(t(test),
        # col = colorRamp2(c(0,1,2,3), c("grey", "#E41A1C", "#377EB8", "#4DAF4A")),
        column_title = "Drugs",
        column_title_gp = gpar(fontsize = 14),
        cluster_rows = F, 
        cluster_columns = F,
        column_dend_height = unit(0,"mm"),
        column_dend_reorder = F,
        show_row_names = T,
        # row_labels = paste("L", c(0:9, "a","b","c"), sep = ""),
        show_column_names = T,
        # column_labels = col_labels,
        column_names_gp = gpar(fontsize = 2),
        row_names_side = "left",
        row_title = "",
        show_heatmap_legend = T,
        name = "Side\neffect")
dev.off()
### point ####
fig3.data.2 = readRDS("../verify_data/GSE161529_HER2/scPharm_drug_rank.rds")
fig3.data.2 = fig3.data.2[c("AH0308", "MH0031", "MH0069", "PM0337")]
show.drugs = c("EGFR, ERBB2", 
               "EGFR, ERBB2, ERBB3", # target 2
               "EGFR",
               "MET, ALK, ROS1",
               "MET",
               "MET, KDR, TIE2, VEGFR3/FLT4, RON, PDGFR, FGFR1, EGFR",
               "IGF1R, IR",
               "IGF1R",
               "PDGFR, KIT, VEGFR",
               "VEGFR, RET, KIT, PDGFR",
               "VEGFR, FLT1, FLT2, FLT3, FLT4, KIT, PDGFRB",
               "FGFR1, FGFR2, FGFR3", # bypass pathway 12
               "MEK1, MEK2",
               "ERK1,ERK2",
               "ERK1, ERK2",
               "KRAS (G12C)",
               "MTOR",
               "JAK1, JAK2",
               "PI3Kalpha",
               "PI3K (beta sparing)",
               "PI3Kalpha, PI3Kdelta, PI3Kbeta, PI3Kgamma",
               "AKT1, AKT2",
               "AKT1, AKT2, AKT3",
               "BRAF", # downstream 24
               "CDK4, CDK6",
               "WEE1, CHEK1",
               "WEE1, PLK1",
               "PARP1, PARP2") # other target
show.drugid = drug_id[which(drug_id$PUTATIVE_TARGET %in% show.drugs),]
rownames(show.drugid) = show.drugid$DRUG_ID
show.drugid$group = "HER2"
show.drugid[which(show.drugid$PUTATIVE_TARGET %in% show.drugs[3:12]),"group"] = "bypass pathway"
show.drugid[which(show.drugid$PUTATIVE_TARGET %in% show.drugs[13:24]),"group"] = "downstream pathway"
show.drugid[which(show.drugid$PUTATIVE_TARGET %in% show.drugs[25:28]),"group"] = "other target"
# top_30_union = c()
# for (i in 1:4) {
#   rank = fig3.data.2[[i]]
#   top_30_union = c(top_30_union, rank$DRUG_ID[1:30])
# }
# top_30_union = unique(top_30_union)

hmp.data = list()
for (i in 1:4) {
  rank = fig3.data.2[[i]]
  rownames(rank) = rank$DRUG_ID
  rank = rank[rownames(show.drugid),c(1,2,5,6)]
  hmp.data[[i]] = rank
}
names(hmp.data) = names(fig3.data.2)
hmp.data = bind_cols(hmp.data)
hmp.data = hmp.data[,-c(5,6,9,10,13,14)]
hmp.score = hmp.data[,seq(3,10,2)]
hmp.score = data.frame(t(hmp.score), check.names = F)
# hmp.score = apply(hmp.score, 2, function(x) {
#   return((x-min(x))/(max(x)-min(x)))
# })
# hmp.score = t(data.frame(hmp.score))
# hmp.score = data.frame(hmp.score)
# hmp.score[5,] = colMeans(hmp.score)
# hmp.score = hmp.score[,order(hmp.score[5,], decreasing = T)]



drug_id = drug_id %>% 
  filter(!if_all(.fns = is.na))
rownames(drug_id) = drug_id$DRUG_ID

column_ha <- HeatmapAnnotation(
  group = show.drugid$group,
  which = "col",
  col = list(
    group = c("HER2" = "#B3E2CD",
              "bypass pathway" = "#FDCDAC",
              "downstream pathway" = "#CBD5E8",
              "other target" = "#F4CAE4")
  )
)
library(circlize)
fig3c = Heatmap(hmp.score,
                col = colorRamp2(seq(min(hmp.score), max(hmp.score), length = 3), c("#8470FF", "#EEEEEE", "red")),
                cluster_rows = F, 
                cluster_columns = F,
                show_row_names = F,
                column_labels = hmp.data$DRUG_NAME...2,
                column_names_gp = gpar(fontsize = 10),
                name = "score",
                top_annotation = column_ha,
                rect_gp = gpar(col = "white", lwd = 1)
                )
  
layout = "
ABCD##
EFGH##
"
fig3 = fig3a_1 + fig3a_2 + fig3a_3 + fig3a_4 + fig3b_1 + fig3b_2 + fig3b_3 + fig3b_4 +
  plot_layout(design = layout, guides = "collect", heights = unit(c(1.6, 1), c('cm', 'null')))

CairoPNG("./scPharm/fig3test.png", units = "in",
         width = 10, height = 3.3, dpi = 300)
print(fig3)
dev.off()

CairoPNG("./scPharm/fig3c_ha.png", units = "in",
         width = 9.2, height = 3, dpi = 300)
print(fig3c)
dev.off()


## Fig3d

fig3.data.3 = readRDS("../verify_data/GSE161529_ER/scPharm_fgsea.rds")
fig3.data.3 = fig3.data.3[c("MH0029-7C", "MH0043-T", "MH0151", "MH0173-T", "PM0360")]
for (i in 1:5) {
  meta.data = fig3.data.3[[i]]@meta.data
  meta.data = meta.data[,c("cell.label", "scPharm_nes_1925_GDC0810")]
  colnames(meta.data) = c("Group", "NES")
  p <- ggplot(meta.data, aes(NES, fill=Group, color=Group, linetype=Group))+
    geom_density(alpha=0.6, linewidth = 0.5)+
    scale_linetype_manual(values = c("dashed","solid"))+ # "dotted",
    scale_color_manual(values = c("#0000FF",'#FF0000'))+ # "#00FFFF",
    geom_rug()+
    # ggtitle(plot_title)+
    theme_bw()+
    # labs(tag = tag.list[i])+
    xlab(NULL)+
    ylab(NULL)+
    # ylim(0,ylimit)+
    # scale_y_continuous(breaks = 0:5, labels = c("0.0","1.0","2.0","3.0","4.0","5.0"))+
    theme(plot.title = element_text(hjust = 0.5,size = 14),
          plot.margin = unit(c(0.2,0.1,0.8,0.4),'cm'),
          # plot.tag = element_text(family = "sans", face = "bold", size = 26),
          # axis.text = element_text(colour = "#000000"),
          axis.text = element_text(family = "sans", face = "bold", size = 10, colour = "#000000"),
          # title = element_text(family = "sans", face = "bold", size = 14),
          axis.title = element_blank(),
          legend.text = element_text(family = "sans", face = "bold", size = 10),
          legend.title = element_text(family = "sans", face = "bold", size = 12),
          legend.position = "right"
    )
  assign(paste("fig3d", i, sep = "_"), p)
}

## Fig3e
kegg.gmt = getGmt("../kegg_pathway.gmt")
for (i in 1:14) {
  sample = fig3.data.3[[i]]
  assay = GetAssayData(sample, slot = "data")
  expr = assay[which(rownames(assay) %in% c("ESR1","ESR2","NCOA1","NCOA3","CCND1","MYC")),]
  ## ER下游MAPKc("ESR1","ESR2","SRC","HRAS","ARAF","MAP2K1","MAPK1","CREB3")
  ## 正常ER通路 c("ESR1","ESR2","NCOA1","NCOA2","NCOA3","BCL2","EBAG9","KRT19","CTSD","TFF1","PGR")
  expr = as.data.frame(t(expr))
  expr$mean_expr = rowMeans(expr)
  sample = AddMetaData(sample, expr[,'mean_expr'], col.name = "target_pathway_activity")
  viol_data = sample@meta.data[,c("cell.label", "scPharm_label_1925_GDC0810", "target_pathway_activity")]
  viol_data = viol_data[viol_data$cell.label == "tumor",c("scPharm_label_1925_GDC0810", "target_pathway_activity")]
  colnames(viol_data) = c("cell.label", "target_pathway_activity")
  # viol_data = viol_data[viol_data$cell.label != "other",]
  viol_data[viol_data$cell.label != "sensitive", "cell.label"] = "no.sensitive"
  # print(table(viol_data[,1]))
  p <- ggplot(viol_data, aes(x=cell.label,
                             y=as.numeric(target_pathway_activity),
                             fill = cell.label))+
    geom_violin(linewidth=0.5)+
    geom_boxplot(width=.2,col="black",fill="white")+
    scale_fill_manual(values = c("#FF0000","#0000FF"))+
    stat_compare_means(aes(group = cell.label), label = "p.format",
                       method = "wilcox.test",
                       label.x = 1.1, label.y = max(viol_data$target_pathway_activity)+0.1,
                       size=4)+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          legend.text = element_text(family = "sans", face = "bold", size = 10),
          legend.title = element_text(family = "sans", face = "bold", size = 12),
          axis.text.y = element_text(family = "sans", face = "bold", size = 10, colour = "#000000"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = 12, family = "sans", face = "bold", colour = "#000000"),
          text = element_text(family = "sans", face = "bold", colour = "#000000"),
          plot.margin = unit(c(0.2,0.1,0.2,0.4),'cm'),
    )+
    ylim(0, max(viol_data$target_pathway_activity)+0.2)+
    xlab(NULL)+
    ylab("Pathway activity")+
    ggtitle("")
  # assign(paste("fig3e",i,sep = "_"), p)
  CairoPNG(paste0("./scPharm/plot/er_pa/", names(fig3.data.3)[i],".png"), units = "in",
           width = 3, height = 2, dpi = 300)
  print(p)
  dev.off()
}
layout = "
ABCDE##
FGHIJ##
"
fig3de = fig3d_1 + fig3d_2 + fig3d_3 + fig3d_4 + fig3d_5 + fig3e_1 + fig3e_2 + fig3e_3 + fig3e_4 + fig3e_5 +
  plot_layout(design = layout, guides = "collect", heights = unit(c(2, 1), c('cm', 'null')))

CairoPNG("./scPharm/fig3de.png", units = "in",
         width = 12, height = 3.3, dpi = 300)
print(fig3de)
dev.off()

## Fig3f
fig3.data.4 = readRDS("../verify_data/GSE161529_ER/scPharm_drug_rank.rds")
fig3.data.4 = fig3.data.4[c("MH0029-7C", "MH0043-T", "MH0151", "MH0173-T", "PM0360")]
show.drugs = c("ESR1, ESR2", # target 
               "ESR",
               "ESR1",
               "AR", # bypass pathway 12
               "MEK1, MEK2",
               "ERK1,ERK2",
               "ERK1, ERK2",
               "KRAS (G12C)",
               "MTOR",
               "JAK1, JAK2",
               "PI3Kalpha",
               "PI3K (beta sparing)",
               "PI3Kalpha, PI3Kdelta, PI3Kbeta, PI3Kgamma",
               "AKT1, AKT2",
               "AKT1, AKT2, AKT3",
               "BRAF", # downstream 24
               "CDK4, CDK6",
               "WEE1, CHEK1",
               "WEE1, PLK1",
               "PARP1, PARP2",
               "EGFR, ERBB2", 
               "EGFR, ERBB2, ERBB3") # other target
show.drugid = drug_id[which(drug_id$PUTATIVE_TARGET %in% show.drugs),]
show.drugid$group = "ER"
show.drugid[which(show.drugid$PUTATIVE_TARGET %in% show.drugs[2:4]),"group"] = "bypass pathway"
show.drugid[which(show.drugid$PUTATIVE_TARGET %in% show.drugs[5:16]),"group"] = "downstream pathway"
show.drugid[which(show.drugid$PUTATIVE_TARGET %in% show.drugs[17:22]),"group"] = "other target"
rownames(show.drugid) = show.drugid$DRUG_ID

hmp.data = list()
for (i in 1:5) {
  rank = fig3.data.4[[i]]
  rownames(rank) = rank$DRUG_ID
  rank = rank[rownames(show.drugid),c(1,2,5,6)]
  hmp.data[[i]] = rank
}
names(hmp.data) = names(fig3.data.4)
hmp.data = bind_cols(hmp.data)
hmp.data = hmp.data[,-c(5,6,9,10,13,14,17,18)]
hmp.score = hmp.data[,seq(3,12,2)]
hmp.score = data.frame(t(hmp.score), check.names = F)
hmp.score = apply(hmp.score, 1, function(x) {
  return((x-min(x))/(max(x)-min(x)))
})
hmp.score = t(hmp.score)
hmp.score = data.frame(hmp.score, check.names = F)



column_ha <- HeatmapAnnotation(
  group = show.drugid$group,
  which = "col",
  col = list(
    group = c("ER" = "#B3E2CD",
              "bypass pathway" = "#FDCDAC",
              "downstream pathway" = "#CBD5E8",
              "other target" = "#F4CAE4")
  )
)
library(circlize)
fig3f = Heatmap(hmp.score,
                col = colorRamp2(seq(min(hmp.score), max(hmp.score), length = 3), c("#8470FF", "#EEEEEE", "red")),
                cluster_rows = F, 
                cluster_columns = F,
                show_row_names = F,
                column_labels = hmp.data$DRUG_NAME...2,
                column_names_gp = gpar(fontsize = 10),
                name = "score",
                top_annotation = column_ha,
                rect_gp = gpar(col = "white", lwd = 1)
)
CairoPNG("./scPharm/fig3f.png", units = "in",
         width = 9.2, height = 3, dpi = 300)
print(fig3f)
dev.off()