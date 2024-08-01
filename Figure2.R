# Title    : Figure 2
# Author   : Wanglab
# Time     : 2024.7

library(Seurat)
library(readxl)
library(ggplot2)
library(ggpubr)
library(Cairo)
library(dplyr)
library(ggsci)

# Figure 2a, 2b----------------------------------------------------------------------

Figure2a_b <- function(data.sc, data.gdsc) {
  # boxplot and line plot depicting the correlation between NES and AUC about docetaxel in LUAD cell lines single cell data
  meta.data = data.sc@meta.data
  drug_id <- data.gdsc[data.gdsc$TCGA_DESC == "LUAD",]
  drug_id <- drug_id[,c(8,9,10,11)]
  drug_id <- drug_id[!duplicated(drug_id$DRUG_ID),]
  drug_id = drug_id %>%
    filter(!if_all(.fns = is.na))
  drug_id = drug_id[drug_id$DRUG_NAME == "Docetaxel",]
  all.corr = c() # all drug correlation
  for (i in 1:nrow(drug_id)) {
    # if (i == 5) {
    #   break
    # }
    # x is the drug data about docetaxel in LUAD cell lines
    x = data.gdsc[(data.gdsc$DRUG_ID == drug_id$DRUG_ID[i] & data.gdsc$TCGA_DESC == "LUAD"), c("CELL_LINE_NAME", "AUC")]
    x = na.omit(x)
    x = x[order(x[,2]),]
    colnames(x) = c("cellline", "AUC")
    x$cellline = sapply(x$cellline, function(name) {return(toupper(gsub("-","",name)))})
    # y is the single cell data associated with cell lines in x
    y = meta.data[, c("orig.ident",paste("scPharm_nes", drug_id$DRUG_ID[i], drug_id$DRUG_NAME[i], sep = "_"))]
    colnames(y) = c("orig.ident", "NES")
    y$cellline = sapply(strsplit(y$orig.ident, split ="_"), function(x) {return(x[[1]][1])})
    show.cl = intersect(unique(x$cellline), unique(y$cellline))
    if (length(show.cl) <=3) {
      next
    }
    x = x[which(x$cellline %in% show.cl),]
    y = y[which(y$cellline %in% show.cl),]
    an.corr = x
    an.corr$median.nes = 0
    # calculate the median NES of cells from each cell line
    for (cl in 1:nrow(an.corr)) {
      cellline = an.corr$cellline[cl]
      an.corr[cl,"median.nes"] = median(y[y$cellline == cellline, "NES"])
    }
    # print(an.corr)
    # calculate the correlation between NES and AUC
    corr = cor(an.corr$AUC, an.corr$median.nes)
    all.corr = c(all.corr, corr)
    x$AUC = x$AUC*8 - 4
    message("1:boxplot ploting")
    p <- ggplot(y) +
      geom_boxplot(aes(x = reorder(cellline, NES, FUN = median), y = NES, fill = "grey"),
                   position=position_dodge(0.4),
                   width=0.4,
                   outlier.color = "#999999",
                   outlier.size = 0.8,
                   linewidth = 0.3) +
      theme_classic()+
      ggtitle(paste(drug_id[drug_id$DRUG_ID == drug_id$DRUG_ID[i], 2], "\nCorr = ", round(corr, 2), sep = "")) +
      scale_y_continuous(limits = c(-4,4),
                         sec.axis = sec_axis(~.,
                                             name = 'AUC',
                                             breaks = c(-4,0,4),
                                             labels = c(0,0.5,1))) +
      geom_point(data = x, aes(x = cellline, y = AUC, colour = "blue")) +
      geom_line(data = x, aes(x = cellline, y = AUC, group = 1, colour = "blue")) +
      scale_fill_manual(values = "grey"
                          #direction = 1
                          ) +
      theme(plot.title = element_text(hjust = 0.5, size = 12),
            axis.text.y = element_blank(), #text(hjust = -1, size = 10, colour = "#000000"),
            axis.title.x = element_blank(),
            # axis.ticks.x = element_blank(),
            axis.text.x = element_text(angle = -90, size = 12, colour = "#000000"),
            axis.title.y = element_blank(), #text(size = 12, colour = "#000000"),
            text = element_text(colour = "#000000"),
            plot.margin = unit(c(0.2,0.2,0.4,0.2),'cm'),
            legend.position = "none",
            legend.title = element_blank(),
            legend.text = element_text(size = 12)
      )
    CairoPDF("./Figure/fig2a.pdf", width = 7.5, height = 4)
    print(p)
    dev.off()
  }
  # Figure 2b boxplot of all correlation
  all.corr = data.frame(corr = all.corr)
  message(nrow(all.corr))
  p3 = ggplot(all.corr, aes(y = corr)) +
    geom_boxplot(width=0.4) +
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.3,size = 12),
          axis.text.y = element_text(size = 12,colour = "#000000"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = 12),
          plot.margin = unit(c(0.3,0.3,0.3,0.3),'cm'),
          legend.position = "none"
    ) +
    ylab(paste0("Corr\nAll drug(n=",nrow(all.corr), ")")) +
    ylim(c(-1,1))
  CairoPDF("./Figure/allcorr.pdf", width = 1.8, height = 4)
  print(p3)
  dev.off()
  # return(all.corr)
}
scp.data = readRDS("../scPharm/result/SCP542_LUAD_scPharm_object_nmcs_50_nfs_200.rds")
gdsc.data = read_xlsx("../scPharm/info/GDSC2_fitted_dose_response_24Jul22.xlsx")
Figure2a_b(scp.data, gdsc.data)

# Figure 2c----------------------------------------------------------------------

calculate_cor = function(data.sc, data.gdsc, cancer_type) {
  meta.data = data.sc@meta.data
  drug_id <- data.gdsc[data.gdsc$TCGA_DESC == cancer_type,]
  drug_id <- drug_id[,c(8,9,10,11)]
  drug_id <- drug_id[!duplicated(drug_id$DRUG_ID),]
  drug_id = drug_id %>%
    filter(!if_all(.fns = is.na))
  all.corr = c()
  for (i in 1:nrow(drug_id)) {
    x = data.gdsc[(data.gdsc$DRUG_ID == drug_id$DRUG_ID[i] & data.gdsc$TCGA_DESC == cancer_type), c("CELL_LINE_NAME", "AUC")]
    x = na.omit(x)
    x = x[order(x[,2]),]
    colnames(x) = c("cellline", "AUC")
    x$cellline = sapply(x$cellline, function(name) {return(toupper(gsub("-","",name)))})
    if (!paste("scPharm_nes", drug_id$DRUG_ID[i], drug_id$DRUG_NAME[i], sep = "_") %in% colnames(meta.data)) {
      print("no scPharm")
      next
    }
    y = meta.data[, c("orig.ident",paste("scPharm_nes", drug_id$DRUG_ID[i], drug_id$DRUG_NAME[i], sep = "_"))]
    colnames(y) = c("orig.ident", "NES")
    y$cellline = sapply(strsplit(y$orig.ident, split ="_"), function(x) {return(x[[1]][1])})
    show.cl = intersect(unique(x$cellline), unique(y$cellline))
    # print(show.cl)
    if (length(show.cl) <=3) {
      print("cell line lack")
      next
    }
    x = x[which(x$cellline %in% show.cl),]
    y = y[which(y$cellline %in% show.cl),]
    an.corr = x
    an.corr$median.nes = 0
    for (cl in 1:nrow(an.corr)) {
      cellline = an.corr$cellline[cl]
      an.corr[cl,"median.nes"] = median(y[y$cellline == cellline, "NES"])
    }
    corr = cor(an.corr$AUC, an.corr$median.nes)
    all.corr = c(all.corr, corr)
  }
  if(length(all.corr) != 0) {
    all.corr = data.frame(corr = all.corr, type = cancer_type)
    return(all.corr)
  }else {
    return(NULL)
  }
}
# Load all cell lines single cell data
SCPFILES = list.files(path = "../scPharm/result", pattern = "SCP542.*nfs_200", full.names = T)
names(SCPFILES) = sapply(strsplit(SCPFILES, split = "_"), function(x) {return(x[2])})
scp.scPharm = lapply(SCPFILES, function(path) {
  data = readRDS(path)
})
# calculate correlation of all drugs in all cancer cell lines
cancer_corr = map2(scp.scPharm, names(scp.scPharm), function(object, name) {
  result = calculate_cor(data.sc = object, data.gdsc = gdsc.data, cancer_type = name)
  print(paste0(name, nrow(result)))
  return(result)
})
cancer_corr = cancer_corr[!sapply(cancer_corr, is.null)]
data = bind_rows(cancer_corr)
# plot
Figure2c <- ggplot(data) +
  stat_boxplot(mapping=aes(x=reorder(type, corr, FUN = median),y=corr),
               geom ="errorbar",                             ##添加箱子的bar为最大、小值
               width=0.15,position=position_dodge(0.4))+
  geom_boxplot(aes(x = type, y = corr, fill = type),
               position=position_dodge(0.4),                 ##因为分组比较，需设组间距
               width=0.6,                                    ##箱子的宽度
               outlier.color = "#999999",
               outlier.size = 0.8,
               linewidth = 0.3) +
  stat_summary(mapping=aes(x=type,y=corr,group=type),                    ##分组计算的变量
               fun="median",                                   ##箱线图添加均值
               geom="point",shape=23,size=2,fill="white",    ##均值图形的设置
               position=position_dodge(0.4))+                ##因为分组比较，需设置两组间距
  scale_fill_manual(values = c("#D4C2AD","#D7A184","#EFDFCC","#BFC7D9","#66796B",
                                "#CDECEF","#DBCDA4","#6B4945","#c1b9b8","#9CBCB7",
                                "#CDD5C6","#B7BCBF")) +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(hjust = -0.5, size = 10, colour = "#000000"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12, colour = "#000000"),
        text = element_text(colour = "#000000"),
        plot.margin = unit(c(0.2,0.2,0.4,0.1),'cm'),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 12)
  ) +
  ylab("corr")
CairoPDF("./Figure/fig2c.pdf", width = 12, height = 3.8)
print(Figure2c)
dev.off()

# Figure 2d ---------------------------------------------------------------------
Figure2d <- function(fig2.data.4, fig2.data.5, cancer_type) {
  # fig2.data4: cancer cell line single cell data
  # fig2.data5: gdsc info
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
    x$cellline = sapply(x$cellline, function(name) {return(toupper(gsub("-","",name)))})
    x$label = "other"
    x1 = x
    # group x by AUC median
    x1[x1$AUC >= median(x1$AUC), 3] = "AUC_high.50%"
    x1[x1$AUC < median(x1$AUC), 3] = "AUC_low.50%"
    message(paste0(drug_id$DRUG_ID[i], ":x done"))
    if (!paste("scPharm_nes", drug_id$DRUG_ID[i], drug_id$DRUG_NAME[i], sep = "_") %in% colnames(meta.data)) {
      print("no scPharm")
      next
    }
    # y is scPharm_nes of cells from cell lines in x1
    y = meta.data[, c("orig.ident",paste("scPharm_nes", drug_id$DRUG_ID[i], drug_id$DRUG_NAME[i], sep = "_"))]
    colnames(y) = c("orig.ident", "response")
    y$cellline = sapply(strsplit(y$orig.ident, split ="_"), function(x) {return(x[[1]][1])})
    y$group = "other"
    y[which(y$cellline %in% x1[x1$label == "AUC_high.50%",]$cellline),"group"] = "AUC_high.50%"
    y[which(y$cellline %in% x1[x1$label == "AUC_low.50%",]$cellline),"group"] = "AUC_low.50%"
    y = y[y$group != "other",]
    message(paste0(drug_id$DRUG_ID[i], "y done"))
    print(paste(sum(y$group == "AUC_low.50%"), sum(y$group == "AUC_high.50%"), sep = "---"))
    if (sum(y$group == "AUC_low.50%") < 3 | sum(y$group == "AUC_high.50%") < 3) {
      print("next")
      next
    }
    # statistics test
    p.val1 = wilcox.test(y[y$group=="AUC_low.50%",]$response, y[y$group=="AUC_high.50%",]$response,
                         alternative = "less") # x compare y
    p.val1 = p.val1$p.value
    scatter.data[nrow(scatter.data)+1,] = c(drug_id$DRUG_ID[i],
                                            drug_id$DRUG_NAME[i],
                                            median(y[y$group == "AUC_low.50%",]$response),
                                            median(y[y$group == "AUC_high.50%",]$response),
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
  # plot
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
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 14),
          axis.title.y = element_text(size = 14, angle = 90),
          legend.text = element_text(size = 12),
          legend.title =  element_blank(),
          legend.position = c(0.76, 0.12),
          legend.background = element_rect(fill = rgb(1,1,1, alpha = 0.001), colour = NA)
    )+
    xlab("NES (AUC_low)")+
    ylab("NES (AUC_high")+
    scale_color_manual(values = c("#6A5ACD", "#FD8D3C")) +
    guides(color=guide_legend(override.aes = list(size=2))) +
    ggtitle(paste0(cancer_type, "\n(# of drug : ",nrow(scatter.data),")"))
  CairoPDF(paste0("./Figure/scatter/fig2d_", cancer_type, ".pdf"), width = 3.4, height = 3.66)
  print(p2)
  dev.off()
  # return(scatter.data)
}

for (i in 1:14) {
  print(names(scp.scPharm)[i])
  Figure2d (fig2.data.4 = scp.scPharm[[i]],
               fig2.data.5 = gdsc.data,
               cancer_type = names(scp.scPharm)[i])
}