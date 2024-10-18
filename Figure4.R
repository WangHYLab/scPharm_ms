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


# figure 4g-i -----------------------------------------------------------------------

source("./functions.R")
#### Fig.4 G ####
setwd("./data/gse158677")
# preprocess sc-data
sc_process2<-function(data,filter=TRUE){
  if(filter){
    data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
    data <- subset(data,nFeature_RNA>200&nFeature_RNA<4000&percent.mt<15)
  }
  data <- NormalizeData(data,verbose = FALSE)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
  data <- ScaleData(data,features=VariableFeatures(data),verbose = FALSE) 
  pc<-20
  data <- RunPCA(data, npcs =pc,verbose = FALSE)
  data <- RunUMAP(data, features = VariableFeatures(data) ,verbose = FALSE)
  data <- FindNeighbors(data, dims = 1:pc, verbose = FALSE)
  data <- FindClusters(data, verbose = FALSE,resolution=0.5)
  return(data)

}

data.list<-list()
for (name in c('T1','T2','T3','T4','T5')) {
  data<-Read10X(data.dir = paste0("./data/gse158677/",name),gene.column = 1)
  data<-CreateSeuratObject(counts = data,project = name,min.cells = 3)
  print(dim(data))
  data<-sc_process2(data)
  data.list[[name]]<-data
  print(dim(data))
}

# cell type identification
ref.se<-ref$ImmGen
matrix <- GetAssayData(data, slot = 'data')
singler.cluster <-SingleR(matrix,ref = ref.se,
                          labels = ref.se$label.main,
                          clusters = data@meta.data$seurat_clusters)

data<-assign_singler(data)

DotPlot(data.list$T5,feature=c("PyMT","Epcam","Krt8","Ptprc","Pdgfra","Ms4a1","Pecam1"))

# merge used celltype
data.list$T1<-subset(data.list$T1,seurat_clusters %in%c(0,1,5))
data.list$T2<-subset(data.list$T2,seurat_clusters %in%c(0,2,3,5))
data.list$T3<-subset(data.list$T3,seurat_clusters %in%c(0,1,2,3,4,6,9))
data.list$T4<-subset(data.list$T4,seurat_clusters %in%c(0,1,7,8))
data.list$T5<-subset(data.list$T5,seurat_clusters %in%c(1,2,10))

## trans homolog gene
datasets<-data.list
for (name in names(datasets)) {
  data<-datasets[[name]]
  mtx<-GetAssayData(data,assay = "RNA",slot = "counts")
  rownames(mtx)<-mus2hsa(rownames(mtx))
  
  seu<-CreateSeuratObject(counts = mtx,project = name,min.cells = 3)
  seu<-sc_process2(seu,filter = FALSE)
  datasets[[name]]<-seu
}

## Run scPharm 
setwd("/home/DATA/zhengjie/scPham/benchmark")
out_gse158677<-benchmark_cal(sc_dataset = datasets,gdsc = gdsc,func = c("scPharm"),
                           scPharm.tissue.or.cellline = "cellline",drug.choose = NULL,process = FALSE,sampling=c("all"),scPharm.type = "BRCA",times=1)
saveRDS(out_gse158677,"./output/gse158677.rank&comb.20231212.rds")
## extract rank info
rank.all<-out_gse158677$rank
rank.table<-c()
for (rank in rank.all) {
  drug<-paste0(rank$DRUG_NAME,"_Id",rank$DRUG_ID)  
  rank.table<-cbind(rank.table,drug)
}
rank.table<-data.frame(rank.table)
colnames(rank.table)<-names(datasets)
write.csv(rank.table,file = "./output/gse158677.ranktable.1212.csv")

## plot
data<-read.csv("./output/gse158677.ranktable.1212.csv",row.names = 1,check.names = F)

drug<-c("ALPELISIB","TASELISIB",  #PI3K
        "TAMOXIFEN","FULVESTRANT",# ESR1
        "IPATASERTIB","MK-2206","UPROSERTIB","AFURESERTIB",# ATK1
        "LAPATINIB","Afatinib","Sapitinib",# ERBB2
        "Gemcitabine", # POLD / POLE
        "Lapatinib","Gefitinib","Erlotinib", #EGFR
        "Palbociclib","Ribociclib", #CDK6 / CDK4
        "Docetaxel","Paclitaxel", #TUBB4B
        "Vinorelbine" # TUBB2B
)
GDSC2<-data$T1
select_drug<-c()
for (d in drug) {
  in_<-GDSC2[grep(d,GDSC2,ignore.case = T,perl = T)]
  select_drug<-c(select_drug,in_)
}
select_drug<-select_drug%>%unique()
d_df<-data.frame(drug=select_drug,index=c(1:length(select_drug)))
rownames(d_df)<-d_df$drug
# write.csv(d_df,"drug-order.gse.csv")
# use all drug
library(paletteer)

df<-c()
for (col in colnames(data)) {
  d<-data[col][,1]
  
  s<-d %in% d_df$drug
  s[s==FALSE]<-0
  s[d %in% d_df$drug]<-d_df[d[d %in% d_df$drug],]$index
  s<-as.numeric(s)
  df<-cbind(df,s)
}
df<-data.frame(df)
colnames(df)<-colnames(data)
library(purrr)
df<-map_df(df, as.numeric)%>%as.matrix()
rownames(df)<-rownames(data)
# sort

pdf("./output/Fig.4g.pdf",width = 10,height = 2)
pheatmap::pheatmap(df%>%t(),cluster_rows = F,cluster_cols = F,
                   show_colnames = T,border_color = "black",
                   
                   color = c("grey80",rainbow(22))
                   # width = 2, 
)
dev.off()

o<-c()
for (d in d_df$drug) {
  index<-c()
  for (l in data) {
    
    index<-c(index,grep(d,l))
  }
  
  o<-rbind(o,c(d,paste0(index,collapse = ","),mean(index)))
}

write.csv(o,file = "./output/order.gse158677.csv")
#### Fig.4 H-I ####
library(ggplot2)
library(tidyverse)
library(openxlsx)

data<-read.xlsx("./data/Mohr &  Yoldi study.xlsx",sheet = "Mohr_1",check.names = F)
data$`concentration[nM]`<-factor(data$`concentration[nM]`)
p1<-ggplot(data,aes(x=`concentration[nM]`,y=`relative.Ki-67.index`))+
  geom_errorbar(aes(ymin=`relative.Ki-67.index`-se,ymax=`relative.Ki-67.index`+se),width=0.2,linewidth=0.2)+
  geom_bar(stat = "identity",width = 0.6)+
  ylim(0,1.2)+
  ylab("Relative Ki-67 index")+
  xlab("Concentration [nM]")+
  ggtitle("Mohr's study",subtitle = "Docetaxel")+
  theme_classic()+
  theme(plot.margin =margin(1,1,1,1,unit = "cm") )

data<-read.xlsx("./data/Mohr &  Yoldi study.xlsx",sheet = "Mohr_2",check.names = F)
data$`concentration[nM]`<-factor(data$`concentration[nM]`)
p2<-ggplot(data,aes(x=`concentration[nM]`,y=`relative.proliferation.rate`))+
  geom_errorbar(aes(ymin=`relative.proliferation.rate`-se,ymax=`relative.proliferation.rate`+se),width=0.2,linewidth=0.2)+
  geom_bar(stat = "identity",width = 0.6)+
  ylim(0,1.2)+
  ylab("Relative proliferation rate")+
  xlab("Concentration [nM]")+
  ggtitle("Mohr's study",subtitle = "Docetaxel")+
  theme_classic()+
  theme(plot.margin =margin(1,1,1,1,unit = "cm") )

data<-read.xlsx("./data/Mohr &  Yoldi study.xlsx",sheet = "Yoldi",check.names = F)
data$Days<-factor(data$Days)
p3<-ggplot(data,aes(x=Days,y=Tumor.volume,color=group,shape=group,group=group))+
  geom_errorbar(aes(ymin=Tumor.volume-se,ymax=Tumor.volume+se),width=0.3,linewidth=0.5)+
  geom_point(size=3)+
  geom_line()+
  # ylim(0,1.4,)+
  ylab("Tumor volume")+
  xlab("Days")+
  ggtitle("Yoldi's study")+
  theme_classic()+
  theme(plot.margin =margin(1,1,1,1,unit = "cm") )+
  expand_limits(y=c(0,1.4))+
  scale_y_continuous(n.breaks = 10 )+
  scale_color_manual(values = c("grey20","darkred"))

ggsave(plot = p1,filename = "./output/Fig.4h-1.docetaxel.pdf",width = 3,height = 3)
ggsave(plot = p2,filename = "./output/Fig.4h-2.docetaxel.pdf",width = 2.5,height = 3)
ggsave(plot = p3,filename = "./output/Fig.4i.docetaxel.pdf",width = 6,height = 3.5)
