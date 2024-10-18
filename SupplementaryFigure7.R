# Title    : Supplementary Figure 7
# Author   : Wanglab
# Time     : 2024.7

library(tidyr)
library(ggplot2)
library(ggalluvial)
library(ggthemes)
library(Cairo)
library(Seurat)

sample = c("MH0001","MH0025","MH0029-7C","MH0029-9C","MH0032","MH0040","MH0042","MH0043-T","MH0064-T",
           "MH0151","MH0163","MH0173-T","PM0360"  )
names(sample) = sample

# SuppleFigure 7a ---------------------------------------------------------------
rank.table<-read.csv("ER.multi-evaluate.ranktable.20231113.csv",row.names = 1,check.names = F)

df<-c()
druglist<-c("Fulvestrant_Id1816","Fulvestrant_Id1200")

for (drug in druglist) {
  for (n in colnames(rank.table)) {
    print(n)
    name<-strsplit(n,split = "_")[[1]]
    sample<-name[1]
    method<-name[2]
    order<-grep(drug,rank.table[,n])
    if(length(order)==1){
      df<-rbind(df,c(drug,method,sample,order))
    }
  }
}
df<-data.frame(df)
colnames(df)<-c('Drug','Method','Sample','Order')
df$Order<-df$Order%>%as.numeric()
df$Method<-factor(df$Method,levels = c("SeuratCCA", "Scissor", "scDEAL", "CaDRReS-Sc", "scPharm"))
df$Score<-1/df$Order


df$Max_conc<-"Low_conc"
df[df$Drug=="Fulvestrant_Id1200",]$Max_conc<-"High_conc"
df$Interaction <- interaction(df$Method, df$Max_conc)
df$Interaction<-factor(df$Interaction,levels =c("SeuratCCA.High_conc","SeuratCCA.Low_conc" ,  "Scissor.High_conc" ,"Scissor.Low_conc" ,     "scDEAL.High_conc",  "scDEAL.Low_conc" ,  
                                                "CaDRReS-Sc.High_conc", "CaDRReS-Sc.Low_conc" , "scPharm.High_conc","scPharm.Low_conc" ) )

## scPharm
df<-df%>%filter(Method=="scPharm")
df$Max_conc<-"High_conc"
df[df$Drug=="Fulvestrant_Id1200",]$Max_conc<-"Low_conc"
df$Max_conc<-factor(df$Max_conc,levels = c("Low_conc","High_conc"))

p<-ggplot(df, aes(x = Max_conc, y = Order, fill = Max_conc)) +
  stat_boxplot(geom="errorbar",width = 0.2, position = position_dodge(width = 0.75)) +
  geom_boxplot(outlier.size = 0,outlier.alpha = 0,width = 0.7, position = position_dodge(width = 0.75)) +
  geom_jitter(size=1.3,alpha=0.5,width = 0.2) +
  
  geom_signif(comparisons = list(c("Low_conc","High_conc")
  ),
  y_position = c(280),
  map_signif_level = T,tip_length = 0.01,vjust = 0.8)+
  
  labs(x = "Drug concentration",
       y = "Order(1-295)") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        plot.margin = unit(c(1,0.5,0.5,0.5),"cm"),title = element_text(size = 10))+
  ggtitle("fulvestrant.ER-positive sample")+
  scale_fill_manual(values = c("#5272A8","#5CAB6F"))

ggsave(filename = "Fig.s7a.pdf",
       plot =p,width =3.5,height = 6)
	   
	   
# SuppleFigure 7b ---------------------------------------------------------------
SuppleFigure7b <- function (sample) {
  data.1 = lapply(sample, function(sam) {
    # Load the single-cell pharmacology result object from the specified path
    rank = readRDS(paste0("./scPharm/result/", sam, "_scPharm_object_nmcs_50_nfs_200.rds"))
    return(rank)
  })
  # Loop through each element in the data.1 list to generate a PDF chart for each
  for (sam in names(data.1)) {
  meta.data = data.1[[sam]]@meta.data
  meta.data = meta.data[meta.data$cell.label == "tumor", c("scPharm_label_1200_Fulvestrant", "scPharm_label_1816_Fulvestrant")]
  meta.data = meta.data[order(meta.data$scPharm_label_1200_Fulvestrant),]
  colnames(meta.data) = c("Fulvestrant_low","Fulvestrant_high")
  meta.data = meta.data %>% gather(key = "dosage", value = "cluster")
  meta.data$cell = c(rep(1:(nrow(meta.data)/2),2))
  CairoPDF(paste("./Figure/supfig7b/alluvial_", sam,"_1112.pdf", sep = ""), width = 2.2, height = 4)
  alluvial = ggplot(meta.data,
                    aes(x = dosage, stratum = cluster, alluvium = cell, y =cell,
                        fill = cluster, label = NULL)) +
    geom_col(width = 0.3,color=NA) +
    scale_x_discrete(limits = c("Fulvestrant_low","Fulvestrant_high"), expand = c(.1, .1)) +
    geom_flow(alpha = 0.6, linewidth = 0) +
    geom_stratum(alpha = 1, width = .3, color = "white", linewidth = 0) +
    scale_fill_manual(values = c("#B5B5B5", "#FF0000", "#0000FF")) +
    # geom_text(stat = "stratum", size = 8) +
    theme_map()+
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 13),
          axis.text.x = element_text(hjust = 0.1, size = 12, angle = -60)) +
    ggtitle(sam)
  print(alluvial)
  dev.off()
  }
}
# statistic test for SuppleFigure 7b
c.switch = function(label) {
  # transform the label to 0,1,2
  if (label == "other") {
    return(0)
  }else if (label == "resistant") {
    return(1)
  }else {
    return(2)
  }
}
STforSF7b <- function (sample) {
  data.list = lapply(sample, function(sam) {
    # Load the single-cell pharmacology result object from the specified path
    rank = readRDS(paste0("./scPharm/result/", sam, "_scPharm_object_nmcs_50_nfs_200.rds"))
    return(rank)
  })
  signif = c()
  for (sam in names(data.list)) {
    meta.data = data.list[[sam]]@meta.data
    meta.data = meta.data[meta.data$cell.label == "tumor",c("scPharm_label_1200_Fulvestrant", "scPharm_label_1816_Fulvestrant")]
    meta.data = meta.data[order(meta.data$scPharm_label_1200_Fulvestrant),]
    meta.data = meta.data[(meta.data$scPharm_label_1200_Fulvestrant != "sensitive" | meta.data$scPharm_label_1816_Fulvestrant != "sensitive"),]
    sample = sample(1:nrow(meta.data), nrow(meta.data))
    meta.data$H0 = meta.data$scPharm_label_1200_Fulvestrant[sample]
    meta.data = meta.data[meta.data$scPharm_label_1200_Fulvestrant == 'other',]
    test.data = data.frame(matrix(0, 1853, 2))
    colnames(test.data) = c("H0","H1")
    for (i in 1:nrow(meta.data)) {
      test.data[i, "H1"] = c.switch(meta.data$scPharm_label_1816_Fulvestrant[i])
      test.data[i, "H0"] = c.switch(meta.data$H0[i])
    }
    pval = wilcox.test(test.data$H1, test.data$H0, alternative = "greater", paired = T)$p.value
    signif = c(signif, pval)
  }
  return(signif)
}

SuppleFigure7b(sample)
