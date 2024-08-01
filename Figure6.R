
library(tidyverse)
library(ggplot2)
library(openxlsx)
 
# ER comb data 
data<-readRDS("data/scPharm_targeted_drug_combo1214.rds")

#### Fig.6b ####
df<-c()
for (i in names(data)) {
  for (j in names(data$AH0319)) {
    d<-data[[i]][[j]]
    d$sample<-i
    df<-rbind(df,d)
  }
}
df.all<-df

## booster effect 
df<-df.all%>%filter(Strategy=="booster effects")
df.freq<-table(df$DRUG_NAME,df$DRUG_FIRST)%>%data.frame()
colnames(df.freq)<-c("Drug-2","Drug-1","Freq")
# sort index 
df1<-df.freq[df.freq$`Drug-1`=="GDC0810",]
df2<-df.freq[df.freq$`Drug-1`=="Fulvestrant",]
df1$sum<-df1$Freq+df2$Freq
df.freq$`Drug-2`<-factor(df.freq$`Drug-2`,levels = df1%>%arrange(desc(sum))%>%pull(`Drug-2`))
# plot
p<-ggplot(df.freq,aes(x=`Drug-2`,y=Freq,fill=`Drug-1`))+
  geom_bar(stat="identity",position =position_dodge(0.75),width = 0.6)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45,vjust = 0.5,hjust = 0),
        plot.margin =margin(1,1,1,1,unit = "cm"))+
  ylim(0,14)+
  ggtitle("Booster effects")+
  ylab("Freq. in ER+ samples")+
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 2,"Paired"))
ggsave(plot = p,filename = "output/Fig.6b.pdf",width = 6.5,height = 3.5)

#### Fig.6c ####
# read data
data<-read.xlsx("GDSC/GDSC2_fitted_dose_response_24Jul22.xlsx")
# select drug
df<-data%>%filter(DRUG_ID %in% c(1816, # Fulvestrant
                                 1925, # GDC0910
                                 1042, # Doramapimod
                                 1632, # Ribociclib
                                 1054, # Palbociclib
                                 1401, #AZD5438
                                 1014, #Refametinib
                                 2096  #VX-11e
                                 ))%>%dplyr::select(SANGER_MODEL_ID,MAX_CONC,TCGA_DESC,DRUG_NAME,LN_IC50,AUC)
# get information
ccl2tcga<-df%>%dplyr::select(SANGER_MODEL_ID,TCGA_DESC)%>%unique.data.frame(.)%>%tibble::remove_rownames()%>%tibble::column_to_rownames(colnames(.)[1])   
drug2maxconc<-df%>%dplyr::select(DRUG_NAME,MAX_CONC)%>%unique.data.frame(.)%>%tibble::remove_rownames()%>%tibble::column_to_rownames(colnames(.)[1])   
drug2maxconc<-log(drug2maxconc)

pdf<-df%>%dplyr::select(SANGER_MODEL_ID,DRUG_NAME,LN_IC50)%>%unique.data.frame()

d<-spread(pdf,"SANGER_MODEL_ID",'LN_IC50')%>%tibble::column_to_rownames(colnames(.)[1])%>%t()  %>%data.frame()
d<-d[,c("Fulvestrant","GDC0810","Doramapimod","Ribociclib","Palbociclib","AZD5438","Refametinib","VX.11e")]
colnames(d)<-c("Fulvestrant","GDC0810","Doramapimod","Ribociclib","Palbociclib","AZD5438","Refametinib","VX.11e")

p.list<-list()
for (d1 in c("Fulvestrant","GDC0810")) {
  for (d2 in c("Doramapimod","Ribociclib","Palbociclib","AZD5438","Refametinib","VX.11e")) {
    
    c<-cor.test(d[,d1],d[,d2])
    p<-c$p.value%>%as.numeric()%>%signif(.,4)
    cor<-c$estimate%>%as.numeric()%>%signif(.,4)
    p<-ggplot(d,aes_string(y=d1,x=d2))+
      geom_point(size=2,color="grey10")+
      geom_vline(xintercept=drug2maxconc[d2,],linetype="dashed",color="grey")+
      geom_hline(yintercept=drug2maxconc[d1,],linetype="dashed",color="grey")+
      annotate("text",x=0,y=drug2maxconc[d1,]-0.1,label="max.conc",color="navy",size=3)+
      annotate("text",x=drug2maxconc[d2,]+0.3,y=0,label="max.conc",color="navy",size=3)+
      theme_classic()+
      ggtitle(label =paste0(d1," - ",d2),subtitle = paste0('Pearson.pval = ',p,"\nPearson.corr = ",cor))+
      theme(plot.margin =margin(1,1,1,1,unit = "cm"))+
      xlab(paste0(d2," (LN_IC50)"))+
      ylab(paste0(d1," (LN_IC50)"))
    p.list[[paste0(d1,"-",d2)]]<-p
  }
}
# function for plot 
integratePlot<-function(fig,ncol){
  library(cowplot)
  text<-"g<-plot_grid("
  for (i in c(1:length(fig))) {
    text<-paste0(text,"fig[[",i,"]],")
  }
  text<-paste0(text,"ncol = ",ncol,")")
  g<-eval(parse(text=text))
  return(g)
}
# save pic
ggsave(plot=integratePlot(p.list,3),filename = "output/Fig.6c.pdf",width = 13,height = 20)

#### Fig.6d ####
# plot drug Fulvestrant+Palbociclib 
data <- read.xlsx("data/Formisano study.xlsx")
data$`Treatment(days)`<-factor(data$`Treatment(days)`)
data$group<-factor(data$group,levels = c("Vehicle","Fulvestrant","Fulvestrant+Palbociclib"))
p<-ggplot(data,aes(x=`Treatment(days)`,y=`Tumor.volume(log2)`,color=group,shape=group,group=group))+
  geom_errorbar(aes(ymin=`Tumor.volume(log2)`-se,ymax=`Tumor.volume(log2)`+se),width=0.3,linewidth=0.2)+
  geom_point(size=3)+
  geom_line()+
  ylab("Tumor volume (log2)")+
  xlab("Treatment (days)")+
  ylim(16,600)+
  theme_classic()+
  theme(plot.margin =margin(1,1,1,1,unit = "cm") )+
  scale_color_manual(values = c("grey20","navy","darkred"))
ggsave(plot = p,filename = "output/Fig.6d.pdf",width = 6,height = 4)

#### Fig.6f ####
df<-df.all%>%filter(Strategy=="compensation effects")
df.freq<-table(df$DRUG_NAME,df$DRUG_FIRST)%>%data.frame()
colnames(df.freq)<-c("Drug-2","Drug-1","Freq")
# sort index
df1<-df.freq[df.freq$`Drug-1`=="GDC0810",]
df2<-df.freq[df.freq$`Drug-1`=="Fulvestrant",]
df1$sum<-df1$Freq+df2$Freq
df.freq$`Drug-2`<-factor(df.freq$`Drug-2`,levels = df1%>%arrange(desc(sum))%>%pull(`Drug-2`))
# plot
p<-ggplot(df.freq,aes(x=`Drug-2`,y=Freq,fill=`Drug-1`))+
  geom_bar(stat="identity",position =position_dodge(0.75),width = 0.6)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45,vjust = 0.5,hjust = 0),
        plot.margin =margin(1,1,1,1,unit = "cm"))+
  ylim(0,14)+
  ggtitle("Compensation effects")+
  ylab("Freq. in ER+ samples")+
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 8,"Paired")[3:4])

ggsave(plot = p,filename = "output/Fig.6f",width = 6.5,height = 3.5)

#### Fig.6g ####
# ER+ sample Compensation effects 
# the data not save in github
datasets.ER<-readRDS("/home/tianpeng/home/project2022/ER_scPharm_result.rds")

library(msigdbr)
library(singleseqgset)
library(dplyr)

kegg.mm <- msigdbr(species="Homo sapiens",category="C2",subcategory = "KEGG")
kegg.names <- unique(kegg.mm$gs_name)
kegg.sets <- vector("list",length=length(kegg.names))
names(kegg.sets) <- kegg.names
for (i in names(kegg.sets)) {
  kegg.sets[[i]] <- pull(kegg.mm[kegg.mm$gs_name==i,"gene_symbol"])
}

kegg_cell_cycle<-kegg.sets$KEGG_CELL_CYCLE
# read drug combination data
data<-readRDS("data/scPharm_targeted_drug_combo1214.rds")

df<-c()
for (i in names(data)) {
  for (j in names(data$AH0319)) {
    d<-data[[i]][[j]]
    d$sample<-i
    df<-rbind(df,d)
  }
}

df.all<-df
d<-df.all%>%filter(DRUG_NAME=='AZD5438')%>%filter(Strategy=="compensation effects")
d$DRUG_FIRST%>%table()
# GDC0810 14   each
# Fulvestrant  12

# add cellcycle score
for (name in names(datasets.ER)) {
  seu<-datasets.ER[[name]]
  ## cellcycle
  seu<-CellCycleScoring(seu,s.features = cc.genes.updated.2019$s.genes,cc.genes.updated.2019$g2m.genes)
  seu$S.Score
  seu$G2M.Score
  seu<-AddModuleScore(seu,features = list(cellcycle=kegg_cell_cycle),name = "KEGG_Cell_Cycle")
  seu$KEGG_Cell_Cycle<-seu$KEGG_Cell_Cycle1
  seu$KEGG_Cell_Cycle1<-NULL
  datasets.ER[[name]]<-seu
}

seu<-datasets.ER$PM0360

seu$drug1<-seu$scPharm_label_1925_GDC0810
seu$drug2<-seu$scPharm_label_1401_AZD5438

p1<-DotPlot(subset(seu,drug1!="other"),features = c("ESR1","CDK2","S.Score","G2M.Score"),scale = T,cols = c("grey","red"),group.by = "drug1")+
  ggtitle("Drug 1: GDC0810")+
  theme(axis.text.x = element_text(angle = -45,vjust = 0.5,hjust = 0),plot.margin =margin(1,1,0.2,1,unit = "cm"),title = element_text(size = 9))+
  xlab("")+ylab("")
p2<-DotPlot(subset(seu,drug2!="other"),features = c("CDK2","S.Score","G2M.Score"),scale = T,cols = c("grey","red"),group.by = "drug2")+
  ggtitle("Drug 2: AZD5438")+
  theme(axis.text.x = element_text(angle = -45,vjust = 0.5,hjust = 0),plot.margin =margin(1,1,0.2,1,unit = "cm"),title = element_text(size = 9))+
  xlab("")+ylab("")

ggsave(plot=cowplot::plot_grid(p1,p2,ncol = 2),filename="output/Fig.6g",width=8,height=4)

#### Fig.6h ####
# PM0360
SEU<-datasets.ER$PM0360
SEU$drug1<-SEU$scPharm_label_1816_Fulvestrant
SEU$drug2<-SEU$scPharm_label_1014_Refametinib

p1<-DotPlot(subset(SEU,drug1!="other"),features = c("ESR1",'MAP2K1','MAP2K2'),scale = T,cols = c("grey","red"),group.by = "drug1")+
  ggtitle("Drug 1: Fulvestrant")+
  # NoLegend()+
  theme(axis.text.x = element_text(angle = -45,vjust = 0.5,hjust = 0),plot.margin =margin(1,1,0.2,1,unit = "cm"),title = element_text(size = 9))+
  xlab("")+ylab("")
p2<-DotPlot(subset(SEU,drug2!="other"),features = c('MAP2K1','MAP2K2'),scale = T,cols = c("grey","red"),group.by = "drug2")+
  ggtitle("Drug 2: Refametinib")+
  # NoLegend()+
  theme(axis.text.x = element_text(angle = -45,vjust = 0.5,hjust = 0),plot.margin =margin(1,1,0.2,1,unit = "cm"),title = element_text(size = 9))+
  xlab("")+ylab("")

ggsave(plot=cowplot::plot_grid(p1,p2,ncol = 2),filename="cc.ER_comb.Fulvestrant.mek.dotplot.PM0360.0801.pdf",width=8,height=4)
