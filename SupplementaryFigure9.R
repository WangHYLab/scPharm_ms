library(tidyverse)
library(ggplot2)
library(openxlsx)

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
ggsave(plot=integratePlot(p.list,3),filename = "output/Fig.s9.pdf",width = 13,height = 20)