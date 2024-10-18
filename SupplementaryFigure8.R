library(Seurat)
library(tidyverse)
library(tictoc)
library(readxl)
library(openxlsx)
library(ggplot2)
library(ggsignif)
library(cowplot)

plot_rank_heatmap<-function(rank.table=data.frame(),drug.list=c(),geo.sample="",title_main="",title_sub=""){
  df<-rank.table
  strings_to_find<-drug.list
  # id mapping
  all_drug<-df[,1]
  strings_to_find_withID<-c()
  for (d in strings_to_find) {
    strings_to_find_withID<-c(strings_to_find_withID,all_drug[grep(d,all_drug,fixed = T)])
  }
  strings_to_find<-strings_to_find_withID %>% unique()
  # Get all unique values in the DataFrame  
  unique_values <- unique(as.vector(df))
  
  # Create a color map that gives the color in strings_to_find, with the rest of the values in gray
  color_map <- setNames(rep("grey", length(unique_values)), unique_values)
  color_map[strings_to_find] <- rainbow(length(strings_to_find))
  
  # Convert DataFrame to long format and add color columns  
  df_long <- df %>%  as.data.frame() %>% 
    rownames_to_column("Row") %>%  
    mutate(Row = as.numeric(Row)) %>%  
    pivot_longer(cols = -Row, names_to = "Column", values_to = "Value") %>%  
    mutate(  
      Label = ifelse(Value %in% strings_to_find, Value, "Other"),
      Color = ifelse(Value %in% strings_to_find, color_map[Value], "grey")
    )
  
  # Ensure that the color mapping matches the unique value in the Label column
  final_color_map <- setNames(c(rainbow(length(strings_to_find)), "grey"), 
                              c(strings_to_find, "Other"))
  
  # plot
  p<-ggplot(df_long, aes(x = Row, y = Column, fill = Label)) +  
    geom_tile() +  
    scale_fill_manual(values = final_color_map, 
                      breaks = c(strings_to_find, "Other"), 
                      labels = c(strings_to_find, "Other")) +
    theme_minimal() +  
    labs(x = "Rank", y = geo.sample, fill = "Drug") +
    theme(legend.position = "right",
          panel.grid = element_blank(),
          axis.text.x = element_blank())+
    ggtitle(label = title_main,subtitle = title_sub)
  return(p)
}

# drug：
BRAF<-c("Dabrafenib","SB590885","PLX-4720")
KRAS<-c("KRAS (G12C) Inhibitor-12")
MEK<-c("Trametinib","Selumetinib","Refametinib","PD0325901")
ERK<-c("Ulixertinib","SCH772984","ERK_2440","ERK_6604","VX-11e")
EGFR<-c("AZD3759","Osimertinib","Gefitinib","Erlotinib")
EGFR_multi<-c("Afatinib","Lapatinib","Sapitinib","Foretinib")
PI3K<-c("GNE-317","Alpelisib","AZD8186","Buparlisib","AZD6482","CZC24832","Taselisib","AMG-319","Pictilisib","Dactolisib")
PI3K_akt<-c("AT13148","Ipatasertib","MK-2206","Uprosertib","Afuresertib","GSK2110183B","Uprosertib")
RAS_mut<-c("THR-101","THR-102","THR-103") #WIMM synthesis	PI3K/MTOR signaling	Mutant RAS
JAK<-c("JAK_8517","GSK2276186C","AZ960","Lestaurtinib","JAK1_8709","Ruxolitinib")
VEGFR<-c("Foretinib","Axitinib","Cediranib","Motesanib","Sorafenib")

# CMML GSE218390 ----
df<-openxlsx::read.xlsx("GSE218390/ranktable.xlsx")
p1<-plot_rank_heatmap(openxlsx::read.xlsx("GSE218390/ranktable.xlsx") %>% dplyr::select(BM.3,BM.4),
                      c(ERK,MEK,"azacytidine"),
                      "GSE218390",
                      title_main = "CMML RAS pathway mutation",
                      title_sub = "Drug: ERK,MEK")+theme(legend.position = "bottom")

p2<-plot_rank_heatmap(openxlsx::read.xlsx("GSE218390/ranktable.xlsx") %>% dplyr::select(-BM.2,-BM.3,-BM.4,-BM.8),
                      c(PI3K_akt,PI3K),
                      "GSE218390",
                      title_main = "CMML RAS pathway mutation",
                      title_sub = "Drug: PI3K")+theme(legend.position = "bottom")

p3<-plot_rank_heatmap(openxlsx::read.xlsx("GSE218390/ranktable.xlsx") %>% dplyr::select(-BM.2,-BM.3,-BM.4,-BM.8),
                      c('MK-2206','Afuresertib',"azacytidine"),
                      "GSE218390",
                      title_main = "CMML RAS pathway mutation",
                      title_sub = "Drug: PI3K")+theme(legend.position = "bottom")
ggsave(plot = cowplot::plot_grid(p1,p2,p3,ncol = 1),filename = "Fig.s8a.pdf",width = 12,height = 10)


p1<-plot_rank_heatmap(openxlsx::read.xlsx("GSE218390/ranktable.xlsx") %>% dplyr::select(-BM.2,-BM.8),
                      c("azacytidine",'MK-2206','Afuresertib',ERK,MEK),
                      "GSE218390",
                      title_main = "CMML RAS pathway mutation",
                      title_sub = "Drug: DNMT耐药测试")+theme(legend.position = "bottom")+scale_color_manual(values = "green")

# SKCM  GSE200218 ----
# SKCM， BRAF mut
p1<-plot_rank_heatmap(openxlsx::read.xlsx("GSE200218/ranktable.xlsx") %>% dplyr::select('GSM6022252','GSM6022253','GSM6022255'),
                      c("Dabrafenib","THR-103","Trametinib"), #BRAF
                      "GSE200218",
                      title_main = "SKCM BRAF/HRAS mutation",
                      title_sub = "Drug: BRAF, RAS")+theme(legend.position = "bottom")
ggsave(plot = p1,filename = "Fig.s8b.pdf",width = 10,height = 4)

# LUAD multi-datasets  ----
# GSE247684
p1<-plot_rank_heatmap(openxlsx::read.xlsx("GSE247684/ranktable.xlsx") %>% dplyr::select("PC9"),
                      c("Ipatasertib","Afuresertib","PD0325901","Ulixertinib","SB590885"), #BRAF
                      "GSE247684",
                      title_main = "LUAD EGFR mutation",
                      title_sub = "Drug: 2*PI3K, 3*MAPK")+theme(legend.position = "bottom")
# GSE148466
p2<-plot_rank_heatmap(openxlsx::read.xlsx("GSE148466/ranktable.xlsx") %>% dplyr::select("GSM4472055"),
                      c("Ipatasertib","Afuresertib","KRAS (G12C) Inhibitor-12"), #BRAF
                      "GSE148466",
                      title_main = "NSCLC KRAS/EGFR mutation",
                      title_sub = "Drug: 2*PI3K, 1*KRAS G12C inhibitor")+theme(legend.position = "bottom")

# GSE241934
p3<-plot_rank_heatmap(openxlsx::read.xlsx("GSE241934/ranktable.xlsx") %>% dplyr::select(-P591,-P547,-P438) ,
                      c("Lapatinib","Selumetinib","Ulixertinib"), #BRAF
                      "GSE241934",
                      title_main = "NSCLC EGFR mutation",
                      title_sub = "Drug: 1*EGFR, 2*MAPK")+theme(legend.position = "bottom")

# GSE223779
p4<-plot_rank_heatmap(openxlsx::read.xlsx("GSE223779/ranktable.xlsx"),
                      c("Sorafenib","Crizo"), #BRAF
                      "GSE223779",
                      title_main = "LUAD EML4-ALK tumor organoid",
                      title_sub = "Drug: 1*TKI")+theme(legend.position = "bottom")
# GSE136246
p5<-plot_rank_heatmap(openxlsx::read.xlsx("GSE136246/ranktable.xlsx"),
                      c("KRAS (G12C) Inhibitor-12"), #BRAF
                      "GSE136246",
                      title_main = "LUAD KRAS mutation",
                      title_sub = "Drug: KRAS (G12C) Inhibitor-12")+theme(legend.position = "bottom")


ggsave(plot = cowplot::plot_grid(p1,p2,p3,p4,p5,ncol = 1),filename = "Fig.s8c.pdf",width = 10,height = 18)



