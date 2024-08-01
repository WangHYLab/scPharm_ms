source("./functions.R")
gdsc<-readGDSC(expr_tpm.file,expr_count.file,GDSC2.file,sample.file,gene.file)

#### Fig.5A ####
datasets.HER2<-readRDS("/home/tianpeng/home/project2022/data2/verify_data/GSE161529_HER2/scPharm_fgsea.rds")


out_HER2<-benchmark_cal(sc_dataset = datasets.HER2,gdsc = gdsc,
                        func = c("Scissor","scDEAL","CaDRReS-Sc","SeuratCCA","scPharm"),
                        drug.choose = NULL,process = FALSE,sampling=c("all"),
                        scPharm.tissue.or.cellline = "tissue",scPharm.type = "BRCA",
                        CaDRReS_Sc.type = "BRCA",times=1)


# ranking plot 
# extract rank table
rank.all<-out_HER2$rank
rank.out<-list()
rank.table<-c()
colname.table<-c()
name<-names(rank.all)
for (n in name) {
  print(strsplit(n,split = "_")[[1]][2])
  me<-strsplit(n,split = "_")[[1]][2]
  da<-strsplit(n,split = "_")[[1]][1]
  if(me=="CaDRReS-Sc"){
    df<-rank.all[[n]]
    drug.info<-read.csv("./CaDRReS-Sc/preprocessed_data/GDSC/drug_stat.csv",row.names = 1,check.names = F)
    df$DrugName<-drug.info[as.character(df$drug_id),]$`Drug Name`
    df$Key<-paste0(df$DrugName,"_Id",df$drug_id)
    df<-df%>%arrange(desc(cell_death))
    df$Rank<-c(1:nrow(df))
    df<-df%>%select("drug_id",  "DrugName" ,   "cell_death",  "Key" ,       "Rank")
    colnames(df)<-c("DRUG_ID", "DRUG_NAME" ,"Score" ,"Key" ,"Rank" )
  }else if(me=="scPharm"){
    df<-rank.all[[n]]
    df$Key<-paste0(df$DRUG_NAME,"_Id",df$DRUG_ID)
  }else{
    df<-rank.all[[n]]
    df<-df[grep(me,rownames(df)),]
    df<-df%>%arrange(desc(Score))
    df$Key<-paste0(df$DRUG_NAME,"_Id",df$DRUG_ID)
    df$Rank<-c(1:nrow(df))
  }
  na<-paste0(da,"_",me)
  rank.out[[na]]<-df
  colname.table<-c(colname.table,na)
  rank.table<-cbind(rank.table,c(df$Key,rep(NA,295-nrow(df))))
}

colnames(rank.table)<-colname.table
write.csv(rank.table,file = "./output/HER2.ranktable.csv")


table<-read.csv("./output/HER2.ranktable.csv",check.names = F)
name<-lapply(colnames(table), function(x){strsplit(x,split = "_")[[1]][2]})%>%as.character()
names(name)<-colnames(table)
name<-sort(name)
order_vector <- c("SeuratCCA", "Scissor", "scDEAL", "CaDRReS-Sc", "scPharm")
sorted_name <- name[order(match(name, order_vector))]

select_name<-setdiff(names(sorted_name),names(sorted_name)[grep('MH0161|MH0176',names(sorted_name))])
data<-table[,select_name]
data<-table[,names(sorted_name)]


data[data=="Lapatinib_Id1558"]<-1
data[data=="Lapatinib_Id119"]<-1 #for some tool use diff GDSC version
data[data=="Afatinib_Id1032"]<-2
data[data=="Sapitinib_Id1549"]<-3
data[is.na(data)]<-0
data[data!=1&data!=2&data!=3]<-0
library(purrr)
df<-map_df(data, as.numeric)%>%as.matrix()
rownames(df)<-rownames(data)

# annotation 
row.anno<-data.frame(sample=lapply(colnames(df),function(x){strsplit(x,"_")[[1]][1]}) %>% as.character(), row.names = colnames(df))
anno.color<- list(sample = c("AH0308"="#EA75A7","MH0031"="#F3C175","MH0069"="#75C799","MH0161"="#D9E675","MH0176"="#BEB2D5","PM0337"="#FFF775"))

pdf("./output/Fig.5a.pdf",width = 10,height = 4)
pheatmap::pheatmap(df%>%t(),cluster_rows = F,cluster_cols = F,
                   show_colnames = F,
                   annotation_row = row.anno,
                   annotation_colors = anno.color,
                   # width = 2, 
                   border_color = NA,gaps_row = c(1,2,3,4,5)*6,
                   legend_labels = c("Other","Lapatinib","Afatinib","Sapitinib"),
                   legend_breaks = c(0.4,1.15,1.9,2.65),
                   color = c("grey80","#4DAF4A","#E41A1C","#377EB8"))
dev.off()


#### Fig.5C ####

rank.table<-read.csv("./output/HER2.ranktable.csv",row.names = 1,check.names = F)
df<-c()
druglist<-c("Afatinib","Lapatinib","Sapitinib")
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
    # index<-grep(me,names(rank.table))
    
  }
}
df<-data.frame(df)
colnames(df)<-c('Drug','Method','Sample','Order')
df$Order<-df$Order%>%as.numeric()
df$Method<-factor(df$Method,levels = c("SeuratCCA","Scissor","scDEAL","CaDRReS-Sc","scPharm"))
df$Score<-1/df$Order

p1<-ggplot(df,aes(y=Order,x=Method,fill=Method))+
  stat_boxplot(geom="errorbar",width=0.3)+
  geom_boxplot(outlier.size = 0,outlier.alpha = 0)+
  geom_jitter(width = 0.25,size=0.5,alpha=0.5)+
  geom_signif(comparisons = list(c('scPharm','CaDRReS-Sc'),
                                 c('scPharm',"scDEAL"),
                                 c('scPharm',"Scissor"),
                                 c('scPharm',"SeuratCCA")
  ),
  y_position = c(220,230,240,250),map_signif_level = T,tip_length = 0.01,vjust = 0.8)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        plot.margin = unit(c(1,0.5,0.5,0.5),"cm"),title = element_text(size = 10))+
  scale_fill_manual(values = c("#936E9F","#58A671","#96991E","#CE7269","#6EAFDE"))

ggsave(filename = "./output/Fig.5c.pdf",plot =p1,width = 4.2,height =5)




#### Fig.5E ####
# running time plot
df<-out_HER2$df
cellnum<-lapply(datasets.HER2, function(x){dim(x)[2]})%>%data.frame()%>%sort()
xlab<-paste0(colnames(cellnum),"\n(",cellnum[1,],")")
df$Data<-factor(df$Data,levels = colnames(cellnum))
p<-ggplot(df,aes(x=Data,y=log10(time),group=Method,color=Method,shape=Method))+
  geom_point(size=3)+
  geom_line(cex=0.5)+
  theme_bw()+
  ylab("log10 Time(s)")+
  xlab("Sample(cells)")+
  scale_x_discrete(labels =str_wrap(xlab, width =6))
ggsave(p,filename = "Fig.5e.pdf",width = 6,height = 5)


#### Fig.5B ####
datasets.ER<-readRDS("/home/tianpeng/home/project2022/data2/verify_data/GSE161529_ER/scPharm_fgsea.rds")
out_ER<-benchmark_cal(sc_dataset = datasets.ER,gdsc = gdsc,
                      func = c("Scissor","scDEAL","CaDRReS-Sc","SeuratCCA","scPharm"),
                      drug.choose = NULL,process = FALSE,sampling=c("all"),
                      scPharm.type = "BRCA",scPharm.tissue.or.cellline = "tissue",
                      CaDRReS_Sc.type = "",CaDRReS_Sc.GDSC2 = TRUE,times=1)

# rank
rank.all<-out_ER$rank

rank.out<-list()
rank.table<-c()
colname.table<-c()
name<-names(rank.all)
for (n in name) {
  
  print(strsplit(n,split = "_")[[1]][2])
  me<-strsplit(n,split = "_")[[1]][2]
  da<-strsplit(n,split = "_")[[1]][1]
  if(me=="CaDRReS-Sc"){
    df<-rank.all[[n]]
    drug.info<-read.csv("./CaDRReS-Sc/preprocessed_data/GDSC/drug_stat2.csv",row.names = 1,check.names = F) # 2为GDSC2
    df$DrugName<-drug.info[as.character(df$drug_id),]$`Drug Name`
    df$Key<-paste0(df$DrugName,"_Id",df$drug_id)
    df<-df%>%arrange(desc(cell_death))
    df$Rank<-c(1:nrow(df))
    df<-df%>%select("drug_id",  "DrugName" ,   "cell_death",  "Key" ,       "Rank")
    colnames(df)<-c("DRUG_ID", "DRUG_NAME" ,"Score" ,"Key" ,"Rank" )
  }else if(me=="scPharm"){
    df<-rank.all[[n]]
    df$Key<-paste0(df$DRUG_NAME,"_Id",df$DRUG_ID)
  }else{
    df<-rank.all[[n]]
    df<-df[grep(me,rownames(df)),]
    df<-df%>%arrange(desc(Score))
    df$Key<-paste0(df$DRUG_NAME,"_Id",df$DRUG_ID)
    df$Rank<-c(1:nrow(df))
  }
  na<-paste0(da,"_",me)
  rank.out[[na]]<-df
  colname.table<-c(colname.table,na)
  rank.table<-cbind(rank.table,c(df$Key,rep(NA,295-nrow(df))))
}

colnames(rank.table)<-colname.table
write.csv(rank.table,file = "ER.ranktable.csv")

# plot 
data<-read.csv("./output/ER.ranktable.csv",row.names = 1,check.names = F)

data[data=="Fulvestrant_Id1816"]<-1  #Fulvestrant_high
data[data=="Fulvestrant_Id1200"]<-2  #Fulvestrant_low 
data[data=="GDC0810_Id1925"]<-3

data[data!=1&data!=2&data!=3]<-0
library(purrr)
df<-map_df(data, as.numeric)%>%as.matrix()
rownames(df)<-rownames(data)
#排序
index<-lapply(colnames(df), function(x){strsplit(x,split = "_")[[1]][2]})%>%as.character()
names(index)<-colnames(df)
index<-index%>%sort()%>%data.frame()
colnames(index)<-'method'
index$name<-rownames(index)
index$method<-factor(index$method,levels = c("SeuratCCA", "Scissor", "scDEAL", "CaDRReS-Sc", "scPharm"))
index <- index[order(index$method), ]
df<-df[,index$name]

row.anno<-data.frame(sample=lapply(colnames(df),function(x){strsplit(x,"_")[[1]][1]}) %>% as.character(), row.names = colnames(df))
anno.color<- list(sample = c('AH0319'=rainbow(14)[1],
                             'MH0001'=rainbow(14)[2],
                             'MH0025'=rainbow(14)[3],
                             'MH0029-7C'=rainbow(14)[4],
                             'MH0029-9C'=rainbow(14)[5],
                             'MH0032'=rainbow(14)[6],
                             'MH0040'=rainbow(14)[7],
                             'MH0042'=rainbow(14)[8],
                             'MH0043-T'=rainbow(14)[9],
                             'MH0064-T'=rainbow(14)[10],
                             'MH0151'=rainbow(14)[11],
                             'MH0163'=rainbow(14)[12],
                             'MH0173-T'=rainbow(14)[13],
                             'PM0360'=rainbow(14)[14]))

pdf("./output/Fig.5b.pdf",width = 12,height = 8)
pheatmap::pheatmap(df%>%t(),cluster_rows = F,cluster_cols = F,
                   show_colnames = F,
                   annotation_row = row.anno,
                   annotation_colors = anno.color,
                   legend_labels = c("Other","Fulvestrant_high","Fulvestrant_low","GDC0810"),
                   border_color = 'grey',gaps_row = c(1,2,3,4,5)*14,
                   color = c("grey80","#4DAF4A","#377EB8","#E41A1C"))
dev.off()

#### Fig.5D ####
df<-c()
druglist<-c("GDC0810_Id1925","Fulvestrant_Id1816","Fulvestrant_Id1200")

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

p1<-ggplot(df,aes(y=Order,x=Method,fill=Method))+
  stat_boxplot(geom="errorbar",width=0.3)+
  geom_boxplot(outlier.size = 0,outlier.alpha = 0)+
  geom_jitter(width = 0.25,size=0.5,alpha=0.5)+
  geom_signif(comparisons = list(c('scPharm','CaDRReS-Sc'),
                                 c('scPharm','scDEAL'),
                                 c('scPharm','Scissor'),
                                 c('scPharm','SeuratCCA')),y_position = c(250,260,270,280),map_signif_level = T,tip_length = 0.01,vjust = 0.8)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        plot.margin = unit(c(1,0.5,0.5,0.5),"cm"),title = element_text(size = 10))+
  ggtitle("ER+ sample")+
  scale_fill_manual(values = c("#60B9EC","#E1746B","#9C9F17","#49AF76","#9E71A6")[5:1])

ggsave(filename = "./output/Fig.5d.pdf",plot =p1,width = 4.2,height = 5)


#### Fig.5F ####
# time
df<-out_ER$df

cellnum<-lapply(datasets.ER, function(x){dim(x)[2]})%>%data.frame(check.names = F)%>%sort()
names(cellnum)
xlab<-paste0(colnames(cellnum),"\n(",cellnum[1,],")")
df$Data<-factor(df$Data,levels = colnames(cellnum))
p<-ggplot(df,aes(x=Data,y=log10(time),group=Method,color=Method,shape=Method))+
  geom_point(size=3)+
  geom_line(cex=0.5)+
  theme_classic()+
  ylab("log10 Time(s)")+
  xlab("Sample(cells)")+
  # scale_x_discrete(labels =str_wrap(xlab, width =6))+
  scale_x_discrete(labels =xlab)+
  theme(axis.text.x = element_text(angle = 270,vjust = 0.6))
ggsave(p,filename = "./output/Fig.5f",width = 8,height = 4)


#### Fig.5G ####

out<-benchmark_cal(datasets.HER2[[1]],drugdata,drug.choose = "Erlotinib",sampling=c(500,1000,1500,2000),scPharm.type = "LUAD",times=3)
ggsave(out$plot,filename = "./output/Fig.5g.pdf",width = 6,height = 4)