library(Scissor)
library(tictoc)
library(readxl)
library(openxlsx)
library(tidyverse)
library(ggplot2)
library(ggsignif)
library(cowplot)

dir=getwd()

#drug source file
{
  expr_count.file <- "GDSC/rnaseq_read_count_20220624.csv"
  expr_tpm.file   <- "GDSC/rnaseq_tpm_20220624.csv"
  GDSC1.file      <- "GDSC/GDSC1_fitted_dose_response_24Jul22.xlsx"
  GDSC2.file      <- "GDSC/GDSC2_fitted_dose_response_24Jul22.xlsx"
  sample.file     <- "GDSC/model_list_20221102.csv"
  gene.file       <- "GDSC/gene_identifiers_20191101.csv"
}

####  functions  ####
getExprAndDrug<-function(expr.file=expr_tpm.file,
                         expr.count.file=expr_count.file,
                         gdsc.version=GDSC2.file,
                         sample.file=sample.file,
                         gene.file=gene.file,
                         tissue.or.cancertype="tissue",
                         drug.choose=drug.choose,
                         sample.choose=sample.choose){
  #expression mtx
  expr<-read.csv(expr.file,row.names = 1,skip = 1,check.names = F)
  expr<-expr[4:nrow(expr),]
  expr<-expr[!duplicated(expr[,1]), ] 
  rownames(expr)<-expr[,1]
  expr<-expr[-1]
  
  #expression count mtx
  expr.count<-read.csv(expr.count.file,row.names = 1,skip = 1,check.names = F)
  expr.count<-expr.count[4:nrow(expr.count),]
  expr.count<-expr.count[!duplicated(expr.count[,1]), ] 
  rownames(expr.count)<-expr.count[,1]
  expr.count<-expr.count[-1]
  
  #drug info
  pham<-read_xlsx(gdsc.version)%>%as.data.frame()
  
  #gene info
  gene<-read.csv(gene.file,row.names = 1)
  
  #sample info
  sample<-read.csv(sample.file,row.names = 1)
  sample$group<-"other"
  
  pham<-pham%>%filter(DRUG_NAME==drug.choose)
  max.conc<-pham[1,'MAX_CONC']%>%as.numeric()%>%log()
  
  if(tissue.or.cancertype=="tissue"){
    pham$tissue<-sample[pham$SANGER_MODEL_ID,]$tissue
  }
  if(tissue.or.cancertype=="cancertype"){
    pham$tissue<-sample[pham$SANGER_MODEL_ID,]$cancer_type_detail
  }
  pham<-pham%>%filter(tissue==sample.choose)
  
  #敏感耐药划分
  pham$drug_group<-""
  pham[pham$LN_IC50<max.conc,]$drug_group<-"sensitive"
  pham[pham$LN_IC50>=max.conc,]$drug_group<-"resistant"
  pham$drug_group_scissor<-0
  pham[pham$drug_group=="sensitive",]$drug_group_scissor<-1
  rownames(pham)<-pham$CELL_LINE_NAME
  
  #overlap
  cell.overlap<-intersect(pham$CELL_LINE_NAME,colnames(expr))
  pham<-pham[cell.overlap,]
  expr<-expr[,cell.overlap]
  expr.count<-expr.count[,cell.overlap]
  #tpm
  out<-as.data.frame(lapply(expr,as.numeric))%>%as.matrix()
  colnames(out)<-colnames(expr)
  rownames(out)<-rownames(expr)
  out<-na.omit(out)
  #count
  out.count<-as.data.frame(lapply(expr.count,as.numeric))%>%as.matrix()
  colnames(out.count)<-colnames(expr.count)
  rownames(out.count)<-rownames(expr.count)
  out.count<-na.omit(out.count)
  
  return(list(expr=out,expr.count=out.count,pham=pham))
}
# read GDSC source files
readGDSC<-function(expr.file,count.file,gdsc.file,sample.file,gene.file){
  
  # expr tpm
  gdsc.expr<-read.csv(expr.file,row.names = 1,skip = 1,check.names = F)
  gdsc.expr<-gdsc.expr[4:nrow(gdsc.expr),]
  gdsc.expr<-gdsc.expr[!duplicated(gdsc.expr[,1]), ] 
  rownames(gdsc.expr)<-gdsc.expr[,1]
  gdsc.expr<-gdsc.expr[-1]
  
  #expression count mtx
  gdsc.count<-read.csv(count.file,row.names = 1,skip = 1,check.names = F)
  gdsc.count<-gdsc.count[4:nrow(gdsc.count),]
  gdsc.count<-gdsc.count[!duplicated(gdsc.count[,1]), ] 
  rownames(gdsc.count)<-gdsc.count[,1]
  gdsc.count<-gdsc.count[-1]
  
  #drug info
  gdsc.drug<-read_xlsx(gdsc.file)%>%as.data.frame()
  print(gdsc.drug$TCGA_DESC%>%unique())
  #gene info
  gdsc.gene<-read.csv(gene.file,row.names = 1)
  
  #sample info
  gdsc.sample<-read.csv(sample.file,row.names = 1)

  return(list(expr=gdsc.expr,count=gdsc.count,gene=gdsc.gene,drug=gdsc.drug,sample=gdsc.sample))
}

# Filter GDSC data based on drug and selected cancer cell line type
filterGDSC<-function(gdsc,drug.choose,sample.choose=sample.choose,by.id=FALSE){
  expr<-gdsc$expr
  pham<-gdsc$drug
  sample<-gdsc$sample
  expr.count<-gdsc$count
  
  if(by.id){
    pham<-pham%>%filter(DRUG_ID==drug.choose)
  }else{
    pham<-pham%>%filter(DRUG_NAME==drug.choose)
  }
  
  max.conc<-pham[1,'MAX_CONC']%>%as.numeric()%>%log()
  
  
  pham<-pham%>%filter(TCGA_DESC==sample.choose)
  
  # split sensitive and resistant
  pham$drug_group<-""
  tryCatch({
    pham[pham$LN_IC50<max.conc,]$drug_group<-"sensitive"
  },error=function(x){
    message("no sensitive CCL")
  })
  tryCatch({
    pham[pham$LN_IC50>=max.conc,]$drug_group<-"resistant"
  },error=function(x){
    message("no rensistant CCL")
  })
  
  pham$drug_group_scissor<-0
  tryCatch({
    pham[pham$drug_group=="sensitive",]$drug_group_scissor<-1
  },error=function(x){
    message("no sensitive label for Scissor")
  })
  
  rownames(pham)<-pham$CELL_LINE_NAME
  
  #overlap
  cell.overlap<-intersect(pham$CELL_LINE_NAME,colnames(expr))
  pham<-pham[cell.overlap,]
  expr<-expr[,cell.overlap]
  expr.count<-expr.count[,cell.overlap]
  #tpm
  out<-as.data.frame(lapply(expr,as.numeric))%>%as.matrix()
  colnames(out)<-colnames(expr)
  rownames(out)<-rownames(expr)
  out<-na.omit(out)
  #count
  out.count<-as.data.frame(lapply(expr.count,as.numeric))%>%as.matrix()
  colnames(out.count)<-colnames(expr.count)
  rownames(out.count)<-rownames(expr.count)
  out.count<-na.omit(out.count)
  message(paste0("Dim of expr : ",paste0(dim(out)%>%as.numeric(),collapse=",")))
  message(paste0("Dim of count: ",paste0(dim(out.count)%>%as.numeric(),collapse=",")))
  message(paste0("Dim of pham : ",paste0(dim(pham)%>%as.numeric(),collapse=",")))
  
  return(list(expr=out,expr.count=out.count,pham=pham))
}

# Single-cell data preprocessing
sc_process<-function(data,...){
  data <- NormalizeData(data,verbose = FALSE)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
  data <- ScaleData(data,features=VariableFeatures(data),verbose = FALSE) 
  # diff pc by cell number
  # Because it doesn't work when the cell number is few but the pc count or other params is high
  cell.num<-dim(data)[2]
  pc<-if( cell.num>50){
    50
  }else if( cell.num>30){
    30
  }else if( cell.num>10){
    10
  }else{
    cell.num-1
  }
  n.n<-if( cell.num>30){
    30
  }else{
    cell.num-1
  }
  data <- RunPCA(data, npcs =pc,verbose = FALSE)
  data <- RunUMAP(data, features = VariableFeatures(data),n.neighbors =n.n ,verbose = FALSE)
  data <- FindNeighbors(data, dims = 1:pc, verbose = FALSE)
  data <- FindClusters(data, verbose = FALSE,resolution=0.5)
  return(data)
}

# Log file saving functions
logging<-function(string="",status="",append=TRUE,file=paste0(dir,"/benchmark.log"),...){
  out<-paste0("[",Sys.time(),"]\t")
  write.table(paste0(out,string,"\t",status),
              file = file,
              row.names = F,col.names = F,quote = F,append=append,...)
}

# Use functions to call Methods for evaluation.
run_scDEAL<-function(sc_dataset,gdsc,drug.choose=NULL,dimreduce="DAE",...){
  library(reticulate)
  library(anndata)
  # use python env 
  use_condaenv("/home/zhengjie/software/anaconda3/envs/stlearn",required = TRUE)
  py_config()
  wd<-getwd()
  setwd("./scDEAL")
  mtx<-GetAssayData(sc_dataset,assay = "RNA",slot = "count")%>%as.matrix()
  write.csv(mtx,file = paste0(dir,"/data/default-data.csv"))
  
  if(is.null(drug.choose)){
    # multi-drugs ----
    message("Run scDEAL in all drugs.")
    time.use<-0
    
    out.metadata<-sc_dataset@meta.data
    drugId.list<-gdsc$drug$DRUG_ID%>%unique()
    pb <- txtProgressBar(style=3)
    len<-length(drugId.list)
    for(d in c(1:len)){
      setTxtProgressBar(pb, d/len)
      drugId<-drugId.list[d]
      drugdata<-filterGDSC(gdsc,drug.choose = drugId,sample.choose,by.id = TRUE)
      drugName<-drugdata$pham$DRUG_NAME[1]
     
      tic()
      tryCatch({
        system(command = paste0("python bulkmodel.py --drug '",drugName,"' --dimreduce '",dimreduce,"' --encoder_h_dims '256,128' --predictor_h_dims '128,64' --bottleneck 512 --data_name 'default-data' --sampling 'upsampling' --dropout 0.1 --lr 0.5 --printgene 'F' -mod 'new' --checkpoint 'False'"))
        system(command = paste0("python scmodel.py --sc_data 'default-data' --drug '",drugName,"' --dimreduce '",dimreduce,"' --bulk_h_dims '256,128'  --predictor_h_dims '128,64' --bottleneck 512  --sampling 'upsampling' --dropout 0.1 --lr 0.5 --printgene 'F' -mod 'new' --checkpoint 'False'"))
        
      },error=function(e){
        message(paste0("Error in drug: ",drugName))
      })
      # record running time
      k<-toc()
      print(as.numeric(k$toc-k$tic))
      
      file<-paste0("./save/adata/default-dataintegrate_data_default-data_drug_",drugName,"_bottle_512_edim_256,128_pdim_128,64_model_",dimreduce,"_dropout_0.1_gene_F_lr_0.5_mod_new_sam_upsampling.h5ad")
      name<-paste0("scDEAL_label_",drugId,"_",drugName)
      out.metadata[[name]]<-"none"
      
      if(file.exists(file)){
        data <- read_h5ad(file)
        # data$obs$sens_label #where 0 represents resistance and 1 represents sensitivity respectively
        outdata<-CreateSeuratObject(t(data$X),meta.data = data$obs) #B ecause of the filtering in scDEAL, there are fewer cells, so the seurat object is reconstructed.
        outdata$scDEAL<-"sensitive"
        outdata@meta.data[outdata$sens_label==0,]$scDEAL<-"resistant"
        
        out.metadata[rownames(outdata@meta.data),name]<-outdata$scDEAL
        
      }
      
      time.use<-time.use+k$toc-k$tic 
    }
    
    close(pb)
    rank<-scPharm_rank_method(out.metadata)
    outdata@meta.data<-out.metadata[outdata@meta.data%>%rownames(),]
    
    setwd(wd)
    return(list(sc_dataset=outdata,time=as.numeric(time.use),rank=rank,comb=""))
    
  }else{
    # single-drug ----
    message(paste0("Run scDEAL in drug: ",drug.choose,"."))
    # record running time
    tic()
    system(command = paste0("python bulkmodel.py --drug '",drug.choose,"' --dimreduce '",dimreduce,"' --encoder_h_dims '256,128' --predictor_h_dims '128,64' --bottleneck 512 --data_name 'default-data' --sampling 'upsampling' --dropout 0.1 --lr 0.5 --printgene 'F' -mod 'new' --checkpoint 'False'"))
    system(command = paste0("python scmodel.py --sc_data 'default-data' --drug '",drug.choose,"' --dimreduce '",dimreduce,"' --bulk_h_dims '256,128'  --predictor_h_dims '128,64' --bottleneck 512  --sampling 'upsampling' --dropout 0.1 --lr 0.5 --printgene 'F' -mod 'new' --checkpoint 'False'"))
    # record running time
    k<-toc()
    print(as.numeric(k$toc-k$tic))
    
    file<-paste0("./save/adata/default-dataintegrate_data_default-data_drug_",drug.choose,"_bottle_512_edim_256,128_pdim_128,64_model_",dimreduce,"_dropout_0.1_gene_F_lr_0.5_mod_new_sam_upsampling.h5ad")
    data <- read_h5ad(file)
    # data$obs$sens_label # where 0 represents resistance and 1 represents sensitivity respectively
    outdata<-CreateSeuratObject(t(data$X),meta.data = data$obs) 
    outdata$scDEAL<-"sensitive"
    outdata@meta.data[outdata$sens_label==0,]$scDEAL<-"resistant"
    
    setwd(wd)

    return(list(sc_dataset=outdata,time=as.numeric(k$toc-k$tic)))
    
  }
  
}
run_Scissor<-function(sc_dataset,gdsc,drug.choose=NULL,sample.choose="BRCA",plot=FALSE,...){
  requireNamespace("Scissor")
  if(is.null(drug.choose)){
    # multi-drugs ----
    message("Run Scissor in all drugs.")
    time.use<-0
    
    out.metadata<-sc_dataset@meta.data
    drugId.list<-gdsc$drug$DRUG_ID%>%unique()
    pb <- txtProgressBar(style=3)
    len<-length(drugId.list)
    for(d in c(1:len)){
      setTxtProgressBar(pb, d/len)
      drugId<-drugId.list[d]
      drugdata<-filterGDSC(gdsc,drug.choose = drugId,sample.choose,by.id = TRUE)
      drugName<-drugdata$pham$DRUG_NAME[1]
      
      bulk_dataset<-drugdata$expr
      #phon 0-1 Data：
      phenotype<-drugdata$pham$drug_group_scissor
      #The different 0-1 encoding in phenotype: `sensitive = 1` and `resistant = 0`
      if(length(unique(phenotype))==1){
        if(unique(phenotype)==0){
          tag <- c('resistant') 
        }else{
          tag <- c('sensitive') 
        }
      }else{
        tag <- c('resistant','sensitive') 
      }
      print(tag)
      
      # record running time
      tic()
      tryCatch({
        out <- Scissor(bulk_dataset, sc_dataset, phenotype, tag = tag, alpha = 0.5, family = "binomial")
      },error=function(e){
        message(paste0("Error in drug: ",drugName))
        out <- list(Scissor_pos=character(0),Scissor_neg=character(0))
      })
     
      # record running time
      k<-toc()
      print(as.numeric(k$toc-k$tic))
      
      sc_dataset$Scissor<-"none"
      tryCatch({
        sc_dataset@meta.data[out$Scissor_pos,]$Scissor<-'sensitive'
      },error=function(e){message("No sensitive cell")})
      
      tryCatch({
        sc_dataset@meta.data[out$Scissor_neg,]$Scissor<-'resistant'
      },error=function(e){message("No resistant cell")})
      
      
      out.metadata[[paste0("Scissor_label_",drugId,"_",drugName)]]<-sc_dataset$Scissor
      time.use<-time.use+k$toc-k$tic
    }
    
    close(pb)
    rank<-scPharm_rank_method(out.metadata)
    sc_dataset@meta.data<-out.metadata[sc_dataset@meta.data%>%rownames(),]
    return(list(sc_dataset=sc_dataset,time=as.numeric(time.use),rank=rank,comb=""))
    
    
  }else{
    # single-cell ----
    message(paste0("Run Scissor in drug: ",drug.choose,"."))
    drugdata<-filterGDSC(gdsc,drug.choose = drug.choose,sample.choose,by.id = FALSE)
    
    
    bulk_dataset<-drugdata$expr
    # 0-1Data：
    phenotype<-drugdata$pham$drug_group_scissor
    #The different 0-1 encoding in phenotype: `sensitive = 1` and `resistant = 0`
    tag <- c('resistant','sensitive') 
    
    # No sensitive cell
    tic()
    out <- Scissor(bulk_dataset, sc_dataset, phenotype, tag = tag, alpha = 0.5, family = "binomial",...)
    # No sensitive cell
    k<-toc()
    print(as.numeric(k$toc-k$tic))
    
    sc_dataset$Scissor<-"none"
    tryCatch({
      sc_dataset@meta.data[out$Scissor_pos,]$Scissor<-'sensitive'
    },error=function(e){message("No sensitive cell")})
    
    tryCatch({
      sc_dataset@meta.data[out$Scissor_neg,]$Scissor<-'resistant'
    },error=function(e){message("No resistant cell")})
    if(plot){
      DimPlot(sc_dataset, reduction = 'umap',group.by = "Scissor", label = T, label.size = 8,repel = T)
    }
    
    return(list(sc_dataset=sc_dataset,time=as.numeric(k$toc-k$tic)))
  }
  
}
run_CaDRReS_Sc<-function(sc_dataset,gdsc,drug.choose=NULL,sample.choose="BRCA",CaDRReS_Sc.GDSC2=FALSE,mode="cell",...){
  # sample.choose="" means no cancer type filtering：TCGA_DESC
  # drug.choose=NULL means return all drugs
  library(reticulate)
  library(anndata)
  library(openxlsx)
  # use python env
  use_condaenv("/home/zhengjie/software/anaconda3/envs/stlearn",required = TRUE)
  py_config()
  # use_python("/home/zhengjie/software/anaconda3/envs/stlearn/bin/python")
  wd<-getwd()
  setwd("./CaDRReS-Sc/notebooks")
  
  if(mode=="cell"){
    mtx<-GetAssayData(sc_dataset,assay = "RNA",slot = "data")
    mtx<-as.matrix(mtx)
    write.table(mtx,file = "../data/sc_input.tsv",sep = "\t")
    prop<-data.frame(patient_id=rep("sc_data",ncol(mtx)),	cluster=colnames(mtx),	percent=rep(1/ncol(mtx),ncol(mtx)))
    write.xlsx(prop,file = "../data/patient/percent_scdata_cluster.xlsx")
    
  }
  
  if(is.null(drug.choose)){
    # multi drug----
    # message("Run CaDRReS_Sc in all drugs.")
    
    tic()
    if(CaDRReS_Sc.GDSC2){
      message(paste0("Run CaDRReS_Sc in all drugs in ",sample.choose," [GDSC2]"))
      system(command = paste0("python train_model_gdsc2.py ",sample.choose))
      system(command = "python predict_sc_comb_gdsc2.py")
      drug.info<-read.csv("../preprocessed_data/GDSC/drug_stat2.csv",row.names = 1,check.names = F)
    }else{
      message(paste0("Run CaDRReS_Sc in all drugs in ",sample.choose," [default]"))
      system(command = paste0("python train_model.py ",sample.choose))
      system(command = "python predict_sc_comb.py")
      drug.info<-read.csv("../preprocessed_data/GDSC/drug_stat.csv",row.names = 1,check.names = F)
    }
    
    # record running time
    k<-toc()
    print(as.numeric(k$toc-k$tic))
    pred.out<-read.csv("../example_result/cadrres-wo-sample-bias_test_pred.csv",check.names = F,row.names = 1)%>%as.data.frame()
    
    rank<-read.csv("../example_result/pred_drug_rank.csv",check.names = F)%>%as.data.frame()%>%arrange(desc(cell_death))
    comb<-read.csv("../example_result/pred_comb.csv",check.names = F)%>%as.data.frame()%>%arrange(desc(cell_death_combi))
    
    # drug.info<-read.csv("../preprocessed_data/GDSC/drug_stat.csv",row.names = 1,check.names = F)
    drug.info.max<-drug.info%>%select("log2_max_conc")
    trans.pred<-c()
    for (id in colnames(pred.out)) {
      trans.pred.col<-pred.out[,id]>=drug.info.max[id,]
      trans.pred<-cbind(trans.pred,trans.pred.col)
    }
    trans.pred<-data.frame(trans.pred)
    trans.pred[trans.pred == 'TRUE'] <- "resistant"
    trans.pred[trans.pred == 'FALSE'] <- "sensitive"
    colnames(trans.pred)<-colnames(pred.out)
    rownames(trans.pred)<-rownames(pred.out)
    
    sc_dataset@meta.data<-cbind(sc_dataset@meta.data,trans.pred[rownames(sc_dataset@meta.data),])
    
    setwd(wd)
    return(list(sc_dataset=sc_dataset,time=as.numeric(k$toc-k$tic),rank=rank,comb=comb))
    
  }else{
    # single-drug----
    # message(paste0("Run CaDRReS_Sc in drug: ",drug.choose,"."))
    drugdata<-filterGDSC(gdsc,drug.choose = drug.choose,sample.choose,by.id = FALSE)
    
    tic()
    if(CaDRReS_Sc.GDSC2){
      message(paste0("Run CaDRReS_Sc in drug: ",drug.choose," in ",sample.choose, "[GDSC2]"))
      system(command = paste0("python train_model_gdsc2.py ",sample.choose))
      system(command = "python predict_sc.py")
      drug.info<-read.csv("../preprocessed_data/GDSC/drug_stat2.csv",row.names = 1,check.names = F)
    }else{
      message(paste0("Run CaDRReS_Sc in drug: ",drug.choose," in ",sample.choose, "[default]"))
      system(command = paste0("python train_model.py ",sample.choose))
      system(command = "python predict_sc.py")
      drug.info<-read.csv("../preprocessed_data/GDSC/drug_stat.csv",row.names = 1,check.names = F)
    }
    
    # record running time
    k<-toc()
    print(as.numeric(k$toc-k$tic))
    pred.out<-read.csv("../example_result/cadrres-wo-sample-bias_test_pred.csv",check.names = F,row.names = 1)%>%as.data.frame()
    
   
    index<-drugdata$pham$DRUG_ID[1]%>%as.character()
    index2<-grep(drug.choose,drug.info$`Drug Name`)%>%as.character()
    if(index2!=index){
      index<-index2
    }
    drug.info.max<-drug.info%>%select("log2_max_conc")
    max_conc<-drug.info.max[index,]
    
    used.drug<-pred.out%>%select(index)
    sc_dataset$CaDRReS_Sc<-"none"
    
    tryCatch({
      sc_dataset@meta.data[pred.out%>%select(index)%>%filter(index>=max_conc)%>%rownames(),]$CaDRReS_Sc<-"resistant"
    }, warning = function(w){
      message("Warning in assignment resistant cells")
    }, error = function(e){
      message("No cell be resistant!")
    })
    
    tryCatch({
      sc_dataset@meta.data[pred.out%>%select(index)%>%filter(index<max_conc)%>%rownames(),]$CaDRReS_Sc<-"sensitive"
    }, warning = function(w){
      message("Warning in assignment sensitive cells")
    }, error = function(e){
      message("No cell be sensitive!")
    })
    
    setwd(wd)
    return(list(sc_dataset=sc_dataset,time=as.numeric(k$toc-k$tic)))
  }
  
}
run_SeuratCCA<-function(sc_dataset,gdsc,drug.choose=NULL,sample.choose="BRCA",...){
  
  
  if(is.null(drug.choose)){
    # multi-drug----
    
    message("Run SeuratCCA in all drugs.")
    time.use<-0
    out.metadata<-sc_dataset@meta.data
    drugId.list<-gdsc$drug$DRUG_ID%>%unique()
    pb <- txtProgressBar(style=3)
    len<-length(drugId.list)
    for(d in c(1:len)){
      setTxtProgressBar(pb, d/len)
      drugId<-drugId.list[d]
      drugdata<-filterGDSC(gdsc,drug.choose = drugId,sample.choose,by.id = TRUE)
      drugName<-drugdata$pham$DRUG_NAME[1]
      metadata<-drugdata$pham%>%dplyr::select(CELL_LINE_NAME, SANGER_MODEL_ID,TCGA_DESC, DRUG_ID,DRUG_NAME,MAX_CONC,LN_IC50,AUC,TCGA_DESC, drug_group)
      sc_dataset.ccl<-CreateSeuratObject(counts = drugdata$expr.count,meta.data = metadata)
      sc_dataset.ccl<-sc_process(sc_dataset.ccl)
      
      # record running time 
      tic()
      anchors <- FindTransferAnchors(
        reference = sc_dataset.ccl,
        query = sc_dataset,
        k.anchor = if(dim(sc_dataset.ccl)[2]>5){5}else{dim(sc_dataset.ccl)[2]-1},
        k.score = if(dim(sc_dataset.ccl)[2]>30){30}else{dim(sc_dataset.ccl)[2]-1},
        # k.filter = if(dim(sc_dataset.ccl)[2]>5){dim(sc_dataset.ccl)[2]-1}else{NA},
        dims = if(dim(sc_dataset.ccl)[2]>30){1:30}else{1:(dim(sc_dataset.ccl)[2]-1)},
        npcs = if(dim(sc_dataset.ccl)[2]>30){30}else{dim(sc_dataset.ccl)[2]-1}
      )
      
      sc_dataset <- MapQuery(
        anchorset = anchors,
        transferdata.args = list(k.weight = if((anchors@anchors%>%nrow())>50){50}else{(anchors@anchors%>%nrow())-1}),
        reference = sc_dataset.ccl,
        query = sc_dataset,
        refdata = c('drug_group')
      )
      # record running time 
      k<-toc()

      
      out.metadata[[paste0("SeuratCCA_label_",drugId,"_",drugName)]]<-sc_dataset$predicted.id
      
      time.use<-time.use+k$toc-k$tic
    }
    close(pb)
    rank<-scPharm_rank_method(out.metadata)
    sc_dataset@meta.data<-out.metadata[sc_dataset@meta.data%>%rownames(),]
    return(list(sc_dataset=sc_dataset,time=as.numeric(time.use),rank=rank,comb=""))
    
    
  }else{
    # single-drug ----
    message(paste0("Run SeuratCCA in drug: ",drug.choose,"."))
    drugdata<-filterGDSC(gdsc,drug.choose = drug.choose,sample.choose,by.id = FALSE)
    
    metadata<-drugdata$pham%>%dplyr::select(CELL_LINE_NAME, SANGER_MODEL_ID,TCGA_DESC, DRUG_ID,DRUG_NAME,MAX_CONC,LN_IC50,AUC,TCGA_DESC, drug_group)
    sc_dataset.ccl<-CreateSeuratObject(counts = drugdata$expr.count,meta.data = metadata)
    sc_dataset.ccl<-sc_process(sc_dataset.ccl)
    
    # record running time 
    tic()
    anchors <- FindTransferAnchors(
      reference = sc_dataset.ccl,
      query = sc_dataset,
      dims = 1:30
    )
    
    sc_dataset <- MapQuery(
      anchorset = anchors,
      reference = sc_dataset.ccl,
      query = sc_dataset,
      refdata = c('drug_group')
    )
    # record running time 
    k<-toc()
    print(as.numeric(k$toc-k$tic))
    sc_dataset$SeuratCCA<-sc_dataset$predicted.id
    return(list(sc_dataset=sc_dataset,time=as.numeric(k$toc-k$tic)))
    
  }
  
}
run_scPharm<-function(sc_dataset,drug.choose=NULL,type="BRCA",tissue.or.cellline="tissue",...){
  wd<-getwd()
  setwd("./scPharm")
  sc.file<-"./data/scdata.rds"
  saveRDS(sc_dataset,file = sc.file)
  if(is.null(drug.choose)){ # NULL mean all drug
    message("Run scPharm return all drug out.")
    # record running time 
    tic()
    system(command = paste0("Rscript scPharm.R -f ","scdata.rds"," -t ",type," -c ",tissue.or.cellline))
    # record running time 
    k<-toc()
    print(as.numeric(k$toc-k$tic))
    

    out.data<-readRDS("./result/scdata_scPharm_object_nmcs_50_nfs_200.rds")
    rank<-read.csv('./result/scdata_scPharm_drug_rank.csv')
    comb<-read.csv('./result/scdata_scPharm_drug_combo.csv')
    
    out.data@meta.data[out.data@meta.data=="other"]<-"none"

    
    setwd(wd)
    return(list(sc_dataset=out.data,time=as.numeric(k$toc-k$tic),rank=rank,comb=comb))
    
    
  }else{
    ## single-drug ----
    message(paste0("Run scPharm return ",drug.choose," out."))
    # record running time 
    tic()
    system(command = paste0("Rscript scPharm.R -f ","scdata.rds"," -t ",type," -c ",tissue.or.cellline," -n ",drug.choose))
    # record running time 
    k<-toc()
    print(as.numeric(k$toc-k$tic))
    
    out.data<-readRDS("./result/scdata_scPharm_object_nmcs_50_nfs_200.rds")
    
    index<-grep(drug.choose,out.data@meta.data%>%colnames())[1]
    out.data$scPharm<-out.data@meta.data[,index]
    out.data@meta.data[out.data$scPharm=="other",]$scPharm<-"none"
    
    setwd(wd)
    return(list(sc_dataset=out.data,time=as.numeric(k$toc-k$tic)))
  }
  
}

# Uniform Calling Interface
benchmark<-function(sc_dataset,gdsc = gdsc,func="Scissor",drug.choose,scPharm.type="BRCA",CaDRReS_Sc.type="",CaDRReS_Sc.GDSC2=FALSE,scPharm.tissue.or.cellline="tissue",plot=FALSE){
  if(func=="Scissor"){
    message(paste0("Run ",func))
    return(run_Scissor(sc_dataset=sc_dataset,gdsc = gdsc,drug.choose=drug.choose,sample.choose=scPharm.type))
  }
  if(func=="scDEAL"){
    message(paste0("Run ",func))
    return(run_scDEAL(sc_dataset=sc_dataset,gdsc = gdsc,drug.choose=drug.choose))
  }
  if(func=="SeuratCCA"){
    message(paste0("Run ",func))
    return(run_SeuratCCA(sc_dataset,gdsc = gdsc,drug.choose=drug.choose,sample.choose=scPharm.type))
  }
  if(func=="scPharm"){
    message(paste0("Run ",func))
    return(run_scPharm(sc_dataset=sc_dataset,drug.choose=drug.choose,type = scPharm.type,tissue.or.cellline=scPharm.tissue.or.cellline))
  }
  if(func=="CaDRReS-Sc"){
    message(paste0("Run ",func))
    return(run_CaDRReS_Sc(sc_dataset=sc_dataset,gdsc = gdsc,sample.choose=CaDRReS_Sc.type,CaDRReS_Sc.GDSC2 = CaDRReS_Sc.GDSC2,drug.choose=drug.choose))
  }
}

# Evaluation functions for sampled data
benchmark_cal<-function(sc_dataset,gdsc = gdsc,drug.choose,func=c("Scissor","scDEAL","CaDRReS-Sc","SeuratCCA","scPharm"),
                        CaDRReS_Sc.GDSC2=FALSE,
                        scPharm.tissue.or.cellline="tissue",
                        CaDRReS_Sc.type="",
                        scPharm.type="BRCA",
                        sampling=c(100,200,400,800,1600,"all"),
                        process=TRUE,times=3,save.data=FALSE){
  
  logging("Start benchmark:","",FALSE)
  metalist<-list()
  dfall<-c()
  ranklist<-list()
  comblist<-list()
  plotlist<-list()
  if(class(sc_dataset)=="list"){  # for list datasets
    message("List data input mode.")
    print(names(sc_dataset))
    
    for (scdata in names(sc_dataset)) {
      message(paste0("Start data: ",scdata))
      
      df<-c()
      data<-sc_dataset[[scdata]]
      print(data)
      for (s in sampling) {
        
        for (f in func) {
          run_times<-c()
        
          #20230629 update
          if(s=="all"){
            data<-data
            s<-"all"
          }else{
            s<-as.numeric(s)
            data<-subset(data,downsample=s)
          }
          
          if(process){
            data<-data%>%sc_process()
          }
          
          for (t in c(1:times)) {
            string<-paste0("data:",scdata,", Method:",f,", Sampling:",s,", Times:",t)
            savename<-paste0(scdata,"_",f,"_",s,"_",t)
            message(string)
            logging(string,"",append = TRUE)
            tryCatch(
              {
                bench_out<-benchmark(data,gdsc = gdsc,func=f,drug.choose,plot=FALSE,scPharm.type = scPharm.type,CaDRReS_Sc.type=CaDRReS_Sc.type, CaDRReS_Sc.GDSC2= CaDRReS_Sc.GDSC2,scPharm.tissue.or.cellline =scPharm.tissue.or.cellline )
                # for test
                # bench_out<-list(sc_dataset=sc_dataset,time=rnorm(1,mean = 100,sd=10))
                if(save.data){
                  write.csv(bench_out$rank,file = paste0(savename,".csv"))
                }
                run_times<-c(run_times,bench_out$time)
                logging(string,paste0(bench_out$time,"s"),append = TRUE)
                
                metalist[[savename]]<-bench_out$sc_dataset@meta.data
                
                if(is.null(drug.choose)){
                  ranklist[[savename]]<-bench_out$rank
                  comblist[[savename]]<-bench_out$comb
                }
              }, error = function(e){
                message("No cell be sensitive!")
                logging(string,"Error!",append = TRUE)
              })
          }
          sd<-sd(run_times)
          mean<-mean(run_times)
          df<-rbind(df,c(scdata,f,s,mean,sd))
          
        }
        
      }
      df<-data.frame(df)
      colnames(df)<-c("Data","Method","cells",'time','sd')
      df$time<-as.numeric(df$time)
      df$sd<-as.numeric(df$sd)
      # df$cells<-as.numeric(df$cells)
      df$cells<-factor(df$cells,levels = unique(df$cells))
      dfall<-rbind(dfall,df)
      
      p<-ggplot(df,aes(cells,time,group=Method,color=Method,shape=Method))+
        geom_point(size=3)+
        geom_line(cex=0.5)+
        geom_errorbar(aes(ymin = time - sd, ymax = time + sd), 
                      width = 0.1,cex=0.5)+
        theme_bw()+
        ylab("Time(s)")+
        xlab("Number of cells")
      
      plotlist[[savename]]<-p
    }
    
    
    
    return(list(metalist=metalist,df=dfall,plot=plotlist,comb=comblist,rank=ranklist))
    
  }else if(class(sc_dataset)=="Seurat"){  # for one data
    
    for (s in sampling) {
      
      for (f in func) {
        run_times<-c()
        
        #20230629 update
        if(s=="all"){
          data<-sc_dataset
        }else{
          data<-subset(sc_dataset,downsample=as.numeric(s))
        }
        
        if(process){
          data<-data%>%sc_process()
        }
        
        for (t in c(1:times)) {
          message(paste0("Method:",f,", Sampling:",s,", Times:",t))
          logging(paste0("Method:",f,", Sampling:",s,", Times:",t),"",append = TRUE)
          tryCatch(
            {
              bench_out<-benchmark(data,drugdata,func=f,drug.choose,plot=FALSE,scPharm.type = scPharm.type)
              # for test
              # bench_out<-list(sc_dataset=sc_dataset,time=rnorm(1,mean = 100,sd=10))
              run_times<-c(run_times,bench_out$time)
              metalist[[paste0(f,"_",s,"_",t)]]<-bench_out$sc_dataset@meta.data
              if(is.null(drug.choose)){
                ranklist[[paste0(f,"_",s,"_",t)]]<-bench_out$rank
                comblist[[paste0(f,"_",s,"_",t)]]<-bench_out$comb
              }
            }, error = function(e){
              message("No cell be sensitive!")
              logging(paste0("Method:",f,", Sampling:",s,", Times:",t),"Error!",append = TRUE)
            })
        }
        sd<-sd(run_times)
        mean<-mean(run_times)
        df<-rbind(df,c(f,s,mean,sd))
      }
      
    }
    
    df<-data.frame(df)
    colnames(df)<-c("Method","cells",'time','sd')
    df$time<-as.numeric(df$time)
    df$sd<-as.numeric(df$sd)
    # df$cells<-as.numeric(df$cells)
    df$cells<-factor(df$cells,levels = unique(df$cells))
    
    p<-ggplot(df,aes(cells,time,group=Method,color=Method,shape=Method))+
      geom_point(size=3)+
      geom_line(cex=0.5)+
      geom_errorbar(aes(ymin = time - sd, ymax = time + sd), 
                    width = 0.1,cex=0.5)+
      theme_bw()+
      ylab("Time(s)")+
      xlab("Number of cells")
    
    return(list(metalist=metalist,df=df,plot=p,comb=comblist,rank=ranklist))
  }

}

# for evaluating 
# need define the $true label of sc_dataset
benchmark_evaluate<-function(benchmark.cal.metalist=out$metalist,sc_dataset,plot.file="evaluate.pdf"){
  skm<-import('sklearn.metrics')
  
  evl_df<-c()
  for(sub.metadata in names(benchmark.cal.metalist)){
    d<-benchmark.cal.metalist[[sub.metadata]]
    
    d$true<-sc_dataset@meta.data[rownames(d),]$true
    
    if(length(grep("scDEAL",sub.metadata))>0){
      d$pred<-d$scDEAL
      me<-"scDEAL"
    }
    if(length(grep("Scissor",sub.metadata))>0){
      d$pred<-d$Scissor
      me<-"Scissor"
    }
    if(length(grep("CaDRReS-Sc",sub.metadata))>0){
      d$pred<-d$CaDRReS_Sc
      me<-"CaDRReS-Sc"
    }
    if(length(grep("SeuratCCA",sub.metadata))>0){
      d$pred<-d$SeuratCCA
      me<-"SeuratCCA"
    }
    if(length(grep("scPharm",sub.metadata))>0){
      d$pred<-d$scPharm
      me<-"scPharm"
    }
    pred<-d$pred
    true<-d$true
    
    if(c(pred,true)%>%unique()%>%length()<3){
      o<-skm$classification_report(y_true = true,y_pred = pred,target_names=c('A', 'B'))
    }else{
      o<-skm$classification_report(y_true = true,y_pred = pred,target_names=c('A', 'B', 'C'))
    }
    o<-strsplit(o,split = "\n")[[1]]
    o<-o[length(o)]
    o<-strsplit(o,split = "      ")[[1]]%>%as.character()%>%as.numeric()
    precision<-o[2]
    recall<-o[3]
    f1<-o[4]
    row<-c(me,strsplit(sub.metadata,split = "_")[[1]][2],precision,recall,f1)
    evl_df<-rbind(evl_df,row)
  }
  evl_df<-evl_df%>%unique.data.frame()%>%as.data.frame()
  colnames(evl_df)<-c("Method","Cells","Precision","Recall","F1-score")
  evl_df$Method<-factor(evl_df$Method)
  evl_df$Cells<-as.numeric(evl_df$Cells)
  evl_df$Precision<-as.numeric(evl_df$Precision)
  evl_df$Recall<-as.numeric(evl_df$Recall)
  evl_df$`F1-score`<-as.numeric(evl_df$`F1-score`)
  
  p1<-ggplot(evl_df,aes(x=Cells,y=Precision,group=Method,fill=Method,shape=Method))+
    geom_point(size=3)+
    geom_line(position = position_dodge(0.1),cex=0.5)+
    theme_classic()+
    xlab("Number of cells")
  p2<-ggplot(evl_df,aes(Cells,Recall,group=Method,fill=Method,shape=Method))+
    geom_point(size=3)+
    geom_line(position = position_dodge(0.1),cex=0.5)+
    theme_classic()+
    xlab("Number of cells")
  p3<-ggplot(evl_df,aes(Cells,`F1-score`,group=Method,fill=Method,shape=Method))+
    geom_point(size=3)+
    geom_line(position = position_dodge(0.1),cex=0.5)+
    theme_classic()+
    scale_fill_manual(values = RColorBrewer::brewer.pal(7,"Set1"))+
    xlab("Number of cells")
  ggsave(cowplot::plot_grid(p1,p2,p3,ncol = 3),filename = plot.file,width = 15,height = 3)
  return(evl_df)
}

# noise benchmark function
addNoise2Data<-function(scdata,
                        ratio=0.1, # noise ratio
                        mode="dropout", # dropout cell expression to 0；
                        do.sc.process=TRUE
){
  if(mode=="dropout"){
    num.gene<-dim(scdata)[1]
    num.cell<-dim(scdata)[2]
    print(paste0("all position: ",num.gene*num.cell))
    o.num<-ceiling(ratio*num.gene*num.cell)
    
    print(paste0("noise position: ",o.num))
    for (g in c(1:o.num)) {
      scdata@assays$RNA@counts[sample(1:num.gene, 1),sample(1:num.cell, 1)]<-0
    }
    if(do.sc.process){
      scdata<-sc_process(scdata)
    }
  }
  
  return(scdata)
}

# Functions for noisy data
benchmark_noise_cal<-function(sc_dataset,drugdata,drug.choose,
                              func=c("Scissor","scDEAL","CaDRReS-Sc","SeuratCCA","scPharm"),
                              noise.ratio=c(0,0.001,0.01,0.1),
                              scPharm.type="LUAD",
                              sampling=c(1600),
                              log.file=paste0(dir,"/benchmark_noise.log"),
                              times=3)
{
  logging("Start benchmark with noise:","",FALSE,file =log.file )
  metalist<-list()
  df<-c()
  for (n in noise.ratio) {
    
    for (s in sampling) {
      
      for (f in func) {
        run_times<-c()
        data<-subset(sc_dataset,downsample=s)%>%addNoise2Data(.,ratio = n,do.sc.process = TRUE)
        
        for (t in c(1:times)) {
          message(paste0("Method:",f,", Sampling:",s,", Times:",t,", Noise:",n))
          logging(paste0("Method:",f,", Sampling:",s,", Times:",t,", Noise:",n),"",append = TRUE,file =log.file)
          tryCatch(
            {
              bench_out<-benchmark(data,drugdata,func=f,drug.choose,plot=FALSE,scPharm.type = scPharm.type)
              logging(paste0("Method:",f,", Sampling:",s,", Times:",t,", Noise:",n),paste0("Time use:",bench_out$time,"s"),append = TRUE,file =log.file)
              # for test
              # bench_out<-list(sc_dataset=sc_dataset,time=rnorm(1,mean = 100,sd=10))
              run_times<-c(run_times,bench_out$time)
              metalist[[paste0(n,"_",f,"_",s,"_",t)]]<-bench_out$sc_dataset@meta.data
            }, error = function(e){
              message("No cell be sensitive!")
              logging(paste0("Method:",f,", Sampling:",s,", Times:",t,", Noise:",n),"Error!",append = TRUE,file =log.file)
            })
        }
        sd<-sd(run_times)
        mean<-mean(run_times)
        df<-rbind(df,c(n,f,s,mean,sd))
      }
      
    }
  }
  df<-data.frame(df)
  colnames(df)<-c("Noise","Method","cells",'time','sd')
  df$time<-as.numeric(df$time)
  df$sd<-as.numeric(df$sd)
  df$cells<-as.numeric(df$cells)
  df$Noise<-factor(df$Noise,levels = unique(df$Noise))
  
  p<-ggplot(df,aes(Noise,time,group=Method,color=Method,shape=Method))+
    geom_point(size=3)+
    geom_line(cex=0.5)+
    geom_errorbar(aes(ymin = time - sd, ymax = time + sd), 
                  width = 0.1,cex=0.5)+
    theme_bw()+
    ylab("Time(s)")+
    xlab("Noise ratio")+
    ggtitle("Dropout Noise benchmark")
  
  return(list(metalist=metalist,df=df,plot=p))
}

# Evaluation functions for noisy data
benchmark_evaluate_noise<-function(benchmark.cal.metalist=out$metalist,sc_dataset,plot.file="evaluate.noise.pdf",xlabel="Noise"){
  skm<-import('sklearn.metrics')
  
  evl_df<-c()
  for(sub.metadata in names(benchmark.cal.metalist)){
    print(sub.metadata)
    d<-benchmark.cal.metalist[[sub.metadata]]
    
    d$true<-sc_dataset@meta.data[rownames(d),]$true
    
    if(length(grep("scDEAL",sub.metadata))>0){
      d$pred<-d$scDEAL
      me<-"scDEAL"
    }
    if(length(grep("Scissor",sub.metadata))>0){
      d$pred<-d$Scissor
      me<-"Scissor"
    }
    if(length(grep("CaDRReS-Sc",sub.metadata))>0){
      d$pred<-d$CaDRReS_Sc
      me<-"CaDRReS-Sc"
    }
    if(length(grep("SeuratCCA",sub.metadata))>0){
      d$pred<-d$SeuratCCA
      me<-"SeuratCCA"
    }
    if(length(grep("scPharm",sub.metadata))>0){
      d$pred<-d$scPharm
      me<-"scPharm"
    }
    pred<-d$pred
    true<-d$true
    if(c(pred,true)%>%unique()%>%length()<2){
      o<-skm$classification_report(y_true = true,y_pred = pred,target_names=c('A'))
    }else if(c(pred,true)%>%unique()%>%length()<3){
      o<-skm$classification_report(y_true = true,y_pred = pred,target_names=c('A', 'B'))
    }else{
      o<-skm$classification_report(y_true = true,y_pred = pred,target_names=c('A', 'B', 'C'))
    }
    
    o<-strsplit(o,split = "\n")[[1]]
    o<-o[length(o)]
    o<-strsplit(o,split = "     ")[[1]]%>%as.character()%>%as.numeric()
    precision<-o[2]
    recall<-o[3]
    f1<-o[4]
    row<-c(me,strsplit(sub.metadata,split = "_")[[1]][1],precision,recall,f1)
    evl_df<-rbind(evl_df,row)
  }
  evl_df<-evl_df%>%unique.data.frame()%>%as.data.frame()
  colnames(evl_df)<-c("Method","Cells","Precision","Recall","F1-score")
  evl_df$Method<-factor(evl_df$Method)
  evl_df$Cells<-factor(evl_df$Cells,levels = unique(evl_df$Cells))
  evl_df$Precision<-as.numeric(evl_df$Precision)
  evl_df$Recall<-as.numeric(evl_df$Recall)
  evl_df$`F1-score`<-as.numeric(evl_df$`F1-score`)
  
  p1<-ggplot(evl_df,aes(x=Cells,y=Precision,group=Method,fill=Method,shape=Method))+
    geom_point(size=3)+
    geom_line(cex=0.5)+
    theme_classic()+
    xlab(xlabel)
  p2<-ggplot(evl_df,aes(Cells,Recall,group=Method,fill=Method,shape=Method))+
    geom_point(size=3)+
    geom_line(cex=0.5)+
    theme_classic()+
    
    xlab(xlabel)
  p3<-ggplot(evl_df,aes(Cells,`F1-score`,group=Method,fill=Method,shape=Method))+
    geom_point(size=3)+
    geom_line(cex=0.5)+
    theme_classic()+
    
    xlab(xlabel)
  ggsave(cowplot::plot_grid(p1,p2,p3,ncol = 3),filename = plot.file,width = 15,height = 3)
  return(evl_df)
}


# Split data by sensitive/resistant ratio
splitDataByLabel<-function(scdata,
                           label="true", #colname
                           ratio=1, #sensitive / rensistant
                           do.sc.process=TRUE)
{
  #origin
  scdata$true<-scdata@meta.data%>%pull(label)
  
  df<-scdata@meta.data$true%>%table()%>%data.frame()%>%tibble::column_to_rownames(colnames(.)[1])
  print(paste0("origin ratio: ",signif(df['sensitive',]/df['resistant',],3)))
  print(paste0("sensitive: ",df['sensitive',]," resistant: ",df['resistant',]))
  o.r<-0
  o.s<-0
  n.r<-df['resistant',]
  n.s<-df['sensitive',]
  
  
  if(ratio==0){
    o.r<-n.r
    o.s<-0
  }else if(ratio==Inf){
    o.r<-0
    o.s<-n.s
  }else{
    if(n.r*ratio<=n.s){
      o.r<-n.r
      o.s<-ceiling(n.r*ratio)
    }else{
      o.r<-ceiling(n.s/ratio)
      o.s<-n.s
    }
  }
  print(paste0("out.sensitive: ",o.s,", out.resistant: ",o.r))
  scdata.r<-subset(scdata,true=="resistant")
  scdata.s<-subset(scdata,true=="sensitive")
  if(scdata$true%>%unique()%>%length()==3){
    n.n<-df['none',]
    scdata.n<-subset(scdata,true=="none")
  }
  
  if(o.r!=0){
    scdata.r.out<-subset(scdata.r,downsample=o.r)
  }
  if(o.s!=0){
    scdata.s.out<-subset(scdata.s,downsample=o.s)
  }
  
  if(scdata$true%>%unique()%>%length()==3){
    if(o.s==0){
      scdata.o<-merge(scdata.n,scdata.r.out)
    }else if(o.r==0){
      scdata.o<-merge(scdata.n,scdata.s.out)
    }else{
      scdata.o<-merge(scdata.n,list(scdata.s.out,scdata.r.out))
    }
  }else{
    if(o.s==0){
      scdata.o<-scdata.r.out
    }else if(o.r==0){
      scdata.o<-scdata.s.out
    }else{
      scdata.o<-merge(scdata.s.out,scdata.r.out)
    }
  }
  
  
  
  if(do.sc.process){
    scdata.o<-sc_process(scdata.o)
  }
  
  return(scdata.o)
}

# Benchmark for sensitive/resistant ratio
benchmark_ratio_cal<-function(sc_dataset,drugdata,drug.choose,
                              func=c("Scissor","scDEAL","CaDRReS-Sc","SeuratCCA","scPharm"),
                              split.ratio=c(0,0.25,0.5,0.75,1,1.33,2,4,Inf),
                              scPharm.type="LUAD",
                              sampling=c(1600),
                              log.file=paste0(dir,"/benchmark_ratio.log"),
                              times=3)
{
  logging("Start benchmark by diff ratio:","",FALSE,file =log.file )
  metalist<-list()
  df<-c()
  for (n in split.ratio) {
    
    for (s in sampling) {
      
      for (f in func) {
        run_times<-c()
        data<-subset(sc_dataset,downsample=s)%>%splitDataByLabel(.,label="true",ratio=n,do.sc.process=TRUE)
        
        for (t in c(1:times)) {
          message(paste0("Method:",f,", Sampling:",s,", Times:",t,", Ratio:",n))
          logging(paste0("Method:",f,", Sampling:",s,", Times:",t,", Ratio:",n),"",append = TRUE,file =log.file)
          tryCatch(
            {
              bench_out<-benchmark(data,drugdata,func=f,drug.choose,plot=FALSE,scPharm.type = scPharm.type)
              logging(paste0("Method:",f,", Sampling:",s,", Times:",t,", Ratio:",n),paste0("Time use:",bench_out$time,"s"),append = TRUE,file =log.file)
              # for test
              # bench_out<-list(sc_dataset=sc_dataset,time=rnorm(1,mean = 100,sd=10))
              run_times<-c(run_times,bench_out$time)
              metalist[[paste0(n,"_",f,"_",s,"_",t)]]<-bench_out$sc_dataset@meta.data
            }, error = function(e){
              message("No cell be sensitive!")
              logging(paste0("Method:",f,", Sampling:",s,", Times:",t,", Ratio:",n),"Error!",append = TRUE,file =log.file)
            })
        }
        sd<-sd(run_times)
        mean<-mean(run_times)
        df<-rbind(df,c(n,f,s,mean,sd))
      }
      
    }
  }
  df<-data.frame(df)
  colnames(df)<-c("Ratio","Method","cells",'time','sd')
  df$time<-as.numeric(df$time)
  df$sd<-as.numeric(df$sd)
  df$cells<-as.numeric(df$cells)
  df$Ratio<-factor(df$Ratio,levels = unique(df$Ratio))
  
  p<-ggplot(df,aes(Ratio,time,group=Method,color=Method,shape=Method))+
    geom_point(size=3)+
    geom_line(cex=0.5)+
    geom_errorbar(aes(ymin = time - sd, ymax = time + sd), 
                  width = 0.1,cex=0.5)+
    theme_bw()+
    ylab("Time(s)")+
    xlab("Sample ratio")+
    ggtitle("Sample ratio benchmark")
  
  return(list(metalist=metalist,df=df,plot=p))
}

# Levenshtein distence 
LevenshteinDistance <- function(list1, list2) {
  len1 <- length(list1)
  len2 <- length(list2)
  

  dp <- matrix(0, nrow = len1 + 1, ncol = len2 + 1)
  

  for (i in 1:(len1 + 1)) {
    dp[i, 1] <- i - 1
  }
  for (j in 1:(len2 + 1)) {
    dp[1, j] <- j - 1
  }
  

  for (i in 2:(len1 + 1)) {
    for (j in 2:(len2 + 1)) {
      if (list1[i - 1] == list2[j - 1]) {
        dp[i, j] <- dp[i - 1, j - 1]
      } else {
        dp[i, j] <- min(dp[i - 1, j], dp[i, j - 1], dp[i - 1, j - 1]) + 1
      }
    }
  }
  

  return(dp[len1 + 1, len2 + 1])
}

## use single to assignment cell 
# ref is MouseRNAseqData()
assign_singler<-function(data){
  library(celldex)
  library(SingleR)
  ref<-MouseRNAseqData()
  
  mtx<-GetAssayData(data,slot = "data")%>%as.matrix()
  sr<-SingleR(test = mtx,ref =ref ,clusters = data$seurat_clusters,labels = ref$label.main)
  # sr$labels
  data$celltype_singler<-""
  for (cl in rownames(sr)) {
    print(cl)
    data@meta.data[data@meta.data$seurat_clusters==cl,]$celltype_singler<-sr[cl,"labels"]
  }
  
}

# homolog gene transfor
mus2hsa<-function(genelist){
  library(homologene)
  mus2hsa<-homologene(genelist, inTax = 10090, outTax = 9606)
  mus2hsa<-mus2hsa%>%distinct(`10090`,.keep_all = T)
  rownames(mus2hsa)<-mus2hsa$`10090`
  colnames(mus2hsa)<-c("mus","hsa","mus.id","hsa.id")
  
  gene<-setdiff(genelist,rownames(mus2hsa))
  mus2hsa_supp<-data.frame(mus=gene,hsa=toupper(gene),`mus.id`=rep(0,length(gene)),`hsa.id`=rep(1,length(gene)))
  mus2hsa<-rbind(mus2hsa,mus2hsa_supp)
  rownames(mus2hsa)<-mus2hsa$mus
  return(mus2hsa[genelist,]$hsa)
}

#### multi drug benchmark supp  ####
library(Scissor)
library(tictoc)
library(readxl)
library(openxlsx)
library(tidyverse)
library(ggplot2)
library(ggsignif)
library(cowplot)


#drug source file
{
  expr_count.file <- "GDSC/rnaseq_read_count_20220624.csv"
  expr_tpm.file   <- "GDSC/rnaseq_tpm_20220624.csv"
  GDSC1.file      <- "GDSC/GDSC1_fitted_dose_response_24Jul22.xlsx"
  GDSC2.file      <- "GDSC/GDSC2_fitted_dose_response_24Jul22.xlsx"
  sample.file     <- "GDSC/model_list_20221102.csv"
  gene.file       <- "GDSC/gene_identifiers_20191101.csv"
  sample.choose   <- "BRCA" # TCGA_DESC
}


# Calculation of fractions in tumour cells(Rd)
scPharm_rank_method<-function(meta.data,cell.orig='tissue',in.tumor=TRUE){
  message("computing drug rank")
  if (cell.orig != "tissue") {
    meta.data$cell.label = "tumor"
  }
  result = meta.data[,c("cell.label", grep("_label",colnames(meta.data),value = T))]
  if(in.tumor){
    result = result[result$cell.label == "tumor", !names(result) %in% "cell.label"]
    
  }else{
    result = result[, !names(result) %in% "cell.label"]
    
  }

  
  sensi.ratio = colMeans(result == "sensitive")
  resis.ratio = colMeans(result == "resistant")
  score = data.frame(DRUG_ID = sapply(strsplit(names(sensi.ratio), split = "_"), function(x) {x[3]}),
                     DRUG_NAME = sapply(strsplit(names(sensi.ratio), split = "_"), function(x) {
                       if (length(x) ==4) {
                         return(x[4])
                       }else {
                         return(paste(x[4:length(x)], collapse = "_"))
                       }
                     }
                     ),
                     SENSI_RATIO = sensi.ratio,
                     RESIS_RATIO = resis.ratio)
  score$Score = score$SENSI_RATIO*(1 - score$RESIS_RATIO)
  score = score[order(-score[,5]),]
  score$Rank = seq(1,nrow(score),1)

  return(score)
}
