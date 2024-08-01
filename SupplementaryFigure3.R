# Title    : Supplementary Figure 3
# Author   : Wanglab
# Time     : 2024.7

library(ggplot2)
library(Cairo)

sample = c("AH0308", "MH0031", "MH0069", "MH0161", "MH0176", "PM0337")
names(sample) = sample
# Load the drug prioritization of HER2 positive BRCA samples from the specified path
fig3.data.2 = lapply(sample, function(sam) {
  rank = read.csv(paste0("./scPharm/result/", sam, "_scPharm_drug_rank.csv"))
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

  CairoPDF(paste0("./Figure/SuppleFigure3_", sam, ".pdf"), width = 2.2, height = 5.8)
  print(p)
  dev.off()
}