# Title    : Supplementary Figure 2
# Author   : Wanglab
# Time     : 2024.9


# Supplemental Figure 2a-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# load results using scPharm in various Cell-ID size
nfeatures = readRDS("./rscript/scPharm/nfeatures.rds")
rank_ct = data.frame(matrix(0, nrow = 295, ncol = 30))
colnames(rank_ct) = paste(names(nfeatures),
                          c(rep("100",6), rep("150", 6), rep("200", 6), rep("250", 6), rep("300", 6)), sep = "_")
names(nfeatures) = colnames(rank_ct)

for (sam in names(nfeatures)) {
  data = nfeatures[[sam]][["rank_ct"]]
  rank_ct[data[data$DRUG_ID == 1032,]$Rank, sam] = 1
  rank_ct[data[data$DRUG_ID == 1549,]$Rank, sam] = 2
  rank_ct[data[data$DRUG_ID == 1558,]$Rank, sam] = 3
}
library(ComplexHeatmap)
library(circlize)
row_info = data.frame(sample = names(nfeatures),
                      group = c(rep("100",6), rep("150", 6), rep("200", 6), rep("250", 6), rep("300", 6)))
rownames(row_info) = row_info$sample
row_anno = HeatmapAnnotation(
  df = row_info$group,
  which = "row",
  col = list(
    df = c("100" = "#ABC6E4",
           "150" = "#C39398",
           "200" = "#FCDABA",
           "250" = "#A7D2BA",
           "300" = "#D0CADE"),
    show_legend = TRUE
  ),
  show_annotation_name = FALSE,
  annotation_label = "Cell-ID size",
  annotation_legend_param = list(
    color_bar = "discrete",
    at = c("100","150","200","250","300"),
    labels = c("100","150", "200", "250","300"),
    legend_gp = gpar(fill = c("#ABC6E4","#C39398", "#FCDABA", "#A7D2BA", "#D0CADE")),
    labels_gp = gpar(col = "#000000", fontsize = 12),
    title_gp = gpar(fontsize = 12)
  )
)

Cairo::CairoPDF("./Rplot/review/copykat_true.pdf", width = 10, height = 5)
Heatmap(t(rank_ct),
        # row_split = 1:4,
        # rect_gp = gpar(col = 'white', lwd = 0.1),
        row_split = row_info$group,
        row_gap = unit(1, "mm"),
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
        row_title = NULL,
        show_heatmap_legend = T,
        name = "Drug",
        heatmap_legend_param = list(
          color_bar = "discrete",
          at = 0:3,
          labels = c("other","Afatinib", "Sapitinib", "Lapatinib"),
          legend_gp = gpar(fill = c("grey","#E41A1C", "#377EB8", "#4DAF4A")), 
          labels_gp = gpar(col = "#000000", fontsize = 12),
          title_gp = gpar(fontsize = 12)
        ),
        left_annotation = row_anno,
        row_labels = sapply(strsplit(names(nfeatures), split = "_"), function(x) x[1])
        # width = unit(15, "cm"), height = unit(1.5, "cm")
)
dev.off()


# Supplemental Figure 2b --------------------------------------------------

## the up panel of Supplemental Figure 2b is from Figure 3d

# Supplemental Figure 2b bottom panel

sample = c("AH0308","MH0031","MH0069","MH0161","MH0176","PM0337")
names(sample) = sample
# load HER-positive sample result
GSE161529_HER2 = lapply(sample, function(sam){
  object = readRDS(paste0("./rscript/scPharm/result/", sam, "_scPharm1208_nmcs_50_nfs_200.rds"))
})

new_cut = lapply(GSE161529_HER2, function(object){
  meta.data = object@meta.data
  for (i in seq(7, ncol(meta.data), 2)) {
    meta.data[,i] = "other"
    meta.data[meta.data[,i+1] < -1.452096, i] = "sensitive"
    meta.data[meta.data[,i+1] > 1.43721, i] = "resistant"
  }
  temp = calcu_rank(meta.data)
})
library(ComplexHeatmap)
library(circlize)
supplefig2b_bot = data.frame(matrix(0, nrow = nrow(new_cut[[1]]), ncol = length(new_cut)))
colnames(supplefig2b_bot) = names(new_cut)
for (sam in names(new_cut)) {
  data = new_cut[[sam]]
  supplefig2b_bot[data[data$DRUG_ID == 1032,]$Rank, sam] = 1
  supplefig2b_bot[data[data$DRUG_ID == 1549,]$Rank, sam] = 2
  supplefig2b_bot[data[data$DRUG_ID == 1558,]$Rank, sam] = 3
}

suppleFig2b = Heatmap(t(supplefig2b_bot),
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
Cairo::CairoPDF("./Rplot/review/figures2_bot.pdf", width = 10, height = 1.6)
print(suppleFig2b)
dev.off()

