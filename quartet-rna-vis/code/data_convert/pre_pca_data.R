library(gmodels)

get_pca_data <- function(file_name_list){
  #rt_pca is a dataframe of which row.names is gene name and colnames is samples 
  load_data <- function(file_name){
    rt_pca <- readRDS(file = paste('./data/rnaseq/expression/',  'Quartet_RNA_', 
    file_name, '_genes_fpkm.csv.rds',  sep = ""))
    return(rt_pca)
  }
  
  rt_pca = file_name_list %>% lapply(load_data) %>% dplyr::bind_cols( ) 
  rt_pca_rm = rt_pca[, c(grep("D5", colnames(rt_pca)),
                         grep("D6", colnames(rt_pca)),
                         grep("F7", colnames(rt_pca)),
                         grep("M8", colnames(rt_pca)))]
  rt_pca_d = rt_pca_rm[apply(rt_pca_rm, 1, function(x){median(x) > 0}), ]
  mat_group = matrix(unlist(strsplit(colnames(rt_pca_d), '_')), byrow = TRUE, ncol = 9)
  group_name = mat_group[, 7]
  set_name = paste(mat_group[, 3], mat_group[, 4], mat_group[, 5], mat_group[, 6], sep = "_")
  
  data.a <- t(as.matrix(rt_pca_d))
  data.pca <- fast.prcomp(data.a, scale = T)  # do PCA
  a <- summary(data.pca)
  tmp <- a$importance  # a include 4 sections which contain importance
  pro1 <- as.numeric(sprintf("%.3f", tmp[2, 1]))*100
  pro2 <- as.numeric(sprintf("%.3f", tmp[2, 2]))*100 
  pro3 <- as.numeric(sprintf("%.3f", tmp[2, 3]))*100 # fetch the proportion of PC1, PC2 and PC3
  pc <- as.data.frame(a$x)  # convert to data.frame
  pc$group = group_name
  pc$dataset_name = set_name
  pc_list <- list(pc, pro1, pro2, pro3)
  return(pc_list)
}