###single file input
get_all_type_ji <- function(file_name, fpkm_det){
  rt <- readRDS(file = paste('./data/rnaseq/expression/', 'Quartet_RNA_', 
  file_name, '_genes_fpkm.csv.rds',  sep = ""))
  
  get_one_type_ji <- function(stype, df, fpkm_det){
    d_mat <- df[, grep(stype, colnames(df))]
    
    count_ji <- function(combn_num, d_mat, fpkm_det){
      d1 = d_mat[, combn_num[1]]
      d2 = d_mat[, combn_num[2]]
      
      names(d1) <- row.names(d_mat)
      names(d2) <- row.names(d_mat)
      
      dd1 <- d1[which(d1 >= fpkm_det)]
      dd2 <- d2[which(d2 >= fpkm_det)]
      
      ji <- length(intersect(names(dd1), names(dd2)))/length(unique(names(dd1), names(dd2)))
      mat_ji <- c(colnames(d_mat)[combn_num[1]], ji)
    }
    
    mat_ji_sim_list <- data.frame(combn(dim(df[, grep(stype, colnames(df))])[2], 2)) %>% lapply(count_ji, d_mat, fpkm_det)
    mat_ji <- do.call(rbind, mat_ji_sim_list)
    colnames(mat_ji) <- c('sample_id', 'ji_value')
    return(mat_ji)
  }
  
  stype_list <- c('D5', 'D6', 'F7', 'M8')
  mat_ji_list <- lapply(stype_list, get_one_type_ji, rt, fpkm_det)
  mat_ji_all <- do.call(rbind, mat_ji_list)
  return(mat_ji_all)
}
