get_cor_data <- function(file_name_list){
  #rt_cor is a dataframe of which row.names is gene name and colnames is samples 
  load_data <- function(file_name){
    df <- readRDS(file = paste('./data/rnaseq/expression/',  'Quartet_RNA_', 
    file_name, '_genes_fpkm.csv.rds',  sep = ""))
    return(df)
  }
  
  df = file_name_list %>% lapply(load_data) %>% dplyr::bind_cols( ) 
  df_rm = df[, c(grep("D5", colnames(df)),
                 grep("D6", colnames(df)),
                 grep("F7", colnames(df)),
                 grep("M8", colnames(df)))]
  df_rm_d = df_rm[apply(df_rm, 1, function(x){median(x) > 0}), ]
  return(df_rm_d)
}