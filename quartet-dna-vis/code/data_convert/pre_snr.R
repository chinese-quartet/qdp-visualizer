library(data.table)
source('./code/data_convert/function_calc_SNR.R')

pre_pca_snr <- function(df_exp_annot_sub, dt_meta, zscore_per_group = c('YES', 'NO')){
  ### df_exp_annot is the dataframe from exprSet_RNA_anno.rds
  ### dt_meta from metadata_RNA.rds
  ### group_pca is one of metadata, such as batch
  ### group_pca_unit is the batch unique element, such as B1, B2 and B3
  ### zscore_per_group is used to determine whether the data is zscore
  
  lst_library_forPlot = unique(df_exp_annot_sub$library)

  # choice z-score or not
  df_exp_annot_sub = df_exp_annot_sub[library%in%lst_library_forPlot]
  dt_meta_batch = unique(df_exp_annot_sub[, c('library', 'batch')])
  batch_list <- unique(dt_meta_batch$batch)
  
  #per batch z-score
  print(Sys.time())
  if(zscore_per_group == 'YES'){
    df_exp_annot_sub[,logfpkm_scale:=(logfpkm-mean(logfpkm))/sd(logfpkm), by=.(batch, gene)]
    exprMat_RNA_FPKM_log_scale = dcast(df_exp_annot_sub, gene~library, value.var='logfpkm_scale') %>% data.frame(row.names = 1) %>% as.matrix()
    exprMat_forSignoise = exprMat_RNA_FPKM_log_scale
  } else if (zscore_per_group == 'NO') {
    exprMat_RNA_FPKM = dcast(df_exp_annot_sub, gene~library, value.var='logfpkm') %>% data.frame(row.names = 1) %>% as.matrix()
    exprMat_forSignoise = exprMat_RNA_FPKM 
  } else {
    stop('please set z-score or not')
  }
  
  print(Sys.time())
  
  # prepare pca list-------------------------------
  expdesign = (dt_meta[lst_library_forPlot,.(library, group = sample)] %>% setkey(.,library))
  
  #get pca list function
  get_pca_list = function(exprMat_forSignoise, expdesign){
    pca_prcomp = prcomp(t(exprMat_forSignoise), scale = F)
    pcs = predict(pca_prcomp) %>% data.frame() 
    pcs$library = row.names(pcs)
    pcs_add_meta = merge(pcs, dt_meta, by = 'library')
    PC1 = round(summary(pca_prcomp)$importance[2,1], digits = 4)
    PC2 = round(summary(pca_prcomp)$importance[2,2], digits = 4)
    PC3 = round(summary(pca_prcomp)$importance[2,3], digits = 4)
    SNR = round(CalcSigNoiseRatio(pca_prcomp, expDesign = expdesign), digits = 2)
    gene_num = dim(exprMat_forSignoise)[1]
    pca_list = list(pcs_add_meta, SNR, gene_num, PC1, PC2, PC3)
    return(pca_list)
  }
  print('signoise')
  print(Sys.time())
    pca_list = get_pca_list(exprMat_forSignoise, expdesign)
  print(Sys.time())
  return(pca_list)
  }

#pca_list_l <- pre_pca_snr(df_exp_annot_sub, dt_meta,
#                        zscore_per_group = 'NO', num_dataset = 2)

