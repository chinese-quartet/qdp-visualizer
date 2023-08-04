#library(plotly)
#library(heatmaply)
### SNV
heatmaply_SNV <- function(snv, Margin=100,Plot_method="ggplot", Key.title="Jaccard Index", Main="SNV"){
  # pre
  snv_lcl5 <- snv[grep('LCL5',rownames(snv)),grep('LCL5',colnames(snv))]
  snv_lcl5_sorted <- snv_lcl5[ order(colnames(snv_lcl5)), order(colnames(snv_lcl5))]
  lcl5_pd <- data.frame(colnames(snv_lcl5_sorted))
  lcl5_pd[,1] <- as.character(lcl5_pd[,1])
  lcl5_pd$platform <- sapply(strsplit(lcl5_pd[,1],"_"),function(x){x[4]})
  lcl5_pd$site <- sapply(strsplit(lcl5_pd[,1],"_"),function(x){x[5]})
  lcl5_pd$library <- c('PCR','PCR','PCR','PCR-free','PCR-free','PCR-free',
                       'PCR-free','PCR-free','PCR-free','PCR-free','PCR-free','PCR-free',
                       'PCR-free','PCR-free','PCR-free','PCR','PCR','PCR',
                       'PCR','PCR','PCR','PCR','PCR','PCR','PCR','PCR','PCR')
  colnames(lcl5_pd)[1] <- 'sample'
  data.table::setnames(lcl5_pd,
                       c("platform","site","library"),
                       c("Sequence platform","Sequence site","Library")
  )
  
  #
  pp <- heatmaply(
    snv_lcl5_sorted,
    scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
      high = "#E41A1C", 
      low = "#377EB8", 
      midpoint = 0.75, 
      limits = c(0.5, 1)
    ),
    #colors = c("#377EB8", "white", "#E41A1C"),
    #seriate = "mean",
    dist_method="euclidean",
    hclust_method="complete",
    dendrogram="both",
    row_dend_left = F,
    plot_method = Plot_method,
    col_side_colors = lcl5_pd[,2:4],
    showticklabels = c(F, F),
    show_dendrogram = c(T, F),
    margins = c(Margin,Margin,Margin,Margin),
    key.title= Key.title,
    main =Main,
    col_side_palette=  c("Nova" = "#7570B3", "SEQ2000" = "#E7298A","T7"="#1B9E77","XTen" = "#D95F02",
                         "ARD"="#4DAF4A","BGI"="#FFFF33",'BRG'='#377EB8','NVG'='#984EA3','WGE'='#E41A1C','WUX'='#FF7F00',
                         'PCR' = '#8e0c24','PCR-free' = '#587aa5'
                         
    )
  )
  pp$x$layout$annotations[[1]]$text <- ''
  pp
}


### SV
heatmaply_SV <- function(ins, Margin=100, Plot_method="ggplot", Key.title="Jaccard Index", Main="INS"){
  # pre
  anno <- data.frame(rownames(ins))
  anno[,1] <- as.character(anno[,1])
  anno$tech <- sapply(strsplit(anno[,1],"_"),function(x){x[1]})
  anno$mapper <- sapply(strsplit(anno[,1],"_"),function(x){x[2]})
  anno$caller <- sapply(strsplit(anno[,1],"_"),function(x){x[3]})
  data.table::setnames(anno,
                       c("tech","mapper","caller"),
                       c("Sequence platform","Mapper","Caller")
  )
  #
  pp <- heatmaply(
    ins,
    scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
      high = "#E41A1C", 
      low = "#377EB8", 
      midpoint = 0.5, 
      limits = c(0, 1)
    ),
    #colors = c("#377EB8", "white", "#E41A1C"),
    #seriate = "mean",
    dist_method="euclidean",
    hclust_method="complete",
    dendrogram="both",
    row_dend_left = F,
    plot_method = Plot_method,
    col_side_colors = anno[,2:4],
    showticklabels = c(F, F),
    show_dendrogram = c(T, F),
    margins = c(Margin,Margin,Margin,Margin),
    key.title= Key.title,
    main =Main,
    col_side_palette=  c("ont" = "#f8766d", "pacbio" = "#7cae00","pacbio2"="#00bfc4",
                         "minimap2"="#E41A1C","ngmlr"="#377EB8",'pbmm2'='#4DAF4A',
                         'cutesv' = '#A6CEE3','nanosv' = '#1F78B4','pbsv'='#4DAF4A','sniffles'="#984EA3",'svim'="#FF7F00"
                         
    )
  )
  pp$x$layout$annotations[[1]]$text <- ''
  pp
}

#heatmaply_SNV(snv,Main = "SNV")
#heatmaply_SNV(indel,Main = "INDEL")

# heatmaply_SV(ins, Main = "Insertion")
# heatmaply_SV(del, Main = "Deletion")
# heatmaply_SV(dup, Main = "Duplication") 
# heatmaply_SV(inv ,Main = "Inversion") 
# heatmaply_SV(bnd, Main = "Breakend")