#sever

shinyServer(function(input, output, session){
  ##########################summary#############################
  output$summary <- renderUI({
    tagList(
      tags$h3("Overview of study design"),
      tags$p(
        style="text-align:justify;",
        "Quartet RNA reference materials were derived from immortalized Epstein-Barr Virus (EBV) infected lymphoblastoid cell lines from a quartet family including monozygotic twin daughters (D5 and D6), and their father (F7) and mother (M8). Multi-batch RNAseq datasets were generated from independent laboratories using different library preparation protocols and sequencing platforms. Intra-batch proficiency and cross-batch reproducibility were then estimated. Based on multi-batch RNAseq data, we constructed reference datasets at relative level, and developed quality metrics."),
      tags$div(
        style="display: flex; justify-content: center; align-items: center;",
        tags$img(src='./assets/design.png', style="width: 80%;")
      )
    )
  })

  ##########################gene expression#############################
  observe({
    if(input$id_type == 'Ensemble ID'){
      gene_list = df_gene[['ensemble_id']] %>% unique
      updateSelectizeInput(session,
                           inputId = "gene_list",
                           choices = gene_list,
                           selected = 'ENSG00000141510',
                           server = TRUE)
     } else {
       gene_list = df_gene[['symbol']] %>% unique()
       updateSelectizeInput(session,
                            inputId = "gene_list",
                            choices = gene_list,
                            selected = 'CFH',
                            server = TRUE)
       }
    })

  ### expression boxplot figure
  output$plot_gene_exp <- renderPlotly({
    ### parameter
    xaxis_ref <- list(
      title = 'group',
      automargin = TRUE,
      tickfont =  list(size = 14),
      titlefont = list(size = 14)
      )

    yaxis_ref <- list(
      title = 'logfpkm',
      automargin = TRUE,
      tickfont =  list(size = 14),
      titlefont = list(size = 14)
    )
    legend_ref <- list(orientation = 'right', font  = list(size = 14))
    margin_ref = list(l = 20, t = 50, r = 20,  b = 20)

    ### plot
    if (input$id_type == 'Ensemble ID') {
      dt_exp_g <- dt_exp[ensemble_id == input$gene_list]
    } else {
      setindex(dt_exp, symbol)
      dt_exp_g <- dt_exp[symbol == input$gene_list]
    }
    
    df_exp_l <- melt(dt_exp_g, id = 1:2)
    df_exp_plot <- df_exp_l[, c('variable', 'value'), with = F]
    colnames(df_exp_plot) <- c('library', 'logfpkm')
    df_exp_annot_plot <- dt_meta_used[df_exp_plot, on = 'library'][, c('logfpkm', 'sample', input$group_a_gene), with = F]
    colnames(df_exp_annot_plot) <- c('logfpkm', 'sample', 'group')
   
     if(input$group_a_gene == 'sample'){
      plot_ly(df_exp_annot_plot,
              x = ~group, 
              y = ~logfpkm, 
              type = 'violin', 
              box = list(visible = T), 
              meanline = list(visible = F),
              color = ~sample,
              colors = pal) %>% 
        layout(title = input$gene_list,
               xaxis = xaxis_ref, 
               yaxis = yaxis_ref,
               legend = legend_ref, 
               margin = margin_ref)
     } else {
       plot_ly(df_exp_annot_plot,
               x = ~group,
               y = ~logfpkm,
               color = ~sample,
               colors = pal,
               type = "box") %>%
         layout(boxmode = "group", 
                title = input$gene_list,
                xaxis = xaxis_ref,
                yaxis = yaxis_ref,
                legend = legend_ref,
                margin = margin_ref)
     }})

  output$plot_gene_exp_legend <- renderUI({
    HTML(
      paste(
        "<br/>In the violin plot, each distinct 'violin' shape represents a separate group. The breadth of each 'violin' at a given point is indicative of the density of samples with an expression level corresponding to that point, effectively illustrating the distribution. Wider sections of the 'violin' are representative of a higher sample density at certain expression levels, suggesting a larger number of samples in the group have those specific expression levels. In contrast, the narrower sections of the 'violin' symbolize fewer samples with that particular level of gene expression.",
        "<br/><br/>The following parameters are used in the plot:",
        "<b>Gene ID Type</b> is the type of Gene ID used. It can be 'ensemblId' or 'GeneSymbol'.",
        "<b>Gene ID</b> is the specific identification number or symbol for the gene of interest.",
        "<b>Group</b> is the category that the sample belongs to. The group can be 'sample', 'batch', 'libraryPrep', 'kit', or 'lab'. These groups represent different stages or aspects of the experimental process:",
        "<br/>sample: This refers to a particular biological specimen (D5, D6, F7, or M8) used for the experiment.",
        "<br/>batch: This refers to a set of samples processed at the same time, under the same conditions.",
        "<br/>libraryPrep: This indicates the method (PolyA or RiboZero) used to prepare the DNA or RNA library for sequencing.",
        "<br/>kit: This refers to the specific commercial kit (BGI, ILM, KAPA, PE, or VAZ) used for library preparation or other experimental steps.",
        "<br/>lab: This indicates the specific laboratory (ABC, ARD, BGI, BRG, FDU, VAZ, WEH or WUX) where the experiment was conducted."
      )
    )
  })
    
  #######################intra batch qc metrics#########################
  df_intra_met <-  reactive(
    fread(file = paste('./data/rnaseq/qc_type/intra_batch/', input$qc_type_intra, '.txt', sep = ''))
  )
  
  
  observe({
    if (input$qc_type_intra == 'Detected_Gene') {
      updateSelectizeInput(session, "intra_metric",
                           choices = "None",
                           selected = "None")
    } else if (input$qc_type_intra == 'Tree_Metrics') {
      metric_list = c('SNR', 'LIR', 'LRR2')
      updateSelectizeInput(session, "intra_metric",
                           choices = metric_list,
                           selected = metric_list[1])
      }
    })

  ### intra batch figures
  output$plot_intra_batch <- renderPlotly({
    ### parameter
    xaxis_ref <- list(
      title = input$group_a_intra,
      automargin = TRUE,
      tickfont =  list(size = 14),
      titlefont = list(size = 14)
      )

    yaxis_ref <- list(
      title = input$intra_metric,
      automargin = TRUE,
      tickfont =  list(size = 14),
      titlefont = list(size = 14)
    )
    legend_ref <- list(orientation = 'right', font  = list(size = 14))
    margin_ref = list(l = 20, t = 50, r = 20,  b = 20)

    ###plot
    if(input$qc_type_intra == 'Tree_Metrics'){
      df_intra = df_intra_met()[, c(input$intra_metric, input$group_a_intra), with = F] %>%
        setnames(input$group_a_intra, "group") %>%
        setnames(input$intra_metric, "value")

      ### reorder x name
      xaxis_ind <- list(
        categoryorder = "array",
        categoryarray =  df_intra[order(value)][['group']],
        title = input$group_a_intra,
        automargin = TRUE,
        tickfont =  list(size = 14),
        titlefont = list(size = 14))

      plot_ly(df_intra, color = I("gray")) %>%
        add_segments(x = ~group, xend = ~group, yend = 0, y = ~value, showlegend = FALSE) %>%
        add_markers(x = ~group, y = ~value, name = input$intra_metric, color = I("red"))  %>%
        layout(title = input$intra_metric,
               xaxis = xaxis_ind,
               yaxis = yaxis_ref,
               legend = legend_ref,
               margin = margin_ref
               )

    } else if (input$qc_type_intra == 'CV'){
      df_intra_cnv = df_intra_met()[, c("CV", "sample", input$group_a_intra), with = F] %>% setnames(input$group_a_intra, "group")
      df_intra_cnv$CV <- as.numeric(df_intra_cnv$CV)
      plot_ly(df_intra_cnv,
              x = ~group,
              y = ~CV,
              color = ~sample,
              colors = pal,
              type = "box") %>%
        layout(boxmode = "group",
               xaxis = xaxis_ref,
               yaxis = yaxis_ref,
               legend = legend_ref,
               margin = margin_ref)

    } else if (input$qc_type_intra == 'Detected_Gene'){
      df_intra_dtg = df_intra_met()[, c("n_rep", "sample", "N", input$group_a_intra), with = F] %>% setnames(input$group_a_intra, "group")
      batch_list = unique(df_intra_met()[[input$group_a_intra]])
      fig_list = lapply(batch_list, make_subbar_plotly, df_intra_dtg)
      subplot(fig_list,
              nrows = 3,
              shareX = FALSE,
              shareY = TRUE,
              titleX = TRUE) 

    }
  })

  output$plot_intra_batch_legend <- renderUI({
    HTML(
      paste(
        "<br/><br/>The following parameters are used in the plot:",
        "<br/><b>Intra QC Type</b> can be 'Detected_Gene' or 'Tree_Metrics'. 'Detected_Gene' refers to the quantity of genes detected in the experiment and 'Tree_Metrics' refers to the metric of tree structure (like hierarchical clustering or phylogenetic tree metrics) used in the study.",
        "<b>Intra Metric</b> outlines the particular metrics within each QC Type used for quality control assessment. This could be SNR (Signal to Noise Ratio), LIR (Log Intensity Ratio), or LRR2 (Log R Ratio squared).",
        "<br/><b>Group</b> is the category that the sample belongs to. The group can be 'sample', 'batch', 'libraryPrep', 'kit', or 'lab'. These groups represent different stages or aspects of the experimental process:",
        "<br/>batch: This refers to a set of samples processed at the same time, under the same conditions.",
        "<br/>libraryPrep: This indicates the method (P - PolyA or R - RiboZero) used to prepare the DNA or RNA library for sequencing.",
        "<br/>kit: This refers to the specific commercial kit (BGI, ILM, KAPA, PE, or VAZ) used for library preparation or other experimental steps.",
        "<br/>lab: This indicates the specific laboratory (ABC, ARD, BGI, BRG, FDU, VAZ, WEH or WUX) where the experiment was conducted.",
        "<br/>quality: This indicates the quality of the library preparation (Low or High)."
      )
    )
  })

  #########################cross batch qc metrics########################
  df_qc_met <-  reactive(
    if(input$figure_type_corss == 'boxplot'){
      fread(file = paste('./data/rnaseq/qc_type/cross_batch/box/', input$qc_type_cross, '.txt', sep = ''))
    } else if (input$figure_type_corss == 'heatmap'){
      fread(file = paste('./data/rnaseq/qc_type/cross_batch/heat/', input$qc_type_cross, '.txt', sep = ''))
    }
  )

  observe({
    if(input$qc_type_cross == 'Relative Correlation'){
      updateSelectizeInput(session,
                           inputId = "heatmap_sample",
                           choices = c("D6/D5", "F7/D5", "M8/D5", "F7/D6", "M8/D6", "M8/F7"),
                           selected = "D6/D5",
                           server = TRUE)
    } else if (input$qc_type_cross == 'Absolute Correlation') {
      updateSelectizeInput(session,
                           inputId = "heatmap_sample",
                           choices = c("D5", "D6", "F7", "M8"),
                           selected = "D5",
                           server = TRUE)
    } else {
      updateSelectizeInput(session,
                           inputId = "heatmap_sample",
                           choices = "None",
                           selected = "None")
    }
  })
  
  ###figure
  output$plot_cross_batch <- renderPlotly({

    xaxis_ref <- list(
      title = 'group',
      automargin = TRUE,
      tickfont =  list(size = 14),
      titlefont = list(size = 14),
      categoryorder = "array",
      categoryarray = c('Intra-batch', 'Cross-batch', 'Cross-lab', 'Cross-platform', 'Cross-protocol')
    )

    yaxis_ref <- list(
      title = input$qc_type_cross,
      automargin = TRUE,
      tickfont =  list(size = 14),
      titlefont = list(size = 14)
    )
    legend_ref <- list(orientation = 'right', font  = list(size = 14))
    margin_ref = list(l = 20, t = 50, r = 20,  b = 20)

    ###plot


    if(input$figure_type_corss == 'boxplot'){
      if(input$qc_type_cross == 'SNR'){
        setnames(df_qc_met(), "SNR", "value")
        
      } else if(input$qc_type_cross == 'Relative Correlation'){
        setnames(df_qc_met(), "Pearson", "value")
        
      } else if(input$qc_type_cross == 'Absolute Correlation'){
        setnames(df_qc_met(), "Cor", "value")
        
      }

      plot_ly(df_qc_met(),
              x = ~batch_type,
              y = ~value,
              color = ~batch_type,
              type = "box") %>%
        layout(
          boxmode = "group",
          title = input$qc_type_cross,
          xaxis = xaxis_ref,
          yaxis = yaxis_ref,
          legend = legend_ref,
          margin = margin_ref
          )

    } else if(input$figure_type_corss == 'heatmap'){
     
       if(input$qc_type_cross == 'SNR'){
        setnames(df_qc_met(), "SNR", "value")
        df_qc_met_s <- df_qc_met()[, c('batchA', 'batchB', 'value')]
        
      } else if(input$qc_type_cross == 'Relative Correlation'){
        df_qc_met_u <- df_qc_met()[compare_group == input$heatmap_sample]
        setnames(df_qc_met_u, "Cor", "value")
        df_qc_met_s <- unique(df_qc_met_u[, c('batchA', 'batchB', 'value')])
        
      } else if(input$qc_type_cross == 'Absolute Correlation'){
        df_qc_met_u <- df_qc_met()[sample == input$heatmap_sample]
        setnames(df_qc_met_u, "Cor", "value")
        df_qc_met_s <- unique(df_qc_met_u[, c('library_A', 'library_B', 'value')])
        colnames(df_qc_met_s) <- c('batchA', 'batchB', 'value')
        
      }
      
      df_qc_met_w <- dcast(df_qc_met_s, batchA~batchB, value.var = 'value')
      cor_het <- as.matrix(df_qc_met_w[, -1])
      
      dd_col <- as.dendrogram(hclust(dist(cor_het)))
      dd_row <- as.dendrogram(hclust(dist(t(cor_het))))
      
      col_ord <- order.dendrogram(dd_col)
      row_ord <- order.dendrogram(dd_row)
      
      x_order <- colnames(cor_het)[col_ord]
      y_order <- df_qc_met_w[[1]][row_ord]
      
      xform <- list(
        title = 'group',
        automargin = TRUE,
        tickfont =  list(size = 14),
        titlefont = list(size = 14),
        categoryorder = "array",
        categoryarray = x_order
      )
      
      yform <- list(
        title = input$qc_type_cross,
        automargin = TRUE,
        tickfont =  list(size = 14),
        titlefont = list(size = 14),
        categoryorder = "array",
        categoryarray = y_order
      )
      
      print(df_qc_met_s)
      print(unique(df_qc_met_s[['batchA']]))
      print(unique(df_qc_met_s[['batchB']]))
  plot_ly(df_qc_met_s,
              x = ~batchA,
              y = ~batchB,
              z = ~value,
              type = "heatmap",
              colors = colorRamp(c("white", "red"))) %>%
        layout(
          xaxis = xform,
          yaxis = yform,
          legend = legend_ref,
          margin = margin_ref
          )
    }
  })

  output$plot_cross_batch_legend <- renderUI({
    HTML(
      paste(
        "<br/><br/>The following parameters are used in the plot:",
        "<br/><b>Cross Batch QC Type</b> can be 'SNR', 'Relative correlation' and 'Absolute correlation'. 'SNR' was established to characterize the ability of a platform, a lab, or a batch to distinguish the intrinsic differences among distinct biological sample groups ('signal') from variations in technical replicates of the same sample group ('noise'). 'Relative correlation' is the Pearson correlation coefficient of ratio-based expression levels of datasets for a given pair of groups. 'Absolute correlation' is the Pearson correlation coefficient of absolute expression levels of datasets for a given sample.",
        "<br/><b>Figure Type</b> can be 'boxplot' and 'heatmap'. 'Boxplot' shows the overall distribution of  the value of SNR, relative correlation and absolute correlation. 'heatmap' shows the value of SNR, relative correlation and absolute correlation between different datasets.",
        "<br/><b>Sample</b> item will change according to your choice of 'Relative correlation' and 'Absolute correlation'. 'Relative correlation' has six groups: D6/D5, F7/D5, M8/D5, F7/D6, M8/D6, M8/F7. 'Absolute correlation' has four samples: D5, D6, F7, M8."
      )
    )
  })

  #########################pca########################
  observe({
    group_pca  = unique(df_exp_annot[[input$group_pca]])
    updatePickerInput(session, "group_pca_unit",
                         choices = group_pca,
                         selected = group_pca[1])
  })

  observe({
    if(input$pca_scale == 'relative'){
      updateSelectizeInput(session, "ref_sample",
                           choices = c('D5', 'D6', 'F7', 'M8'),
                           selected = 'D5')
    } else {
      updateSelectizeInput(session, "ref_sample",
                           choices = 'None',
                           selected = 'None')
    }
  })
  
  df_exp_annot_sub <-  reactive(
    df_exp_annot[input$group_pca_unit, on = input$group_pca][, c(input$group_pca, "logfpkm", "library", "ensemble_id"), with = F] %>%
      setnames(input$group_pca, "batch")
  )

  ###figure oca
  output$plot_pca <- renderPlotly({
    df_pca <- df_exp_annot_sub()[.(pca_gene_list$nam), on = .(ensemble_id)]
      pca_list = pre_pca_snr(df_exp_annot_sub = df_pca, dt_meta = dt_meta,
                             ref_sample = input$ref_sample, scale_per_group = input$pca_scale)
      mat_pca = pca_list[[1]]
      snr = pca_list[[2]]
      gene_num = pca_list[[3]]
      pro1 = pca_list[[4]]
      pro2 = pca_list[[5]]
      pro3 = pca_list[[6]]

    ###parameter
      xaxis_ref <- list(
        title = paste("PC1 (", pro1*100,  "%)", sep = ""),
        automargin = TRUE,
        zeroline = FALSE, 
        tickfont =  list(size = 14),
        titlefont = list(size = 14)
      )

      yaxis_ref <- list(
        title = paste("PC2 (", pro2*100,  "%)", sep = ""),
        automargin = TRUE,
        zeroline = FALSE, 
        tickfont =  list(size = 14),
        titlefont = list(size = 14)
      )

      zaxis_ref <- list(
      title = paste("PC3 (", pro3*100,  "%)", sep = ""),
      automargin = TRUE,
      tickfont =  list(size = 14),
      titlefont = list(size = 14)
      )

      legend_ref <- list(orientation = 'right', font  = list(size = 14))
      margin_ref = list(l = 20, t = 50, r = 20,  b = 20)
      pal <- c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745")
      title_ref <-  paste('SNR = ', snr, '\n', '(N = ', gene_num, ")", sep = "")

    ###plot
    if(input$d2_d3_switch == TRUE){
      plot_ly(data = mat_pca,
              x = ~PC1,
              y = ~PC2,
              z = ~PC3,
              color = ~sample,
              colors = pal,
              marker = list(size = 10),
              alpha = 1,
              symbol = ~symbol.batch,
              symbols = c('circle','triangle-up','square'),
              text = ~library) %>%
        layout(
          title = title_ref,
          xaxis = xaxis_ref,
          yaxis = yaxis_ref,
          legend = legend_ref,
          margin = margin_ref)
    } else {
      plot_ly(data = mat_pca,
              x = ~PC1,
              y = ~PC2,
              color = ~sample,
              colors = pal,
              marker = list(size = 10),
              alpha = 1,
              symbol = ~symbol.batch,
              symbols = c('circle','triangle-up','square'),
              text = ~library) %>%
        layout(
          title = title_ref,
          xaxis = xaxis_ref,
          yaxis = yaxis_ref,
          legend = legend_ref,
          margin = margin_ref)
    }
  })

  output$plot_pca_legend <- renderUI({
    HTML(
      paste(
        "<br/><br/>The following parameters are used in the plot:",
        "<br/><b>Group</b> can be 'batch', 'libraryPrep', 'lab', 'sample', 'kit' and 'sequcingSite'. batch: This refers to a set of samples processed at the same time, under the same conditions.",
        "<br/>libraryPrep: This indicates the method (PolyA or RiboZero) used to prepare the DNA or RNA library for sequencing.",
        "<br/>lab: This indicates the specific laboratory (ABC, ARD, BGI, BRG, FDU, VAZ, WEH or WUX) where the experiment was conducted.",
        "<br/>sample: This refers to a particular biological specimen (D5, D6, F7, or M8) used for the experiment.",
        "<br/>kit: This refers to the specific commercial kit (BGI, ILM, KAPA, PE, or VAZ) used for library preparation or other experimental steps.",
        "<br/>sequcingSite: This refers to the Sequencing Platforms (BGI and ILM).",
        "<br/><b>PCA Scale</b> can be 'no-scale', 'z-score' and 'relative'. 'no-scale' represents the absolute expression level matrix. 'z-score' represents the z-score corrected expression matrix. 'relative' represents the relative expression level matrix.",
        "<br/><b>Element of group</b> item will change according to your choice of Group: 'batch', 'libraryPrep', 'lab', 'sample', 'kit' and 'sequcingSite'. For example, when Group selects 'batch', you can select one or more datasets from all batches.",
        "<br/><b>Reference Sample</b> is only available when PCA Scale selects 'relative', and any sample from D5, D6, D7 and M8 can be selected as the base for calculating the relative expression matrix."
      )
    )
  })

  ##############################reference dataset##########################
  print('RD')

  observe({
    if(input$data_type %in% c('Consensus of DEGs', 'Volcano for Ref DEGs')){
      group_list = unique(df_dataset()[['compare']])
      updateSelectizeInput(session, "rd_compare_group",
                           choices = group_list,
                           selected = tail(group_list, 1))

    } else if ( input$data_type == 'Distribution of Gene Type'){
      group_list = unique(df_dataset()[['sample']])
      updateSelectizeInput(session, "rd_compare_group",
                           choices = group_list,
                           selected = tail(group_list, 1))

    } else if ( input$data_type == 'Performance Summary'){
      group_list = c('summary', 'correlation', 'library protocols')
      updateSelectizeInput(session, "rd_compare_group",
                           choices = group_list,
                           selected = tail(group_list, 1))

    } else {
      updateSelectizeInput(session, "rd_compare_group",
                           choices = "None",
                           selected = "None")
    }
  })

  df_dataset <-  reactive(
    fread(file = paste('./data/rnaseq/qc_type/reference_dataset/', input$data_type, '.txt', sep = ''))
  )

  df_compare <- reactive(
    df_dataset()[compare == input$rd_compare_group]
  )

  output$plot_RD <- renderPlotly({
    xaxis_ref <- list(
      title = 'group',
      automargin = TRUE,
      tickfont =  list(size = 14),
      titlefont = list(size = 14)
    )

    yaxis_ref <- list(
      title = input$gene_list,
      automargin = TRUE,
      tickfont =  list(size = 14),
      titlefont = list(size = 14)
    )
    legend_ref <- list(orientation = 'right', font  = list(size = 14))
    margin_ref = list(l = 20, t = 50, r = 20,  b = 20)

    #sample color
    pal <- c('#4CC3D9', '#7BC8A4', '#FFC65D', '#F16745')
    pal <- setNames(pal, c("D5", "D6", "F7", "M8"))

    ###plot
    if (input$data_type == 'Consensus of Detected Genes'){
      con_det_gene <- df_dataset()
      plot_ly(con_det_gene,
              x = ~Freq,
              y = ~Var1,
              type = 'funnel',
              name = ~sample,
              color = ~sample,
              colors = pal) %>%
        layout(
          font = list(size = 14),
          margin = list(t = 100),
          yaxis = list(title = 'Detection Agreement', automargin = TRUE, size = 14),
          xaxis = list(title = 'Detected Genes Value', automargin = TRUE, size = 14))

    } else if (input$data_type == 'Distribution of Gene Type'){
      dis_gene_type <-  df_dataset()[sample == input$rd_compare_group]
      plot_ly(dis_gene_type,
              x = ~genereps,
              y = ~ratio,
              color = ~gene_type,
              name = ~gene_type,
              type = 'bar') %>%
        layout(barmode = 'stack',
               font = list(size = 14),
               margin = list(t = 100),
               title = input$rd_compare_group,
               yaxis = list(title = 'Percentage', automargin = TRUE, size = 14),
               xaxis = list(title = 'Detection Agreement', automargin = TRUE, size = 14))

    } else if (input$data_type == 'Consensus of DEGs'){
      plot_ly(df_compare(), x = ~value, y = ~Var1, type = 'funnel', color = ~groupcol,
                      colors = c('blue', 'red', 'grey')) %>%
        layout(title = input$rd_compare_group,
               yaxis = list(title = 'Agreement', automargin = TRUE, size = 14),
               xaxis = list(title = 'Number of DEGs', automargin = TRUE, size = 14),
               font = list(size = 14),
               margin = list(t = 100)
               )

    } else if (input$data_type == 'Volcano for Ref DEGs'){
      plot_ly(df_compare(),
              x = ~meanlogFC,
              y = ~lg_medianP,
              color = ~Final,
              colors = c('blue', 'grey', 'red'),
              marker = list(size = 5),
              text = ~Gene_Name) %>%
        layout(title = input$rd_compare_group,
               font = list(size = 14),
               margin = list(t = 100),
               yaxis = list(title = '-log10(median p)', automargin = TRUE, size = 14),
               xaxis = list(title = paste('mean log (', input$rd_compare_group, ")", sep = '')), automargin = TRUE, size = 14)

    } else if (input$data_type == 'Performance Evaluation at RE'){
      per_summary <- df_dataset()
      fig1 <- plot_ly(data = per_summary,
                      x = ~corr,
                      y = ~MCC,
                      marker = list(size = 10),
                      color = ~batch,
                      showlegend = FALSE,
                      symbol = ~libprep_qual,
                      symbols = c('square','triangle-up', 'circle', 'diamond'),
                      text = ~batch) %>%
        layout(font = list(size = 14))
      fig2 <- plot_ly(data = per_summary,
                      x = ~libprep_qual,
                      y = ~MCC,
                      type = "box",
                      color = ~libprep_qual,
                      showlegend = FALSE) %>%
        layout(font = list(size = 14))
      fig3 <- plot_ly(data = per_summary,
                      x = ~corr,
                      y = ~libprep_qual,
                      type = "box",
                      color = ~libprep_qual) %>%
        layout(font = list(size = 14))
      s1 <- subplot(fig3,
                    fig1,
                    shareX = TRUE,
                    nrows = 2,
                    heights = c(0.2, 0.8))
      subplot(s1,
              fig2,
              widths = c(0.8, 0.2),
              titleX = TRUE,
              titleY = TRUE)

    } else if (input$data_type == 'Performance Summary') {
      if(input$rd_compare_group == 'summary'){
        xform <- list(categoryorder = "array",
                      categoryarray =  unique(df_dataset()[['metric']]))
        yform <- list(categoryorder = "array",
                      categoryarray =  df_dataset()[metric == 'total score'][order(value)][['batch']])

        plot_ly(df_dataset(),
                x = ~metric,
                y = ~batch,
                z = ~value,
                type = 'scatter',
                color = ~log10(value),
                size = ~log10(value),
                colors = c('darkviolet', 'gold') ,
                text = ~value,
                mode = 'markers',
                marker = list(symbol = 'circle', sizemode = 'diameter',
                              line = list(width = 2, color = '#FFFFFF'),
                              sizeref = 3.5)) %>%
          layout(xaxis = xform,
                 yaxis = yform,
                 font = list(size = 10))

      } else if(input$rd_compare_group == 'library protocols'){
        plot_ly(df_dataset(),
                x = ~metric,
                y = ~value,
                color = ~libprep,
                type = "box") %>%
          layout(boxmode = "group",
                 font = list(size = 14))

      } else if(input$rd_compare_group == 'correlation'){
        df_cor_scatter <- dcast(df_dataset(), batch~metric)[, !"batch"]
        df_cor_map <- map2(df_cor_scatter, names(df_cor_scatter), ~list(values = .x, label = .y))
        plot_ly(type = "splom",
                dimensions = setNames(df_cor_map, NULL),
                showupperhalf = FALSE, diagonal = list(visible = FALSE)) %>%
          layout(font = list(size = 10),
                 margin = list(l = 100, b = 50))

      }
      }
  })

  output$plot_RD_legend <- renderUI({
    HTML(
      paste(
        "<br/><br/>The following parameters are used in the plot:",
        "<br/><b>Data Type</b> can be 'Consensus of Detected Genes', 'Distribution of Gene Type', 'Consensus of DEGs', 'Volcano for Ref DEGs', 'Performance Evaluation of RE' and 'Performance Summary'",
        "<br/>Consensus of Detected Genes: shows the distribution of the number of genes in D5, D6, F7 and M8 samples that can be detected simultaneously in 0-15 datasets.",
        "<br/>Distribution of Gene Type: shows the distribution of the proportions of different gene types in D5, D6, F7 or M8 samples from different datasets.",
        "<br/>Consensus of DEGs: shows the distribution of the DEGs in six groups (D6/D5, F7/D5, M8/D5, F7/D6, M8/D6, M8/F7) that can be detected simultaneously in 0-15 datasets.",
        "<br/>Volcano for Ref DEGs: is a volcano plot showing the distribution of fold change and p-value in six groups (D6/D5, F7/D5, M8/D5, F7/D6, M8/D6, M8/F7)",
        "<br/>Performance Evaluation of RE: shows the values of relative correlation and MCC for each sample in each dataset through a scatter plot.",
        "<br/>Performance Summary: is a summary of the SNR, relative correlation, absolute correlation, and relative correlation for the reference datasets. total score is calculated to measure the overall quality of all datasets.",
        "<br/><b>Group</b> item will change according to your choice of <b>Data Type</b>: Distribution of Geene Type', 'Consenesus of DEGs', 'Volcano for Ref DEGs' and 'Performance Summary'. For example, when Group selects 'Distribution of Geene Type', you can select any sample, such as D5, D6, F7 or M8."
      )
    )
  })

  session$onSessionEnded(stopApp)
})
