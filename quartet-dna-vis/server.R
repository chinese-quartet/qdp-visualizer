#sever

shinyServer(function(input, output, session){
  ##########################summary#############################
  output$summary <- renderUI({
    tagList(
      tags$h3("Overview of study design"),
      tags$p(
        style="text-align:justify;",
        "(a) Overview of study design. Briefly, DNA reference materials were constructed from a Quartet family with father (F7), mother (M8) and monozygotic twin daughters (D5 and D6). They were sequenced by four short- and three long-read platforms at seven labs. Small variant and structural variant benchmark calls were integrated from massive sequencing datasets. Performance of a test dataset can be evaluated by comparing with benchmark calls or genetic built-in truth within the Quartet family. (b) Schematic overview of short-read sequencing datasets. Three replicates for each of the four Quartet DNA reference materials were sequenced in nine batches, by both PCR and PCR-free libraries on four sequencing platforms at six labs, resulting in 108 WGS libraries. (c) Schematic overview of long-read sequencing datasets. One replicate for each of the four Quartet DNA reference materials were sequenced per batch by PacBio Sequel, PacBio Sequel II and ONT. Eleven combinations of three mappers and five callers were used to detect structural variants, resulting in 120 variants calling datasets."),
      tags$div(
        style="display: flex; justify-content: center; align-items: center;",
        tags$img(src='./design.png', style="width: 80%;")
      )
    )
  })

  ############################## Mendelian Violations##########################
  ########## add 4 figures
  
  df_mendelian <-  reactive(
    
    if(input$Mendelian_Violations =="Mendelian Violations Ratio"){
      if(input$Mendelian_Violations_type3=="SNV"){
        readRDS('./data/DNAseq/Mendelian_Violations/Mendelian_Violations_ratio/mendelian_ratio_SNV.rds')
      }else if(input$Mendelian_Violations_type3=="Indel"){
        readRDS('./data/DNAseq/Mendelian_Violations/Mendelian_Violations_ratio/mendelian_ratio_indel.rds')
      }else if(input$Mendelian_Violations_type3=="Deletion"){
        readRDS('./data/DNAseq/Mendelian_Violations/Mendelian_Violations_ratio/SV_mendelian_ratio_DEL.rds')
      }else if(input$Mendelian_Violations_type3=="Insertion"){
        readRDS('./data/DNAseq/Mendelian_Violations/Mendelian_Violations_ratio/SV_mendelian_ratio_INS.rds')
      }
    }else if(input$Mendelian_Violations =="Mendelian Violation Quality"){
      if(input$Mendelian_Violations_Quality == "Low Quality"){
        if(input$Mendelian_Violations_type1 == "All Quartet Sample"){
          readRDS('./data/DNAseq/Mendelian_Violations/Low_Quality/conflict.annotated.v.rds')
        }else if(input$Mendelian_Violations_type1 == "Spanning Deletions"){
          readRDS('./data/DNAseq/Mendelian_Violations/Low_Quality/spanning_deletions.anno.v.rds')
        }else if(input$Mendelian_Violations_type1 == "Detectable Mendelian Violations"){
          readRDS('./data/DNAseq/Mendelian_Violations/Low_Quality/mendelian.violation.annotated.v.rds')
        }
      }else if(input$Mendelian_Violations_Quality == "High Quality"){
        if(input$Mendelian_Violations_type2 == "Small Variants"){
          readRDS('./data/DNAseq/Mendelian_Violations/High_Quality/High_quality_Mendelian_Violations_SNV.rds')
        }else if(input$Mendelian_Violations_type2 == "Structural Variants"){
          readRDS('./data/DNAseq/Mendelian_Violations/High_Quality/High_quality_Mendelian_Violations_SV.rds')
        }
      }
      
      
      
    }
  )

  output$plot_mendelian_legend <- renderUI({
    HTML(
      paste(
        "<br/>When the '<b>Mendelian Violation Ratio</b>' section is selected, the other parameter, Mendelian Violations Type, can be used, which includes 'SNV', 'Indel', 'Deletion', or 'Insertion'. In this section, the line plot displays the Mendelian violation ratio of variants with different numbers of times being detected. The bar plot shows the number of variants detected at different numbers of times. The data were calculated based on Quartet_D5.",
        "<br/>When the '<b>Mendelian Violation Low Quality</b>' section is selected, the other parameter, Mendelian Violations Type, can be used, which includes 'All Quartet Sample', 'Spanning Deletions', or 'Detectable Mendelian Violations'. In this section, the bar plots display the distribution of three types of Mendelian violated small variants on genomic repeat regions, large deletions from reference SVs, and 500 bp surrounding regions of SVs (INS and DELs).",
        "<br/>When the '<b>Mendelian Violation High Quality</b>' section is selected, the other parameter, Mendelian Violations Type, can be used, which includes 'Small Variant' or 'Structural Variant'. In this section, the bars represent the percentage of small and structural variants on difficult-to-map regions of each call set."
      )
    )
  })
  
  output$mendelian_plot <- renderPlotly({
    
    if(input$Mendelian_Violations =="Mendelian Violations Ratio"){
      if(input$Mendelian_Violations_type3=="SNV"){
        mendel_ratio_SV_SNV_plotly(df_mendelian(),
                                   titlesize=input$mendelian_title_size,  #15
                                   ytitle=input$mendelian_plot_xy_size,
                                   Title="SNV",
                                   Margin =input$mendelian_plot_margin,
                                   figHeight=input$mendelian_plot_figHeight,#0.5
                                   figMargin=input$mendelian_plot_figMargin, #0.01
                                   Showlegend=as.logical(input$mendelian_plot_Showlegend),# T
                                   variantType="SNV",
                                   Group = "Samples",
                                   Color = "library"
                                   )
      }else if(input$Mendelian_Violations_type3=="Indel"){
        mendel_ratio_SV_SNV_plotly(df_mendelian(),
                                   titlesize=input$mendelian_title_size,  #15
                                   ytitle=input$mendelian_plot_xy_size,
                                   Title="Indel",
                                   Margin =input$mendelian_plot_margin,
                                   figHeight=input$mendelian_plot_figHeight,#0.5
                                   figMargin=input$mendelian_plot_figMargin, #0.01
                                   Showlegend=as.logical(input$mendelian_plot_Showlegend),# T
                                   variantType="SNV",
                                   Group = "Samples",
                                   Color = "library"
        )

      }else if(input$Mendelian_Violations_type3=="Deletion"){
        mendel_ratio_SV_SNV_plotly(df_mendelian(),
                                   titlesize=input$mendelian_title_size,  #15
                                   ytitle=input$mendelian_plot_xy_size,
                                   Title="Deletion",
                                   Margin =input$mendelian_plot_margin,
                                   figHeight=input$mendelian_plot_figHeight,#0.5
                                   figMargin=input$mendelian_plot_figMargin, #0.01
                                   Showlegend=as.logical(input$mendelian_plot_Showlegend),# T
                                   variantType="SV",
                                   Group = "D5_pipelines",
                                   Color = "caller"
        )
        
      }else if(input$Mendelian_Violations_type3=="Insertion"){
        mendel_ratio_SV_SNV_plotly(df_mendelian(),
                                   titlesize=input$mendelian_title_size,  #15
                                   ytitle=input$mendelian_plot_xy_size,
                                   Title="Insertion",
                                   Margin =input$mendelian_plot_margin,
                                   figHeight=input$mendelian_plot_figHeight,#0.5
                                   figMargin=input$mendelian_plot_figMargin, #0.01
                                   Showlegend=as.logical(input$mendelian_plot_Showlegend),# T
                                   variantType="SV",
                                   Group = "D5_pipelines",
                                   Color = "caller"
        )
      }
      
    }else if(input$Mendelian_Violations =="Mendelian Violation Quality"){
      if(input$Mendelian_Violations_Quality == "Low Quality"){
        p <-  ggplot(data=df_mendelian(), aes(x=region, y=percentage)) +
          geom_bar(stat="identity", fill="steelblue")+
          xlim(c("Simple Repeat","Satellite","SINE","LINE","LTR","Low Complexity","Segmental duplication","Centromere",
                 "Tier1 DEL","Tier1 DEL breakpoint","Tier1 INS breakpoint"))+
          xlab("")+
          ylab("Percentage of Variants")+
          geom_text(aes(label=percentage), vjust= -0.5, color="black", size=input$mendelian_plot_annotation_size)+
          theme_classic()+theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11))+ # text on X axis
          theme(axis.text.y = element_text(size = input$mendelian_plot_xy_size))+
          theme(axis.title.y = element_text(size = input$mendelian_plot_xy_size))+
          theme(axis.text.y = element_text(size = input$mendelian_plot_xy_size))
        ggplotly(p) %>% layout(margin=list(l = input$mendelian_plot_margin, t = input$mendelian_plot_margin, 
                                           r = input$mendelian_plot_margin, b = input$mendelian_plot_margin))
        
      }else if(input$Mendelian_Violations_Quality == "High Quality"){
        if(input$Mendelian_Violations_type2 == "Small Variants"){
          p <- ggplot(df_mendelian(), aes(x=variable, y=value)) + 
            geom_bar(stat="identity", color="black", 
                     position=position_dodge()) +
            geom_errorbar(aes(ymin=value-sd, ymax=value+sd, width=.2)) +
            theme_bw()+
            theme(strip.text.x = element_text(size= input$mendelian_plot_xy_size, face="bold"),
                  strip.text.y = element_text(size= input$mendelian_plot_xy_size, face="bold"))+
            xlab('') + 
            #xlim(c('satellite.1','simple_repeat.1','SINE.1','LINE.1','LTR.1','Low_complexity.1','SD.1','Telomere_5M.1','tier1.del.1','tier1.del.breakpoint.1','tier1.ins.breakpoint.1'))+
            ylab('Precentage of variants (%)')+
            theme(axis.text.x = element_text( color="black", 
                                              size= input$mendelian_plot_xy_size, angle=45,hjust=1,vjust=1),
                  axis.text.y = element_text( color="black", 
                                              size= input$mendelian_plot_xy_size))+
            theme(axis.title.y = element_text(size = input$mendelian_plot_xy_size))+
            scale_x_discrete(limits=c("Satellite", "Simple Repeat", "SINE","LINE","LTR",
                                      "Low Complexity","Segmental Duplication","Telomere (5M)",
                                      "Tier1 DEL","Tier1 DEL breakpoint","Tier1 INS breakpoint"
            ))
          ggplotly(p) %>% layout(margin=list(l = input$mendelian_plot_margin, t = input$mendelian_plot_margin, 
                                             r = input$mendelian_plot_margin, b = input$mendelian_plot_margin))
          
        }else if(input$Mendelian_Violations_type2 == "Structural Variants"){
          p <- ggplot(df_mendelian(), aes(x=variable, y=value)) + 
            geom_bar(stat="identity", color="black", 
                     position=position_dodge()) +
            geom_errorbar(aes(ymin=value-sd, ymax=value+sd, width=.2)) +
            theme_bw()+
            theme(strip.text.x = element_text(size= input$mendelian_plot_xy_size, face="bold"),
                  strip.text.y = element_text(size= input$mendelian_plot_xy_size, face="bold"))+
            xlab('') + 
            #xlim(c('satellite.1','simple_repeat.1','trf100bp.1','SINE.1','LINE.1','Low_complexity.1','SD.1','telomere.1'))+
            ylab('Precentage of variants (%)')+
            theme(axis.title.y = element_text(size = input$mendelian_plot_xy_size))+
            theme(axis.text.x = element_text( color="black", 
                                              size= input$mendelian_plot_xy_size, angle=45,hjust=1,vjust=1),
                  axis.text.y = element_text( color="black", 
                                              size= input$mendelian_plot_xy_size))+
            scale_x_discrete(limits=c("Satellite", "Simple Repeat","Tandem Repeat ≥100bp",
                                      "SINE","LINE",
                                      "Low Complexity","Segmental Duplication","Telomere (5M)"
            ))
          ggplotly(p) %>% layout(margin=list(l = input$mendelian_plot_margin, t = input$mendelian_plot_margin, 
                                             r = input$mendelian_plot_margin, b = input$mendelian_plot_margin))
        }
      } 
    }
  })
  
  

  ############################## variant quality ##########################
  df_Variant_Quality <-  reactive(
    #readRDS(file = paste('./data/DNAseq/Variant_Quality/', input$Variant_Quality_type, '.rds', sep = ''))
    if(input$Variant_Quality_type=="Small Variants"){
      read_feather('./data/DNAseq/Variant_Quality/Small_variants_mapping_quality.feather')
    }else{
      read_feather('./data/DNAseq/Variant_Quality/Structural_variants_RE.feather')
    }
    
  )

  output$plot_Variant_Quality_legend <- renderUI({
    HTML(
      paste(
        "<br/>The '<b>Structural Variants</b>' section represents the allele frequency and alternative allele depth distribution of filtered SVs and reference SVs.",
        "<br/>The '<b>Small Variants</b>' section represents the depth, allele frequency, genotype quality, and mapping quality distribution of filtered small variants (green) and reference small variants."
      )
    )
  })

  output$Variant_Quality_plot <- renderPlotly({
    if(input$Variant_Quality_type=="Small Variants"){
      big_df<- df_Variant_Quality()
     # data.table::setnames(big_df,c("dp","af","gq","mq"),c("Depth","Allele Frequency","Genotype Quality","Mapping Quality"))
      colnames(big_df)[1:4] <-  c("Depth","Allele Frequency","Genotype Quality","Mapping Quality")
       p <- small_vaiant_quality(big_df,xfactor=input$Variant_Quality_distribution_type1,
                           axisxysize=input$Variant_Quality_plot_xy_size,
                           Xlab=input$Variant_Quality_distribution_type,
                           Bins=input$Variant_Quality_plot_his_height,
                           Position=input$Variant_Quality_plot_his_type
                           )
      ggplotly(p) %>% layout(margin=list(l = input$Variant_Quality_plot_margin, t = input$Variant_Quality_plot_margin, 
                                         r = input$Variant_Quality_plot_margin, b = input$Variant_Quality_plot_margin))
      
      
    }else{
      big_df<- df_Variant_Quality()
      #data.table::setnames(big_df,c("af","alt"),c("Allele Frequency","RE"))
      colnames(big_df)[c(10,8)] <- c("Allele Frequency","RE")
      p <- SV_quality(big_df,xfactor=input$Variant_Quality_distribution_type2,
                                axisxysize=input$Variant_Quality_plot_xy_size,
                                Xlab=input$Variant_Quality_distribution_type,
                                Bins=input$Variant_Quality_plot_his_height,
                                Position=input$Variant_Quality_plot_his_type
      )
      ggplotly(p) %>% layout(margin=list(l = input$Variant_Quality_plot_margin, t = input$Variant_Quality_plot_margin, 
                                         r = input$Variant_Quality_plot_margin, b = input$Variant_Quality_plot_margin))
    }

  })
  
  #######################Difficult Genomic Regions###############################
  
  df_diff_region <-  reactive(
    if(input$diff_region_variant_type == "Small Variants"){
      readRDS('./data/DNAseq/Difficult_Genomic_Regions/diff_region_snv.rds')
    }else{
      readRDS('./data/DNAseq/Difficult_Genomic_Regions/diff_region_sv.rds')
    }
  )

  output$plot_diff_region_legend <- renderUI({
    HTML(
      paste(
        "<br/>The '<b>Structural Variants</b>' section represents the Mendelian consistent ratio of SNVs and INDELs on genomic repeat regions, 500 bp surrounding regions of Tier 1 reference SVs, and high-confidence bed regions.",
        "<br/>The '<b>Small Variants</b>' section represents the Mendelian consistent ratio of DELs and INS on genomic repeat regions and high-confidence bed regions that represent an additional 1000 bp surrounding regions of Tier 1 and Tier 2 reference SVs.",
        "<br/>The library preparation includes PCR and PCR-free methods. The callers used are cutesv, nanosv, pbsv, sniffle, and svim."
      )
    )
  })

  output$diff_region_plot <- renderPlotly({
    if(input$diff_region_variant_type == "Small Variants"){
      
      p <- ggplot(df_diff_region(), aes(x=variable, y=(1-value) * 100, fill=Library)) +
        geom_boxplot()+
        scale_fill_manual(values=c("#8e0c24", "#587aa5"))+
        facet_grid(. ~ type)+
        xlim(c('satellite','simple_repeat','SINE','LINE','LTR','Low_complexity','SD','Telomere_5M','tier1.del.breakpoint','tier1.ins.breakpoint','high.confidence'))+
        theme_bw()+
        theme(strip.text.x = element_text(size=input$diff_region_plot_subtitle_size, face="bold"),
              strip.text.y = element_text(size=input$diff_region_plot_subtitle_size, face="bold"))+
        xlab('') + 
        ylab('Precentage of Mendelian Consistent')+
        theme(axis.text.x = element_text( color="black", 
                                          size=input$diff_region_plot_xy_size, angle=45,hjust=1,vjust=1),
              axis.text.y = element_text( color="black", 
                                          size=input$diff_region_plot_xy_size))
      ggplotly(p)%>%layout(boxmode = "group")%>% 
        layout(margin=list(l = input$diff_region_plot_margin, t = input$diff_region_plot_margin, 
                                         r = input$diff_region_plot_margin, b = input$diff_region_plot_margin))
      
    }else if(input$diff_region_variant_type == "Structural Variants"){
      p <- ggplot(df_diff_region(), aes(x=Region, y=Mendelian_rate * 100, fill=Caller)) +
        geom_boxplot()+
        facet_grid(. ~ SV)+
        xlim(c('Satellite','Simple_repeat','trf100bp','SINE','LINE','Low_complexity','SD','telomere','confidence'))+
        theme_bw()+
        theme(strip.text.x = element_text(size=input$diff_region_plot_subtitle_size, face="bold"),
              strip.text.y = element_text(size=input$diff_region_plot_subtitle_size, face="bold"))+
        xlab('') + 
        ylab('Precentage of Mendelian consistent (%)')+
        theme(axis.text.x = element_text( color="black", 
                                          size=input$diff_region_plot_xy_size, angle=45,hjust=1,vjust=1),
              axis.text.y = element_text( color="black", 
                                          size=input$diff_region_plot_xy_size))
      ggplotly(p)%>%layout(boxmode = "group")%>% 
        layout(margin=list(l = input$diff_region_plot_margin, t = input$diff_region_plot_margin, 
                           r = input$diff_region_plot_margin, b = input$diff_region_plot_margin))
      
    }
    
  })
  
  #######################Variant Validation###############################
  
  df_variant_validation <-  reactive(
    if(input$variant_validation_type == "Small Variant Validation"){
      readRDS('./data/DNAseq/Variant_Validation/snv_validation.rds')
    }else{
      readRDS('./data/DNAseq/Variant_Validation/sv_validation.rds')
    }
  )

  output$plot_variant_validation_legend <- renderUI({
    HTML(
      paste(
        "<br/>The '<b>SV Validation</b>' section represents the percentage of Mendelian consistent small variants validated by PMRA inside and outside the high-confidence bed regions.",
        "<br/>The '<b>Small Variants Validation</b>' section represents the percentage of Mendelian consistent SVs validated by other orthogonal technologies (Illumina short-read sequencing, Illumina linked-read sequencing, Bionano, and long-read assembly-based SVs calling) inside and outside the high-confidence bed regions.",
        "<br/>The library preparation includes PCR and PCR-free methods. The callers used are cutesv, nanosv, pbsv, sniffle, and svim."
      )
    )
  })

  output$variant_validation_plot <- renderPlotly({
    if(input$variant_validation_type == "Small Variant Validation"){
  
      p <- ggplot(df_variant_validation(), aes(x=type, y=value * 100, fill=Library)) + 
        geom_bar(stat='identity', color="black", 
                 position=position_dodge()) +
        geom_errorbar(aes(ymin=(value-sd)*100, ymax=(value+sd)*100), width=.2,
                      position=position_dodge(.9)) +
        facet_grid(. ~ Region)+
        theme_bw()+
        theme(strip.text.x = element_text(size=input$variant_validation_subtitle_size, face="bold"),
              strip.text.y = element_text(size=input$variant_validation_subtitle_size, face="bold"))+
        xlab('') + 
        ylab('Precentage of variants validated by PMRA (%)')+
        scale_fill_manual(values=c("#8e0c24", "#587aa5"))+
        theme(axis.text.x = element_text( color="black", 
                                          size=input$variant_validation_xy_size, angle=45,hjust=1,vjust=1),
              axis.text.y = element_text( color="black", 
                                          size=input$variant_validation_xy_size))

      ggplotly(p)%>%
        layout(margin=list(l = input$variant_validation_margin, t = input$variant_validation_margin, 
                           r = input$variant_validation_margin, b = input$dvariant_validation_margin))
      
    }else if(input$variant_validation_type == "SV Validation"){
      
      p <- ggplot(df_variant_validation(), aes(x=SVTYPE, y=Validated_rate * 100, fill=Caller)) + 
        geom_bar(stat="identity", color="black", 
                 position=position_dodge()) +
        geom_errorbar(aes(ymin=(Validated_rate-sd)*100, ymax=(Validated_rate+sd)*100), width=.2,
                      position=position_dodge(.9)) +
        facet_grid(. ~ Region)+
        theme_bw()+
        theme(strip.text.x = element_text(size=input$variant_validation_subtitle_size, face="bold"),
              strip.text.y = element_text(size=input$variant_validation_subtitle_size, face="bold"))+
        xlab('') + 
        xlim(c('DEL','INS'))+
        ylab('Variants validated by other technologies (%)')+
        #  scale_fill_manual(values=c("#8e0c24", "#587aa5"))+
        theme(axis.text.x = element_text( color="black", 
                                          size=input$variant_validation_xy_size, angle=45,hjust=1,vjust=1),
              axis.text.y = element_text( color="black", 
                                          size=input$variant_validation_xy_size))
      ggplotly(p)%>%
        layout(margin=list(l = input$variant_validation_margin, t = input$variant_validation_margin, 
                           r = input$variant_validation_margin, b = input$variant_validation_margin))
      
    }
    
  })
  
  
  ####################### SV_Reference###############################

  df_SV_Reference <-  reactive(
    if(input$SV_Reference_type=="Summary"){
      readRDS(file = paste('./data/DNAseq/SV_Reference/Validation_Summary.rds')) 
    }else{
      readRDS(file = paste('./data/DNAseq/SV_Reference/Validation_detail.rds')) 
    }
  )
  
  output$plot_SV_Reference_legend <- renderUI({
    HTML(
      paste(
        "<br/>In the '<b>Summary</b>' section, the bar plot shows the fraction of reference SVs validated by other technologies (Illumina short-read sequencing, Illumina linked-read sequencing, Bionano, and long-read assembly-based SVs calling). Yellow indicates SVs supported by at least one of the other technologies, while grey indicates SVs not supported by any of the other technologies.",
        "<br/>In the '<b>Detail</b>' section, the bar plot shows the percentage of reference SVs validated by short-reads, linked-reads, Bionano, and long-read assembly-based SVs call sets. Each technology is represented by a different color for clarity and distinction."
      )
    )
  })

  output$SV_Reference_plot <- renderPlotly({
    if(input$SV_Reference_type=="Summary"){
      p <- ggplot(data=df_SV_Reference(), aes(x=type, y=value, fill=variable)) +
        geom_bar(stat="identity")+
        scale_fill_manual(values=c('#E69F00','#999999'))+
        theme_classic()+theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = input$SV_Reference_xy_size))+ # text on X axis
        theme(axis.text.y = element_text(size =  input$SV_Reference_xy_size))+
        theme(axis.title.y = element_text(size =  input$SV_Reference_xy_size))+
        theme( legend.title = element_blank())+
        theme(axis.title.x = element_blank())+
        theme(plot.title = element_text( face = "bold", #字体("plain", "italic", "bold", "bold.italic")
                                         colour = "black", #字体颜色
                                         size = 14,#字体大小
                                         hjust = .5, #调整轴标题1：纵轴靠上，横轴靠右；0.5居中；0：纵轴靠下，横轴靠左
                                         vjust = .5, #1：靠图边框；0靠近轴线；.5居中
                                         angle = 0 #为什么只对横轴标题有作用？
        ))
        
      ggplotly(p) %>%
        add_annotations(text=input$SV_Reference_legend_title, xref="paper", yref="paper",
                         x=1.03, xanchor="left",
                         y=0.8, yanchor="bottom",    # Same y as legend below
                         legendtitle=TRUE, showarrow=F ) %>%
        layout(legend=list(y=0.8, yanchor="top" ) )%>%
        layout(margin=list(l = input$SV_Reference_margin, t = input$SV_Reference_margin, 
                           r = input$SV_Reference_margin, b = input$SV_Reference_margin))
      
    }else {
      p <- ggplot(data=df_SV_Reference(), aes(x=type, y=value, fill=variable)) +
        geom_bar(stat="identity", position=position_dodge())+
        scale_fill_manual(values=c('#f15a57','#f2c46e','#83a4d7','#877dbb'))+
        theme_classic()+
        theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = input$SV_Reference_xy_size))+ # text on X axis
        theme(axis.text.y = element_text(size =  input$SV_Reference_xy_size))+
        theme(axis.title.y = element_text(size =  input$SV_Reference_xy_size))+
        theme(legend.title = element_blank())+
        theme(axis.title.x = element_blank())+
        theme(plot.title = element_text( face = "bold", #字体("plain", "italic", "bold", "bold.italic")
                                         colour = "black", #字体颜色
                                         size = 14,#字体大小
                                         hjust = .5, #调整轴标题1：纵轴靠上，横轴靠右；0.5居中；0：纵轴靠下，横轴靠左
                                         vjust = .5, #1：靠图边框；0靠近轴线；.5居中
                                         angle = 0 #为什么只对横轴标题有作用？
        ))
      
      ggplotly(p) %>%
        add_annotations(text=input$SV_Reference_legend_title, xref="paper", yref="paper",
                        x=1.03, xanchor="left",
                        y=0.8, yanchor="bottom",    # Same y as legend below
                        legendtitle=TRUE, showarrow=F ) %>%
        layout(legend=list(y=0.8, yanchor="top" ) )%>%
        layout(margin=list(l = input$SV_Reference_margin, t = input$SV_Reference_margin, 
                           r = input$SV_Reference_margin, b = input$SV_Reference_margin))
    }
    
    
  })
  
  
  ####################### Performance_Assessment###############################
  df_Performance_Assessment <-  reactive(
    if(input$PA_variant_type=="SNV"& input$PA_Performance=="Precision"){
      readRDS(file = paste('./data/DNAseq/Performance_Assessment/snv_precision_recall.rds')) 
    }else if(input$PA_variant_type=="SNV"& input$PA_Performance=="Reproducibility"){
      readRDS(file = paste('./data/DNAseq/Performance_Assessment/snv_mendelian_jaccard.rds')) 
    }else if(input$PA_variant_type=="SV"& input$PA_Performance=="Precision"){
      readRDS(file = paste('./data/DNAseq/Performance_Assessment/sv_precision_recall.rds')) 
    }else{
      readRDS(file = paste('./data/DNAseq/Performance_Assessment/sv_mendelian_jaccard.rds')) 
    }
  )

  output$plot_Performance_Assessment_legend <- renderUI({
    HTML(
      paste(
        "<br/>The '<b>SNV</b>' section represents the precision and recall of short-read sequencing datasets compared with reference small variants and Mendelian concordance within Quartet families. It also shows reproducibility between twins of short-read sequencing datasets. Boxplots on the edges display the distribution of precision and recall of PCR (red) and PCR-free (blue) datasets.",
        "<br/>The '<b>SV</b>' section represents the precision and recall of long-read sequencing datasets compared with Tier 1 reference SVs, considering genotypes. It includes Mendelian concordance within Quartet families and reproducibility between twins of long-read sequencing datasets. Boxplots on the edges show precision and recall of different callers. The callers used are cutesv, nanosv, pbsv, sniffle, and svim. The mappers used are ngmlr, minimaps, and pbmm2.",
        "<br/>Comparing with reference datasets can only estimate the performance for easy-to-detect variants. However, calculating Mendelian concordance within the four Quartet samples can estimate the performance on challenging variants and regions not included in the reference datasets."
      )
    )
  })

  output$Performance_Assessment_plot <- renderPlotly({
    subp <- function(){ 
      subplot(
        pp2, 
        plotly_empty(), 
        ggplotly(pmain),
        pp1,
        nrows = 2, heights = c(input$PA_porportion,1-input$PA_porportion),
        widths = c(1-input$PA_porportion, input$PA_porportion), 
        margin = 0,
        shareX = T, shareY = T, titleX =T, titleY = T
      ) %>%  
        layout(margin=list(l = input$PA_margin, t = input$PA_margin,
                           r = input$PA_margin, b = input$PA_margin))
    }
    dat <- df_Performance_Assessment()
    if(input$PA_variant_type=="SNV"& input$PA_Performance=="Precision"){
      pmain <- ggplot(dat, aes(x = METRIC.Precision, y = METRIC.Recall, color = library))+
        geom_point(aes(shape=Type))+
        scale_color_manual(values=c("#8e0c24", "#587aa5"))+
        scale_shape_manual(values=c(1,2,3,4,5,6,7,8,9,10,11))+ 
        theme_bw()+
        theme(legend.title = element_blank())+
        labs(x="Precision",y="Recall")
      # Marginal densities along x axis
      pp1 <- plot_ly(dat,
                     x = ~library,
                     y = ~ METRIC.Recall,
                     color = ~library,
                     colors =c("#8e0c24", "#587aa5"),
                     type='box',showlegend = FALSE) %>% 
        layout(xaxis=list(title = ""),yaxis=list(title = "Recall"))
      
      pp2 <- plot_ly(dat,
                     x = ~METRIC.Precision,
                     y = ~ library,
                     color = ~library,
                     colors =c("#8e0c24", "#587aa5"),
                     type='box' ) %>% 
        layout(xaxis=list(title="Precision"),yaxis=list(title = ""))
      subp()
      
    }else if(input$PA_variant_type=="SNV"& input$PA_Performance=="Reproducibility"){
      pmain <- ggplot(dat, aes(x = reproducibility, y = mendelian, color = library))+
        geom_point(aes(shape=type))+
        scale_color_manual(values=c("#8e0c24", "#587aa5"))+
        scale_shape_manual(values=c(1,2,3,4,5,6,7,8,9,10,11))+ 
        theme_bw()+
        theme(legend.title = element_blank())+
        labs(x="Reproducibility between twins",y="Mendelian concordant ratio")
      pp1 <- plot_ly(dat,
                     x = ~type,
                     y = ~ mendelian,
                     color = ~type,
                     colors =c("#8e0c24", "#587aa5"),
                     type='box' ) %>% 
        layout(xaxis=list(title=""),yaxis=list(title = "Mendelian concordant ratio"))
      pp2 <- plot_ly(dat,
                     x = ~reproducibility,
                     y = ~type,
                     color = ~type,
                     colors =c("#8e0c24", "#587aa5"),
                     type='box',showlegend = FALSE) %>% 
        layout(xaxis=list(title = "Reproducibility between twins"),yaxis=list(title = ""))
       subp()
      
    }else if(input$PA_variant_type=="SV"& input$PA_Performance=="Precision"){
      pmain <-  ggplot(dat, aes(x = Precision, y = Recall, color = caller))+
        geom_point(aes(shape=Pipelines))+
        scale_shape_manual(values=c(1,2,3,4,5,6,7,8,9,10,11))+ 
        theme_bw()+
        theme(legend.title = element_blank())+
        labs(x="Precision",y="Recall")
      # Marginal densities along x axis
      pp2 <- plot_ly(dat,
                     x = ~Precision,
                     y = ~ caller,
                     color = ~caller,
                     #colors =c("#8e0c24", "#587aa5"),
                     type='box',showlegend = FALSE) %>% 
        layout(xaxis=list(title = "Precision"),yaxis=list(title = ""))
      pp1 <- plot_ly(dat,
                     x = ~caller,
                     y = ~ Recall,
                     color = ~caller,
                     # colors =c("#8e0c24", "#587aa5"),
                     type='box' ) %>% 
        layout(xaxis=list(title=" "),yaxis=list(title = "Recall"))
      subp()
    }else{
      pmain <-  ggplot(dat, aes(x = Jaccard.index, y = Mendelian.Rate, color = caller))+
        geom_point(aes(shape=Pipelines))+
        scale_shape_manual(values=c(1,2,3,4,5,6,7,8,9,10,11))+ 
        theme_bw()+
        theme(legend.title = element_blank())+
        labs(x="Reproducibility between twins",y="Mendelian concordant ratio")
      pp2 <- plot_ly(dat,
                     x = ~Jaccard.index,
                     y = ~caller,
                     color = ~caller,
                     #colors =c("#8e0c24", "#587aa5"),
                     type='box',showlegend = FALSE) %>% 
        layout(xaxis=list(title = "Reproducibility between twins"),yaxis=list(title = ""))
      pp1 <- plot_ly(dat,
                     x = ~caller,
                     y = ~ Mendelian.Rate,
                     color = ~caller,
                     #colors =c("#8e0c24", "#587aa5"),
                     type='box' ) %>% 
        layout(xaxis=list(title=""),yaxis=list(title = "Mendelian concordant ratio"))
      subp()
      }
  })
  
  
  ####################### Variant_Statistics###############################
  
  df_Variant_Statistics <- reactive(
    if(input$Variant_Statistics_type=="SNV"){
      readRDS('./data/DNAseq/Variant_Statistics/snv_statistics.rds')
    }else{
      readRDS('./data/DNAseq/Variant_Statistics/sv_statistics.rds')
    }
  )

  output$plot_Variant_Statistics_legend <- renderUI({
    HTML(
      paste(
        "<br/>The bar colors represent how many times a variant was detected among all call sets. Green to red indicates an increase in the number of times a variant was detected among all replicates.",
        "<br/>The sequencing sites include ARD, BGI, BRG, NVG, WGE, and WUX. The sequencing platforms consist of DNB-SEQT7, Illumina Hiseq XTen, Illumina Novaseq6000, and MGI-SEQ2000. The library preparation methods include PCR and PCR-free. The callers used are cutesv, nanosv, pbsv, sniffle, and svim. The platforms for SV calling include ont (Nanopore), pb (Pacbio Sequel), and pb2 (Pacbio Sequel II)."
      )
    )
  })

  output$Variant_Statistics_plot <- renderPlotly({
    if(input$Variant_Statistics_type=="SNV"){
      dat_long <-  df_Variant_Statistics()
      colfunc <- colorRampPalette(c("#c6ffdd","#fbd786","#f7797d"))
      h1 <- ggplot(data= df_Variant_Statistics(), 
                   aes(x=plot_order, y=value, fill=variable)) +
        geom_bar(stat="identity")+
        scale_fill_manual(values=colfunc(30))+
        theme_bw()+
        theme(axis.text.x = element_blank(),
              axis.title.x=element_blank(),
              axis.text.y = element_text(size = 13,face = "bold",color = 'black'))+
        theme(legend.position = "none")+
        labs(y="Small Variants Number(million)")
      
      h2 <- ggplot(dat_long)+
        geom_bar(mapping = aes(x = plot_order, y = 1, fill = Library), 
                 stat = "identity", 
                 width = 1)+
        scale_fill_manual(values=c('#8e0c24','#587aa5'))+
        theme_void()+
        #facet_grid(~Quartet.Sample, scales = 'free_x', space = 'free_x')+
        theme(strip.background = element_blank(), strip.text = element_blank())+
        theme(panel.spacing.x = unit(1, "mm"))+
        theme(legend.position = "none")
    
      h3 <- ggplot(dat_long)+
        geom_bar(mapping = aes(x = plot_order, y = 1, fill = Sequence.Platform), 
                 stat = "identity", 
                 width = 1)+
        theme_void()+
        scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3","#E7298A"))+
        #facet_grid(~Quartet.Sample, scales = 'free_x', space = 'free_x')+
        theme_void()+
        theme(strip.background = element_blank(), strip.text = element_blank())+
        theme(panel.spacing.x = unit(1, "mm"))+
        theme(legend.position = "none")
  
      h4 <- ggplot(dat_long)+
        geom_bar(mapping = aes(x = plot_order, y = 1, fill = Sequence.Site), 
                 stat = "identity", 
                 width = 1)+
        scale_fill_manual(values=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33" ))+
        theme_void()+
        #facet_grid(~Quartet.Sample, scales = 'free_x', space = 'free_x')+
        theme(strip.background = element_blank(), strip.text = element_blank())+
        theme(panel.spacing.x = unit(1, "mm"))+
        theme(legend.position = "none")+
        theme(legend.position = "none")
   
      h6 <- ggplot(dat_long)+
        geom_bar(mapping = aes(x = plot_order, y = 1, fill = Mean.Deduped.Coverage), 
                 stat = "identity", 
                 width = 1)+
        theme_void()+
        theme(strip.background = element_blank(), strip.text = element_blank())+
        #facet_grid(~Quartet.Sample, scales = 'free_x', space = 'free_x')+
        scale_fill_gradient2(low = "green", 
                             mid = "yellow",
                             high = "red", 
                             midpoint = median(dat_long$Mean.Deduped.Coverage)) + 
        theme(panel.spacing.x = unit(1, "mm"))+
        theme(legend.position = "none")
      
      figH <- input$Variant_Statistics_figH
      figH1 <- (1-figH)/4
      
      subplot(ggplotly(h1),ggplotly(h6),ggplotly(h4),ggplotly(h3),ggplotly(h2),
              nrows = 5,
              heights = c(figH,figH1,figH1,figH1,figH1),
              margin=input$Variant_Statistics_fig_margin)%>% 
        layout(margin=list(l = input$Variant_Statistics_margin, 
                           t = input$Variant_Statistics_margin, 
                           r =input$Variant_Statistics_margin, 
                           b = input$Variant_Statistics_margin)) %>% 
        layout(showlegend = as.logical(input$Variant_Statistics_showlegend)) %>% 
        layout(yaxis=list(title="Small Variants Number(million)",
                          titlefont = list(size = input$Variant_Statistics_yaxisfont)
        )) 

    }else{
      dat_long <-  df_Variant_Statistics()
      colfunc <- colorRampPalette(c("#c6ffdd","#fbd786","#f7797d"))
      h1 <- ggplot(data=dat_long, aes(x=plot_order, y=value, fill=variable)) +
        geom_bar(stat="identity",color="grey")+
        scale_fill_manual(values=colfunc(30))+
        theme_bw()+
        theme(axis.text.x = element_blank(),
              axis.title.x=element_blank(),
              axis.text.y = element_text(size = 13,face = "bold",color = 'black'))+
        theme(legend.position = "none")

      h4 <- ggplot(dat_long)+
        geom_bar(mapping = aes(x = plot_order, y = 1, fill = caller), 
                 stat = "identity", 
                 width = 1)+
        scale_fill_brewer(palette="Paired")+
        theme_void()+
        theme(legend.position = "none")+
        theme(panel.spacing.x = unit(1, "mm"))
      
      h3 <- ggplot(dat_long)+
        geom_bar(mapping = aes(x = plot_order, y = 1, fill = mapper), 
                 stat = "identity", 
                 width = 1)+
        scale_fill_brewer(palette="Set1")+
        theme_void()+
        theme(legend.position = "none")+
        theme(panel.spacing.x = unit(1, "mm"))
      
      h2 <- ggplot(dat_long)+
        geom_bar(mapping = aes(x = plot_order, y = 1, fill = platform), 
                 stat = "identity", 
                 width = 1)+
        theme_void()+
        theme(legend.position = "none")+
        theme(panel.spacing.x = unit(1, "mm"))
      
      
      figH <- input$Variant_Statistics_figH
      figH1 <- (1-figH)/3
      subplot(ggplotly(h1),ggplotly(h4),ggplotly(h3),ggplotly(h2),
              nrows = 4,
              heights = c(figH,figH1,figH1,figH1),
              margin=input$Variant_Statistics_fig_margin)%>% 
        layout(margin=list(l = input$Variant_Statistics_margin, 
                           t = input$Variant_Statistics_margin, 
                           r =input$Variant_Statistics_margin, 
                           b = input$Variant_Statistics_margin)) %>% 
        layout(showlegend = as.logical(input$Variant_Statistics_showlegend)) %>% 
        layout(yaxis=list(title="Structural Variants Number",
                          titlefont = list(size = input$Variant_Statistics_yaxisfont)
        ))
    }
  
  })
  
  
  
  ####################### Jaccard_Index ###############################
  df_Jaccard_Index <- reactive(
    if(input$Jaccard_Index_type=="SNV"){
      read.table('./data/DNAseq/Jaccard_Index/small_variants_snv_jaccard_index.txt')
    }else if(input$Jaccard_Index_type=="INDEL"){
      read.table('./data/DNAseq/Jaccard_Index/small_variants_indel_jaccard_index.txt')
    }else if(input$Jaccard_Index_type=="Insertion"){
      read.table('./data/DNAseq/Jaccard_Index/sv_ins.txt')
    }else if(input$Jaccard_Index_type=="Deletion"){
      read.table('./data/DNAseq/Jaccard_Index/sv_del.txt')
    }else if(input$Jaccard_Index_type=="Duplication"){
      read.table('./data/DNAseq/Jaccard_Index/sv_dup.txt')
    }else if(input$Jaccard_Index_type=="Inversion"){
      read.table('./data/DNAseq/Jaccard_Index/sv_inv.txt')
    }else if(input$Jaccard_Index_type=="Breakend"){
      read.table('./data/DNAseq/Jaccard_Index/sv_bnd.txt')
    }
  )

  output$plot_Jaccard_Index_legend <- renderUI({
    HTML(
      paste(
        "<br/>For pairwise comparison of SNVs and INDELs between 27 small variant call sets of Quartet_D5, the color of each cell corresponds to the reproducibility of any two replicates. This is calculated as the number of shared small variants divided by the union of two call sets, with high similarity shown in red and low similarity in blue.",
        "<br/>For pairwise comparison of Insertions and Deletions between 30 SVs (Structural Variants) call sets of Quartet_D5, the Jaccard index is calculated based on the number of SVs shared divided by the union of two call sets, ignoring genotypes.",
        "<br/>For pairwise comparison of duplications, breakends, and inversions between 30 SVs call sets of Quartet_D5, the color of each cell corresponds to the Jaccard index, which is calculated as the number of variants intersection divided by the union of two call sets, with high similarity shown in red and low similarity in blue."
      )
    )
  })

  output$Jaccard_Index_plot <- renderPlotly({
    if(input$Jaccard_Index_type=="SNV"){
    heatmaply_SNV(df_Jaccard_Index(),
                  #Plot_method=input$Jaccard_Index_Plot_method, #ggplot,plotly
                  Margin= input$Jaccard_Index_margin,
                  Key.title = input$Jaccard_Index_legend_title,
                  Main = "SNV")
    }else if(input$Jaccard_Index_type=="INDEL"){
      heatmaply_SNV(df_Jaccard_Index(),
                    Margin=input$Jaccard_Index_margin,
                    Key.title = input$Jaccard_Index_legend_title,
                    Main = "INDEL")
    }else if(input$Jaccard_Index_type=="Insertion"){
      heatmaply_SV(df_Jaccard_Index(),
                   Margin= input$Jaccard_Index_margin,
                   Key.title = input$Jaccard_Index_legend_title,
                    Main = "Insertion")
    }else if(input$Jaccard_Index_type=="Deletion"){
      heatmaply_SV(df_Jaccard_Index(),
                   Key.title = input$Jaccard_Index_legend_title,
                   Margin= input$Jaccard_Index_margin,
                   Main = "Deletion")
    }else if(input$Jaccard_Index_type=="Duplication"){
      heatmaply_SV(df_Jaccard_Index(),
                   Key.title = input$Jaccard_Index_legend_title,
                   Margin= input$Jaccard_Index_margin,
                   Main = "Duplication")
    }else if(input$Jaccard_Index_type=="Inversion"){
      heatmaply_SV(df_Jaccard_Index(),
                   Key.title = input$Jaccard_Index_legend_title,
                   Margin= input$Jaccard_Index_margin,
                   Main = "Inversion")
    }else if(input$Jaccard_Index_type=="Breakend"){
      heatmaply_SV(df_Jaccard_Index(),
                   Key.title = input$Jaccard_Index_legend_title,
                   Margin= input$Jaccard_Index_margin,
                   Main = "Breakend")
    }
  })
  
  
  ####################### Quartet Advantage###############################
  df_Quartet_Advantage<-  reactive(
    if(input$Quartet_Advantage_type=="Small variants"){
      readRDS('./data/DNAseq/Quartet_Advantage/quartet_advantage_SNV.rds')
    }else{
      readRDS('./data/DNAseq/Quartet_Advantage/quartet_advantage_SV.rds')
    }
  )
  
  output$Quartet_Advantage_plot <- renderPlotly({
    if(input$Quartet_Advantage_type=="Small variants"){
      p <- ggplot(data=df_Quartet_Advantage(), aes(x=variable, y=value, fill=variable)) +
        geom_boxplot(fill="gray")+
        theme_classic()+
        theme(axis.text.x = element_text(size = input$Quartet_Advantage_xy_size))+ # text on X axis
        theme(axis.text.y = element_text(size = input$Quartet_Advantage_xy_size))+
        theme(axis.title.y = element_text(size = input$Quartet_Advantage_xy_size))+
        theme(axis.title.x = element_text(size = input$Quartet_Advantage_xy_size))+
        labs(title="Small variants",x=" ",y="Potential sequence errors (%)")+
        #theme(axis.title.y = element_text(size = 13))+
        theme(plot.title = element_text( face = "bold", #字体("plain", "italic", "bold", "bold.italic")
                                         colour = "black", #字体颜色
                                         size = input$Quartet_Advantage_title_size,#字体大小
                                         hjust = .5, #调整轴标题1：纵轴靠上，横轴靠右；0.5居中；0：纵轴靠下，横轴靠左
                                         vjust = .5, #1：靠图边框；0靠近轴线；.5居中
                                         angle = 0 #为什么只对横轴标题有作用？
        ))
      ggplotly(p)%>% layout(margin=list(l = input$Quartet_Advantage_margin, t = input$Quartet_Advantage_margin, 
                                         r = input$Quartet_Advantage_margin, b = input$Quartet_Advantage_margin))
      
   }else {
     p <- ggplot(data=df_Quartet_Advantage(), aes(x=variable, y=value, fill=variable)) +
       geom_boxplot(fill="gray")+
       theme_classic()+
       theme(axis.text.x = element_text(size = input$Quartet_Advantage_xy_size))+ # text on X axis
       theme(axis.text.y = element_text(size = input$Quartet_Advantage_xy_size))+
       theme(axis.title.y = element_text(size = input$Quartet_Advantage_xy_size))+
       theme(axis.title.x = element_text(size = input$Quartet_Advantage_xy_size))+
       labs(title="Structural variants",x=" ",y="Potential sequence errors (%)")+
       #theme(axis.title.y = element_text(size = 13))+
       theme(plot.title = element_text( face = "bold", #字体("plain", "italic", "bold", "bold.italic")
                                        colour = "black", #字体颜色
                                        size = input$Quartet_Advantage_title_size,#字体大小
                                        hjust = .5, #调整轴标题1：纵轴靠上，横轴靠右；0.5居中；0：纵轴靠下，横轴靠左
                                        vjust = .5, #1：靠图边框；0靠近轴线；.5居中
                                        angle = 0 #为什么只对横轴标题有作用？
       ))
     ggplotly(p)%>% layout(margin=list(l = input$Quartet_Advantage_margin, t = input$Quartet_Advantage_margin, 
                                       r = input$Quartet_Advantage_margin, b = input$Quartet_Advantage_margin))
    }
    
    
  })
  
  
  ####################### Small_Variants_Distribution###############################
  df_Small_Variants_Distribution<- reactive(
  read_feather("./data/DNAseq/Small_Variants_Distribution/Small_Variants_Distribution.feather")
  )

  output$plot_Small_Variants_Distribution_legend <- renderUI({
    HTML(
      paste(
        "<br/>The density plots display the sequencing depth, allele frequency, genotype quality, and mapping quality distributions for three types of Mendelian violated small variants and Mendelian consistent small variants. Pink indicates no call or conflicting variants in at least one of the four Quartet samples. Green indicates Mendelian consistent variants. Blue indicates variants overlapped with spanning deletions. Purple indicates detectable variants."
      )
    )
  })
  
  output$Small_Variants_Distribution_plot <- renderPlotly({
    snv_distri_plotly(df_Small_Variants_Distribution(),
                      xfactor=input$Small_Variants_Distribution_type,
                      Xlab=input$Small_Variants_Distribution_type,
                      axisxysize=input$Small_Variants_Distribution_xy_size,
                      Margin=input$Small_Variants_Distribution_margin
                      )
  })
  
  
  #######################Reference Datasets Summary###############################
  
  df_Reference_Datasets_Summary <- reactive(
    if(input$Reference_Datasets_Summary_type=="Small Variants"){
      readRDS('./data/DNAseq/Reference_Datasets_Summary/small_variant_reference_summary.rds')
    }else{
      readRDS('./data/DNAseq/Reference_Datasets_Summary/structral_variant_reference_summary.rds')
    }
  )

  output$plot_Reference_Datasets_Summary_legend <- renderUI({
    HTML(
      paste(
        "<br/>In the '<b>Small Variants</b>' section, the plot represents the count of small variants in 108 short-read whole-genome sequencing (WGS) libraries. The bar colors indicate the respective types of small variants.",
        "<br/>In the '<b>Structural Variants</b>' section, the plot displays the count of structural variants in 120 long-read sequencing call sets. The bar colors represent distinct SV (Structural Variant) types such as insertions (INS), deletions (DEL), duplications (DUP), inversions (INV), and breakends (BND).",
        "<br/>The sample consists of a father (F7), mother (M8), and monozygotic twin daughters (D5 and D6). The sequence sites include ARD, BGI, BRG, NVG, WGE, and WUX. The sequence platforms consist of DNB-SEQT7, Illumina Hiseq XTen, Illumina Novaseq6000, and MGI-SEQ2000. The library preparation methods include PCR and PCR-free. The callers used are cutesv, nanosv, pbsv, sniffle, and svim. The platforms for SV calling include ont (Nanopore), pb (Pacbio Sequel), and pb2 (Pacbio Sequel II)."
      )
    )
  })
  
  output$Reference_Datasets_Summary_plot <- renderPlotly({
    if(input$Reference_Datasets_Summary_type=="Small Variants"){
      dat_long <- df_Reference_Datasets_Summary()
      h1 <- ggplot(data=dat_long, aes(x=plot_order, y=value, fill=variable)) +
        geom_bar(stat="identity",color="black")+
        scale_fill_brewer(palette="Paired")+
        facet_grid(~Quartet.Sample, scales = 'free_x', space = 'free_x')+
        theme(strip.background = element_blank(),strip.text.x = element_blank())+
        theme_tufte()+
        labs(y="Variants Number(million)")+
        theme(legend.position = "none")+
        theme(axis.text.x = element_blank(),
              axis.title.x=element_blank(),
              axis.text.y = element_text(size = 13,face = "bold",color = 'black'))
    
      h2 <- ggplot(dat_long)+
        geom_bar(mapping = aes(x = plot_order, y = 1, fill = Library), 
                 stat = "identity", 
                 width = 1)+
        scale_fill_manual(values=c('#8e0c24','#587aa5'))+
        theme_void()+
        facet_grid(~Quartet.Sample, scales = 'free_x', space = 'free_x')+
        theme(strip.background = element_blank(),strip.text.x = element_blank())+
        theme(panel.spacing.x = unit(1, "mm"))+
        theme(legend.position = "none")

      h3 <- ggplot(dat_long)+
        geom_bar(mapping = aes(x = plot_order, y = 1, fill = Sequence.Platform), 
                 stat = "identity", 
                 width = 1)+
        theme_void()+
        scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3","#E7298A"))+
        facet_grid(~Quartet.Sample, scales = 'free_x', space = 'free_x')+
        theme_void()+
        theme(strip.background = element_blank(),strip.text.x = element_blank())+
        theme(panel.spacing.x = unit(1, "mm"))+
        theme(legend.position = "none")

      h4 <- ggplot(dat_long)+
        geom_bar(mapping = aes(x = plot_order, y = 1, fill = Sequence.Site), 
                 stat = "identity", 
                 width = 1)+
        scale_fill_manual(values=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33" ))+
        theme_void()+
        facet_grid(~Quartet.Sample, scales = 'free_x', space = 'free_x')+
        theme(strip.background = element_blank(),strip.text.x = element_blank())+
        theme(panel.spacing.x = unit(1, "mm"))+
        theme(legend.position = "none")

      h5 <- ggplot(dat_long)+
        geom_bar(mapping = aes(x = plot_order, y = 1, fill = Quartet.Sample), 
                 stat = "identity", 
                 width = 1)+
        theme_void()+
        scale_fill_manual(values=c('#4cc3d9','#7bc8a4','#ffc65d','#f16745'))+
        theme(strip.background = element_blank(),strip.text.x = element_blank())+
        facet_grid(~Quartet.Sample, scales = 'free_x', space = 'free_x')+
        theme(panel.spacing.x = unit(1, "mm"))+
        theme(legend.position = "none")

      h6 <- ggplot(dat_long)+
        geom_bar(mapping = aes(x = plot_order, y = 1, fill = Mean.Deduped.Coverage), 
                 stat = "identity", 
                 width = 1)+
        theme_void()+
        theme(strip.background = element_blank(), strip.text.x = element_blank())+
        facet_grid(~Quartet.Sample, scales = 'free_x', space = 'free_x')+
        scale_fill_gradient2(low = "green", 
                             mid = "yellow",
                             high = "red", 
                             midpoint = median(dat_long$Mean.Deduped.Coverage)) + 
        theme(panel.spacing.x = unit(1, "mm"))+
        theme(legend.position = "none")
 
      figH <- input$Reference_Datasets_Summary_figH
      figH1 <- (1-input$Reference_Datasets_Summary_figH)/5
      
      subplot(ggplotly(h1),ggplotly(h6),ggplotly(h5),ggplotly(h4),ggplotly(h3),ggplotly(h2),
              nrows = 6,heights = c(figH,figH1,figH1,figH1,figH1,figH1),
              margin=input$Reference_Datasets_Summary_fig_margin )%>% 
        layout(margin=list(l = input$Reference_Datasets_Summary_margin, 
                           t = input$Reference_Datasets_Summary_margin,
                           r = input$Reference_Datasets_Summary_margin,
                           b = input$Reference_Datasets_Summary_margin)) %>% 
        layout(showlegend = as.logical(input$Reference_Datasets_Summary_showlegend)) %>% 
        layout(yaxis=list(title="Variants Number",
                          titlefont = list(size = input$Reference_Datasets_Summary_yaxisfont)
        ))
    }else{
      dat_long <- df_Reference_Datasets_Summary()
      h1 <- ggplot(data=dat_long, aes(x=plot_order, y=value, fill=variable)) +
        geom_bar(stat="identity",color="black")+
        scale_fill_brewer(palette="Paired")+
        facet_grid(~Sample, scales = 'free_x', space = 'free_x')+
        theme_tufte()+
        theme(strip.background = element_blank(), strip.text.x = element_blank())+
        labs(y="Variants Number")+
        theme(axis.text.x =  element_blank(),
              axis.title.x = element_blank(),
              axis.text.y =  element_text(size = 13,face = "bold",color = 'black'),
              axis.title.y = element_text(size = 13,face = "bold",color = 'black'))+
        theme(legend.position = "none")

      h4 <- ggplot(dat_long)+
        geom_bar(mapping = aes(x = plot_order, y = 1, fill = Caller), 
                 stat = "identity", 
                 width = 1)+
        scale_fill_brewer(palette="Paired")+
        facet_grid(~Sample, scales = 'free_x', space = 'free_x')+
        theme_void()+
        theme(strip.background = element_blank(), strip.text.x = element_blank())+
        theme(panel.spacing.x = unit(1, "mm"))+
        theme(legend.position = "none")

      h3 <- ggplot(dat_long)+
        geom_bar(mapping = aes(x = plot_order, y = 1, fill = Mapper), 
                 stat = "identity", 
                 width = 1)+
        facet_grid(~Sample, scales = 'free_x', space = 'free_x')+
        theme_void()+
        theme(strip.background = element_blank(), strip.text.x = element_blank())+
        scale_fill_brewer(palette="Set1")+
        theme(panel.spacing.x = unit(1, "mm"))+
        theme(legend.position = "none")
 
      h5 <- ggplot(dat_long)+
        geom_bar(mapping = aes(x = plot_order, y = 1, fill = Sample), 
                 stat = "identity", 
                 width = 1)+
        facet_grid(~Sample, scales = 'free_x', space = 'free_x')+
        theme_void()+
        theme(strip.background = element_blank(), strip.text.x = element_blank())+
        scale_fill_manual(values=c('#4cc3d9','#7bc8a4','#ffc65d','#f16745'))+
        theme(panel.spacing.x = unit(1, "mm"))+
        theme(legend.position = "none")
  
      h2 <- ggplot(dat_long)+
        geom_bar(mapping = aes(x = plot_order, y = 1, fill = Platform), 
                 stat = "identity", 
                 width = 1)+
        facet_grid(~Sample, scales = 'free_x', space = 'free_x')+
        theme_void()+
        theme(strip.background = element_blank(), strip.text.x = element_blank())+
        theme(panel.spacing.x = unit(1, "mm"))+
        theme(legend.position = "none")

      figH <- input$Reference_Datasets_Summary_figH
      figH1 <- (1-input$Reference_Datasets_Summary_figH)/4
      
      subplot(ggplotly(h1),ggplotly(h5),ggplotly(h4),ggplotly(h3),ggplotly(h2),
              nrows =5,heights = c(figH,figH1,figH1,figH1,figH1),
              margin=input$Reference_Datasets_Summary_fig_margin
      )%>% 
        layout(margin=list(l = input$Reference_Datasets_Summary_margin,
                           t = input$Reference_Datasets_Summary_margin,
                           r = input$Reference_Datasets_Summary_margin, 
                           b = input$Reference_Datasets_Summary_margin)) %>% 
        layout(showlegend = as.logical(input$Reference_Datasets_Summary_showlegend)) %>% 
        layout(yaxis=list(title="Variants Number",
                          titlefont = list(size = input$Reference_Datasets_Summary_yaxisfont)
        ))

    }
    
  })
  
  
  
  ####################### Large Deletion###############################
  
  
  df_Large_Deletion<- reactive(
    if(input$Large_Deletion_Variant_type=="Insertion"){
      if(input$Large_Deletion_show_type=="Variant Number"){
      readRDS("./data/DNAseq/Large_Deletion/Insertion_Variant_Number.rds")
      }else if(input$Large_Deletion_show_type=="Filtered Variant"){
        readRDS("./data/DNAseq/Large_Deletion/Insertion_filtered_variant.rds")
      }
    }else{
      if(input$Large_Deletion_show_type=="Variant Number"){
       readRDS("./data/DNAseq/Large_Deletion/Deletion_Variant_Number.rds")
      }else if(input$Large_Deletion_show_type=="Filtered Variant"){
        readRDS("./data/DNAseq/Large_Deletion/Deletion_filtered_variant.rds")
      }
    }
  )

  output$plot_Large_Deletion_legend <- renderUI({
    HTML(
      paste(
        "<br/>In the '<b>Insertion</b>' section, the plot shows the total number of small variants aggregated over Tier 1 and Tier 2 INSs and DELs, divided by distance to an SV breakpoint in 500-bp bins. Red indicates variants that passed the voting process, while green indicates irreproducible and filtered variants among replicates.",
        "<br/>In the '<b>Deletion</b>' section, the plot displays the percentage of filtered variants aggregated over Tier 1 and Tier 2 INSs and DELs, divided by the distance to an SV breakpoint in 500-bp bins. The color and shape of lines represent different libraries and sequencing platforms.",
        "<br/>The sequence platforms consist of DNB-SEQT7, Illumina Hiseq XTen, Illumina Novaseq6000, and MGI-SEQ2000. The library preparation methods include PCR and PCR-free."
      )
    )
  })
  
  output$Large_Deletion_plot <- renderPlotly({
    if(input$Large_Deletion_show_type=="Filtered Variant"){
       p<-ggplot( df_Large_Deletion(), aes(x=gap, y=value, group=variable)) +
        geom_line(aes(color=SequencePlatform,linetype=Library))+
        geom_point(aes(color=SequencePlatform))+
        scale_color_manual(values=c('#c64b2b','#f9ce8c','#c3e3e5','#589bad'))+
        ylab(c('Percentage of Filtered Variants'))+
        xlab(c('Distance from SV breakpoints'))+
        theme_classic()+
        theme(axis.text.x = element_text(angle=90, 
                                         size=input$Large_Deletion_xy_size),
              axis.text.y = element_text( angle=90,
                                          size=input$Large_Deletion_xy_size))+
        theme(
          axis.title.x = element_text( size=input$Large_Deletion_xy_size),
          axis.title.y = element_text( size=input$Large_Deletion_xy_size)
        )+
        scale_x_continuous(breaks=seq(-10000,10000,1000))
      p <- ggplotly(p)
      p$x$layout$annotations[[1]]$text <- ""
      p %>% 
        layout(margin=list(l = input$Large_Deletion_margin,
                           t = input$Large_Deletion_margin,
                           r = input$Large_Deletion_margin, 
                           b = input$Large_Deletion_margin))
      
    }else if(input$Large_Deletion_show_type=="Variant Number"){
      p <- ggplot(data=df_Large_Deletion(), aes(x=gap, y=value, fill=type)) +
        geom_bar(stat="identity")+
        ylab(c('Number of Variants'))+
        xlab(c('Distance from SV breakpoints'))+
        theme_classic()+
        theme(axis.text.x = element_text(angle=90, 
                                         size=input$Large_Deletion_xy_size),
              axis.text.y = element_text( angle=90,
                                          size=input$Large_Deletion_xy_size))+
        theme(
          axis.title.x = element_text( size=input$Large_Deletion_xy_size),
          axis.title.y = element_text( size=input$Large_Deletion_xy_size)
        )+
        scale_x_continuous(breaks=seq(-10000,10000,1000))
      p <- ggplotly(p)
      p$x$layout$annotations[[1]]$text <- ""
      p %>% layout(margin=list(l = input$Large_Deletion_margin,
                           t = input$Large_Deletion_margin,
                           r = input$Large_Deletion_margin, 
                           b = input$Large_Deletion_margin))

    }
    
  })
  session$onSessionEnded(stopApp)
})






