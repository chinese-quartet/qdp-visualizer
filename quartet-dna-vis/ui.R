#ui.R

sidebar <- dashboardSidebar(
  tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
  tags$style(".left-side, .main-sidebar {padding-top: 0px}"),
  width = 300,
  sidebarMenu(id="tabs",
              menuItem("Summary", 
                       tabName = "summary", 
                       icon = icon("angle-double-right"),
                       selected = TRUE), 
              menuItem("Mendelian Violations", 
                       tabName = "mendelian", 
                       icon = icon("angle-double-right")), 
              fluidRow(
                column(11,
                       conditionalPanel(
                         "input.tabs === 'mendelian'",
                         div(
                           #style = "font-size:10px; font-family:serif; font-weight:normal",
                           column(12, 
                                  
                                  selectInput("Mendelian_Violations",
                                              label = 'Mendelian Violations',
                                              choices = c("Mendelian Violations Ratio", "Mendelian Violation Quality"),
                                              selected = "Mendelian Violations Ratio"),
                                  conditionalPanel("input.Mendelian_Violations == 'Mendelian Violations Ratio'",
                                                   selectInput("Mendelian_Violations_type3",
                                                               label = 'Mendelian Violations Type',
                                                               choices = c("SNV","Indel","Deletion","Insertion"),
                                                               selected = "SNV"),
                                                   sliderInput("mendelian_title_size", "Title Size",
                                                               min = 0, max = 30, value = 15,  step = 1),
                                                   sliderInput("mendelian_plot_figHeight", "Main Figure Height",
                                                               min = 0, max = 1, value = 0.5,  step = 0.1),
                                                   sliderInput("mendelian_plot_figMargin", "Main Figure Margin",
                                                               min = 0, max = 1, value = 0.01,  step = 0.01),
                                                   #mendelian_plot_Showlegend
                                                   selectInput("mendelian_plot_Showlegend",
                                                               label = 'Show Legend',
                                                               choices=c("YES"=TRUE, "NO"=FALSE),
                                                               selected =TRUE) #,
                                                   
                                  ),
                                  
                                  conditionalPanel("input.Mendelian_Violations == 'Mendelian Violation Quality'",
                                                   selectInput("Mendelian_Violations_Quality",
                                                               label = 'Mendelian Violations Quality',
                                                               choices = c('Low Quality', 'High Quality'),
                                                               selected = 'Low Quality'),
                                                   conditionalPanel("input.Mendelian_Violations_Quality == 'Low Quality'",
                                                                    selectInput("Mendelian_Violations_type1",
                                                                                label = 'Mendelian Violations Type',
                                                                                choices = c("All Quartet Sample","Spanning Deletions","Detectable Mendelian Violations"),
                                                                                selected = "All Quartet Sample"),
                                                                    sliderInput("mendelian_plot_annotation_size", "Annotation Size",
                                                                                min = 0, max = 10, value = 4,  step = 0.1)
                                                   ),
                                                   conditionalPanel("input.Mendelian_Violations_Quality == 'High Quality'",
                                                                    selectInput("Mendelian_Violations_type2",
                                                                                label = 'Mendelian Violations Type',
                                                                                choices = c("Small Variants","Structural Variants"),
                                                                                selected = "Small Variants")
                                                   ) #,
                                                   
                                  ),
                                  
                                  sliderInput("mendelian_plot_xy_size", "X & Y Tick Size",
                                              min = 1, max = 30, value = 11, step = 1), 
                                  # sliderInput("mendelian_plot_annotation_size", "Annotation Size",
                                  #             min = 0, max = 5, value = 2.5,  step = 0.1),
                                  sliderInput("mendelian_plot_margin", "Margin",
                                              min = 0, max = 200, value = 50 , step = 5)
                                  
                                  
                           )
                         )
                       )
                )
              ),
              
              menuItem("Variant Quality", 
                       tabName = 'variant_quality', 
                       icon = icon("angle-double-right")),
              fluidRow(
                column(11,
                       conditionalPanel(
                         "input.tabs === 'variant_quality'",
                         div(
                           style = "font-size:12px; font-family:serif; font-weight:normal",
                           column(12,
                                  selectInput("Variant_Quality_type",
                                              label = 'Variant Quality Type',
                                              choices = c('Small Variants', 'Structural Variants'),
                                              selected = 'Structural Variants'),
                                  
                                  conditionalPanel("input.Variant_Quality_type == 'Small Variants'",
                                                   selectInput("Variant_Quality_distribution_type1",
                                                               label = 'Distribution Type',
                                                               choices = c("Depth","Allele Frequency","Genotype Quality","Mapping Quality"),
                                                               selected = "Depth")
                                  ),
                                  conditionalPanel("input.Variant_Quality_type == 'Structural Variants'",
                                                   selectInput("Variant_Quality_distribution_type2",
                                                               label = 'Distribution Type',
                                                               choices = c("Allele Frequency","RE"),
                                                               selected = "Allele Frequency")
                                  ),
                                  selectInput("Variant_Quality_plot_his_type",
                                              label = 'Histogram Type',
                                              choices = c('identity','stack','dodge'),
                                              selected = 'identity'),
                                  sliderInput("Variant_Quality_plot_his_height", "Histogram Height",
                                              min = 10, max = 200, value = 50, step = 10), 
                                  sliderInput("Variant_Quality_plot_xy_size", "X & Y Tick Size",
                                              min = 1, max = 30, value = 13, step = 1), 
                                  sliderInput("Variant_Quality_plot_margin", "Margin",
                                              min = 0, max = 200, value = 50 , step = 5)
                                  
                                  
                                  
                           )
                         )
                       )
                )
              ),
              
              menuItem("Difficult Genomic Regions", 
                       tabName = 'diff_region', 
                       icon = icon("angle-double-right")),
              fluidRow(
                column(11,
                       conditionalPanel(
                         "input.tabs === 'diff_region'",
                         div(
                           style = "font-size:12px; font-family:serif; font-weight:normal",
                           column(12,      
                                  selectInput("diff_region_variant_type",
                                              label = 'Variant Type',
                                              choices = c('Small Variants', 'Structural Variants'),
                                              selected = 'Small Variants'),
                                  # selectInput("Variant_Quality_plot_his_type",
                                  #             label = 'Histogram type',
                                  #             choices = c('identity','stack','dodge'),
                                  #             selected = 'identity'),
                                  sliderInput("diff_region_plot_subtitle_size", "Subtitle Size",
                                              min = 1, max = 30, value = 11, step = 1), 
                                  sliderInput("diff_region_plot_xy_size", "X & Y Tick Size",
                                              min = 1, max = 30, value = 13, step = 1), 
                                  sliderInput("diff_region_plot_margin", "Margin",
                                              min = 0, max = 200, value = 50 , step = 5)
                                  
                                  
                           )
                         )
                       )
                )
              ),
              
              menuItem("Variant Validation", 
                       tabName = 'variant_validation', 
                       icon = icon("angle-double-right")),
              fluidRow(
                column(11,
                       conditionalPanel(
                         "input.tabs === 'variant_validation'",
                         div(
                           style = "font-size:12px; font-family:serif; font-weight:normal",
                           column(12,
                                  selectInput("variant_validation_type",
                                              label = 'Variant type',
                                              choices = c('Small Variant Validation', 'SV Validation'),
                                              selected = 'Small Variant Validation'),
                                  sliderInput("variant_validation_subtitle_size", "Subtitle Size",
                                              min = 1, max = 30, value = 11, step = 1), 
                                  sliderInput("variant_validation_xy_size", "X & Y Tick Size",
                                              min = 1, max = 30, value = 13, step = 1), 
                                  sliderInput("variant_validation_margin", "Margin",
                                              min = 0, max = 200, value = 50 , step = 5)
                           )
                         )
                       )
                )
              ),
              
              menuItem("SV Reference", 
                       tabName = 'SV_Reference', 
                       icon = icon("angle-double-right")),
              fluidRow(
                column(11,
                       conditionalPanel(
                         "input.tabs === 'SV_Reference'",
                         div(
                           style = "font-size:12px; font-family:serif; font-weight:normal",
                           column(12,
                                  selectInput("SV_Reference_type",
                                              label = 'SV Reference Type',
                                              choices = c('Summary', 'Detail'),
                                              selected = 'Summary'),
                                  textInput("SV_Reference_legend_title", "Legend Title", value = "Type"),
                                  sliderInput("SV_Reference_xy_size", "X & Y Tick Size",
                                              min = 1, max = 30, value = 13, step = 1), 
                                  sliderInput("SV_Reference_margin", "Margin",
                                              min = 0, max = 200, value = 50 , step = 5)
                                  
                                  
                                  
                                  
                                  
                           )
                         )
                       )
                )
              ),
              
              # add new
              
              menuItem("Performance Assessment", 
                       tabName = 'Performance_Assessment', 
                       icon = icon("angle-double-right")),
              fluidRow(
                column(11,
                       conditionalPanel(
                         "input.tabs === 'Performance_Assessment'",
                         div(
                           style = "font-size:12px; font-family:serif; font-weight:normal",
                           column(12,  
                                  selectInput("PA_variant_type",
                                              label = 'Variant Type',
                                              choices = c('SNV', 'SV'),
                                              selected = 'SNV'),
                                  selectInput("PA_Performance",
                                              label = 'Performance',
                                              choices = c('Precision', 'Reproducibility'),
                                              selected = 'Precision'),
                                  #textInput("PA_Performance", "Legend_Title", value = "Type"),
                                  sliderInput("PA_porportion", "Figure Height",
                                              min = 0, max = 1, value = 0.2, step = 0.1), 
                                  sliderInput("PA_margin", "Margin",
                                              min = 0, max = 200, value = 50 , step = 5)
                                  
                                  
                           )
                         )
                       )
                )
              ),
              
              
              menuItem("Variant Statistics", 
                       tabName = 'Variant_Statistics', 
                       icon = icon("angle-double-right")),
              fluidRow(
                column(11,
                       conditionalPanel(
                         "input.tabs === 'Variant_Statistics'",
                         div(
                           style = "font-size:12px; font-family:serif; font-weight:normal",
                           column(12,
                                  selectInput("Variant_Statistics_type",
                                              label = 'Variant Type',
                                              choices = c('SNV', 'SV'),
                                              selected = 'SNV'),
                                  selectInput("Variant_Statistics_showlegend",
                                              label = 'Show Legend',
                                              choices=c("YES"=TRUE, "NO"=FALSE),
                                              selected =FALSE),
                                  sliderInput("Variant_Statistics_figH", "Figure Height",
                                              min = 0, max = 1, value = 0.8, step = 0.1), 
                                  sliderInput("Variant_Statistics_yaxisfont", "Yaxis Font",
                                              min = 0, max = 30, value = 15, step = 1),
                                  sliderInput("Variant_Statistics_fig_margin", "Figure Margin",
                                              min = 0, max = 0.1, value = 0, step = 0.01),
                                  sliderInput("Variant_Statistics_margin", "Margin",
                                              min = 0, max = 200, value = 50 , step = 5)
                                  
                           )
                         )
                       )
                )
              ),
              
              menuItem("Jaccard Index", 
                       tabName = 'Jaccard_Index', 
                       icon = icon("angle-double-right")),
              fluidRow(
                column(11,
                       conditionalPanel(
                         "input.tabs === 'Jaccard_Index'",
                         div(
                           style = "font-size:12px; font-family:serif; font-weight:normal",
                           column(12,
                                  selectInput("Jaccard_Index_type",
                                              label = 'Variant Type',
                                              choices = c('SNV', 'INDEL',
                                                          'Insertion','Deletion','Duplication','Inversion','Breakend'),
                                              selected = 'SNV'),
                                  textInput("Jaccard_Index_legend_title", "Legend Title", value = "Jaccard Index"),
                                  sliderInput("Jaccard_Index_margin", "Margin",
                                              min = 0, max = 200, value = 100 , step = 5)
                                  
                           )
                         )
                       )
                )
              ),
              
              menuItem("Quartet Advantage", 
                       tabName = 'Quartet_Advantage', 
                       icon = icon("angle-double-right")),
              fluidRow(
                column(11,
                       conditionalPanel(
                         "input.tabs === 'Quartet_Advantage'",
                         div(
                           style = "font-size:12px; font-family:serif; font-weight:normal",
                           column(12,
                                  selectInput("Quartet_Advantage_type",
                                              label = 'Variant Type',
                                              choices = c('Small variants', 'Structural variants'),
                                              selected = 'Small variants'),
                                  textInput("Jaccard_Index_legend_title", "Legend Title", value = "Jaccard Index"),
                                  sliderInput("Quartet_Advantage_xy_size", "X & Y Tick Size",
                                              min = 0, max = 20, value = 11 , step = 1),
                                  sliderInput("Quartet_Advantage_title_size", "Title Size",
                                              min = 0, max = 20, value = 14 , step = 1),
                                  sliderInput("Quartet_Advantage_margin", "Margin",
                                              min = 0, max = 200, value = 120 , step = 5)
                                  
                                  
                           )
                         )
                       )
                )
              ), 
              menuItem("Small Variants Distribution", 
                       tabName = 'Small_Variants_Distribution', 
                       icon = icon("angle-double-right")),
              fluidRow(
                column(11,
                       conditionalPanel(
                         "input.tabs === 'Small_Variants_Distribution'",
                         div(
                           style = "font-size:12px; font-family:serif; font-weight:normal",
                           column(12,
                                  selectInput("Small_Variants_Distribution_type",
                                              label = 'Distribution Type',
                                              choices = c("Depth","Allele Frequency","Genotype Quality","Mapping Quality"),
                                              selected = "Depth"),
                                  sliderInput("Small_Variants_Distribution_xy_size", "X & Y Tick Size",
                                              min = 0, max = 20, value = 11 , step = 1),
                                  sliderInput("Small_Variants_Distribution_margin", "Margin",
                                              min = 0, max = 300, value = 100 , step = 5)
                                  
                           )
                         )
                       )
                )
              ), 
              menuItem("Reference Datasets Summary", 
                       tabName = 'Reference_Datasets_Summary', 
                       icon = icon("angle-double-right")),
              fluidRow(
                column(11,
                       conditionalPanel(
                         "input.tabs === 'Reference_Datasets_Summary'",
                         div(
                           style = "font-size:12px; font-family:serif; font-weight:normal",
                           column(12,
                                  selectInput("Reference_Datasets_Summary_type",
                                              label = 'Variant Type',
                                              choices = c('Small Variants', 'Structural Variants'),
                                              selected = 'Small Variants'),
                                  selectInput("Reference_Datasets_Summary_showlegend",
                                              label = 'Show Legend',
                                              choices=c("YES"=TRUE, "NO"=FALSE),
                                              selected =FALSE),
                                  sliderInput("Reference_Datasets_Summary_figH", "Figure Height",
                                              min = 0, max = 1, value = 0.8, step = 0.1), 
                                  sliderInput("Reference_Datasets_Summary_yaxisfont", "Yaxis Font",
                                              min = 0, max = 30, value = 15, step = 1),
                                  sliderInput("Reference_Datasets_Summary_fig_margin", "Figure Margin",
                                              min = 0, max = 0.1, value = 0, step = 0.01),
                                  sliderInput("Reference_Datasets_Summary_margin", "Margin",
                                              min = 0, max = 200, value = 50 , step = 5)
                                  
                           )
                         )
                       )
                )
              ), 
              menuItem("Large Deletion", 
                       tabName = 'Large_Deletion', 
                       icon = icon("angle-double-right")),
              fluidRow(
                column(11,
                       conditionalPanel(
                         "input.tabs === 'Large_Deletion'",
                         div(
                           style = "font-size:12px; font-family:serif; font-weight:normal",
                           column(12,
                                  selectInput("Large_Deletion_Variant_type",
                                              label = 'Variant Type',
                                              choices = c('Insertion', 'Deletion'),
                                              selected = 'Insertion'),
                                  selectInput("Large_Deletion_show_type",
                                              label = 'Show Type',
                                              choices=c("Variant Number", "Filtered Variant"),
                                              selected ="Filtered Variant"),
                                  sliderInput("Large_Deletion_xy_size", "X & Y Tick Size",
                                              min = 0, max = 30, value = 15, step = 1),
                                  sliderInput("Large_Deletion_margin", "Margin",
                                              min = 0, max = 200, value = 50 , step = 5)
                                  
                           )
                         )
                       )
                )
              )
              
              # add new
              
              
              
  )
)

body <- dashboardBody(
  theme_onenote,
  #shinyDashboardThemes(
    #theme = "onenote",
    #theme = "poor_mans_flatly"
    
 # ),

  tabItems(
    # summary
    tabItem(
      tabName = "summary",
      selected = TRUE,
      
      fluidRow(
        column(width = 12,
                box(width = NULL, 
                    collapsible = FALSE,
                    solidHeader = TRUE,
                    uiOutput("summary", height = "500px", class = "content-box"))
                )
        )
    ),

    # Mendelian Violations
    tabItem(
      tabName = "mendelian",
      
      fluidRow(
        column(width = 12,
               box(width = NULL, 
                   withSpinner(plotlyOutput("mendelian_plot", height = "600px")), 
                   collapsible = FALSE,
                  #  title = "Mendelian Violations", 
                   solidHeader = TRUE)
        )
      )
    ),
    
    #Variant Quality
    tabItem(
      tabName = "variant_quality",
      
      fluidRow(
        column(width = 12,
               box(width = NULL, 
                   withSpinner(plotlyOutput("Variant_Quality_plot", height = "600px")), 
                   collapsible = FALSE,
                   title = "Variant Quality",
                   solidHeader = TRUE)
        )
      )
    ),
    
    #Difficult Genomic Regions
    tabItem(
      tabName = "diff_region",
      
      fluidRow(
        column(width = 12,
               box(width = NULL, 
                   withSpinner(plotlyOutput("diff_region_plot", height = "600px")), 
                   collapsible = FALSE,
                   title = "Difficult Genomic Regions",  
                   solidHeader = TRUE)
        )
      )
    ),
    
    #Variant Validation
    tabItem(
      tabName = "variant_validation",
      
      fluidRow(
        column(width = 12,
               box(width = NULL, 
                   withSpinner(plotlyOutput("variant_validation_plot", height = "600px")), 
                   collapsible = FALSE,
                   title = "Variant Validation",  
                   solidHeader = TRUE)
        )
      )
    ), 
    
    #SV Reference
    tabItem(
      tabName = "SV_Reference",
      
      fluidRow(
        column(width = 12,
               box(width = NULL, 
                   withSpinner(plotlyOutput("SV_Reference_plot", height = "600px")), 
                   collapsible = FALSE,
                   title = "SV Reference", 
                   solidHeader = TRUE)
        )
      )
    ),
    
    #Performance_Assessment
    tabItem(
      tabName = "Performance_Assessment",
      
      fluidRow(
        column(width = 12,
               box(width = NULL, 
                   withSpinner(plotlyOutput("Performance_Assessment_plot", height = "600px")), 
                   collapsible = FALSE,
                   title = "Performance Assessment",  
                   solidHeader = TRUE)
        )
      )
    ), 
    #Variant_Statistics
    tabItem(
      tabName = "Variant_Statistics",
      
      fluidRow(
        column(width = 12,
               box(width = NULL, 
                   withSpinner(plotlyOutput("Variant_Statistics_plot", height = "600px")), 
                   collapsible = FALSE,
                   title = "Variant Statistics",  
                   solidHeader = TRUE)
        )
      )
    ), 
    #Jaccard_Index
    tabItem(
      tabName = "Jaccard_Index",
      
      fluidRow(
        column(width = 12,
               box(width = NULL, 
                   withSpinner(plotlyOutput("Jaccard_Index_plot", height = "600px")), 
                   collapsible = FALSE,
                   title = "Jaccard Index",  
                   solidHeader = TRUE)
        )
      )
    ), 
    #Quartet_Advantage
    tabItem(
      tabName = "Quartet_Advantage",
      
      fluidRow(
        column(width = 12,
               box(width = NULL, 
                   withSpinner(plotlyOutput("Quartet_Advantage_plot", height = "600px")), 
                   collapsible = FALSE,
                   title = "Quartet Advantage",  
                   solidHeader = TRUE)
        )
      )
    ), 
    #Small_Variants_Distribution
    tabItem(
      tabName = "Small_Variants_Distribution",
      
      fluidRow(
        column(width = 12,
               box(width = NULL, 
                   withSpinner(plotlyOutput("Small_Variants_Distribution_plot", height = "600px")), 
                   collapsible = FALSE,
                   title = "Small Variants Distribution",  
                   solidHeader = TRUE)
        )
      )
    ), 
    #Reference_Datasets_Summary
    tabItem(
      tabName = "Reference_Datasets_Summary",
      
      fluidRow(
        column(width = 12,
               box(width = NULL, 
                   withSpinner(plotlyOutput("Reference_Datasets_Summary_plot", height = "600px")), 
                   collapsible = FALSE,
                   title = "Reference Datasets Summary",  
                   solidHeader = TRUE)
        )
      )
    ), 
    #Large_Deletion
    tabItem(
      tabName = "Large_Deletion",
      
      fluidRow(
        column(width = 12,
               box(width = NULL, 
                   withSpinner(plotlyOutput("Large_Deletion_plot", height = "600px")), 
                   collapsible = FALSE,
                   title = "Large Deletion",  
                   solidHeader = TRUE)
        )
      )
    )
    
    #add new
    
  )
)

dashboardPage(
  dashboardHeader(title = "Quartet Vis DNA",titleWidth = 300, disable = TRUE),
  sidebar,
  body
)
