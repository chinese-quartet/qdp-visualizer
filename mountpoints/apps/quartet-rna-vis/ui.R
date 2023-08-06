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
              menuItem("Gene Expression", 
                       tabName = "gene_exp", 
                       icon = icon("angle-double-right")), 
              fluidRow(
                column(11,
                       conditionalPanel(
                         "input.tabs === 'gene_exp'",
                         div(
                           style = "font-size:12px; font-family:sans-serif", 
                           column(12,
                                  selectizeInput("id_type",
                                                 label = 'Gene ID Type',
                                                 choices = c('Ensemble ID', 'Gene Symbol'),
                                                 selected = 'Ensemble ID'), 
                                  selectizeInput("gene_list",
                                                 label = 'Gene ID',
                                                 choices = NULL), 
                                  selectizeInput("group_a_gene",
                                                 label = 'Group',
                                                 choices = c('sample', 'batch', 'libraryPrep', 'kit', 'lab'),
                                                 selected = 'sample',
                                                 multiple = FALSE)
                                  )
                           )
                         )
                       )
                ),
              
              menuItem("Intra Batch", 
                       tabName = 'intra_batch', 
                       icon = icon("angle-double-right")),
              fluidRow(
                column(11,
                       conditionalPanel(
                         "input.tabs === 'intra_batch'",
                         div(
                           style = "font-size:12px; font-family:serif; font-weight:normal",
                           column(12,
                                  selectizeInput("qc_type_intra",
                                                 label = 'Intra QC Type',
                                                 choices = c('Detected_Gene', 'Tree_Metrics'),
                                                 selected = 'Tree_Metrics'), 
                                  selectizeInput("intra_metric",
                                                 label = 'Intra Metric',
                                                 choices = NULL),
                                  selectizeInput("group_a_intra",
                                                 label = 'Group',
                                                 choices = c('batch', 'libraryPrep', 'kit', 'lab', 'quality'),
                                                 selected = 'batch')
                                 )
                           )
                         )
                       )
                ),
              
              menuItem("Cross Batch", 
                       tabName = 'cross_batch', 
                       icon = icon("angle-double-right")),
              fluidRow(
                column(11,
                       conditionalPanel(
                         "input.tabs === 'cross_batch'",
                         div(
                           style = "font-size:12px; font-family:serif; font-weight:normal",
                           column(12,
                                  selectizeInput("qc_type_cross",
                                                 label = 'Cross Batch QC Type',
                                                 choices = c('SNR', 'Relative Correlation', 'Absolute Correlation'),
                                                 selected = 'SNR'),
                                  selectizeInput("figure_type_corss",
                                                 label = 'Figure Type',
                                                 choices = c('boxplot', 'heatmap'),
                                                 selected = 'boxplot'),
                                  selectizeInput("heatmap_sample",
                                              label = 'Sample',
                                              choices = NULL)
                                  )
                           )
                         )
                       )
                ),
              
              menuItem("PCA", 
                       tabName = 'pca', 
                       icon = icon("angle-double-right")),
              fluidRow(
                column(11,
                       conditionalPanel(
                         "input.tabs === 'pca'",
                         div(
                           style = "font-size:12px; font-family:serif; font-weight:normal",
                           column(12,
                                  materialSwitch(inputId = "d2_d3_switch", 
                                                 label = 'D2/D3', 
                                                 status = 'success'),
                                  selectizeInput("group_pca",
                                                 label = 'Group',
                                                 choices = pca_group,
                                                 selected = 'batch'),
                                  selectizeInput("pca_scale",
                                                 label = 'PCA Scale',
                                                 choices = c('no-scale', 'z-score', 'relative'),
                                                 selected = 'no-scale'),
                                  shinyWidgets::pickerInput("group_pca_unit",
                                                            label = 'Element of Group',
                                                            choices = NULL,
                                                            options = list('actions-box' = TRUE),
                                                            multiple = TRUE),
                                  selectizeInput("ref_sample",
                                                 label = 'Reference Sample',
                                                 choices = NULL)
                                  )
                           )
                         )
                       )
                ),
              
              menuItem("Reference Dataset", 
                       tabName = 'RD', 
                       icon = icon("angle-double-right")),
              fluidRow(
                column(11,
                       conditionalPanel(
                         "input.tabs === 'RD'",
                         div(
                           style = "font-size:12px; font-family:serif; font-weight:normal",
                           column(12,
                                  selectInput("data_type",
                                              label = 'Data Type',
                                              choices = c('Consensus of Detected Genes', 'Distribution of Gene Type', 'Consensus of DEGs', 
                                                          'Volcano for Ref DEGs', 'Performance Evaluation at RE', 'Performance Summary'),
                                              selected = 'Consensus of Detected Genes'),
                                  selectInput("rd_compare_group",
                                              label = 'Group',
                                              choices = NULL)
                                  )
                           )
                         )
                       )
                )
              )
  )

body <- dashboardBody(
  
  theme_onenote,

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

    # gene expression
    tabItem(
      tabName = "gene_exp",
      
      fluidRow(
        column(width = 12,
               box(width = NULL, 
                   tags$p("The plot is used to visualize the distribution of gene expression values across different groups of samples. If you want to know the expression of a specific gene, please select the gene ID type and gene ID in the left panel."),
                   withSpinner(plotlyOutput("plot_gene_exp", height = "500px")), 
                   uiOutput("plot_gene_exp_legend"),
                   collapsible = FALSE,
                  #  title = "Gene Expression", 
                   solidHeader = TRUE)
               )
        )
      ),
    
    #intra batch qc metrics
    tabItem(
      tabName = "intra_batch",
      
      fluidRow(
        column(width = 12,
               box(width = NULL, 
                   conditionalPanel(
                     "input.qc_type_intra === 'Tree_Metrics'",
                     tags$p("The lollipop plot serves to visually display the distribution of quality metrics across a variety of groups, where the metrics could include SNR, LIR, or LRR2. This distribution gives you insight into the level of quality as indicated by these metrics.")
                   ),
                   withSpinner(plotlyOutput("plot_intra_batch", height = "500px")), 
                   uiOutput("plot_intra_batch_legend"),
                   collapsible = FALSE,
                  #  title = "Intra-Batch",
                   solidHeader = TRUE)
               )
        )
      ),
    
    #cross batch qc metrics
    tabItem(
      tabName = "cross_batch",
      
      fluidRow(
        column(width = 12,
               box(width = NULL, 
                   tags$p("The plot is used to visualize the distribution of the value of SNR, relative correlation and absolute correlation of different datasets."),
                   withSpinner(plotlyOutput("plot_cross_batch", height = "500px")), 
                   uiOutput("plot_cross_batch_legend"),
                   collapsible = FALSE,
                  #  title = "Cross-Batch",  
                   solidHeader = TRUE)
               )
        )
      ),
    
    #PCA SNR
    tabItem(
      tabName = "pca",

      fluidRow(
        column(width = 12,
               box(width = NULL, 
                   tags$p("The plot is used to visualize the PCA distribution for any single or multiple datasets at different scales."),
                   withSpinner(plotlyOutput("plot_pca", height = "500px")), 
                   uiOutput("plot_pca_legend"),
                   collapsible = FALSE,
                  #  title = "PCA",  
                   solidHeader = TRUE)
               )
        )
      ), 
    
    #reference dataset
    tabItem(
      tabName = "RD",
      
      fluidRow(
        column(width = 12,
               box(width = NULL, 
                   tags$p("The plot summarises and presents the reference dataset performance metrics."),
                   withSpinner(plotlyOutput("plot_RD", height = "500px")), 
                   uiOutput("plot_RD_legend"),
                   collapsible = FALSE,
                  #  title = "Reference Dataset", 
                   solidHeader = TRUE)
               )
        )
      )
    )
  )

fillPage(
  dashboardPage(
    dashboardHeader(title = "Quartet Vis", titleWidth = 300, disable = TRUE),
    sidebar,
    body
  )
)


