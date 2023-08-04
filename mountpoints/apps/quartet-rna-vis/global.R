library('shiny')
library('shinyjs')
library('shinydashboard')
library('shinycssloaders')
library('shinyWidgets')
library('magrittr')
library('markdown')
library('ggplot2')
library('plotly')
library('DT')
library('rlang')
library('RColorBrewer')
library('lattice')
library('scales')
library('plyr')
library('dplyr')
library('purrr')
library('dashboardthemes')

source('./code/data_convert/function_calc_SNR.R')
source('./code/data_convert/pre_snr.R')
source('./code/data_convert/theme_vis.R')
##########################function####################
make_subbar_plotly <- function(batch_name, df_detect_gene){
  dt_c = df_detect_gene[n_rep == 1][, N1:=df_detect_gene[n_rep == 1][['N']]][, N2:=df_detect_gene[n_rep == 2][['N']]][, N3:=df_detect_gene[n_rep == 3][['N']]]
  print(dt_c)
  dt_b = dt_c[group == batch_name]
  fig = plot_ly(dt_b, x = ~sample, y = ~N3, type = 'bar', name = '3', color = '#1f77b4', alpha = 1) %>% 
    add_trace(x = ~sample, y = ~N2, name = '2', color = '#17becf', alpha = 0.5) %>%
    add_trace(x = ~sample, y = ~N1, name = '1', color = 'blue', alpha = 0.1) %>%
    layout(yaxis = list(title = 'Count'), barmode = 'stack', xaxis = list(title = batch_name))
  return(fig)
}

############################data#######################
print(Sys.time())
# meta data
dt_meta <-  readRDS(file = './data/rnaseq/qc_type/metadata_RNA.rds')
dt_meta_used <- dt_meta[, c('library', 'libraryPrep', 'lab', 'sample', 'kit', 'sequencingSite', 'batch')]

# expression data
dt_exp <- fread('./data/rnaseq/qc_type/exprSet_RNA_wide.txt')
df_gene <- dt_exp[, 1:2]
dt_exp_long <- data.table::melt(dt_exp, id = 1:2)
setnames(dt_exp_long, c("variable", "value"), c("library","logfpkm"))
df_exp_annot <- dt_meta_used[dt_exp_long, on = 'library']

# pca gene list
pca_group <- c('batch', 'libraryPrep', 'lab', 'sample', 'kit', 'sequencingSite')
pca_gene_list <- fread(file = './data/rnaseq/qc_type/pca_gene_list.txt')
pal <- c("#4CC3D9","#7BC8A4", "#FFC65D", "#F16745")
print(Sys.time())
