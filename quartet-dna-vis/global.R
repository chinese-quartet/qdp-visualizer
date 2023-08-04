library(shiny)
library(shinyjs)
library(shinydashboard)
library(shinycssloaders)
library(shinyWidgets)
library(magrittr)
library(markdown)
library(ggplot2)
library(plotly)
library(DT)
library(rlang)
library(RColorBrewer)
library(lattice)
library(scales)
library(plyr)
library(dplyr)

# new
library(feather)
library(heatmaply)
library(ggthemes)
library(dashboardthemes)

#source('./code/data_convert/function_calc_SNR.R')
#source('./code/data_convert/pre_snr.R')

# new
source('./code/data_convert/plotly_heatmap.R')
source('./code/data_convert/SNV_SV_quality.R')
source('./code/data_convert/mendel_ratio_SV_SNV_plotly.R')
source('./code/data_convert/theme_vis.R')

####################### snv_distri_plotly #################

snv_distri_plotly <- function(all_df,xfactor="DP",Xlab='Depth',axisxysize=13,
                              Margin=100
){
  p <- ggplot(all_df, aes(x=.data[[xfactor]],  color=tag, fill=tag)) +
    geom_histogram( aes(y=..density..),
                    alpha=0.5, bins=50,
                    position="identity")+
    #  scale_color_manual(values=c("#00bfc4","#f8766d"))+
    #  scale_fill_manual(values=c("#00bfc4","#f8766d"))+
    theme_light()+
    ylab('Density') +
    xlab(Xlab)+
    labs(fill = " ")+
    theme(axis.text.x = element_text(size = axisxysize),
          axis.text.y = element_text(size = axisxysize),
          axis.title.x = element_text(size = axisxysize),
          axis.title.y = element_text(size = axisxysize))
  p <- ggplotly(p) %>% layout(margin=list(l=Margin,t=Margin,r=Margin,b=Margin))
  p$x$layout$annotations[[1]]$text <- ""
  p
}










