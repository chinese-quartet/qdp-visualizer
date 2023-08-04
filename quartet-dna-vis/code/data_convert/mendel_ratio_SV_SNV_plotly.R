mendel_ratio_SV_SNV_plotly <- function(df_mendel,titlesize=15,ytitle=13,
                                       Title="SNV",figHeight=0.5,figMargin=0.01,
                                       Showlegend=T,Margin=100,variantType="SNV",
                                       Group = "D5_pipelines", #Samples
                                       Color = "caller"   #library
                                       
){
  if(variantType=="SV"){
    VV <- c("#f15a57",'#f2c46e','#b4c6e6','#83a4d7','#877dbb')
  }else{
    VV <- c("#8e0c24",'#587aa5')
  }
  
  ins_h1 <- ggplot(data=df_mendel[[1]], aes(x=vote, y=measurement, group=.data[[Group]],color=.data[[Color]])) +
    geom_line(aes())+
    geom_point()+
    scale_color_manual(values=VV)+
    ylim(c(0,100))+
    theme(plot.title = element_text(size=titlesize, color="black", face= "bold", vjust=0.5, hjust=0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=ytitle, color="black", vjust=0.5, hjust=0.5),
          axis.text.x = element_blank()
    )+ labs(title = Title,y="Mendalian Violation (%)")
  
  ins_h2 <- ggplot(df_mendel[[2]], aes(x=variable, y=value))+
    geom_bar(stat="identity", color="black",fill="#d1d0c8")+
    geom_errorbar(aes(ymin=value-se, ymax=value+se),width=.2,position=position_dodge(.9))+
    #scale_y_sqrt()+
    #ylim(c(0,3000))+
    theme(axis.title.x=element_blank(),
          axis.title.y = element_text(size=ytitle, color="black", vjust=0.5, hjust=0.5),
    )+labs(y="Variants Number")
  
  if(variantType=="SV"){
    ins_h2 <- ins_h2 +ylim(c(0,3000))
  }else{
    ins_h2 <- ins_h2+scale_y_sqrt()
  }
  
  pp <- subplot(
    ggplotly(ins_h1),ggplotly(ins_h2),
    nrows = 2, heights = c(figHeight, 1-figHeight), margin = figMargin,
    shareX = TRUE, shareY = TRUE, titleX =T, titleY = T
  ) %>%   layout(showlegend = Showlegend) %>% 
    layout(margin=list(l = Margin, t =Margin, r = Margin, b = Margin))
  pp$x$layout$annotations[[1]]$text <- ''
  
  pp
  
}

# df_mendel <- readRDS("SV_mendelian_ratio_DEL.rds")
# mendel_ratio_SV_SNV_plotly(df_mendel,variantType="SV")
# df_mendel <- readRDS("SV_mendelian_ratio_INS.rds")
# mendel_ratio_SV_SNV_plotly(df_mendel,variantType="SV")
# 
# df_mendel <- readRDS("mendelian_ratio_SNV.rds")
# mendel_ratio_SV_SNV_plotly(df_mendel,variantType="SNV",Group= "Samples",Color = "library")
# df_mendel <- readRDS("mendelian_ratio_SNV.rds")
# mendel_ratio_SV_SNV_plotly(df_mendel,variantType="SNV",Group= "Samples",Color = "library")
