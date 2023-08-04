# small_vaiant_quality 
# SV_quality

small_vaiant_quality <- function(big_df,xfactor="dp",axisxysize=13,
                                 Xlab="Mapping Quality",Bins=50,
                                 Position="identity"){
  ggplot(big_df, aes(x=.data[[xfactor]],color=tag, fill=tag)) + 
    geom_histogram(aes(y=..density..), alpha=0.5, bins=Bins,
                   position=Position)+
    #  geom_density(alpha=.2)+
    scale_color_manual(values=c("#00bfc4","#f8766d"))+
    scale_fill_manual(values=c("#00bfc4","#f8766d"))+
    theme_light()+
    ylab('Density') + 
    xlab(Xlab)+
    theme(axis.text.x = element_text(size = axisxysize),
          axis.text.y = element_text(size = axisxysize),
          axis.title.x = element_text(size = axisxysize),
          axis.title.y = element_text(size = axisxysize))
}

SV_quality <- function(big_df,xfactor="af",axisxysize=13,
                       Xlab="Allele Frequence",Bins=50,
                       Position="identity" ){
  ggplot(big_df, aes(x=.data[[xfactor]],color=tag, fill=tag)) + 
    geom_histogram(aes(y=..density..), alpha=0.5, bins=Bins,
                   position=Position)+
    #  geom_density(alpha=.2)+
    # scale_color_manual(values=c("#00bfc4","#f8766d"))+
    # scale_fill_manual(values=c("#00bfc4","#f8766d"))+
    theme_light()+
    ylab('Density') + 
    xlab(Xlab)+
    theme(axis.text.x = element_text(size = axisxysize),
          axis.text.y = element_text(size = axisxysize),
          axis.title.x = element_text(size = axisxysize),
          axis.title.y = element_text(size = axisxysize))
}


