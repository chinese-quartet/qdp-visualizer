
CalcSigNoiseRatio <- function(pca_prcomp, expDesign){
  
  pcs <- as.data.frame(predict(pca_prcomp))
  dt.perc.pcs <- data.table(PCX=1:nrow(pcs),
                            Percent=summary(pca_prcomp)$importance[2,],
                            AccumPercent=summary(pca_prcomp)$importance[3,])
  
  dt.dist <- data.table(ID.A = rep(rownames(pcs),each=nrow(pcs)),
                        ID.B = rep(rownames(pcs),time=nrow(pcs)))
  
  dt.dist$Group.A <- expDesign[dt.dist$ID.A]$group
  dt.dist$Group.B <- expDesign[dt.dist$ID.B]$group
  
  dt.dist[,Type:=ifelse(ID.A==ID.B,'Same',
                        ifelse(Group.A==Group.B,'Intra','Inter'))]
  dt.dist[,Dist:=sqrt(dt.perc.pcs[1]$Percent*(pcs[ID.A,1]-pcs[ID.B,1])^2+dt.perc.pcs[2]$Percent*(pcs[ID.A,2]-pcs[ID.B,2])^2)]
  
  dt.dist.stats <- dt.dist[,.(Avg.Dist=mean(Dist)),by=.(Type)]
  setkey(dt.dist.stats,Type)
  signoise <- dt.dist.stats['Inter']$Avg.Dist/dt.dist.stats['Intra']$Avg.Dist  
  return(signoise)
}

