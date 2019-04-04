#' 2) Runs QC analyses and normalizes data.
#' 
#' @export
#' @param raw ExpressionFeatureSet from either processGEOfiles or processCELfiles function
#' @param path Directory for plots
#' @return Normalized ExpressionFeatureSet object with normalized data and phenotype information
#' @examples 
#' norm_celfiles = QCnorm(celfiles,'Users/name/folderName/plots') 
#' @references See packages rgl, Biobase, heatmaply, oligo
#' @note Normalizes using rma, see oligo package
#' @note Outputs 3D PCA and similarity heatmap
#' @note Outputs pre-normalization plots, QC plots, post-normalization plots, see oligo package for details


QCnorm = function(raw,path) {
  library(rgl)
  library(Biobase)
  library(heatmaply)
  library(reshape2)
  library(plotly)
  library(reshape2)
  
  QCnorm_ERR = file(paste0(path,'/QCnorm.err'),open='wt')
  sink(QCnorm_ERR,type='message',append=TRUE)
  
  dat = hist(raw,which="all")  
  dat_x = melt(dat[[1]]$x)
  dat_y = melt(dat[[1]]$y)
  dat = data.frame(x=dat_x$value,y=dat_y$value,sample=dat_x$Var2)
 
  HistplotBN = plot_ly(dat, x=~x, y=~y, color=~sample, text = ~sample) %>%
    add_lines() %>%
    layout(title='Distribution before Normalization')
  htmlwidgets::saveWidget(HistplotBN, paste0(path,"/histBeforeNorm.html"))


  nbfacs=nrow(pData(raw))
  MAplotBN<-List()
  for (i in 1:nbfacs) {
    MAplotBN<-c(MAplotBN,paste0("/MAplotsBeforeNorm",i,".jpg"))
    jpeg(paste0(path,"/MAplotsBeforeNorm",i,".jpg"),width=5, height=5,units = "in", res = 300)
    MAplot(raw,which=i,plotFun=smoothScatter,refSamples=c(1:nbfacs), main='', cex=2)  #Raw MAplots
    dev.off() 
  }
  
  #Export boxplotBN, RLE and NUSE data for plots
  boxplotDataBN = oligo::boxplot(raw, which="all", plot=FALSE)                       
  colnames(boxplotDataBN[[1]]$stats) = boxplotDataBN[[1]]$names
  boxplotDataBN = boxplotDataBN[[1]]$stats
  
  qc = fitProbeLevelModel(raw)                                                        
  RLEdata = RLE(qc,type='values')
  NUSEdata = NUSE(qc, type='values')
  
  
  #Normalize data
  if (raw@annotation=="pd.hg.u133.plus.2" | raw@annotation=="pd.clariom.s.human.ht" | raw@annotation=="pd.clariom.s.human" | raw@annotation=="pd.clariom.s.mouse.ht" | raw@annotation=="pd.clariom.s.mouse" | raw@annotation=='pd.mouse430.2' | raw@annotation=='pd.hg.u133a' | raw@annotation=='pd.hg.u133a.2' | raw@annotation=='pd.hg.u219' | raw@annotation=='pd.mg.u74av2' | raw@annotation=='pd.mouse430a.2' | raw@annotation=='pd.moe430a' | raw@annotation=='pd.hg.u95av2' | raw@annotation=='pd.hg.u133b') {
    norm =rma(raw, background=TRUE, normalize=TRUE, subset=NULL)
  } else {
    norm =rma(raw, background=TRUE, normalize=TRUE, subset=NULL, target="core")
  }
  
  
  dat = hist(norm)  
  dat_x = melt(dat[['x']])
  dat_y = melt(dat[['y']])
  dat = data.frame(x=dat_x$value,y=dat_y$value,sample=dat_x$Var2)
  
  HistplotAN = plot_ly(dat, x=~x, y=~y, color=~sample, text = ~sample) %>%
    add_lines() %>%
    layout(title='Distribution after Normalization')
  htmlwidgets::saveWidget(HistplotAN, paste0(path,"/histAfterNorm.html"))
  
  
  MAplotAN<-List()
  for (i in 1:nbfacs) {
    MAplotAN<-c(MAplotAN,paste0("/MAplotsAfterNorm",i,".jpg"))
    jpeg(paste0(path,"/MAplotsAfterNorm",i,".jpg"),width=5, height=5,units = "in", res = 300)
    MAplot(norm,which=i,plotFun=smoothScatter,refSamples=c(1:nbfacs), main='', cex=2) #Normalized MAplots
    dev.off() 
  }
  
  # Output boxplotAN data for plots
  boxplotDataAN = oligo::boxplot(norm, which="all", plot=FALSE)                       
  colnames(boxplotDataAN$stats) = boxplotDataAN$names
  boxplotDataAN = boxplotDataAN$stats
  
  # Output data for 3D PCA #                                                                         
  tedf= t(exprs(norm))
  if (length(which(apply(tedf, 2, var)==0)) >= 0){
    tedf = tedf[ , apply(tedf, 2, var) != 0]
  }
  pca=prcomp(tedf, scale. = T)
  #pcaData = pca$x[,1:3]
  

  mat=as.matrix(dist(t(exprs(norm))))
  rownames(mat)=pData(norm)$title
  colnames(mat)=rownames(mat)
  
  heatmaply(
    mat,margins = c(80,120,60,40),
    colorRampPalette(colors = c("red", "yellow")),
    file = paste0(path,"/heatmapAfterNorm.html")
  )
  Heatmapolt<-"/heatmapAfterNorm.html"
  print("+++QCnorm+++")
  #return (List(HistplotBN,MAplotBN,boxplotDataBN,RLEdata,NUSEdata,HistplotAN,MAplotAN,boxplotDataAN,pca,Heatmapolt,norm))
  return (List(MAplotBN,boxplotDataBN,RLEdata,NUSEdata,MAplotAN,boxplotDataAN,pca,Heatmapolt,norm))
  sink(type='message')
}