#' 2) Runs QC analyses and normalizes data.
#' 
#' @export
#' @param raw ExpressionFeatureSet from either processGEOfiles or processCELfiles function
#' @param path Directory for plots
#' @param contrast Contrast for subsetting QC plots by samples in chosen groups
#' @return Normalized ExpressionFeatureSet object with normalized data and phenotype information
#' @examples 
#' norm_celfiles = RMA_QCnorm(celfiles,'Users/name/folderName/plots',c("RNA_1-Ctl","RNA_2-Ctl")) 
#' @references See packages rgl, Biobase, heatmaply, oligo
#' @note Normalizes using rma, see oligo package
#' @note Outputs 3D PCA and similarity heatmap
#' @note Outputs pre-normalization plots, QC plots, post-normalization plots, see oligo package for details


RMA_QCnorm = function(raw,path,contrast) {
  library(rgl)
  library(Biobase)
  library(heatmaply)
  library(reshape2)
  library(plotly)
  library(reshape2)
  library(amap)
  library(gplots)
  library(limma)
  
  QCnorm_ERR = file(paste0(path,'/RMA_QCnorm.err'),open='wt')
  sink(QCnorm_ERR,type='message',append=TRUE)
  
  # histogram before normalization
  svg(paste0(path,"/temp.svg"),width=8, height=8)
  dat = hist(raw,which="all")
  dev.off()
  
  dat_x = melt(dat[[1]]$x)
  dat_y = melt(dat[[1]]$y)
  dat = data.frame(x=dat_x$value,y=dat_y$value,sample=dat_x$Var2)
  dat$sample = gsub('\n','',dat$sample)
 
  HistplotBN = plot_ly(dat, x=~x, y=~y, color=~sample, text = ~sample) %>%
    add_lines() %>%
    layout(title='Distribution before Normalization',xaxis=list(title='log-intensity'),yaxis=list(title='density'),legend=list(x=13,y=0.1))
  htmlwidgets::saveWidget(HistplotBN, paste0(path,"/histBeforeNorm.html"))

  # MA plots before normalization
  nbfacs=nrow(pData(raw))
  MAplotBN<-List()
  for (i in 1:nbfacs) {
    MAplotBN<-c(MAplotBN,paste0("/MAplotsBeforeNorm",i,".jpg"))
    jpeg(paste0(path,"/MAplotsBeforeNorm",i,".jpg"),width=5, height=5,units = "in", res = 300)
    MAplot(raw,which=i,plotFun=smoothScatter,refSamples=c(1:nbfacs), main='', cex=2)  #Raw MAplots
    dev.off() 
  }
  
  #Export boxplotBN, RLE and NUSE data for plots
  boxplotDataBN = list()
  temp = oligo::boxplot(raw, which="all", plot=FALSE)                       
  colnames(temp[[1]]$stats) = temp[[1]]$names
  boxplotDataBN[[1]] = temp[[1]]$stats
  boxplotDataBN[[2]] = 'log-intensity'
  
  RLEdata = list()
  NUSEdata = list()
  qc = fitProbeLevelModel(raw)                                                        
  RLEdata[[1]] = RLE(qc,type='values')
  RLEdata[[2]] = 'RLE'
  NUSEdata[[1]] = NUSE(qc, type='values')
  NUSEdata[[2]] = 'NUSE'
  
  
  #Normalize data
  if (raw@annotation=="pd.hg.u133.plus.2" | raw@annotation=="pd.clariom.s.human.ht" | raw@annotation=="pd.clariom.s.human" | raw@annotation=="pd.clariom.s.mouse.ht" | raw@annotation=="pd.clariom.s.mouse" | raw@annotation=='pd.mouse430.2' | raw@annotation=='pd.hg.u133a' | raw@annotation=='pd.hg.u133a.2' | raw@annotation=='pd.hg.u219' | raw@annotation=='pd.mg.u74av2' | raw@annotation=='pd.mouse430a.2' | raw@annotation=='pd.moe430a' | raw@annotation=='pd.hg.u95av2' | raw@annotation=='pd.hg.u133b') {
    norm =rma(raw, background=TRUE, normalize=TRUE, subset=NULL)
  } else {
    norm =rma(raw, background=TRUE, normalize=TRUE, subset=NULL, target="core")
  }
  
  # histogram after normalization
  svg(paste0(path,"/temp2.svg"),width=8, height=8)
  dat = hist(norm)
  dev.off()
  
  dat_x = melt(dat[['x']])
  dat_y = melt(dat[['y']])
  dat = data.frame(x=dat_x$value,y=dat_y$value,sample=dat_x$Var2)
  dat$sample = gsub('\n','',dat$sample)
  
  HistplotAN = plot_ly(dat, x=~x, y=~y, color=~sample, text = ~sample) %>%
    add_lines() %>%
    layout(title='Distribution after RMA Normalization',xaxis=list(title='log-intensity'),yaxis=list(title='density'),legend=list(x=13,y=0.1))
  htmlwidgets::saveWidget(HistplotAN, paste0(path,"/histAfterRMAnorm.html"))
  
  # MAplot after normalization
  MAplotAN<-List()
  for (i in 1:nbfacs) {
    MAplotAN<-c(MAplotAN,paste0("/MAplotsAfterRMAnorm",i,".jpg"))
    jpeg(paste0(path,"/MAplotsAfterRMAnorm",i,".jpg"),width=5, height=5,units = "in", res = 300)
    MAplot(norm,which=i,plotFun=smoothScatter,refSamples=c(1:nbfacs), main='', cex=2) #Normalized MAplots
    dev.off() 
  }
  
  # boxplot after normalization
  boxplotDataAN = list()
  temp = oligo::boxplot(norm, which="all", plot=FALSE)                       
  colnames(temp$stats) = temp$names
  boxplotDataAN[[1]] = temp$stats
  boxplotDataAN[[2]] = 'log-intensity'
  
  # Output data for 3D PCA #                                                                         
  tedf= t(exprs(norm))
  if (length(which(apply(tedf, 2, var)==0)) >= 0){
    tedf = tedf[ , apply(tedf, 2, var) != 0]
  }
  pca=prcomp(tedf, scale. = T)

  # similarity heatmap
  # Only output samples in selected contrast
  sampleColumns = c(which(raw@phenoData@data$groups==gsub("-.*$","",contrast)),which(raw@phenoData@data$groups==gsub("^.*-","",contrast)))
  mat=as.data.frame(as.matrix(Dist(t(exprs(norm)[,sampleColumns]),method = 'pearson',diag = TRUE)))
  mat = 1 - mat
  #sample color palette for heatmap
  x = col2hex(raw@phenoData@data$colors)[sampleColumns]
  sampleColors = x            # to return
  mat$annotation = x
  mat$Groups = raw@phenoData@data$groups[sampleColumns]
  x = unique(x)
  names(x) = unique(raw@phenoData@data$groups[sampleColumns])
  
  heatmaply(mat[,1:(ncol(mat)-2)],
            plot_method = 'plotly',
            margins = c(80,120,60,40),
            # row_side_colors = mat[,'Group', drop=FALSE],
            # row_side_palette = x,
            col_side_palette = x,
            col_side_colors = mat[,'Groups', drop=FALSE],
            colorRampPalette(colors = c("yellow", "red")),
            file = paste0(path,"/heatmapAfterRMAnorm.html"),
            symm = TRUE  )
  Heatmapolt<-"/heatmapAfterRMAnorm.html"
  
  norm_all = exprs(norm)
  norm = norm[,sampleNames(norm)[sampleColumns]]
  
  print("+++QCnorm+++")
  return (List(MAplotBN,boxplotDataBN,RLEdata,NUSEdata,HistplotAN,MAplotAN,boxplotDataAN,pca,Heatmapolt,norm,sampleColors,norm_all))
  sink(type='message')
}