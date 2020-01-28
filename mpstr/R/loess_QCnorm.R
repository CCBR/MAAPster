#' 2) Runs QC analyses and normalizes data.
#' 
#' @export
#' @param raw ExpressionFeatureSet from either processGEOfiles or processCELfiles function
#' @param path Directory for plots
#' @param contrast Contrast for subsetting QC plots by samples in chosen groups
#' @param listBatches Optional list of batches for each sample, follow alphabetical order of samples
#' @return Normalized ExpressionFeatureSet object with normalized data and phenotype information
#' @examples 
#' norm_celfiles = loess_QCnorm(celfiles,'Users/name/folderName/plots',c("RNA_1-Ctl")) 
#' @references See packages rgl, Biobase, heatmaply, oligo
#' @note Normalizes using rma and cyclic loess, see oligo  and limma packages
#' @note Outputs 3D PCA and similarity heatmap
#' @note Outputs pre-normalization plots, QC plots, post-normalization plots, see oligo package for details


loess_QCnorm = function(raw,path,contrast,listBatches=NULL) {
  library(rgl)
  library(Biobase)
  library(heatmaply)
  library(reshape2)
  library(plotly)
  library(reshape2)
  library(amap)
  library(gplots)
  library(limma)
  library(oligo)
  library(sva)
  
  QCnorm_ERR = file(paste0(path,'/loess_QCnorm.err'),open='wt')
  sink(QCnorm_ERR,type='message',append=TRUE)
  
  # subset to samples in contrast
  sampleColumns = c(which(raw@phenoData@data$groups==gsub("-.*$","",contrast)),which(raw@phenoData@data$groups==gsub("^.*-","",contrast)))
  raw = raw[,sampleNames(raw)[sampleColumns]]
  
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
  
  loess = normalizeCyclicLoess(norm, weights = NULL,  span=0.7, iterations = 3, method = "fast")
  pheno = new("AnnotatedDataFrame",data=norm@phenoData@data)
  loess = ExpressionSet(assayData = loess,phenoData = pheno)
  loess@annotation = norm@annotation
  norm = loess

  
  # add optional batch correction for plots
  if (is.null(listBatches)) {
    norm_plots = norm
  } else {
    # set up variables
    temp = norm
    modCombat = model.matrix(~1,data=pData(temp))
    batch = pData(temp)$batch
    # run combat
    norm_plots = ComBat(dat=exprs(temp),batch=batch,mod=modCombat,par.prior=TRUE,prior.plots=FALSE)
    # turn back into expression set
    pheno = new("AnnotatedDataFrame",data=temp@phenoData@data)
    norm_plots = ExpressionSet(assayData=norm_plots,phenoData=pheno)
    norm_plots@annotation = temp@annotation
  }
  
  
  # histogram after normalization
  svg(paste0(path,"/temp2.svg"),width=8, height=8)
  dat = hist(norm_plots)
  dev.off()
  
  dat_x = melt(dat[['x']])
  dat_y = melt(dat[['y']])
  dat = data.frame(x=dat_x$value,y=dat_y$value,sample=dat_x$Var2)
  dat$sample = gsub('\n','',dat$sample)
  
  HistplotAN = plot_ly(dat, x=~x, y=~y, color=~sample, text = ~sample) %>%
    add_lines() %>%
    layout(title='Distribution after RMA and Loess Normalization',xaxis=list(title='log-intensity'),yaxis=list(title='density'),legend=list(x=13,y=0.1))
  htmlwidgets::saveWidget(HistplotAN, paste0(path,"/histAfterLoessNorm.html"))
  
  # MAplot after normalization
  MAplotAN<-List()
  for (i in 1:nbfacs) {
    MAplotAN<-c(MAplotAN,paste0("/MAplotsAfterLoessNorm",i,".jpg"))
    jpeg(paste0(path,"/MAplotsAfterLoessNorm",i,".jpg"),width=5, height=5,units = "in", res = 300)
    MAplot(norm_plots,which=i,plotFun=smoothScatter,refSamples=c(1:nbfacs), main='', cex=2) #Normalized MAplots
    dev.off() 
  }
  
  # boxplot after normalization
  boxplotDataAN = list()
  temp = oligo::boxplot(norm_plots, which="all", plot=FALSE)                       
  colnames(temp$stats) = temp$names
  boxplotDataAN[[1]] = temp$stats
  boxplotDataAN[[2]] = 'log-intensity'
  
  # Output data for 3D PCA #                                                                         
  tedf= t(exprs(norm_plots))
  if (length(which(apply(tedf, 2, var)==0)) >= 0){
    tedf = tedf[ , apply(tedf, 2, var) != 0]
  }
  pca=prcomp(tedf, scale. = T)
  
  # similarity heatmap
  # Only output samples in selected contrast
  mat=as.data.frame(as.matrix(Dist(t(exprs(norm_plots)),method = 'pearson',diag = TRUE)))
  mat = 1 - mat
  #sample color palette for heatmap
  x = col2hex(raw@phenoData@data$colors)
  sampleColors = x            # to return
  mat$annotation = x
  mat$Groups = raw@phenoData@data$groups
  x = unique(x)
  names(x) = unique(raw@phenoData@data$groups)
  
 
  heatmaply(mat[,1:(ncol(mat)-2)],
            plot_method = 'plotly',
            margins = c(80,120,60,40),
            # row_side_colors = mat[,'Group', drop=FALSE],
            # row_side_palette = x,
            col_side_palette = x,
            col_side_colors = mat[,'Groups', drop=FALSE],
            colorRampPalette(colors = c("yellow", "red")),
            file = paste0(path,"/heatmapAfterLoessNorm.html"),
            symm = TRUE  )
  Heatmapolt<-"/heatmapAfterLoessNorm.html"
  
  norm_all = exprs(norm)
  # norm = norm[,sampleNames(norm)[sampleColumns]]
  
  print("+++QCnorm+++")
  return (List(MAplotBN,boxplotDataBN,RLEdata,NUSEdata,HistplotAN,MAplotAN,boxplotDataAN,pca,Heatmapolt,norm,sampleColors,norm_all,norm_plots))
  sink(type='message')
}