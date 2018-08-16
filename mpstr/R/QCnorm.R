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
  
  HistplotBN<-"/histBeforeNorm.svg"
  svg(paste0(path,"/histBeforeNorm.svg"),width=8, height=8)
  hist(raw,which="all", main =" Raw Samples distribution")          
  dev.off()  
  #Raw histogram
  nbfacs=nrow(pData(raw))
  MAplotBN<-List()
  for (i in 1:nbfacs) {
    MAplotBN<-c(MAplotBN,paste0("/MAplotsBeforeNorm",i,".jpg"))
    jpeg(paste0(path,"/MAplotsBeforeNorm",i,".jpg"),width=5, height=5,units = "in", res = 300)
    MAplot(raw,which=i,plotFun=smoothScatter,refSamples=c(1:nbfacs), main='', cex=2)  #Raw MAplots
    dev.off() 
  }
  
  BoxplotBN<-"/boxplotBeforeNorm.svg"
  svg(paste0(path,"/boxplotBeforeNorm.svg"),width=6, height=6) 
  boxplot(raw, which="all", main="Boxplots before normalization",
          las=2,names=pData(raw)$title, col=pData(raw)$colors)                        #Raw boxplot
  dev.off()
  
  
  qc = fitProbeLevelModel(raw)                                                        #Calculate QC
  RLEplotBN<-"/RLEBeforeNorm.svg"
  svg(paste0(path,"/RLEBeforeNorm.svg"),width=6, height=6) 
  RLE(qc, main="RLE plot",names=pData(raw)$title, las=2, col=pData(raw)$colors)       #RLE
  dev.off() 
  
  NUSEplotBN<-"/NUSEBeforeNorm.svg"
  svg(paste0(path,"/NUSEBeforeNorm.svg"),width=6, height=6) 
  NUSE(qc, main="NUSE plot",names=pData(raw)$title, las=2, col=pData(raw)$colors)     #NUSE
  dev.off()
  
  #Normalize data
  if (raw@annotation=="pd.hg.u133.plus.2" | raw@annotation=="pd.clariom.s.human.ht" | raw@annotation=="pd.clariom.s.human" | raw@annotation=="pd.clariom.s.mouse.ht" | raw@annotation=="pd.clariom.s.mouse" | raw@annotation=='pd.mouse430.2' | raw@annotation=='pd.hg.u133a' | raw@annotation=='pd.hg.u133a.2' | raw@annotation=='pd.hg.u219' | raw@annotation=='pd.mg.u74av2' | raw@annotation=='pd.mouse430a.2' | raw@annotation=='pd.moe430a' | raw@annotation=='pd.hg.u95av2' | raw@annotation=='pd.hg.u133b') {
    norm =rma(raw, background=TRUE, normalize=TRUE, subset=NULL)
  } else {
    norm =rma(raw, background=TRUE, normalize=TRUE, subset=NULL, target="core")
  }
  
  HistplotAN<-"/histAfterNorm.svg"
  svg(paste0(path,"/histAfterNorm.svg"),width=6, height=6)
  hist(norm, main ="Distribution after Normalization")                             #Normalized histogram
  dev.off()
  
  MAplotAN<-List()
  for (i in 1:nbfacs) {
    MAplotAN<-c(MAplotAN,paste0("/MAplotsAfterNorm",i,".jpg"))
    jpeg(paste0(path,"/MAplotsAfterNorm",i,".jpg"),width=5, height=5,units = "in", res = 300)
    MAplot(norm,which=i,plotFun=smoothScatter,refSamples=c(1:nbfacs), main='', cex=2) #Normalized MAplots
    dev.off() 
  }
  
  BoxplotAN<-"/boxplotsAfterNorm.svg"
  svg(paste0(path,"/boxplotsAfterNorm.svg"),width=6, height=6)
  boxplot(norm, main="Boxplots after RMA normalization",las=2,
          names=pData(raw)$title, col=pData(raw)$colors)                              #Normalized boxplot
  dev.off() 
  
  # 3D PCA #                                                                          #3D PCA
  tedf= t(exprs(norm))
  if (length(which(apply(tedf, 2, var)==0)) >= 0){
    tedf = tedf[ , apply(tedf, 2, var) != 0]
  }
  pca=prcomp(tedf, scale. = T)
  open3d()
  bg3d('white')
  plot3d(pca$x[,1:3], type='s',size=2, col=pData(raw)$colors)
  group.v=as.vector(pData(raw)$title)
  text3d(pca$x, pca$y, pca$z, group.v, cex=0.6, adj=1.5)
  par3d(mouseMode = "trackball")
  # END 3D PCA / BEGIN HEATMAP #                                                      #Heatmap
  writeWebGL(dir = file.path(path, "webGL"),width = 650, reuse = TRUE)
  
  PCA<-"/webGL/index.html"
  
  mat=as.matrix(dist(t(exprs(norm))))
  rownames(mat)=pData(norm)$title
  colnames(mat)=rownames(mat)
  
  heatmaply(
    mat,margins = c(80,120,60,40),
    colorRampPalette(colors = c("red", "yellow")),
    file = paste0(path,"/heatmapAfterNorm.html")
  )
  Heatmapolt<-"/heatmapAfterNorm.html"
  print("+++cal+++")
  return (List(HistplotBN,MAplotBN,BoxplotBN,RLEplotBN,NUSEplotBN,HistplotAN,MAplotAN,BoxplotAN,PCA,Heatmapolt,norm))
}