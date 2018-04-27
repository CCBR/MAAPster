projectId = 'testing'
#### 1) Process GEO files function takes gseid and returns ExpressionFeatureSet object  ####
processGEOfiles <- function(id,listGroups){
  library(GEOquery)
  library(oligo)
  library(Biobase)
  id = gsub(" ","",id,fixed=TRUE) 
  system(paste0('rm *.[cC][eE][lL].gz'))                           #removes previous CEL files if run consecutively
  getGEOSuppFiles(id, makeDirectory = T, baseDir = getwd())
  fileID = paste0(id, '_RAW.tar')
  untar(paste0(getwd(),'/',id,'/',fileID))
  SampleName = list.files(pattern = '/*CEL.gz', ignore.case = T)    #list contents of new directory with zipped CEL files
  celfiles = read.celfiles(SampleName)
  gds <- getGEO(id, GSEMatrix = F,getGPL=T,AnnotGPL=T)              #get meta data 
  tableNames=c("gsm","title","description","groups")
  pData(celfiles)[tableNames] = NA
  for (k in 1:length(GSMList(gds)))                                 #fill table with meta data
  {
    if (is.null(Meta(GSMList(gds)[[k]])$description)) {    
      pData(celfiles)[k,2:4] <-c(Meta(GSMList(gds)[[k]])$geo_accession[1], Meta(GSMList(gds)[[k]])$title[1], 'No data available')
    } else {
      pData(celfiles)[k,2:4] <-c(Meta(GSMList(gds)[[k]])$geo_accession[1], Meta(GSMList(gds)[[k]])$title[1], Meta(GSMList(gds)[[k]])$description[1])
    }
  }
  pData(celfiles)$groups = listGroups                               #assign groups to samples
  ####creates a list of colors specific to each group
  fs = factor(pData(celfiles)$groups)
  lFs=levels(fs)
  numFs=length(lFs)
  colors = list()
  for (i in 1:numFs){
    colors[which(fs==lFs[i])] = i*5
  }
  colors = unlist(colors)
  pData(celfiles)$colors = colors
  ####end colors
  return(celfiles)
}
#if user selects 'ANALYZE GEO FILES', call this function, input GSE ID  (length of group assignments must match number of files for testing purposes):
celfiles = processGEOfiles('GSE61989', c('Ctl','Ctl','Ctl','Ctl','RNA_1','RNA_1','RNA_1','RNA_1','RNA_2','RNA_2','RNA_2','RNA_2'))   #human example     
celfiles = processGEOfiles('GSE37874', c('Ctl','Ctl','Ctl','Ctl','RNA_1','RNA_1','RNA_1','RNA_1','RNA_2','RNA_2','RNA_2','RNA_2'))    #mouse example    

#### 1) Process files function takes path to celfiles and returns ExpressionFeatureSet object ####
processCELfiles <- function(path,listGroups) {
  library(tools)
  library(Biobase)
  library(oligo)
  SampleName = list.files(path = path, pattern = '/*CEL*', ignore.case = T, full.names=T)
  celfiles = read.celfiles(SampleName)
  pData(celfiles)$title = basename(file_path_sans_ext(SampleName))  #add sample name to pheno
  pData(celfiles)$groups = listGroups                               #add groups to pheno
  ####creates a list of colors specific to each group
  fs = factor(pData(celfiles)$groups)
  lFs=levels(fs)
  numFs=length(lFs)
  colors = list()
  for (i in 1:numFs){
    colors[which(fs==lFs[i])] = i*5
  }
  colors = unlist(colors)
  pData(celfiles)$colors = colors
  ####end
  return(celfiles)
}
#If user selects 'ANALYZE CEL FILES', call this function, input path of files (length of group assignments must match number of files for testing purposes):
celfiles = processCELfiles('/Users/valdezkm/Documents/2___Combined',c('Ctl_1','Ctl_1','Ctl_1','KO_1','KO_1','KO_1','Ctl_2','Ctl_2','Ctl_2','KO_2','KO_2','KO_2'))

#### 2) QC / Normalize data function takes ExpressionFeatureSet from above and prints pre-normalization plots, QC plots, post-normalization plots.  Returns normalized data ExpressionFeatureSet ####
calc = function(raw) {
  library(rgl)
  library(Biobase)
  library(heatmaply)
  hist(raw,which="all", main =" Raw Samples distribution")                            #Raw histogram
  nbfacs=nrow(pData(raw))
  for (i in 1:nbfacs) {
    MAplot(raw,which=i,plotFun=smoothScatter,refSamples=c(1:nbfacs), main='', cex=2)  #Raw MAplots
  }
  boxplot(raw, which="all", main="Boxplots before normalization",
          las=2,names=pData(raw)$title, col=pData(raw)$colors)                        #Raw boxplot
  qc = fitProbeLevelModel(raw)                                                        #Calculate QC
  RLE(qc, main="RLE plot",names=pData(raw)$title, las=2, col=pData(raw)$colors)       #RLE
  NUSE(qc, main="NUSE plot",names=pData(raw)$title, las=2, col=pData(raw)$colors)     #NUSE
  #Normalize data
  if (raw@annotation=="pd.hg.u133.plus.2" | raw@annotation=="pd.clariom.s.human.ht" | raw@annotation=="pd.clariom.s.human" | raw@annotation=="pd.clariom.s.mouse.ht" | raw@annotation=="pd.clariom.s.mouse" | raw@annotation=='pd.mouse430.2' | raw@annotation=='pd.hg.u133a' | raw@annotation=='pd.hg.u133a.2' | raw@annotation=='pd.hg.u219' | raw@annotation=='pd.mg.u74av2' | raw@annotation=='pd.mouse430a.2' | raw@annotation=='pd.moe430a' | raw@annotation=='pd.hg.u95av2' | raw@annotation=='pd.hg.u133b') {
    norm =rma(raw, background=TRUE, normalize=TRUE, subset=NULL)
  } else {
    norm =rma(raw, background=TRUE, normalize=TRUE, subset=NULL, target="core")
  }
  hist(norm, main ="Distribution after Normalization")                                #Normalized histogram
  for (i in 1:nbfacs) {
    MAplot(norm,which=i,plotFun=smoothScatter,refSamples=c(1:nbfacs), main='', cex=2) #Normalized MAplots
  }
  boxplot(norm, main="Boxplots after RMA normalization",las=2,
          names=pData(raw)$title, col=pData(raw)$colors)                              #Normalized boxplot
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
  mat=as.matrix(dist(t(exprs(norm))))
  rownames(mat)=pData(norm)$title
  colnames(mat)=rownames(mat)
  print(heatmaply(mat,margins = c(80,120,60,40),colorRampPalette(colors = c("red", "yellow"))))
  return(norm)
}
norm_celfiles = calc(celfiles)                                                        #Call function

#### 3) Differentially Expressed Genes function takes files, group and contrast data. Returns list of DEGs for each contrast, annotated normalized data, and pheno data ####
# Output should dynamically respond to user-selected contrast
deg = function(norm, listContrasts) {
  library(limma)
  library(pd.mogene.2.0.st)
  library(mogene20sttranscriptcluster.db)
  library(pd.hg.u133.plus.2)
  library(hgu133plus2.db)
  library(pd.hugene.2.0.st)
  library(hugene20sttranscriptcluster.db)
  library(pd.clariom.s.human.ht)
  library(clariomshumanhttranscriptcluster.db)
  library(pd.clariom.s.human)
  library(clariomshumantranscriptcluster.db)
  library(pd.clariom.s.mouse.ht)
  library(clariomsmousehttranscriptcluster.db)
  library(pd.clariom.s.mouse)
  library(clariomsmousetranscriptcluster.db)
  library(pd.mouse430.2)
  library(mouse4302.db)
  library(pd.hg.u133a)
  library(hgu133a.db)
  library(pd.hugene.1.0.st.v1)
  library(hugene10sttranscriptcluster.db)
  library(pd.mogene.1.0.st.v1)
  library(mogene10sttranscriptcluster.db)
  library(pd.hg.u133a.2)
  library(hgu133a2.db)
  library(pd.huex.1.0.st.v2)
  library(huex10sttranscriptcluster.db)   
  library(pd.hg.u219)
  library(hgu219.db)
  library(pd.mg.u74av2)
  library(mgu74av2.db)
  library(pd.mouse430a.2)
  library(mouse430a2.db)
  library(pd.moe430a)
  library(moe430a.db)
  library(pd.hg.u95av2)
  library(hgu95av2.db)
  library(pd.hta.2.0)
  library(hta20transcriptcluster.db)
  library(pd.moex.1.0.st.v1)
  library(moex10sttranscriptcluster.db)
  library(pd.hg.u133b)
  library(hgu133b.db)
  library(pd.hugene.1.1.st.v1)
  library(hugene11sttranscriptcluster.db)
  library(pd.mogene.1.1.st.v1)
  library(mogene11sttranscriptcluster.db)
  myfactor <- factor(pData(norm)$groups)
  design1 <- model.matrix(~0+myfactor)
  colnames(design1) <- levels(myfactor)
  fit1 <- lmFit(norm,design1)
  contrast.matrix <- makeContrasts(contrasts=cons,levels=design1)
  fit2 <- contrasts.fit(fit1, contrast.matrix)
  ebayes.fit2=eBayes(fit2) 
  if (norm@annotation=="pd.mogene.2.0.st") {  
    Annot <- data.frame(ACCNUM=sapply(contents(mogene20sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mogene20sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(mogene20sttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(mogene20sttranscriptclusterENTREZID), paste, collapse=", "))     
  } else {
    # if (input$Platform=="h133p2") {
    if (norm@annotation=="pd.hg.u133.plus.2") {
      Annot <- data.frame(ACCNUM=sapply(contents(hgu133plus2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133plus2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133plus2GENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hgu133plus2ENTREZID), paste, collapse=", "))
    } else {
      if (norm@annotation=="pd.hugene.2.0.st") {
        Annot <- data.frame(ACCNUM=sapply(contents(hugene20sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hugene20sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hugene20sttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hugene20sttranscriptclusterENTREZID), paste, collapse=", "))
      } else {
        if (norm@annotation=="pd.clariom.s.human.ht") {
          Annot <- data.frame(ACCNUM=sapply(contents(clariomshumanhttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(clariomshumanhttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(clariomshumanhttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(clariomshumanhttranscriptclusterENTREZID), paste, collapse=", "))
        } else {
          if (norm@annotation=="pd.clariom.s.mouse.ht") {
            Annot <- data.frame(ACCNUM=sapply(contents(clariomsmousehttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(clariomsmousehttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(clariomsmousehttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(clariomsmousehttranscriptclusterENTREZID), paste, collapse=", "))
          } else {
            if (norm@annotation=="pd.clariom.s.mouse") {
              Annot <- data.frame(ACCNUM=sapply(contents(clariomsmousetranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(clariomsmousetranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(clariomsmousetranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(clariomsmousetranscriptclusterENTREZID), paste, collapse=", "))
            } else {
              if (norm@annotation=="pd.clariom.s.human") {
                Annot <- data.frame(ACCNUM=sapply(contents(clariomshumantranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(clariomshumantranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(clariomshumantranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(clariomshumantranscriptclusterENTREZID), paste, collapse=", "))
              } else {
                if (norm@annotation=="pd.mouse430.2") {
                  Annot <- data.frame(ACCNUM=sapply(contents(mouse4302ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mouse4302SYMBOL), paste, collapse=", "), DESC=sapply(contents(mouse4302GENENAME), paste, collapse=", "), ENTREZ=sapply(contents(mouse4302ENTREZID), paste, collapse=", "))
                } else {
                  if (norm@annotation=='pd.hg.u133a') {
                    Annot <- data.frame(ACCNUM=sapply(contents(hgu133aACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133aSYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133aGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hgu133aENTREZID), paste, collapse=", "))
                  } else {
                    if (norm@annotation=='pd.hugene.1.0.st.v1') {
                      Annot <- data.frame(ACCNUM=sapply(contents(hugene10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hugene10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hugene10sttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hugene10sttranscriptclusterENTREZID), paste, collapse=", "))
                    } else {
                      if (norm@annotation=='pd.mogene.1.0.st.v1') {
                        Annot <- data.frame(ACCNUM=sapply(contents(mogene10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mogene10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(mogene10sttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(mogene10sttranscriptclusterENTREZID), paste, collapse=", "))
                      } else {
                        if (norm@annotation=='pd.hg.u133a.2') {
                          Annot <- data.frame(ACCNUM=sapply(contents(hgu133a2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133a2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133a2GENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hgu133a2ENTREZID), paste, collapse=", "))
                        } else {
                          if (norm@annotation=='pd.huex.1.0.st.v2') {
                            Annot <- data.frame(ACCNUM=sapply(contents(huex10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(huex10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(huex10sttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(huex10sttranscriptclusterENTREZID), paste, collapse=", "))
                          } else {
                            if (norm@annotation=='pd.hg.u219') {
                              Annot <- data.frame(ACCNUM=sapply(contents(hgu219ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu219SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu219GENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hgu219ENTREZID), paste, collapse=", "))
                            } else {
                              if (norm@annotation=='pd.ht.hg.u133.plus.pm') {
                                Annot <- data.frame(ACCNUM=sapply(contents(hgu133plus2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133plus2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133plus2GENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hgu133plus2ENTREZID), paste, collapse=", "))
                              } else {
                                if (norm@annotation=='pd.mg.u74av2') {
                                  Annot <- data.frame(ACCNUM=sapply(contents(mgu74av2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mgu74av2SYMBOL), paste, collapse=", "), DESC=sapply(contents(mgu74av2GENENAME), paste, collapse=", "), ENTREZ=sapply(contents(mgu74av2ENTREZID), paste, collapse=", "))
                                } else {
                                  if (norm@annotation=='pd.mouse430a.2') {
                                    Annot <- data.frame(ACCNUM=sapply(contents(mouse430a2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mouse430a2SYMBOL), paste, collapse=", "), DESC=sapply(contents(mouse430a2GENENAME), paste, collapse=", "), ENTREZ=sapply(contents(mouse430a2ENTREZID), paste, collapse=", "))
                                  } else {
                                    if (norm@annotation=='pd.moe430a') {
                                      Annot <- data.frame(ACCNUM=sapply(contents(moe430aACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(moe430aSYMBOL), paste, collapse=", "), DESC=sapply(contents(moe430aGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(moe430aENTREZID), paste, collapse=", "))
                                    } else {
                                      if (norm@annotation=='pd.hg.u95av2') {
                                        Annot <- data.frame(ACCNUM=sapply(contents(hgu95av2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu95av2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu95av2GENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hgu95av2ENTREZID), paste, collapse=", "))
                                      } else {
                                        if (norm@annotation=='pd.hta.2.0') {
                                          Annot <- data.frame(ACCNUM=sapply(contents(hta20transcriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hta20transcriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hta20transcriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hta20transcriptclusterENTREZID), paste, collapse=", "))
                                        } else {
                                          if (norm@annotation=='pd.moex.1.0.st.v1') {
                                            Annot <- data.frame(ACCNUM=sapply(contents(moex10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(moex10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(moex10sttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(moex10sttranscriptclusterENTREZID), paste, collapse=", "))
                                          } else {
                                            if (norm@annotation=='pd.hg.u133b') {
                                              Annot <- data.frame(ACCNUM=sapply(contents(hgu133bACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133bSYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133bGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hgu133bENTREZID), paste, collapse=", "))
                                            } else {
                                              if (norm@annotation=='pd.hugene.1.1.st.v1') {
                                                Annot <- data.frame(ACCNUM=sapply(contents(hugene11sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hugene11sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hugene11sttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hugene11sttranscriptclusterENTREZID), paste, collapse=", "))
                                              } else {
                                                if (norm@annotation=='pd.mogene.1.1.st.v1') {
                                                  Annot <- data.frame(ACCNUM=sapply(contents(mogene11sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mogene11sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(mogene11sttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(mogene11sttranscriptclusterENTREZID), paste, collapse=", "))
                                                }
                                              }
                                            }
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  numContrasts = length(cons)
  listDEGs = vector("list",numContrasts)                                            #initialize output list for each contrast
  for (i in 1:numContrasts)
  {
    all.genes.con = topTable(ebayes.fit2, coef = i, number=nrow(ebayes.fit2))
    all <- merge(all.genes.con, Annot,by.x=0, by.y=0, all.x=T)                      #annotate
    all=all[order(all$P.Value),]
    colnames(all)[1]="probsetID"
    all$FC = ifelse(all$logFC<0, -1/(2^all$logFC), 2^all$logFC)                     #add fold change and rearrange columns
    all = all[,c(9,12,2,5,6,3,8,10,11,1,4,7)]
    # Write out to a file
    write.table(all,file=paste(projectId,"_",cons[i],"_all_genes.txt",sep=""),sep="\t",row.names=F)
    listDEGs[[i]]=all
  }
  names(listDEGs)=cons
  norm_annotated <- merge(exprs(norm), Annot,by.x=0, by.y=0, all.x=T)               #write out normalized annotated data
  y<-paste("_",projectId, sep="")
  tNorm = tempfile(pattern = "normalized_data_", tmpdir = getwd(), fileext = paste0(y,'.txt'))
  write.table(norm_annotated,file=tNorm,sep="\t",row.names=F)  
  for (i in 1:length(listDEGs)) {                                                   #Volcano plots
    dat=listDEGs[[i]]
    dat = dat[dat$SYMBOL!='NA',]
    log_FC=dat$logFC
    log_pval=-log10(dat$P.Value)
    Significant=rep("NotSignificant",length(log_FC))
    Significant[which(dat$P.Value<0.05 & abs(dat$logFC)>=1)]="AbsLogFoldChange>1 & PValue<0.05"
    Significant[which(dat$P.Value<0.05 & abs(dat$logFC)<1)]="PValue<0.05"
    Significant[which(dat$P.Value>=0.05 & abs(dat$logFC)>=1)]="AbsLogFoldChange>1"
    gene=dat$SYMBOL
    volcano_data=as.data.frame(cbind(gene,log_FC,log_pval,Significant))
    print(plot_ly(type='scatter', data = volcano_data, x = log_FC, y = log_pval, text = gene, mode = "markers", color = Significant) %>% layout(title=paste0('Volcano plot for: ',names(listDEGs)[i]),xaxis=list(title="Fold Change",range =c(-5,5),tickvals=c(-5,-4,-3,-2,-1,0,1,2,3,4,5),ticktext=c('-32','-16','-8','-4','-2','1','2','4','8','16','32')),yaxis=list(title="-Log10 pvalue",range =c(0,15))))
  }
  return(list(listDEGs=listDEGs, norm_annotated=norm_annotated, pheno=pData(norm)))
}
# if using processGEOfiles() function for test example, create this contrasts variable:
cons = c("RNA_1-Ctl","RNA_2-Ctl")
# or if using processCELfiles() function for test example, create this contrasts variable:
cons = c("KO_1-Ctl_1","KO_2-Ctl_2")
diff_expr_genes = deg(norm_celfiles,cons)                                           #Call function


#### 4) l2p pathway analysis function, takes DEGs and species as input, returns list of up and downregulated pathways for each contrast ####
# Output should dynamically respond to user-selected contrast
pathways = function(degs,species) {
  listPathways = vector("list",length(degs$listDEGs))
  for (i in 1:length(degs$listDEGs)) {
    up_down = vector("list",2)
    all = degs$listDEGs[[i]]
    iup=which(all$P.Value<0.05 & all$logFC>=0)
    idw=which(all$P.Value<0.05 & all$logFC<0)
    fin.up=all[iup,]
    if (length(iup) > 500)
    {
      fin.up=fin.up[order(fin.up$P.Value),]
      fin.up=fin.up[1:500,]
    }
    fin.dw=all[idw,]
    if (length(idw) > 500)
    {
      fin.dw=fin.dw[order(fin.dw$P.Value),]
      fin.dw=fin.dw[1:500,]
    }
    fin.up$SYMBOL = as.character(fin.up$SYMBOL)
    fin.dw$SYMBOL = as.character(fin.dw$SYMBOL)
    if (species=='human') 
    {
      cat(fin.up$SYMBOL,file=(paste0(projectId,'_',names(degs$listDEGs[i]),'_Top500_Up.txt')), sep='\n')
      cat(fin.dw$SYMBOL,file=(paste0(projectId,'_',names(degs$listDEGs[i]),'_Top500_Down.txt')),sep='\n')
    }
    else
    {
      cat(fin.up$SYMBOL,file=paste0(names(degs$listDEGs[i]),"_Top500temp_Up.txt"),sep='\n')
      cat(fin.dw$SYMBOL,file=paste0(names(degs$listDEGs[i]),"_Top500temp_Dw.txt"),sep='\n')
      
      system(paste0("cat ",names(degs$listDEGs[i]),"_Top500temp_Up.txt | grep -v \"^NA\" | ./m2h | grep -v XXXX | cut -f2 -d\" \" >",projectId,'_',names(degs$listDEGs[i]),"_Top500_Up.txt"))
      system(paste0("cat ",names(degs$listDEGs[i]),"_Top500temp_Dw.txt | grep -v \"^NA\" | ./m2h | grep -v XXXX | cut -f2 -d\" \" >",projectId,'_',names(degs$listDEGs[i]),"_Top500_Down.txt"))
    }
    system(paste0("cat ",projectId,'_',names(degs$listDEGs[i]),"_Top500_Up.txt |sort | uniq | ./l2p >",projectId,'_',names(degs$listDEGs[i]),"_Pathways_Up.txt"))
    system(paste0("cat ",projectId,'_',names(degs$listDEGs[i]),"_Top500_Down.txt |sort | uniq | ./l2p >",projectId,'_',names(degs$listDEGs[i]),"_Pathways_Down.txt"))
    
    addUpCol = read.delim(paste0(projectId,'_',names(degs$listDEGs[i]),"_Pathways_Up.txt"), sep = '\t')
    addDwCol = read.delim(paste0(projectId,'_',names(degs$listDEGs[i]),"_Pathways_Down.txt"), sep = '\t')
    
    colnames(addUpCol)=c("pval","fdr","ratio","nb.hits","nb.genes.path","nb.user.genes","tot.back.genes","path_id","source","description","type","gene.list")
    colnames(addDwCol)=c("pval","fdr","ratio","nb.hits","nb.genes.path","nb.user.genes","tot.back.genes","path_id","source","description","type","gene.list")
    addUpCol = addUpCol[order(addUpCol$pval),]
    addDwCol = addDwCol[order(addDwCol$pval),]
    addUpCol = addUpCol[,c(8,9,10,11,1,2,3,12,4,5,6,7)]
    addDwCol = addDwCol[,c(8,9,10,11,1,2,3,12,4,5,6,7)]
    write.table(addUpCol, file = paste0(projectId,'_',names(degs$listDEGs[i]),"_Pathways_Up.txt"), sep = '\t', row.names = F)
    write.table(addDwCol, file = paste0(projectId,'_',names(degs$listDEGs[i]),"_Pathways_Down.txt"), sep = '\t', row.names = F)
    up_down[[1]]=addUpCol
    up_down[[2]]=addDwCol
    names(up_down) = c("upregulated_pathways","downregulated_pathways")
    listPathways[[i]] = up_down
  }
  names(listPathways) = names(degs$listDEGs)  
  return(listPathways)
}
l2p_pathways = pathways(diff_expr_genes,'human')


#### 5) Function for gene heatmap from l2p pathways.  Input: deg function output, l2p pathways, contrast, choice of up or down pathways, and name of pathway. ####
#Output should change dynamically with user-selected contrast and l2p pathway (from either up OR downregulated pathways) 
geneHeatmap = function(degs, paths, contrast, upOrDown, pathway_name) {
  library(pheatmap)
  human2mouse = read.delim('human2mouse.csv',sep=',')
  paths = paths[[contrast]][[upOrDown]]
  genes = paths$gene.list[paths$description==pathway_name]              #select user input pathway, extract genes
  genes = strsplit(as.character(genes),' ')
  genes = unlist(genes)
  exp = degs$norm_annotated                                             #extract normalized expression, subset by genes, aggregate duplicate symbols by mean
  exp = exp[exp$SYMBOL %in% genes,]
  if (nrow(exp)==0) {
    genes = human2mouse$mouse[human2mouse$human %in% genes]
    genes = as.character(genes)
    exp = degs$norm_annotated                                            
    exp = exp[exp$SYMBOL %in% genes,]
  }
  exp = subset(exp, select = -c(ACCNUM,DESC,ENTREZ,Row.names))
  exp = aggregate(.~SYMBOL,data=exp,mean)
  sampleColumns = c(which(degs$pheno$groups==gsub("-.*$","",contrast)),which(degs$pheno$groups==gsub("^.*-","",contrast)))
  rownames(exp) = exp$SYMBOL
  exp = subset(exp, select = -c(SYMBOL))
  exp = exp[,sampleColumns]
  if(nrow(exp)>100){                                                    #limit to 100 genes
    exp = exp[1:100,]
  }
  matCol = data.frame(group=degs$pheno$groups[sampleColumns])           #set heatmap parameters
  rownames(matCol) = rownames(degs$pheno)[sampleColumns]
  matColors = list(group = unique(degs$pheno$colors[sampleColumns]))
  names(matColors$group) = unique(degs$pheno$groups[sampleColumns])
  path_name = pathway_name
  exp = t(scale(t(exp)))                                                #get z-scores by row
  if (nrow(exp) > 30){
    pheatmap(exp, main=path_name, annotation_col=matCol, annotation_colors=matColors, drop_levels=TRUE, fontsize_row = 6)
  } else {
    pheatmap(exp, main=path_name, annotation_col=matCol, annotation_colors=matColors, drop_levels=TRUE, fontsize_row = 10)
  }
}
geneHeatmap(diff_expr_genes, l2p_pathways, 'RNA_1-Ctl', 'upregulated_pathways','oxidation-reduction process')   #if GEO
geneHeatmap(diff_expr_genes, l2p_pathways, 'KO_1-Ctl_1', 'upregulated_pathways','oxidation-reduction process')   #if CEL file upload


#### 6) ssGSEA function, takes as input: output from deg function, species, and gene set modules(.gmt). Outputs one table of enrichment scores and tables of diff expr pathways per contrast. Prints ssGSEA heatmap ####
# Output should dynamically respond to user-selected contrast
ss = function(deg_normAnnot, species, geneSet) {
  library(GSEABase)
  library(GSVA)
  normAnnot = deg_normAnnot$norm_annotated
  ssgs = normAnnot[normAnnot$SYMBOL!='NA',]
  #if human or mouse, prepare data for gsva
  if (species=='human') {
    ssgs = subset(ssgs, select=-c(ACCNUM,DESC,Row.names,ENTREZ))
    ssgs = aggregate(.~SYMBOL,data=ssgs,mean)                               #aggregate duplicate probes by mean
    rownames(ssgs) = ssgs$SYMBOL
    ssgs = subset(ssgs, select=-c(SYMBOL))
    ssgs = as.matrix(ssgs)
    getSet = switch(geneSet, "H: Hallmark Gene Sets"="h.all.v6.1.symbols.gmt", "C1: Positional Gene Sets"="c1.all.v6.1.symbols.gmt", "C2: Curated Gene Sets"="c2.all.v6.1.symbols.gmt", 
                    "C3: Motif Gene Sets"="c3.all.v6.1.symbols.gmt", "C4: Computational Gene Sets"="c4.all.v6.1.symbols.gmt","C5: GO gene sets"="c5.all.v6.1.symbols.gmt", 
                    "C6: Oncogenic Signatures"="c6.all.v6.1.symbols.gmt", "C7: Immunologic Signatures"="c7.all.v6.1.symbols.gmt")
  } else {
    ssgs = subset(ssgs, select=-c(ACCNUM,DESC,Row.names,SYMBOL))
    ssgs = aggregate(.~ENTREZ,data=ssgs,mean)                               #aggregate duplicate probes by mean
    rownames(ssgs) = ssgs$ENTREZ
    ssgs = subset(ssgs, select=-c(ENTREZ))
    ssgs = as.matrix(ssgs)
    getSet = switch(geneSet, "H: Hallmark Gene Sets"="mouse_H_v5p2.gmt", "C2: Curated Gene Sets"="mouse_C2_v5p2.gmt", "C3: Motif Gene Sets"="mouse_C3_v5p2.gmt", "C4: Computational Gene Sets"="mouse_C4_v5p2.gmt",
                    "C5: GO gene sets"="mouse_C5_v5p2.gmt", "C6: Oncogenic Signatures"="mouse_C6_v5p2.gmt", "C7: Immunologic Signatures"="mouse_C7_v5p2.gmt")
  }
  gset = getGmt(getSet) 
  ssgsResults = gsva(ssgs, gset, method='ssgsea')                           #run ssGSEA
  y<-paste("_",projectId, sep="")                                           #write out results
  tSS = tempfile(pattern = "ssGSEA_enrichmentScores_", tmpdir = getwd(), fileext = paste0(y,'.txt'))
  write.table(ssgsResults,file=tSS,sep="\t",col.names=NA)
  myfactor <- factor(deg_normAnnot$pheno$groups)
  design1 <- model.matrix(~0+myfactor)
  colnames(design1) <- levels(myfactor)
  fit1 = lmFit(ssgsResults,design1)                                                                                                      #DE analysis of ssGSEA enrichment scores
  cons = names(deg_normAnnot$listDEGs)
  contrast.matrix = makeContrasts(contrasts=cons,levels=design1)
  fit2 = contrasts.fit(fit1,contrast.matrix)
  ebayes.fit2 = eBayes(fit2)
  DEss=vector("list",length(deg_normAnnot$listDEGs))                        
  for (i in 1:length(deg_normAnnot$listDEGs))                           
  {
    all.pathways = topTable(ebayes.fit2, coef=i, number=nrow(ebayes.fit2))                                                               #Determine DE pathways
    all.pathways = all.pathways[order(abs(all.pathways$P.Value)),]
    colnames(all.pathways)[2] = 'Avg.Enrichment.Score'
    write.table(all.pathways,file=paste0(projectId,"_",cons[i],"_ssGSEA_pathways.txt"),sep="\t",row.names=T,col.names=NA)
    DEss[[i]] = all.pathways
  }
  names(DEss)=cons
  for (i in 1:length(DEss)){                                                                                                             #Heatmap
    sampleColumns = c(which(deg_normAnnot$pheno$groups==gsub("-.*$","",cons[i])),which(deg_normAnnot$pheno$groups==gsub("^.*-","",cons[i])))   #Subset columns (samples)  
    paths = ssgsResults[rownames(ssgsResults) %in% rownames(DEss[[i]])[1:50],]                                                           #Subset rows (pathways)
    paths = paths[,sampleColumns]
    matCol = data.frame(group=deg_normAnnot$pheno$groups[sampleColumns])
    rownames(matCol) = rownames(deg_normAnnot$pheno)[sampleColumns]
    matColors = list(group = unique(deg_normAnnot$pheno$colors[sampleColumns]))
    names(matColors$group) = unique(deg_normAnnot$pheno$groups[sampleColumns])
    paths = t(scale(t(paths))) 
    pheatmap(paths,annotation_col=matCol,annotation_colors=matColors,drop_levels=TRUE,fontsize=7, main='Enrichment Scores for Top 50 Differentially Expressed ssGSEA Pathways')
  }
  return(list(ssgsResults=ssgsResults, DEss=DEss))
}
ssGSEA_results = ss(diff_expr_genes,'human','C2: Curated Gene Sets')






