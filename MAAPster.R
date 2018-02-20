library(shiny)
library(shinyjs)
library(GEOquery)
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
#library(primeviewprobe)
#library(GSEA)
library(limma)
library(oligo)
library(gplots)
library(geneplotter)
library(multtest)
library(rgl)
library(rglwidget)
library(DT)
library(getopt)
library(annotate)
library(knitr)
library(reshape)
library(RColorBrewer)
library(mixOmics)
library(calibrate)
library(rmarkdown)
library(ggplot2)
library(ggfortify)
library(shinyRGL)
library(plotly)
library(htmltools)
library(heatmaply)
library(Biobase)
library(GSVA)
library(GSEABase)
library(pheatmap)
library(viridis)
library(dendsort)

setwd('/Users/valdezkm/Documents/MicroarrayPipeline/CodeInProgress/MicroArrayPipeline')

#user input: Project ID (cannot be NA)
projectId = 'test'

#user input: GSE number
id = 'GSE44392'
id = gsub(" ","",id,fixed=TRUE)  

#gets meta data 
gds <- getGEO(id, GSEMatrix = F,getGPL=T,AnnotGPL=T)

#creates empty table
mytable=matrix("",length(GSMList(gds)),3)
colnames(mytable)=c("gsm","title","description")

#Populates table with meta data
# (For CEL file input, this table will have one column for CEL file name, 
# an option would be to allow the user to select unique ids or upload a file assigning the samples to ids)
for (k in 1:length(GSMList(gds)))
{
  if (is.null(Meta(GSMList(gds)[[k]])$description)) {    
    mytable[k,] <-c(Meta(GSMList(gds)[[k]])$geo_accession[1], Meta(GSMList(gds)[[k]])$title[1], 'No data available')
  } else {
    mytable[k,] <-c(Meta(GSMList(gds)[[k]])$geo_accession[1], Meta(GSMList(gds)[[k]])$title[1], Meta(GSMList(gds)[[k]])$description[1])
  }
}

#example of a way to process CEL files once they're uploaded:
# cels = upload_cel_files_here
# for (k in 1:length(cels))
# {
#   mytable[k,]<-c(cels[k])
# }
# mytable <- data.frame(mytableCEL)


# a group column will be populated by the user by clicking on the samples and assigning them to groups
mytable = data.frame(mytable)
mytable$group = '...'

#user input: group names and indices chosen by user (min 2, no max)
mytable$group[1:4] = make.names('Control')
mytable$group[13:16] = make.names('Stim')
mytable$group[17:20] = make.names('Stim+3.3')


system(paste0('rm *.[cC][eE][lL].gz'))        #removes previous CEL files
getGEOSuppFiles(id, makeDirectory = T, baseDir = getwd())
fileID = paste0(id, '_RAW.tar')
untar(paste0(getwd(),'/',id,'/',fileID))

SampleName = list.files(pattern = '/*CEL.gz', ignore.case = T)    #list contents of new directory with zipped CEL files
if (length(grep('*CEL*',SampleName,ignore.case = T)) == 0) {
  info("Raw files must be CEL files")
}
rownames(mytable) = mytable$title
cels = SampleName

pd = AnnotatedDataFrame(mytable)
celfiles = read.celfiles(cels, phenoData = pd)
colnames(pData(celfiles))[2] = 'SampleID' 

y<-paste("_",projectId, sep="")
tAnnot = tempfile(pattern = "annotation_", tmpdir = getwd(), fileext = paste0(y,'.txt'))
cat(celfiles@annotation,file=tAnnot)    

#if platform is not supported display error message
if (celfiles@annotation!="pd.hg.u133.plus.2" & celfiles@annotation!="pd.mogene.2.0.st" & celfiles@annotation!="pd.hugene.2.0.st" & celfiles@annotation!="pd.clariom.s.human.ht" & celfiles@annotation!="pd.clariom.s.human" & celfiles@annotation!="pd.clariom.s.mouse.ht" & celfiles@annotation!="pd.clariom.s.mouse" & celfiles@annotation!='pd.mouse430.2' & celfiles@annotation!='pd.hg.u133a' & celfiles@annotation!='pd.hugene.1.0.st.v1' & celfiles@annotation!='pd.mogene.1.0.st.v1' & celfiles@annotation!='pd.hg.u133a.2' & celfiles@annotation!='pd.huex.1.0.st.v2' & celfiles@annotation!='pd.hg.u219' & celfiles@annotation!='pd.mg.u74av2' & celfiles@annotation!='pd.mouse430a.2' & celfiles@annotation!='pd.moe430a' & celfiles@annotation!='pd.hg.u95av2' & celfiles@annotation!='pd.hta.2.0' & celfiles@annotation!='pd.moex.1.0.st.v1' & celfiles@annotation!='pd.hg.u133b' & celfiles@annotation!='pd.hugene.1.1.st.v1' & celfiles@annotation!='pd.mogene.1.1.st.v1') {
  #cat("Please sort your phenotype on sample name and upload it again. \n")
  info(paste0("Affymetrix platform: ",celfiles@annotation," NOT supported. Leaving..."))
  stopApp(-1)
}

#normalization
if (celfiles@annotation=="pd.hg.u133.plus.2" | celfiles@annotation=="pd.clariom.s.human.ht" | celfiles@annotation=="pd.clariom.s.human" | celfiles@annotation=="pd.clariom.s.mouse.ht" | celfiles@annotation=="pd.clariom.s.mouse" | celfiles@annotation=='pd.mouse430.2' | celfiles@annotation=='pd.hg.u133a' | celfiles@annotation=='pd.hg.u133a.2' | celfiles@annotation=='pd.hg.u219' | celfiles@annotation=='pd.mg.u74av2' | celfiles@annotation=='pd.mouse430a.2' | celfiles@annotation=='pd.moe430a' | celfiles@annotation=='pd.hg.u95av2' | celfiles@annotation=='pd.hg.u133b') {
  celfiles.rma =rma(celfiles, background=TRUE, normalize=TRUE, subset=NULL)
} else {
  celfiles.rma =rma(celfiles, background=TRUE, normalize=TRUE, subset=NULL, target="core")
}

#QC 
celfiles.qc=fitProbeLevelModel(celfiles)


#Contrast (user input, many contrasts may be added as rows, contrasts must be chosen from groups assigned above)
contra = data.frame(k1 = unique(mytable$group[1:4]), k2 = unique(mytable$group[13:16]))
contra = rbind(contra, data.frame(k1 = unique(mytable$group[1:4]), k2 = unique(mytable$group[17:20])))

#####  DEG  #####
nb=dim(contra)[1]
cons=c()

#order: experimental samples first, control/baseline samples second
for (k in 1:nb) {
  cons=c(cons,paste(contra[k,2],"-",contra[k,1],sep=""))
}

myfactor <- factor(celfiles$group)
design1 <- model.matrix(~0+myfactor)
colnames(design1) <- levels(myfactor)

fit1 <- lmFit(celfiles.rma,design1)
contrast.matrix <- makeContrasts(contrasts=cons,levels=design1)

fit2 <- contrasts.fit(fit1, contrast.matrix)
ebayes.fit2=eBayes(fit2) # smooths the std error

if (celfiles@annotation=="pd.mogene.2.0.st") {  
  Annot <- data.frame(ACCNUM=sapply(contents(mogene20sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mogene20sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(mogene20sttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(mogene20sttranscriptclusterENTREZID), paste, collapse=", "))     
  #Annot <- data.frame(ACCNUM=sapply(contents(mogene20sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mogene20sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(mogene20sttranscriptclusterGENENAME), paste, collapse=", "))   
} else {
  # if (input$Platform=="h133p2") {
  if (celfiles@annotation=="pd.hg.u133.plus.2") {
    #Annot <- data.frame(ACCNUM=sapply(contents(hgu133plus2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133plus2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133plus2GENENAME), paste, collapse=", "))
    Annot <- data.frame(ACCNUM=sapply(contents(hgu133plus2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133plus2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133plus2GENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hgu133plus2ENTREZID), paste, collapse=", "))
  } else {
    if (celfiles@annotation=="pd.hugene.2.0.st") {
      #Annot <- data.frame(ACCNUM=sapply(contents(hugene20sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hugene20sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hugene20sttranscriptclusterGENENAME), paste, collapse=", "))
      Annot <- data.frame(ACCNUM=sapply(contents(hugene20sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hugene20sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hugene20sttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hugene20sttranscriptclusterENTREZID), paste, collapse=", "))
    } else {
      if (celfiles@annotation=="pd.clariom.s.human.ht") {
        #Annot <- data.frame(ACCNUM=sapply(contents(clariomshumanhttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(clariomshumanhttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(clariomshumanhttranscriptclusterGENENAME), paste, collapse=", "))
        Annot <- data.frame(ACCNUM=sapply(contents(clariomshumanhttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(clariomshumanhttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(clariomshumanhttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(clariomshumanhttranscriptclusterENTREZID), paste, collapse=", "))
      } else {
        if (celfiles@annotation=="pd.clariom.s.mouse.ht") {
          #Annot <- data.frame(ACCNUM=sapply(contents(clariomsmousehttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(clariomsmousehttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(clariomsmousehttranscriptclusterGENENAME), paste, collapse=", "))
          Annot <- data.frame(ACCNUM=sapply(contents(clariomsmousehttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(clariomsmousehttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(clariomsmousehttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(clariomsmousehttranscriptclusterENTREZID), paste, collapse=", "))
        } else {
          if (celfiles@annotation=="pd.clariom.s.mouse") {
            #Annot <- data.frame(ACCNUM=sapply(contents(clariomsmousetranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(clariomsmousetranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(clariomsmousetranscriptclusterGENENAME), paste, collapse=", "))
            Annot <- data.frame(ACCNUM=sapply(contents(clariomsmousetranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(clariomsmousetranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(clariomsmousetranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(clariomsmousetranscriptclusterENTREZID), paste, collapse=", "))
          } else {
            if (celfiles@annotation=="pd.clariom.s.human") {
              #Annot <- data.frame(ACCNUM=sapply(contents(clariomshumantranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(clariomshumantranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(clariomshumantranscriptclusterGENENAME), paste, collapse=", "))
              Annot <- data.frame(ACCNUM=sapply(contents(clariomshumantranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(clariomshumantranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(clariomshumantranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(clariomshumantranscriptclusterENTREZID), paste, collapse=", "))
            } else {
              if (celfiles@annotation=="pd.mouse430.2") {
                #Annot <- data.frame(ACCNUM=sapply(contents(mouse4302ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mouse4302SYMBOL), paste, collapse=", "), DESC=sapply(contents(mouse4302GENENAME), paste, collapse=", "))
                Annot <- data.frame(ACCNUM=sapply(contents(mouse4302ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mouse4302SYMBOL), paste, collapse=", "), DESC=sapply(contents(mouse4302GENENAME), paste, collapse=", "), ENTREZ=sapply(contents(mouse4302ENTREZID), paste, collapse=", "))
              } else {
                if (celfiles@annotation=='pd.hg.u133a') {
                  #Annot <- data.frame(ACCNUM=sapply(contents(hgu133aACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133aSYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133aGENENAME), paste, collapse=", "))
                  Annot <- data.frame(ACCNUM=sapply(contents(hgu133aACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133aSYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133aGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hgu133aENTREZID), paste, collapse=", "))
                } else {
                  if (celfiles@annotation=='pd.hugene.1.0.st.v1') {
                    #Annot <- data.frame(ACCNUM=sapply(contents(hugene10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hugene10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hugene10sttranscriptclusterGENENAME), paste, collapse=", "))
                    Annot <- data.frame(ACCNUM=sapply(contents(hugene10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hugene10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hugene10sttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hugene10sttranscriptclusterENTREZID), paste, collapse=", "))
                  } else {
                    if (celfiles@annotation=='pd.mogene.1.0.st.v1') {
                      #Annot <- data.frame(ACCNUM=sapply(contents(mogene10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mogene10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(mogene10sttranscriptclusterGENENAME), paste, collapse=", "))
                      Annot <- data.frame(ACCNUM=sapply(contents(mogene10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mogene10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(mogene10sttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(mogene10sttranscriptclusterENTREZID), paste, collapse=", "))
                    } else {
                      if (celfiles@annotation=='pd.hg.u133a.2') {
                        #Annot <- data.frame(ACCNUM=sapply(contents(hgu133a2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133a2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133a2GENENAME), paste, collapse=", "))
                        Annot <- data.frame(ACCNUM=sapply(contents(hgu133a2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133a2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133a2GENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hgu133a2ENTREZID), paste, collapse=", "))
                      } else {
                        if (celfiles@annotation=='pd.huex.1.0.st.v2') {
                          #Annot <- data.frame(ACCNUM=sapply(contents(huex10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(huex10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(huex10sttranscriptclusterGENENAME), paste, collapse=", "))
                          Annot <- data.frame(ACCNUM=sapply(contents(huex10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(huex10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(huex10sttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(huex10sttranscriptclusterENTREZID), paste, collapse=", "))
                        } else {
                          if (celfiles@annotation=='pd.hg.u219') {
                            #Annot <- data.frame(ACCNUM=sapply(contents(hgu219ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu219SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu219GENENAME), paste, collapse=", "))
                            Annot <- data.frame(ACCNUM=sapply(contents(hgu219ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu219SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu219GENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hgu219ENTREZID), paste, collapse=", "))
                          } else {
                            if (celfiles@annotation=='pd.ht.hg.u133.plus.pm') {
                              #Annot <- data.frame(ACCNUM=sapply(contents(hgu133plus2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133plus2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133plus2GENENAME), paste, collapse=", "))
                              Annot <- data.frame(ACCNUM=sapply(contents(hgu133plus2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133plus2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133plus2GENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hgu133plus2ENTREZID), paste, collapse=", "))
                            } else {
                              if (celfiles@annotation=='pd.mg.u74av2') {
                                #Annot <- data.frame(ACCNUM=sapply(contents(mgu74av2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mgu74av2SYMBOL), paste, collapse=", "), DESC=sapply(contents(mgu74av2GENENAME), paste, collapse=", "))
                                Annot <- data.frame(ACCNUM=sapply(contents(mgu74av2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mgu74av2SYMBOL), paste, collapse=", "), DESC=sapply(contents(mgu74av2GENENAME), paste, collapse=", "), ENTREZ=sapply(contents(mgu74av2ENTREZID), paste, collapse=", "))
                              } else {
                                if (celfiles@annotation=='pd.mouse430a.2') {
                                  #Annot <- data.frame(ACCNUM=sapply(contents(mouse430a2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mouse430a2SYMBOL), paste, collapse=", "), DESC=sapply(contents(mouse430a2GENENAME), paste, collapse=", "))
                                  Annot <- data.frame(ACCNUM=sapply(contents(mouse430a2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mouse430a2SYMBOL), paste, collapse=", "), DESC=sapply(contents(mouse430a2GENENAME), paste, collapse=", "), ENTREZ=sapply(contents(mouse430a2ENTREZID), paste, collapse=", "))
                                } else {
                                  if (celfiles@annotation=='pd.moe430a') {
                                    #Annot <- data.frame(ACCNUM=sapply(contents(moe430aACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(moe430aSYMBOL), paste, collapse=", "), DESC=sapply(contents(moe430aGENENAME), paste, collapse=", "))
                                    Annot <- data.frame(ACCNUM=sapply(contents(moe430aACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(moe430aSYMBOL), paste, collapse=", "), DESC=sapply(contents(moe430aGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(moe430aENTREZID), paste, collapse=", "))
                                  } else {
                                    if (celfiles@annotation=='pd.hg.u95av2') {
                                      #Annot <- data.frame(ACCNUM=sapply(contents(hgu95av2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu95av2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu95av2GENENAME), paste, collapse=", "))
                                      Annot <- data.frame(ACCNUM=sapply(contents(hgu95av2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu95av2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu95av2GENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hgu95av2ENTREZID), paste, collapse=", "))
                                    } else {
                                      if (celfiles@annotation=='pd.hta.2.0') {
                                        #Annot <- data.frame(ACCNUM=sapply(contents(hta20transcriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hta20transcriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hta20transcriptclusterGENENAME), paste, collapse=", "))
                                        Annot <- data.frame(ACCNUM=sapply(contents(hta20transcriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hta20transcriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hta20transcriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hta20transcriptclusterENTREZID), paste, collapse=", "))
                                      } else {
                                        if (celfiles@annotation=='pd.moex.1.0.st.v1') {
                                          #Annot <- data.frame(ACCNUM=sapply(contents(moex10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(moex10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(moex10sttranscriptclusterGENENAME), paste, collapse=", "))
                                          Annot <- data.frame(ACCNUM=sapply(contents(moex10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(moex10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(moex10sttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(moex10sttranscriptclusterENTREZID), paste, collapse=", "))
                                        } else {
                                          if (celfiles@annotation=='pd.hg.u133b') {
                                            #Annot <- data.frame(ACCNUM=sapply(contents(hgu133bACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133bSYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133bGENENAME), paste, collapse=", "))
                                            Annot <- data.frame(ACCNUM=sapply(contents(hgu133bACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133bSYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133bGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hgu133bENTREZID), paste, collapse=", "))
                                          } else {
                                            if (celfiles@annotation=='pd.hugene.1.1.st.v1') {
                                              #Annot <- data.frame(ACCNUM=sapply(contents(hugene11sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hugene11sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hugene11sttranscriptclusterGENENAME), paste, collapse=", "))
                                              Annot <- data.frame(ACCNUM=sapply(contents(hugene11sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hugene11sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hugene11sttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hugene11sttranscriptclusterENTREZID), paste, collapse=", "))
                                            } else {
                                              if (celfiles@annotation=='pd.mogene.1.1.st.v1') {
                                                #Annot <- data.frame(ACCNUM=sapply(contents(mogene11sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mogene11sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(mogene11sttranscriptclusterGENENAME), paste, collapse=", "))
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

mylist=vector("list",nb)

for (i in 1:nb)
{
  all.genes.con = topTable(ebayes.fit2, coef = i, number=nrow(ebayes.fit2))
  
  all <- merge(all.genes.con, Annot,by.x=0, by.y=0, all.x=T)
  all=all[order(all$P.Value),]
  colnames(all)[1]="probsetID"
  
  #add fold change and rearrange columns
  all$FC = ifelse(all$logFC<0, -1/(2^all$logFC), 2^all$logFC)
  all = all[,c(9,12,2,5,6,3,8,10,11,1,4,7)]
  
  # Write out to a file
  write.table(all,file=paste(projectId,"_",cons[i],"_all_genes.txt",sep=""),sep="\t",row.names=F)

  #GUI displays table for users selected contrast
  mylist[[i]]=all
}
nAll <- merge(exprs(celfiles.rma), Annot,by.x=0, by.y=0, all.x=T)

y<-paste("_",projectId, sep="")
tNorm = tempfile(pattern = "normalized_data_", tmpdir = getwd(), fileext = paste0(y,'.txt'))
write.table(nAll,file=tNorm,sep="\t",row.names=F)  #not in GUI
names(mylist)=cons
#####  END DEG  #####


##### PATHWAYS #####
#user inputs contrast, in this case the first constrast, output will dynamically change
userInput = 'Stim-Control'
num = which(userInput==cons)

pAll = mylist[[num]]

iup=which(pAll$P.Value<0.05 & pAll$logFC>=0)
idw=which(pAll$P.Value<0.05 & pAll$logFC<0)
fin.up=pAll[iup,]

if (length(iup) > 500)
{
  fin.up=fin.up[order(fin.up$P.Value),]
  fin.up=fin.up[1:500,]
}

fin.dw=pAll[idw,]
if (length(idw) > 500)
{
  fin.dw=fin.dw[order(fin.dw$P.Value),]
  fin.dw=fin.dw[1:500,]
}

fin.up$SYMBOL = as.character(fin.up$SYMBOL)
fin.dw$SYMBOL = as.character(fin.dw$SYMBOL)

if (celfiles@annotation=="pd.hg.u133.plus.2" | celfiles@annotation=="pd.hugene.2.0.st" | celfiles@annotation=="pd.clariom.s.human.ht" | celfiles@annotation=="pd.clariom.s.human" | celfiles@annotation=='pd.hg.u133a' | celfiles@annotation=='pd.hugene.1.0.st.v1' | celfiles@annotation=='pd.hg.u133a.2' | celfiles@annotation=='pd.huex.1.0.st.v2' | celfiles@annotation=='pd.hg.u219' | celfiles@annotation=='pd.ht.hg.u133.plus.pm' | celfiles@annotation=='pd.hg.u95av2' | celfiles@annotation=='pd.hta.2.0' | celfiles@annotation=='pd.hg.u133b' | celfiles@annotation=='pd.hugene.1.1.st.v1') 
{
  cat(fin.up$SYMBOL,file=(paste0(projectId,'_',names(mylist[num]),'_Top500_Up.txt')), sep='\n')
  cat(fin.dw$SYMBOL,file=(paste0(projectId,'_',names(mylist[num]),'_Top500_Down.txt')),sep='\n')
} else {
  cat(fin.up$SYMBOL,file=paste0(names(mylist[num]),"_Top500temp_Up.txt"),sep='\n')
  cat(fin.dw$SYMBOL,file=paste0(names(mylist[num]),"_Top500temp_Dw.txt"),sep='\n')
  
  system(paste0("cat ",names(mylist[num]),"_Top500temp_Up.txt | grep -v \"^NA\" | ./m2h | grep -v XXXX | cut -f2 -d\" \" >",projectId,'_',names(mylist[num]),"_Top500_Up.txt"))
  system(paste0("cat ",names(mylist[num]),"_Top500temp_Dw.txt | grep -v \"^NA\" | ./m2h | grep -v XXXX | cut -f2 -d\" \" >",projectId,'_',names(mylist[num]),"_Top500_Down.txt"))
}
system(paste0("cat ",projectId,'_',names(mylist[num]),"_Top500_Up.txt |sort | uniq | ./l2p >",projectId,'_',names(mylist[num]),"_Pathways_Up.txt"))
system(paste0("cat ",projectId,'_',names(mylist[num]),"_Top500_Down.txt |sort | uniq | ./l2p >",projectId,'_',names(mylist[num]),"_Pathways_Down.txt"))

addUpCol = read.delim(paste0(projectId,'_',names(mylist[num]),"_Pathways_Up.txt"), sep = '\t')
addDwCol = read.delim(paste0(projectId,'_',names(mylist[num]),"_Pathways_Down.txt"), sep = '\t')

colnames(addUpCol)=c("pval","fdr","ratio","nb.hits","nb.genes.path","nb.user.genes","tot.back.genes","path_id","source","description","type","gene.list")
colnames(addDwCol)=c("pval","fdr","ratio","nb.hits","nb.genes.path","nb.user.genes","tot.back.genes","path_id","source","description","type","gene.list")
addUpCol = addUpCol[order(addUpCol$pval),]
addDwCol = addDwCol[order(addDwCol$pval),]
addUpCol = addUpCol[,c(8,9,10,11,1,2,3,12,4,5,6,7)]
addDwCol = addDwCol[,c(8,9,10,11,1,2,3,12,4,5,6,7)]
write.table(addUpCol, file = paste0(projectId,'_',names(mylist[num]),"_Pathways_Up.txt"), sep = '\t', row.names = F)
write.table(addDwCol, file = paste0(projectId,'_',names(mylist[num]),"_Pathways_Down.txt"), sep = '\t', row.names = F)
##### END PATHWAYS #####


##### ssGSEA #####
ssgs =  nAll
ssgs = ssgs[ssgs$SYMBOL!='NA',]

#if human or mouse, prepare data for gsva
if (celfiles@annotation=="pd.hg.u133.plus.2" | celfiles@annotation=="pd.hugene.2.0.st" | celfiles@annotation=="pd.clariom.s.human.ht" | celfiles@annotation=="pd.clariom.s.human" | celfiles@annotation=='pd.hg.u133a' | celfiles@annotation=='pd.hugene.1.0.st.v1' | celfiles@annotation=='pd.hg.u133a.2' | celfiles@annotation=='pd.huex.1.0.st.v2' | celfiles@annotation=='pd.hg.u219' | celfiles@annotation=='pd.ht.hg.u133.plus.pm' | celfiles@annotation=='pd.hg.u95av2' | celfiles@annotation=='pd.hta.2.0' | celfiles@annotation=='pd.hg.u133b' | celfiles@annotation=='pd.hugene.1.1.st.v1') {
  ssgs = subset(ssgs, select=-c(ACCNUM,DESC,Row.names,ENTREZ))
  ssgs = aggregate(.~SYMBOL,data=ssgs,mean)                               #aggregate duplicate probes by mean
  rownames(ssgs) = ssgs$SYMBOL
  ssgs = subset(ssgs, select=-c(SYMBOL))
  ssgs = as.matrix(ssgs)
} else {
  ssgs = subset(ssgs, select=-c(ACCNUM,DESC,Row.names,SYMBOL))
  ssgs = aggregate(.~ENTREZ,data=ssgs,mean)                               #aggregate duplicate probes by mean
  rownames(ssgs) = ssgs$ENTREZ
  ssgs = subset(ssgs, select=-c(ENTREZ))
  ssgs = as.matrix(ssgs)
}

#user input gene set (mouse or human depending on dataset)
userInputGeneSet = 'h.all.v6.1.symbols.gmt'
if (userInputGeneSet=="") {
  info("Please select human or mouse gene set from list")
} else {
  gset = getGmt(userInputGeneSet) 
  ssgsResults = gsva(ssgs, gset, method='ssgsea')                           #run GSVA
}

fit1 = lmFit(ssgsResults,design1)
contrast.matrix = makeContrasts(contrasts=cons,levels=design1)
fit2 = contrasts.fit(fit1,contrast.matrix)
ebayes.fit2 = eBayes(fit2)

myPathways=vector("list",nb)
for (i in 1:nb)
{
  all.pathways = topTable(ebayes.fit2, coef=i, number=nrow(ebayes.fit2))
  all.pathways = all.pathways[order(abs(all.pathways$P.Value)),]
  colnames(all.pathways)[2] = 'Avg.Enrichment.Score'
  write.table(all.pathways,file=paste0(projectId,"_",cons[i],"_ssGSEA_pathways.txt"),sep="\t",row.names=T,col.names=NA)
  myPathways[[i]] = all.pathways
}
names(myPathways)=cons

####creates a list of colors specific to each group
fs = factor(celfiles$group)
lFs=levels(fs)
numFs=length(lFs)
colors = list()
for (i in 1:numFs){
  colors[which(fs==lFs[i])] = i*5
}
colors = unlist(colors)
####end




##################   All outputs below: ##################  

#### Raw data histogram ####
hist(celfiles,which="all", main =" Raw Samples distribution")

#### Raw data MAplot ####
facs <- pData(celfiles)$SampleID
nbfacs=length(facs)
for (i in 1:nbfacs) {
    MAplot(celfiles,which=i,plotFun=smoothScatter,refSamples=c(1:nbfacs), main='', cex=2)
}

#### Raw data boxplot ####
boxplot(celfiles, col=colors, which="all", main="Boxplots before normalization",las=2,names=celfiles$SampleID)

#### RLE plot ####
RLE(celfiles.qc, main="RLE plot",names=celfiles$SampleID, col=colors, las=2)

#### NUSE plot ####
NUSE(celfiles.qc, main="NUSE plot",names=celfiles$SampleID, col=colors, las=2)

#### Normalized data histogram ####
hist(celfiles.rma, main ="Distribution after Normalization")

#### Normalized data MA plot ####
facs <- pData(celfiles)$SampleID
nbfacs=length(facs)
for (i in 1:nbfacs) {
    MAplot(celfiles.rma,which=i,plotFun=smoothScatter,refSamples=c(1:nbfacs), main='', cex=2)
}

#### Normalized box plot ####
boxplot(celfiles.rma,col=colors, main="Boxplots after RMA normalization",las=2,names=celfiles$SampleID)

#### 3D PCA ####
tedf= t(exprs(celfiles.rma))
if (length(which(apply(tedf, 2, var)==0)) >= 0){
  tedf = tedf[ , apply(tedf, 2, var) != 0]
}
pca=prcomp(tedf, scale. = T)
rgl.open(useNULL=T)
bg3d('white')
plot3d(pca$x[,1:3],col=colors, type='s',size=2)
group.v=as.vector(celfiles.rma$SampleID)
text3d(pca$x, pca$y, pca$z, group.v, cex=0.6, adj=1.5)
par3d(mouseMode = "trackball")
rglwidget()

#### Sample heatmap ####
mat=as.matrix(dist(t(exprs(celfiles.rma))))
rownames(mat)=celfiles.rma$SampleID
colnames(mat)=rownames(mat)
heatmaply(mat,margins = c(120,120,120,120),colorRampPalette(colors = c("red", "yellow")))

#### DEG table output ####
#userInput is user-selected contrast, will dynamically change
num = which(userInput==cons)
dat = mylist[[num]]
dat = dat[,-3]
dat[,3:4] = format(dat[,3:4], scientific = TRUE)
#Show table on GUI based on user parameters for pvalue and fold change
userinput_pval = 0.05
userinput_FC = 2
if (is.na(userinput_pval) & is.na(userinput_FC)) {
  dat
} else if (is.na(userinput_pval))  {
  dat = dat[(abs(as.numeric(dat[,2])) >= userinput_FC),]
} else if (is.na(userinput_FC)) {
  dat = dat[(as.numeric(dat[,3]) <= userinput_pval),]
} else {
  dat = dat[(as.numeric(dat[,3]) <= userinput_pval & abs(as.numeric(dat[,2])) >= userinput_FC),]
  dat
}
#add hyperlink
dat$ENTREZ <- sapply(dat$ENTREZ, function(x)
  toString(tags$a(href=paste0("https://www.ncbi.nlm.nih.gov/gene/", x), x)))

#### Pathways from top 500 upregulated genes ####
#Table on GUI changes based on user parameters for pathway pvalue
userInput_pathwayPv = 0.05
topUp = addUpCol
if(is.na(userInput_pathwayPv)) {
  topUp
} else {
  topUp = topUp[(as.numeric(topUp[,5]) <= userInput_pathwayPv),]
}
topUp

#### Gene heatmap for upregulated pathway ####
#Heatmap dynamically changes as user selects one pathway at a time from list above
userInput_upPath = 'ko03030'
genes = addUpCol[addUpCol$path_id==userInput_upPath, 'gene.list']
genes = strsplit(as.character(genes),' ')
genes = unlist(genes)

#extract normalized expression, subset by genes, filter columns, aggregate duplicate symbols with means
exp = nAll
exp = exp[exp$SYMBOL %in% genes,]
exp = subset(exp, select = -c(ACCNUM,DESC,ENTREZ,Row.names))
exp = aggregate(.~SYMBOL,data=exp,mean)

#userInput is user-selected contrast from above, will dynamically change
num = which(userInput==cons)
sampleColumns = c(which(celfiles.rma$group==contra$k2[num]),which(celfiles.rma$group==contra$k1[num])) 
rownames(exp) = exp$SYMBOL
exp = subset(exp, select = -c(SYMBOL))
exp = exp[,sampleColumns]
#limit to 100 genes
if(nrow(exp)>100){
  exp = exp[1:100,]
}
#set heatmap parameters
matCol = data.frame(group=celfiles.rma$group[sampleColumns])
rownames(matCol) = celfiles.rma$SampleID[sampleColumns]
matColors = list(group = unique(colors[sampleColumns]))
names(matColors$group) = unique(celfiles.rma$group[sampleColumns])
path_name = paste0(toupper(addUpCol[addUpCol$path_id==userInput_upPath, 'description']),' PATHWAY (max 100 genes)')
#get z-scores by row
exp = t(scale(t(exp)))  

if (nrow(exp) > 30){
  pheatmap(exp, color=inferno(10), main=path_name, annotation_col=matCol, annotation_colors=matColors, drop_levels=TRUE, fontsize_row = 6)
} else {
  pheatmap(exp, color=inferno(10), main=path_name, annotation_col=matCol, annotation_colors=matColors, drop_levels=TRUE, fontsize_row = 10)
}

#### Pathways from top 500 downregulated genes ####
#Table on GUI changes based on user parameters for pathway pvalue
userInput_pathwayPv = 0.05
topDw = addDwCol
if (is.na(userInput_pathwayPv)) {
  topDw
} else {
  topDw = topDw[(as.numeric(topDw[,5]) <= userInput_pathwayPv),]
}
topDw

#### Gene heatmap for downregulated pathway ####
#Heatmap dynamically changes as user selects one pathway at a time from list above
userInput_upPath = 'GO:0006351'
genes = addDwCol[addDwCol$path_id==userInput_upPath, 'gene.list']
genes = strsplit(as.character(genes),' ')
genes = unlist(genes)

#extract normalized expression, subset by genes, filter columns, aggregate duplicate symbols with means
exp = nAll
exp = exp[exp$SYMBOL %in% genes,]
exp = subset(exp, select = -c(ACCNUM,DESC,ENTREZ,Row.names))
exp = aggregate(.~SYMBOL,data=exp,mean)

#userInput is user-selected contrast from above, will dynamically change
num = which(userInput==cons)
sampleColumns = c(which(celfiles.rma$group==contra$k2[num]),which(celfiles.rma$group==contra$k1[num])) 
rownames(exp) = exp$SYMBOL
exp = subset(exp, select = -c(SYMBOL))
exp = exp[,sampleColumns]
#limit to 100 genes
if(nrow(exp)>100){
  exp = exp[1:100,]
}
#set heatmap parameters
matCol = data.frame(group=celfiles.rma$group[sampleColumns])
rownames(matCol) = celfiles.rma$SampleID[sampleColumns]
matColors = list(group = unique(colors[sampleColumns]))
names(matColors$group) = unique(celfiles.rma$group[sampleColumns])
path_name = paste0(toupper(addUpCol[addUpCol$path_id==userInput_upPath, 'description']),' PATHWAY (max 100 genes)')
#get z-scores by row
exp = t(scale(t(exp)))  

if (nrow(exp) > 30){
  pheatmap(exp, color=inferno(10), main=path_name, annotation_col=matCol, annotation_colors=matColors, drop_levels=TRUE, fontsize_row = 6)
} else {
  pheatmap(exp, color=inferno(10), main=path_name, annotation_col=matCol, annotation_colors=matColors, drop_levels=TRUE, fontsize_row = 10)
}

#### Volcano plot ####
#dynamically changes based on user-input contrast
num = which(userInput==cons)
dat=mylist[[num]]

dat = dat[dat$SYMBOL!='NA',]
log_FC=dat$logFC
log_pval=-log10(dat$P.Value)
Significant=rep("NotSignificant",length(log_FC))
Significant[which(dat$P.Value<0.05 & abs(dat$logFC)>=1)]="AbsLogFoldChange>1 & PValue<0.05"
Significant[which(dat$P.Value<0.05 & abs(dat$logFC)<1)]="PValue<0.05"
Significant[which(dat$P.Value>=0.05 & abs(dat$logFC)>=1)]="AbsLogFoldChange>1"
gene=dat$SYMBOL
volcano_data=as.data.frame(cbind(gene,log_FC,log_pval,Significant))
plot_ly(type='scatter', data = volcano_data, x = log_FC, y = log_pval, text = gene, mode = "markers", color = Significant) %>% layout(title=paste0('Volcano plot for: ',names(mylist)[num]),xaxis=list(title="Fold Change",range =c(-5,5),tickvals=c(-5,-4,-3,-2,-1,0,1,2,3,4,5),ticktext=c('-32','-16','-8','-4','-2','1','2','4','8','16','32')),yaxis=list(title="-Log10 pvalue",range =c(0,15)))

#### ssGSEA table ####
num = which(userInput==cons)
dat = myPathways[[num]]
#user defines p-val and fc, dynamically changes output
userInput_ss_pval = 0.05
userInput_ss_fc = 0.1

if (is.na(userInput_ss_pval) & is.na(userInput_ss_fc)) {   
  dat
} else if (is.na(userInput_ss_pval))  {
  dat = dat[(abs(as.numeric(dat[,1])) >= userInput_ss_fc),]
} else if (is.na(userInput_ss_fc)) {
  dat = dat[(as.numeric(dat[,4]) <= userInput_ss_pval),]
} else {
  dat = dat[(as.numeric(dat[,4]) <= userInput_ss_pval & abs(as.numeric(dat[,1])) >= userInput_ss_fc),]
  dat
}

#### ssGSEA heatmap ####
num = which(userInput==cons)
each = myPathways[[num]]
userInput_ss_pval = 0.05
userInput_ss_fc = 0.1

if (is.na(userInput_ss_pval) & is.na(userInput_ss_fc)) {   
  each
} else if (is.na(userInput_ss_pval))  {
  each = each[(abs(as.numeric(each[,1])) >= userInput_ss_fc),]
} else if (is.na(userInput_ss_fc)) {
  each = each[(as.numeric(each[,4]) <= userInput_ss_pval),]
} else {
  each = each[(as.numeric(each[,4]) <= userInput_ss_pval & abs(as.numeric(each[,1])) >= userInput_ss_fc),]
}

sampleColumns = c(which(celfiles.rma$group==contra$k2[num]),which(celfiles.rma$group==contra$k1[num]))  #subset columns (samples) for user input contrast
paths = ssgsResults[rownames(ssgsResults) %in% rownames(each)[1:50],]   #subset diff exprs pathways for user input contrast
paths = paths[,sampleColumns]

matCol = data.frame(group=celfiles.rma$group[sampleColumns])
rownames(matCol) = celfiles.rma$SampleID[sampleColumns]
matColors = list(group = unique(colors[sampleColumns]))
names(matColors$group) = unique(celfiles.rma$group[sampleColumns])

#calculate z scores for rows
paths = t(scale(t(paths))) 

pheatmap(paths,color=inferno(10),annotation_col=matCol,annotation_colors=matColors,drop_levels=TRUE,fontsize=7, main='Enrichment Scores for Top 50 Differentially Expressed ssGSEA Pathways')


#### Download report ####
file = rmarkdown::render('MAAPster_report.Rmd','html_document',paste(projectId,"_","report.html",sep=""))
mytables=list.files(pattern=paste(projectId,".*.html",sep=""))




    
