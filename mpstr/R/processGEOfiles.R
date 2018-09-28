#' 1) Process GEO mouse or human data from an experiment with raw Affymetrix CEL files, see https://www.ncbi.nlm.nih.gov/geo/
#' 
#' @export
#' @param projectId A unique identifier for the project
#' @param id A GSE id number
#' @param listGroups Group assignments for each sample, follow order of samples on GEO series page
#' @param workspace Working directory
#' @return ExpressionFeatureSet object with raw data and phenotype information
#' @examples 
#' celfiles = processGEOfiles('NCI_Project_1','GSE61989',c('Ctl','Ctl','Ctl','Ctl','RNA_1','RNA_1','RNA_1','RNA_1','RNA_2','RNA_2','RNA_2','RNA_2'),'/Users/name/folderName')      
#' celfiles = processGEOfiles('NCI_Project_2','GSE106988', c('Ctl_1','Ctl_1','Ctl_1','KO_1','KO_1','KO_1'),'/Users/name/folderName')
#' @references See packages oligo, GEOquery, Biobase


processGEOfiles <- function(projectId,id,listGroups,workspace){
  library(GEOquery)
  library(oligo)
  library(Biobase)
  
  if(dir.exists(workspace)) {
    unlink(workspace,recursive = TRUE)                                      #Delete the directory and files in that dir 
  }
  
  if(!dir.exists(workspace)) {
    dir.create(workspace)                                                   #Create a directory 
  }
  processGEOfiles_ERR = file(paste0(workspace,"/processGEOfiles.err"),open='wt')    #Open file to write error messages
  sink(processGEOfiles_ERR,type='message',append=TRUE)                      #Save error messages to file
  
  id = gsub(" ","",id,fixed=TRUE) 
  #system(paste0('rm *.[cC][eE][lL].gz'))                                   #removes previous CEL files if run consecutively
  
  saveMessage = ''
  readID = function(id) {                                                   #error handling: wrong GSE id
    out = tryCatch(
      {
        getGEOSuppFiles(id, makeDirectory = T, baseDir = workspace)
        fileID = paste0(id, '_RAW.tar')
        untar(paste0(workspace,'/',id,'/',fileID),exdir=workspace)
        SampleName = list.files(path = workspace, pattern = '/*CEL*.gz', ignore.case = T, full.names=T)
        if (length(grep('*CEL*',SampleName,ignore.case = T)) == 0) {
          saveMessage = "Raw files must be CEL files"
        }
        celfiles = read.celfiles(SampleName)
        gds <- getGEO(id, GSEMatrix = F,getGPL=T,AnnotGPL=T)             #get meta data 
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
        print("+++loadGSE+++")
        celfiles
      },
      error=function(cond) {
        return(saveMessage)
      }
    )
  }
  celfiles = readID(id)
  
  #Error handling, return data if wrong id
  returnData = function(cels) {
    tableNames=c("gsm","title","description","groups")
    out = tryCatch(
      {
        return (list(
          files=as.data.frame(apply(pData(celfiles),c(1,2),utils::URLencode)),
          tableOrder=tableNames
        ))
      },
      error=function(cond) {
        message(paste0(cels,"\nPlease enter a correct GSE id, this is the number after 'Series' on the GEO series webpage (ex: GSE106988)"))
        return(paste0(cels,"\nPlease enter a correct GSE id, this is the number after 'Series' on the GEO series webpage (ex: GSE106988)"))
      }
    )
  }
    returnData(celfiles)
}
sink(type='message')      #Return messages to console

