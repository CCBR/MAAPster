#' 1) Process GEO mouse or human data from an experiment with raw Affymetrix CEL files, see https://www.ncbi.nlm.nih.gov/geo/
#' 
#' @export
#' @param projectId A unique identifier for the project
#' @param id A GSE id number
#' @param listGroups Group assignments for each sample, follow order of samples on GEO series page
#' @param workspace Working directory
#' @return ExpressionFeatureSet object with raw data and phenotype information
#' @examples 
#' celfiles = getLocalGEOfiles('NCI_Project_1','GSE61989',c('Ctl','Ctl','Ctl','Ctl','RNA_1','RNA_1','RNA_1','RNA_1','RNA_2','RNA_2','RNA_2','RNA_2'),'/Users/name/folderName')      
#' celfiles = getLocalGEOfiles('NCI_Project_2','GSE106988', c('Ctl_1','Ctl_1','Ctl_1','KO_1','KO_1','KO_1'),'/Users/name/folderName')
#' @references See packages oligo, GEOquery, Biobase

getLocalGEOfiles <- function(projectId,id,listGroups,workspace){
  library(GEOquery)
  library(oligo)
  library(Biobase)
  
  id = gsub(" ","",id,fixed=TRUE) 
  
  #list contents of new directory with zipped CEL files
  SampleName = list.files(path = workspace, pattern = '/*CEL*.gz', ignore.case = T, full.names=T)
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
  ####end colors
  
  
  return(celfiles)
}