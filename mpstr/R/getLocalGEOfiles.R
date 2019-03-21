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
  
  getLocalGEOfiles_ERR = file(paste0(workspace,'/getLocalGEOfiles.err'),open='wt')
  sink(getLocalGEOfiles_ERR,type='message',append=TRUE)
  
  id = gsub(" ","",id,fixed=TRUE) 
  
  #list contents of new directory with zipped CEL files
  SampleName = list.files(path = workspace, pattern = '/*CEL.gz|/*CEL$', ignore.case = T, full.names=T)
  gds <- getGEO(id, GSEMatrix = F,getGPL=T,AnnotGPL=T)             #get meta data 
  mytable=matrix("",length(GSMList(gds)),4)
  colnames(mytable) = c("gsm","title","description","groups")
  for (k in 1:length(GSMList(gds)))                                 #fill table with meta data
  {
    if (is.null(Meta(GSMList(gds)[[k]])$description)) {    
      mytable[k,1:3] <-c(Meta(GSMList(gds)[[k]])$geo_accession[1], Meta(GSMList(gds)[[k]])$title[1], 'No data available')
    } else {
      mytable[k,1:3] <-c(Meta(GSMList(gds)[[k]])$geo_accession[1], Meta(GSMList(gds)[[k]])$title[1], Meta(GSMList(gds)[[k]])$description[1])
    }
  }
  mytable = as.data.frame(mytable)
  mytable$groups = listGroups                               #assign groups to samples
  
  ####creates a list of colors specific to each group
  fs = factor(mytable$groups)
  lFs=levels(fs)
  numFs=length(lFs)
  colors = list()
  for (i in 1:numFs){
    colors[which(fs==lFs[i])] = i*5
  }
  colors = unlist(colors)
  mytable$colors = colors
  ####end colors
  
  # Restructure pheno data to change rownames from filenames to sample titles (for MAplot sample titles)
  pd = mytable
  rownames(pd) = mytable$title
  pd = AnnotatedDataFrame(pd)
  celfiles = read.celfiles(SampleName,phenoData = pd)
  
  
  return(celfiles)
}
sink(type='message')
