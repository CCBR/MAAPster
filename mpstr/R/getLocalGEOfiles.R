#' 1) Process GEO mouse or human data from an experiment with raw Affymetrix CEL files, see https://www.ncbi.nlm.nih.gov/geo/
#' 
#' @export
#' @param projectId A unique identifier for the project
#' @param id A GSE id number
#' @param listGroups Group assignments for each sample, follow order of samples on GEO series page
#' @param listBatches Optional list of batches for each sample, follow alphabetical order of samples
#' @param workspace Working directory
#' @return ExpressionFeatureSet object with raw data and phenotype information
#' @examples 
#' celfiles = getLocalGEOfiles(projectId='NCI_Project_1',id='GSE61989',listGroups=c('Ctl','Ctl','Ctl','Ctl','RNA_1','RNA_1','RNA_1','RNA_1','RNA_2','RNA_2','RNA_2','RNA_2'),workspace='/Users/name/folderName')      
#' celfiles = getLocalGEOfiles(projectId='NCI_Project_2',id='GSE106988',listGroups=c('Ctl_1','Ctl_1','Ctl_1','KO_1','KO_1','KO_1'),listBatches=c(rep('A',6),rep('B',6)),workspace='/Users/name/folderName')
#' @references See packages oligo, GEOquery, Biobase

getLocalGEOfiles <- function(projectId,id,listGroups,listBatches=NULL,workspace){
  library(GEOquery)
  library(oligo)
  library(Biobase)
  
  getLocalGEOfiles_ERR = file(paste0(workspace,'/getLocalGEOfiles.err'),open='wt')
  sink(getLocalGEOfiles_ERR,type='message',append=TRUE)
  
  id = gsub(" ","",id,fixed=TRUE) 
  
  #list contents of new directory with zipped CEL files
  SampleName = list.files(path = workspace, pattern = '/*CEL.gz|/*CEL$', ignore.case = T, full.names=T)
  gds <- getGEO(id, GSEMatrix = F,getGPL=T,AnnotGPL=T)             #get meta data 
  
  if(is.null(listBatches)) {
    mytable=matrix("",length(GSMList(gds)),4)
    colnames(mytable) = c("gsm","title","description","groups")
  } else {
    mytable=matrix("",length(GSMList(gds)),5)
    colnames(mytable) = c("gsm","title","description","groups","batch")
  }
  
  for (k in 1:length(GSMList(gds)))                                 #fill table with meta data
  {
    if (is.null(Meta(GSMList(gds)[[k]])$description)) {    
      mytable[k,1:3] <-c(Meta(GSMList(gds)[[k]])$geo_accession[1], Meta(GSMList(gds)[[k]])$title[1], 'No data available')
    } else {
      mytable[k,1:3] <-c(Meta(GSMList(gds)[[k]])$geo_accession[1], Meta(GSMList(gds)[[k]])$title[1], Meta(GSMList(gds)[[k]])$description[1])
    }
  }
  mytable = as.data.frame(mytable)
  mytable$groups = listGroups                       #assign groups to samples
  if(!is.null(listBatches)) {
    mytable$batch = listBatches                     # assign optional batch to samples
  }
  
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
  
  #parse sample names
  colNames = as.character(mytable$title)
  parsedNames = vector("list",length(colNames))
  for (i in 1:length(colNames)) {
    temp = substring(colNames[i], seq(1, nchar(colNames[i])-1, nchar(colNames[i])/2), seq(nchar(colNames[i])/2, nchar(colNames[i]), nchar(colNames[i])-nchar(colNames[i])/2))
    temp = paste(temp, collapse = '\n')
    parsedNames[[i]] = temp
  } #end parse
  
  rownames(pd) = parsedNames
  pd = AnnotatedDataFrame(pd)
  celfiles = read.celfiles(SampleName,phenoData = pd)
  
  
  return(celfiles)
}
sink(type='message')
