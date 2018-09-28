#' 1) Process local Affymetrix CEL files for mouse or human data 
#' 
#' @export
#' @param projectId A unique identifier for the project
#' @param listGroups Group assignments for each sample, follow alphabetical order of samples
#' @param workspace Working directory
#' @return ExpressionFeatureSet object with raw data and phenotype information
#' @examples 
#' celfiles = processCELfiles('NCI_Project_1',c('Ctl_1','Ctl_1','Ctl_1','KO_1','KO_1','KO_1'),'/Users/name/folderName')
#' celfiles = processCELfiles('NCI_Project_2',c('Ctl','Ctl','Ctl','Ctl','RNA_1','RNA_1','RNA_1','RNA_1','RNA_2','RNA_2','RNA_2','RNA_2'),'/Users/name/folderName')      
#' @references See packages tools, Biobase, oligo

processCELfiles <- function(projectId,listGroups,workspace) {
  library(tools)
  library(Biobase)
  library(oligo)
  
  processCELfiles_ERR = file(paste0(workspace,'/processCELfiles.err'),open='wt')
  sink(processCELfiles_ERR,type='message',append=TRUE)
  
  print(workspace)
  SampleName = list.files(path = workspace, pattern = '/*CEL*.gz', ignore.case = T, full.names=T)
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
  print("+++getCELfiles+++")
  return (list(
    files=as.data.frame(apply(pData(celfiles),c(1,2),utils::URLencode))
  ))
}
sink(type='message')
