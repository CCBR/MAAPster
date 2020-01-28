#' 1) Process local Affymetrix CEL files for mouse or human data 
#' 
#' @export
#' @param projectId A unique identifier for the project
#' @param listGroups Group assignments for each sample, follow alphabetical order of samples
#' @param listBatches Optional list of batches for each sample, follow alphabetical order of samples
#' @param workspace Working directory
#' @return ExpressionFeatureSet object with raw data and phenotype information
#' @examples 
#' celfiles = processCELfiles(projectId='NCI_Project_1',listGroups=c('Ctl_1','Ctl_1','Ctl_1','KO_1','KO_1','KO_1'),workspace='/Users/name/folderName')
#' celfiles = processCELfiles(projectId='NCI_Project_2',listGroups=c('Ctl','Ctl','Ctl','Ctl','RNA_1','RNA_1','RNA_1','RNA_1','RNA_2','RNA_2','RNA_2','RNA_2'),listBatches=c(rep('A',6),rep('B',6)),workspace='/Users/name/folderName')      
#' @references See packages tools, Biobase, oligo

processCELfiles <- function(projectId,listGroups,listBatches=NULL,workspace) {
  library(tools)
  library(Biobase)
  library(oligo)
  library(stringr)
  
  processCELfiles_ERR = file(paste0(workspace,'/processCELfiles.err'),open='wt')
  sink(processCELfiles_ERR,type='message',append=TRUE)
  
  print(workspace)
  SampleName = list.files(path = workspace, pattern = '/*CEL.gz|/*CEL$', ignore.case = T, full.names=T)
  celfiles = read.celfiles(SampleName)
  pData(celfiles)$title = basename(file_path_sans_ext(SampleName))  #add sample name to pheno
  
  pData(celfiles)$groups = listGroups                               #add groups to pheno
  if(!is.null(listBatches)) {
    pData(celfiles)$batch = listBatches                             # add optional batch to pheno
  }
  
  
  # testing - split up multiple groups per sample - delete later
  # maxNgroups = str_count(pData(celfiles)$groups,',')
  # maxNgroups = max(maxNgroups)
  # newGroups = str_split_fixed(pData(celfiles)$groups,',',maxNgroups+1)
  # colnames(newGroups) = paste0('Group_',1:(maxNgroups+1))
  # pData(celfiles) = data.frame(pData(celfiles),newGroups)
  
  
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
