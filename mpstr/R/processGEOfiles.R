#' 1) Process GEO mouse or human data from an experiment with raw Affymetrix CEL files, see https://www.ncbi.nlm.nih.gov/geo/
#' 
#' @export
#' @param projectId A unique identifier for the project
#' @param id A GSE id number
#' @param listGroups Group assignments for each sample, follow order of samples on GEO series page
#' @param listBatches Optional list of batches for each sample, follow alphabetical order of samples
#' @param workspace Working directory
#' @param chip Optional argument specifying platform (ie GPL571) in multi chip GEO datasets
#' @return ExpressionFeatureSet object with raw data and phenotype information
#' @examples 
#' celfiles = processGEOfiles(projectId='NCI_Project_1',id='GSE61989',listGroups=c('Ctl','Ctl','Ctl','Ctl','RNA_1','RNA_1','RNA_1','RNA_1','RNA_2','RNA_2','RNA_2','RNA_2'),workspace='/Users/name/folderName')      
#' celfiles = processGEOfiles(projectId='NCI_Project_2',id='GSE106988',listGroups=c('Ctl_1','Ctl_1','Ctl_1','KO_1','KO_1','KO_1'),listBatches=c(rep('A',6),rep('B',6)),workspace='/Users/name/folderName',chip='GPL571')
#' @references See packages oligo, GEOquery, Biobase


processGEOfiles <- function(projectId,id,listGroups,listBatches=NULL,workspace,chip=NULL){
  library(GEOquery)
  library(oligo)
  library(Biobase)
  
  listPlatforms = c('pd.mogene.2.0.st', 'pd.hg.u133.plus.2', 'pd.hugene.2.0.st', 'pd.clariom.s.human.ht', 'pd.clariom.s.human', 'pd.clariom.s.mouse.ht', 'pd.clariom.s.mouse', 'pd.mouse430.2', 'pd.hg.u133a', 'pd.hugene.1.0.st.v1', 'pd.mogene.1.0.st.v1', 'pd.hg.u133a.2', 'pd.huex.1.0.st.v2', 'pd.hg.u219', 'pd.mg.u74av2', 'pd.mouse430a.2', 'pd.moe430a', 'pd.hg.u95av2', 'pd.hta.2.0', 'pd.moex.1.0.st.v1', 'pd.hg.u133b', 'pd.hugene.1.1.st.v1', 'pd.mogene.1.1.st.v1', 'pd.hugene.2.1.st', 'pd.ht.hg.u133a', 'pd.clariom.d.human')
  
  if(dir.exists(workspace)) {
    unlink(workspace,recursive = TRUE)                                      #Delete the directory and files in that dir 
  }
  
  if(!dir.exists(workspace)) {
    dir.create(workspace)                                                   #Create a directory 
  }
  processGEOfiles_ERR = file(paste0(workspace,"/processGEOfiles.err"),open='wt')    #Open file to write error messages
  sink(processGEOfiles_ERR,type='message',append=TRUE)                      #Save error messages to file
  
  processGEOfiles_OUT = file(paste0(workspace,'/processGEOfiles.out'),open='wt')
  sink(processGEOfiles_OUT,type='output',append=TRUE)
  
  id = gsub(" ","",id,fixed=TRUE) 
  #system(paste0('rm *.[cC][eE][lL].gz'))                                   #removes previous CEL files if run consecutively
  
  # error handling: invalid GSE 
  check_GSE = function(id) {
    tryCatch(
      {
        getGEOSuppFiles(id, makeDirectory = T, baseDir = workspace)
        return(NULL)
      },
    error = function(cond) return("Invalid GSE id, check https://www.ncbi.nlm.nih.gov/geo/ and enter a valid id, ex: GSE118295")
    )
  }
  isGSE = check_GSE(id)
  if(!is.null(isGSE)) return(isGSE)
  
  # error handling: valid GSE but not Affymetrix
  check_GSE_2 = function(id) {
    tryCatch(
      {
        fileID = paste0(id, '_RAW.tar')
        untar(paste0(workspace,'/',id,'/',fileID),exdir=workspace)
        tryCels = list.files(path = workspace, pattern = '/*CEL.gz|/*CEL$', ignore.case = T, full.names=T)
        if (length(grep('*CEL*',tryCels,ignore.case = T)) == 0) {
          stop()
        }
        gds <- getGEO(id, GSEMatrix = F,getGPL=T,AnnotGPL=T)             #get meta data 
        return(gds)
      },
      error = function(cond) {
        gds <- getGEO(id)  
        platform = c()
        for (i in 1:length(gds)) {
          platform = c(platform,gds[[i]]@annotation)
        }
        platform = unlist(platform)
        platform = paste0(platform,collapse = ', ')
        return(paste0("This GSE id does not correspond to an Affymetrix dataset. Your platform id(s): ",platform,". Check https://www.ncbi.nlm.nih.gov/geo/ for details."))
      }  
    )
  }
  gds = check_GSE_2(id)
  if(class(gds)=='character') return(gds)
  
  SampleName = list.files(path = workspace, pattern = '/*CEL.gz|/*CEL$', ignore.case = T, full.names=T)
  
  # check for multiple chips
  if (length(gds@gpls) > 1) {                   # if multi chips present
    if (!is.null(chip)) {                       # if chip is specified

      # subset by chip
      shortList = Filter(function(x) {Meta(x)$platform_id==chip},GSMList(gds))
      
      # create fillable table, add batch column if necessary
      if(is.null(listBatches)) {
        mytable=matrix("",length(shortList),4)
        colnames(mytable) = c("gsm","title","description","groups")
      } else {
        mytable=matrix("",length(shortList),5)
        colnames(mytable) = c("gsm","title","description","groups","batch")
      }
      
      #fill table with meta data
      for (k in 1:length(shortList))                                
      {
        if (is.null(shortList[[k]]@header$description)) {    
          mytable[k,1:3] <-c(shortList[[k]]@header$geo_accession[1], shortList[[k]]@header$title[1], 'No data available')
        } else {
          mytable[k,1:3] <-c(shortList[[k]]@header$geo_accession[1], shortList[[k]]@header$title[1], shortList[[k]]@header$description[1])
        }
      }
      
      # if multi chips detected, but no chip specified, return list with chip names and samples
    } else {
      chipChoices = list()
      for (i in 1:length(names(GPLList(gds)))) {
        each = Filter(function(x) {Meta(x)$platform_id==names(GPLList(gds))[i]},GSMList(gds))
        chipChoices[[i]] = data.frame(matrix(ncol=5,nrow=length(each)))
        colnames(chipChoices[[i]]) = c('gsm','title','description','groups','color')
        for (j in 1:length(each)) {
          if(is.null(each[[j]]@header$description)) {
            chipChoices[[i]][j,1:3] = c(each[[j]]@header$geo_accession[1], each[[j]]@header$title[1], 'No data available')
          } else {
            chipChoices[[i]][j,1:3] = c(each[[j]]@header$geo_accession[1], each[[j]]@header$title[1], each[[j]]@header$description[1])
          }
        }
        names(chipChoices)[i] = names(GPLList(gds))[i]
      }
      return(chipChoices)
    }
    
    # if only 1 chip detected, proceed
  } else {
    
    # create fillable table, add batch column if nec
    if(is.null(listBatches)) {
      mytable=matrix("",length(GSMList(gds)),4)
      colnames(mytable) = c("gsm","title","description","groups")
    } else {
      mytable=matrix("",length(GSMList(gds)),5)
      colnames(mytable) = c("gsm","title","description","groups","batch")
    }
    
    #fill table with meta data
    for (k in 1:length(GSMList(gds)))                                 
    {
      if (is.null(Meta(GSMList(gds)[[k]])$description)) {    
        mytable[k,1:3] <-c(Meta(GSMList(gds)[[k]])$geo_accession[1], Meta(GSMList(gds)[[k]])$title[1], 'No data available')
      } else {
        mytable[k,1:3] <-c(Meta(GSMList(gds)[[k]])$geo_accession[1], Meta(GSMList(gds)[[k]])$title[1], Meta(GSMList(gds)[[k]])$description[1])
      }
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
  print("+++loadGSE+++")
  
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
  
  if (length(gds@gpls) == 1) {
    celfiles = read.celfiles(SampleName,phenoData = pd, verbose=FALSE)
  } else {
    celfiles = read.celfiles(SampleName[gsub('(_|\\.).*','',basename(SampleName)) %in% names(shortList)], phenoData = pd, verbose = FALSE)
  }
  
  # check if supported Affymetrix chip
  if (!celfiles@annotation%in%listPlatforms) {
    return(paste0('Your Affymetrix platform ',celfiles@annotation,' is not yet supported.  You may request to have it added by contacting NCIMicroArrayWebAdmin@mail.nih.gov.'))
  }
  
  if (is.null(listBatches)) {
    tableNames=c("gsm","title","description","groups")
  } else {
    tableNames=c('gsm','title','description','groups','batch')
  }
  
  return (list(
    files=as.data.frame(apply(pData(celfiles),c(1,2),utils::URLencode)),
    tableOrder=tableNames
  ))
  
}
sink(type='message')  
