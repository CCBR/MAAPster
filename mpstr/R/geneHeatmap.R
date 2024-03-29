#' 5) Gene heatmap with selected genes from l2p pathways
#' 
#' @export
#' @param degs Differentially expressed genes object from deg function
#' @param paths Pathways object from pathways function
#' @param contrast User-selected contrast to display, must be a contrast from listContrasts argument of deg function
#' @param upOrDown Is user-selected pathway from up or down list? (matches names of up and down pathways from output of pathways function)
#' @param pathway_name What is the pathway name? Found in description column of pathway table
#' @param saveImageFileName Name of heatmap
#' @param path Path to app configuration files
#' @param workspace Path to output plots
#' @param species Either human or mouse
#' @examples 
#' geneHeatmap(diff_expr_genes, l2p_pathways, 'RNA_1-Ctl', 'upregulated_pathways','oxidation-reduction process','Heatmap_Redox','/Users/name/folderName','/Users/name/config','human') 
#' geneHeatmap(diff_expr_genes, l2p_pathways, 'KO_1-Ctl_1', 'downregulated_pathways','G-protein coupled receptor activity','Heatmap_GPCR','/Users/name/folderName','/Users/name/config','human') 
#' @note Nothing to return, outputs heatmap
#' @references See pheatmap package, mouse/human homologs extracted from http://www.informatics.jax.org/homology.shtml

geneHeatmap = function(degs, paths, contrast, upOrDown, pathway_name,saveImageFileName,path,workspace,species) {
  library(pheatmap)
  
  geneHeatmap_ERR = file(paste0(workspace,'/geneHeatmap.err'),open='wt')
  sink(geneHeatmap_ERR,type='message',append=TRUE)
  
  human2mouse = read.delim(paste0(path,'/human2mouse.csv', sep = ''),sep=',')
  paths = paths[[contrast]][[upOrDown]]
  genes = paths$Gene_List[paths$Pathway_Name==pathway_name]              #select user input pathway, extract genes
  genes = strsplit(as.character(genes),' ')
  genes = unlist(genes)
  exp = degs$norm_plots_annotated                                             #extract normalized expression, subset by genes, aggregate duplicate symbols by mean
  exp = exp[exp$SYMBOL %in% genes,]
  if (tolower(species)=='mouse') {
    genes = human2mouse$mouse[human2mouse$human %in% genes]
    genes = as.character(genes)
    exp = degs$norm_plots_annotated                                            
    exp = exp[exp$SYMBOL %in% genes,]
  }
  exp = subset(exp, select = -c(ACCNUM,DESC,ENTREZ,Row.names))
  exp = aggregate(.~SYMBOL,data=exp,mean)
  sampleColumns = c(which(degs$pheno$groups==gsub("-.*$","",contrast)),which(degs$pheno$groups==gsub("^.*-","",contrast)))
  rownames(exp) = exp$SYMBOL
  exp = subset(exp, select = -c(SYMBOL))
  exp = exp[,sampleColumns]
  matCol = data.frame(group=degs$pheno$groups[sampleColumns])           #set heatmap parameters
  rownames(matCol) = rownames(degs$pheno)[sampleColumns]
  matColors = list(group = unique(degs$pheno$colors[sampleColumns]))
  names(matColors$group) = unique(degs$pheno$groups[sampleColumns])
  path_name = pathway_name
  exp = t(scale(t(exp)))                                                #get z-scores by row
  
  # break up sample names
  # colNames = rownames(matCol)
  # parsedNames = vector("list",length(colNames))
  # for (i in 1:length(colNames)) {
  #   temp = substring(colNames[i], seq(1, nchar(colNames[i])-1, nchar(colNames[i])/2), seq(nchar(colNames[i])/2, nchar(colNames[i]), nchar(colNames[i])-nchar(colNames[i])/2))
  #   temp = paste(temp, collapse = '\n')
  #   parsedNames[[i]] = temp
  # }
  
  # if (nrow(exp) > 30){
  #   pheatmap(exp, main=paste0(path_name,'\n(Row Z-Scores)'),annotation_col=matCol, annotation_colors=matColors, drop_levels=TRUE, fontsize_row = 4,labels_col = parsedNames,filename=paste0(workspace,'/',saveImageFileName))
  # } else {
  #   pheatmap(exp, main=paste0(path_name,'\n(Row Z-Scores)'),annotation_col=matCol, annotation_colors=matColors, drop_levels=TRUE, fontsize_row = 10,labels_col = parsedNames,filename=paste0(workspace,'/',saveImageFileName))
  # }
  
  if (nrow(exp) > 30){
    pheatmap(exp, main=paste0(path_name,'\n(Row Z-Scores)'),annotation_col=matCol, annotation_colors=matColors, drop_levels=TRUE, fontsize_row = 4,labels_col = as.character(rownames(matCol)),filename=paste0(workspace,'/',saveImageFileName))
  } else {
    pheatmap(exp, main=paste0(path_name,'\n(Row Z-Scores)'),annotation_col=matCol, annotation_colors=matColors, drop_levels=TRUE, fontsize_row = 10,labels_col = as.character(rownames(matCol)),filename=paste0(workspace,'/',saveImageFileName))
  }
  
  sink(type='message')
}