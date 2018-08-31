#' 5) Gene heatmap with selected genes from l2p pathways
#' 
#' @export
#' @param degs Differentially expressed genes object from deg function
#' @param paths Pathways object from pathways function
#' @param contrast User-selected contrast to display, must be a contrast from listContrasts argument of deg function
#' @param upOrDown Is user-selected pathway from up or down list? (matches names of up and down pathways from output of pathways function)
#' @param pathway_name What is the pathway name? Found in description column of pathway table
#' @param saveImageFileName Name of heatmap
#' @examples 
#' geneHeatmap(diff_expr_genes, l2p_pathways, 'RNA_1-Ctl', 'upregulated_pathways','oxidation-reduction process','Heatmap_Redox','/Users/name/folderName') 
#' geneHeatmap(diff_expr_genes, l2p_pathways, 'KO_1-Ctl_1', 'downregulated_pathways','G-protein coupled receptor activity','Heatmap_GPCR','/Users/name/folderName') 
#' @note Nothing to return, outputs heatmap
#' @references See pheatmap package, mouse/human homologs extracted from http://www.informatics.jax.org/homology.shtml

geneHeatmap = function(degs, paths, contrast, upOrDown, pathway_name,saveImageFileName,path) {
  library(pheatmap)
  human2mouse = read.delim(paste0(path,'/human2mouse.csv', sep = ''),sep=',')
  paths = paths[[contrast]][[upOrDown]]
  genes = paths$Gene_List[paths$Description==pathway_name]              #select user input pathway, extract genes
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
    pheatmap(exp, main=path_name, annotation_col=matCol, annotation_colors=matColors, drop_levels=TRUE, fontsize_row = 4,filename=saveImageFileName)
  } else {
    pheatmap(exp, main=path_name, annotation_col=matCol, annotation_colors=matColors, drop_levels=TRUE, fontsize_row = 10,filename=saveImageFileName)
  }
}