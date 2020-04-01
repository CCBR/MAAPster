#' l2p Pathway Analysis
#'
#' @export
#' @param degs Differentially expressed genes object from deg function
#' @param species Species of organism (human or mouse)
#' @param workspace Working directory
#' @param projectId A unique identifier for the project
#' @param configuration_path Path to configuration files
#' @return List of up and downregulated pathways for each contrast
#' @examples
#' l2p_pathways = l2pPathways(diff_expr_genes,'human','/Users/name/folderName','NCI_Project_1','/Users/name/config')
#' l2p_pathways = l2pPathways(diff_expr_genes,'mouse','/Users/name/folderName','NCI_Project_2','/Users/name/config')
#' @references l2p and m2h from CCBR/CCR/NCI/NIH

l2pPathways = function(degs,species,workspace,projectId,configuration_path) {
  library(l2p)

  l2pPathways_ERR = file(paste0(workspace,'/l2pPathways.err'),open='wt')
  sink(l2pPathways_ERR,type='message',append=TRUE)
  
  listPathways = vector("list",length(degs$listDEGs))
  for (i in 1:length(degs$listDEGs)) {
    up_down = vector("list",2)
    all = degs$listDEGs[[i]]
    all$SYMBOL = as.character(all$SYMBOL)
    all = all[which(all$SYMBOL!='NA'),]
    iup=which(all$P.Value<0.05 & all$logFC>=0)
    idw=which(all$P.Value<0.05 & all$logFC<0)
    fin.up=all[iup,]
    if (length(iup) > 500)
    {
      fin.up=fin.up[order(fin.up$P.Value),]
      fin.up=fin.up[1:500,]
    }
    fin.dw=all[idw,]
    if (length(idw) > 500)
    {
      fin.dw=fin.dw[order(fin.dw$P.Value),]
      fin.dw=fin.dw[1:500,]
    }
    fin.up = unique(fin.up$SYMBOL)
    fin.dw = unique(fin.dw$SYMBOL)
    univ = unique(all$SYMBOL)
    if (tolower(species)=='mouse') {
      fin.up = m2h(fin.up)
      fin.dw = m2h(fin.dw)
      univ = m2h(univ)
    }
    
    addUpCol = l2p(fin.up,universe=univ)
    addDwCol = l2p(fin.dw,universe=univ)
    
    addUpCol = addUpCol[order(addUpCol$pval),]
    addDwCol = addDwCol[order(addDwCol$pval),]
    colnames(addUpCol)=c("P_Value","FDR","Ratio","Number_Hits","Number_Misses","Number_User_Genes","Total_Genes_Minus_Input","Pathway_ID","Source","Description","Gene_List")
    colnames(addDwCol)=c("P_Value","FDR","Ratio","Number_Hits","Number_Misses","Number_User_Genes","Total_Genes_Minus_Input","Pathway_ID","Source","Description","Gene_List")
    addUpCol = addUpCol[,c(8,9,10,1:7,11)]
    addDwCol = addDwCol[,c(8,9,10,1:7,11)]
    write.table(addUpCol, file = paste0(workspace,'/',projectId,'_',names(degs$listDEGs[i]),"_Pathways_Up.txt"), sep = '\t', row.names = F)
    write.table(addDwCol, file = paste0(workspace,'/',projectId,'_',names(degs$listDEGs[i]),"_Pathways_Down.txt"), sep = '\t', row.names = F)
    up_down[[1]]=addUpCol
    up_down[[2]]=addDwCol
    names(up_down) = c("upregulated_pathways","downregulated_pathways")
    listPathways[[i]] = up_down
  }
  names(listPathways) = names(degs$listDEGs)
  print("+++pathways+++")
  return(listPathways)
  sink(type='message')
}
