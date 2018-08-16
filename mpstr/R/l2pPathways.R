#' l2p Pathway Analysis
#' 
#' @param degs Differentially expressed genes object from deg function
#' @param species Species of organism (human or mouse)
#' @param workspace Working directory
#' @param projectId A unique identifier for the project
#' @return List of up and downregulated pathways for each contrast
#' @examples 
#' l2p_pathways = l2pPathways(diff_expr_genes,'human','/Users/name/folderName','NCI_Project_1')
#' l2p_pathways = l2pPathways(diff_expr_genes,'mouse','/Users/name/folderName','NCI_Project_2')
#' @references l2p and m2h from CCBR/CCR/NCI/NIH

l2pPathways = function(degs,species,workspace,projectId) {
  
  listPathways = vector("list",length(degs$listDEGs))
  for (i in 1:length(degs$listDEGs)) {
    up_down = vector("list",2)
    all = degs$listDEGs[[i]]
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
    fin.up$SYMBOL = as.character(fin.up$SYMBOL)
    fin.dw$SYMBOL = as.character(fin.dw$SYMBOL)
    
    if (species=='human') 
    {
      cat(fin.up$SYMBOL,file=(paste0(workspace,projectId,'_',names(degs$listDEGs[i]),'_Top500_Up.txt')), sep='\n')
      cat(fin.dw$SYMBOL,file=(paste0(workspace,projectId,'_',names(degs$listDEGs[i]),'_Top500_Down.txt')),sep='\n')
    }
    else
    {
      
      cat(fin.up$SYMBOL,file=paste0(workspace,names(degs$listDEGs[i]),"_Top500temp_Up.txt"),sep='\n')
      cat(fin.dw$SYMBOL,file=paste0(workspace,names(degs$listDEGs[i]),"_Top500temp_Dw.txt"),sep='\n')
      
      system(paste0("cat ",workspace,names(degs$listDEGs[i]),"_Top500temp_Up.txt | grep -v \"^NA\" |  ./m2h | grep -v XXXX | cut -f2 -d\" \" >",workspace,projectId,'_',names(degs$listDEGs[i]),"_Top500_Up.txt"))
      system(paste0("cat ",workspace,names(degs$listDEGs[i]),"_Top500temp_Dw.txt | grep -v \"^NA\" |  ./m2h | grep -v XXXX | cut -f2 -d\" \" >",workspace,projectId,'_',names(degs$listDEGs[i]),"_Top500_Down.txt"))
    }
    
    
    system(paste0("cat ",workspace,projectId,'_',names(degs$listDEGs[i]),"_Top500_Up.txt |sort | uniq |   ./l2p>",workspace,projectId,'_',names(degs$listDEGs[i]),"_Pathways_Up.txt"))
    #print(paste0("cat ",workspace,projectId,'_',names(degs$listDEGs[i]),"_Top500_Down.txt |sort | uniq |  ./l2p>",workspace,projectId,'_',names(degs$listDEGs[i]),"_Pathways_Down.txt"))
    system(paste0("cat ",workspace,projectId,'_',names(degs$listDEGs[i]),"_Top500_Down.txt |sort | uniq | ./l2p >",workspace,projectId,'_',names(degs$listDEGs[i]),"_Pathways_Down.txt"))
    
    
    addUpCol = read.delim(paste0(workspace,projectId,'_',names(degs$listDEGs[i]),"_Pathways_Up.txt"), sep = '\t')
    addDwCol = read.delim(paste0(workspace,projectId,'_',names(degs$listDEGs[i]),"_Pathways_Down.txt"), sep = '\t')
    
    colnames(addUpCol)=c("P_Value","FDR","Ratio","Number_Hits","Number_Genes_Pathway","Number_User_Genes","Total_Number_Genes","Pathway_ID","Source","Description","Type","Gene_List")
    colnames(addDwCol)=c("P_Value","FDR","Ratio","Number_Hits","Number_Genes_Pathway","Number_User_Genes","Total_Number_Genes","Pathway_ID","Source","Description","Type","Gene_List")
    addUpCol = addUpCol[order(addUpCol$P_Value),]
    addDwCol = addDwCol[order(addDwCol$P_Value),]
    addUpCol = addUpCol[,c(8,9,10,11,1,2,3,12,4,5,6,7)]
    addDwCol = addDwCol[,c(8,9,10,11,1,2,3,12,4,5,6,7)]
    write.table(addUpCol, file = paste0(workspace,projectId,'_',names(degs$listDEGs[i]),"_Pathways_Up.txt"), sep = '\t', row.names = F)
    write.table(addDwCol, file = paste0(workspace,projectId,'_',names(degs$listDEGs[i]),"_Pathways_Down.txt"), sep = '\t', row.names = F)
    up_down[[1]]=addUpCol
    up_down[[2]]=addDwCol
    names(up_down) = c("upregulated_pathways","downregulated_pathways")
    listPathways[[i]] = up_down
  }
  names(listPathways) = names(degs$listDEGs)  
  print("+++pathways+++")
  return(listPathways)
}