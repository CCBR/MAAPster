#' 3) Differentially Expressed Genes (DEGs)
#' 
#' @export
#' @param norm ExpressionFeatureSet object with normalized data from QCnorm function
#' @param cons List groups to compare, groups must match assignments in listGroups param of previous functions
#' @param projectId A unique identifier for the project
#' @param workspace Working directory
#' @return List of DEGs for each contrast, annotated normalized data, and pheno data
#' @examples 
#' diff_expr_genes = diffExprGenes(norm_celfiles,c("RNA_1-Ctl","RNA_2-Ctl"),'NCI_Project_1','/Users/name/folderName')
#' diff_expr_genes = diffExprGenes(norm_celfiles,c("KO_1-Ctl_1"),'NCI_Project_2','/Users/name/folderName')
#' @note Baseline/control group should be listed second in contrast assignment: c("experimentalGroup-ControlGroup","diseasedGroup-ControlGroup")
#' @note Outputs volcano plot (using plotly) of differentially expressed genes
#' @references See plotly, limma and Bioconductor platform design/annotation documentation

diffExprGenes = function(norm,cons,projectId,workspace) {
  library(limma)
  library(pd.mogene.2.0.st)
  library(mogene20sttranscriptcluster.db)
  library(pd.hg.u133.plus.2)
  library(hgu133plus2.db)
  library(pd.hugene.2.0.st)
  library(hugene20sttranscriptcluster.db)
  library(pd.clariom.s.human.ht)
  library(clariomshumanhttranscriptcluster.db)
  library(pd.clariom.s.human)
  library(clariomshumantranscriptcluster.db)
  library(pd.clariom.s.mouse.ht)
  library(clariomsmousehttranscriptcluster.db)
  library(pd.clariom.s.mouse)
  library(clariomsmousetranscriptcluster.db)
  library(pd.mouse430.2)
  library(mouse4302.db)
  library(pd.hg.u133a)
  library(hgu133a.db)
  library(pd.hugene.1.0.st.v1)
  library(hugene10sttranscriptcluster.db)
  library(pd.mogene.1.0.st.v1)
  library(mogene10sttranscriptcluster.db)
  library(pd.hg.u133a.2)
  library(hgu133a2.db)
  library(pd.huex.1.0.st.v2)
  library(huex10sttranscriptcluster.db)   
  library(pd.hg.u219)
  library(hgu219.db)
  library(pd.mg.u74av2)
  library(mgu74av2.db)
  library(pd.mouse430a.2)
  library(mouse430a2.db)
  library(pd.moe430a)
  library(moe430a.db)
  library(pd.hg.u95av2)
  library(hgu95av2.db)
  library(pd.hta.2.0)
  library(hta20transcriptcluster.db)
  library(pd.moex.1.0.st.v1)
  library(moex10sttranscriptcluster.db)
  library(pd.hg.u133b)
  library(hgu133b.db)
  library(pd.hugene.1.1.st.v1)
  library(hugene11sttranscriptcluster.db)
  library(pd.mogene.1.1.st.v1)
  library(mogene11sttranscriptcluster.db)
  library(pd.hugene.2.1.st)
  library(hugene21sttranscriptcluster.db)
  
  diffExprGenes_ERR = file(paste0(workspace,'/diffExprGenes.err'),open='wt')
  sink(diffExprGenes_ERR,type='message',append=TRUE)
  
  myfactor <- factor(pData(norm)$groups)
  design1 <- model.matrix(~0+myfactor)
  colnames(design1) <- levels(myfactor)
  fit1 <- lmFit(norm,design1)
  contrast.matrix <- makeContrasts(contrasts=cons,levels=design1)
  fit2 <- contrasts.fit(fit1, contrast.matrix)
  
  # try/catch for 1 sample/group
  attemptEbayes = function(theFit) {
    out = tryCatch(
      {
        ebayes.fit2=eBayes(fit2) 
      },
      error=function(cond){
        return('At least 2 samples per group required for differential expression analysis.')
      }
    ) 
    return(out)
  } 
  outMsg = attemptEbayes(fit2)
  
  if (norm@annotation=="pd.mogene.2.0.st") {  
    Annot <- data.frame(ACCNUM=sapply(contents(mogene20sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mogene20sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(mogene20sttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(mogene20sttranscriptclusterENTREZID), paste, collapse=", "))     
  } else {
    # if (input$Platform=="h133p2") {
    if (norm@annotation=="pd.hg.u133.plus.2") {
      Annot <- data.frame(ACCNUM=sapply(contents(hgu133plus2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133plus2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133plus2GENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hgu133plus2ENTREZID), paste, collapse=", "))
    } else {
      if (norm@annotation=="pd.hugene.2.0.st") {
        Annot <- data.frame(ACCNUM=sapply(contents(hugene20sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hugene20sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hugene20sttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hugene20sttranscriptclusterENTREZID), paste, collapse=", "))
      } else {
        if (norm@annotation=="pd.clariom.s.human.ht") {
          Annot <- data.frame(ACCNUM=sapply(contents(clariomshumanhttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(clariomshumanhttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(clariomshumanhttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(clariomshumanhttranscriptclusterENTREZID), paste, collapse=", "))
        } else {
          if (norm@annotation=="pd.clariom.s.mouse.ht") {
            Annot <- data.frame(ACCNUM=sapply(contents(clariomsmousehttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(clariomsmousehttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(clariomsmousehttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(clariomsmousehttranscriptclusterENTREZID), paste, collapse=", "))
          } else {
            if (norm@annotation=="pd.clariom.s.mouse") {
              Annot <- data.frame(ACCNUM=sapply(contents(clariomsmousetranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(clariomsmousetranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(clariomsmousetranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(clariomsmousetranscriptclusterENTREZID), paste, collapse=", "))
            } else {
              if (norm@annotation=="pd.clariom.s.human") {
                Annot <- data.frame(ACCNUM=sapply(contents(clariomshumantranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(clariomshumantranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(clariomshumantranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(clariomshumantranscriptclusterENTREZID), paste, collapse=", "))
              } else {
                if (norm@annotation=="pd.mouse430.2") {
                  Annot <- data.frame(ACCNUM=sapply(contents(mouse4302ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mouse4302SYMBOL), paste, collapse=", "), DESC=sapply(contents(mouse4302GENENAME), paste, collapse=", "), ENTREZ=sapply(contents(mouse4302ENTREZID), paste, collapse=", "))
                } else {
                  if (norm@annotation=='pd.hg.u133a') {
                    Annot <- data.frame(ACCNUM=sapply(contents(hgu133aACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133aSYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133aGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hgu133aENTREZID), paste, collapse=", "))
                  } else {
                    if (norm@annotation=='pd.hugene.1.0.st.v1') {
                      Annot <- data.frame(ACCNUM=sapply(contents(hugene10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hugene10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hugene10sttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hugene10sttranscriptclusterENTREZID), paste, collapse=", "))
                    } else {
                      if (norm@annotation=='pd.mogene.1.0.st.v1') {
                        Annot <- data.frame(ACCNUM=sapply(contents(mogene10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mogene10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(mogene10sttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(mogene10sttranscriptclusterENTREZID), paste, collapse=", "))
                      } else {
                        if (norm@annotation=='pd.hg.u133a.2') {
                          Annot <- data.frame(ACCNUM=sapply(contents(hgu133a2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133a2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133a2GENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hgu133a2ENTREZID), paste, collapse=", "))
                        } else {
                          if (norm@annotation=='pd.huex.1.0.st.v2') {
                            Annot <- data.frame(ACCNUM=sapply(contents(huex10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(huex10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(huex10sttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(huex10sttranscriptclusterENTREZID), paste, collapse=", "))
                          } else {
                            if (norm@annotation=='pd.hg.u219') {
                              Annot <- data.frame(ACCNUM=sapply(contents(hgu219ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu219SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu219GENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hgu219ENTREZID), paste, collapse=", "))
                            } else {
                              if (norm@annotation=='pd.ht.hg.u133.plus.pm') {
                                Annot <- data.frame(ACCNUM=sapply(contents(hgu133plus2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133plus2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133plus2GENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hgu133plus2ENTREZID), paste, collapse=", "))
                              } else {
                                if (norm@annotation=='pd.mg.u74av2') {
                                  Annot <- data.frame(ACCNUM=sapply(contents(mgu74av2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mgu74av2SYMBOL), paste, collapse=", "), DESC=sapply(contents(mgu74av2GENENAME), paste, collapse=", "), ENTREZ=sapply(contents(mgu74av2ENTREZID), paste, collapse=", "))
                                } else {
                                  if (norm@annotation=='pd.mouse430a.2') {
                                    Annot <- data.frame(ACCNUM=sapply(contents(mouse430a2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mouse430a2SYMBOL), paste, collapse=", "), DESC=sapply(contents(mouse430a2GENENAME), paste, collapse=", "), ENTREZ=sapply(contents(mouse430a2ENTREZID), paste, collapse=", "))
                                  } else {
                                    if (norm@annotation=='pd.moe430a') {
                                      Annot <- data.frame(ACCNUM=sapply(contents(moe430aACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(moe430aSYMBOL), paste, collapse=", "), DESC=sapply(contents(moe430aGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(moe430aENTREZID), paste, collapse=", "))
                                    } else {
                                      if (norm@annotation=='pd.hg.u95av2') {
                                        Annot <- data.frame(ACCNUM=sapply(contents(hgu95av2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu95av2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu95av2GENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hgu95av2ENTREZID), paste, collapse=", "))
                                      } else {
                                        if (norm@annotation=='pd.hta.2.0') {
                                          Annot <- data.frame(ACCNUM=sapply(contents(hta20transcriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hta20transcriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hta20transcriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hta20transcriptclusterENTREZID), paste, collapse=", "))
                                        } else {
                                          if (norm@annotation=='pd.moex.1.0.st.v1') {
                                            Annot <- data.frame(ACCNUM=sapply(contents(moex10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(moex10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(moex10sttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(moex10sttranscriptclusterENTREZID), paste, collapse=", "))
                                          } else {
                                            if (norm@annotation=='pd.hg.u133b') {
                                              Annot <- data.frame(ACCNUM=sapply(contents(hgu133bACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133bSYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133bGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hgu133bENTREZID), paste, collapse=", "))
                                            } else {
                                              if (norm@annotation=='pd.hugene.1.1.st.v1') {
                                                Annot <- data.frame(ACCNUM=sapply(contents(hugene11sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hugene11sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hugene11sttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hugene11sttranscriptclusterENTREZID), paste, collapse=", "))
                                              } else {
                                                if (norm@annotation=='pd.mogene.1.1.st.v1') {
                                                  Annot <- data.frame(ACCNUM=sapply(contents(mogene11sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mogene11sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(mogene11sttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(mogene11sttranscriptclusterENTREZID), paste, collapse=", "))
                                                } else {
                                                  if (raw()@annotation=='pd.hugene.2.1.st') {
                                                    Annot <- data.frame(ACCNUM=sapply(contents(hugene21sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hugene21sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hugene21sttranscriptclusterGENENAME), paste, collapse=", "), ENTREZ=sapply(contents(hugene21sttranscriptclusterENTREZID), paste, collapse=", "))
                                                  }
                                                }
                                              }
                                            }
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  catchAndReturn = function(outMsg) {
    out = tryCatch(
      {
        ebayes.fit2 = outMsg
        numContrasts = length(cons)
        listDEGs = vector("list",numContrasts)                                            #initialize output list for each contrast
        for (i in 1:numContrasts)
        {
          all.genes.con = topTable(ebayes.fit2, coef = i, number=nrow(ebayes.fit2))
          all <- merge(all.genes.con, Annot,by.x=0, by.y=0, all.x=T)                      #annotate
          all=all[order(all$P.Value),]
          colnames(all)[1]="probsetID"
          all$FC = ifelse(all$logFC<0, -1/(2^all$logFC), 2^all$logFC)                     #add fold change and rearrange columns
          all = all[,c(9,12,2,5,6,3,8,10,11,1,4,7)]
          # Write out to a file
          write.table(all,file=paste(workspace,'/',projectId,"_",cons[i],"_all_genes.txt",sep=""),sep="\t",row.names=F)
          listDEGs[[i]]=all
        }
        names(listDEGs)=cons
        norm_annotated <- merge(exprs(norm), Annot,by.x=0, by.y=0, all.x=T)               #write out normalized annotated data
        y<-paste("_",projectId, sep="")
        tNorm = tempfile(pattern = "normalized_data_", tmpdir =workspace, fileext = paste0(y,'.txt'))
        write.table(norm_annotated,file=tNorm,sep="\t",row.names=F)  
        for (i in 1:length(listDEGs)) {                                                   #Volcano plots
          dat=listDEGs[[i]]
          dat = dat[dat$SYMBOL!='NA',]
          log_FC=dat$logFC
          log_pval=-log10(dat$P.Value)
          Significant=rep("NotSignificant",length(log_FC))
          Significant[which(dat$P.Value<0.05 & abs(dat$logFC)>=1)]="AbsLogFoldChange>1 & PValue<0.05"
          Significant[which(dat$P.Value<0.05 & abs(dat$logFC)<1)]="PValue<0.05"
          Significant[which(dat$P.Value>=0.05 & abs(dat$logFC)>=1)]="AbsLogFoldChange>1"
          gene=dat$SYMBOL
          volcano_data=as.data.frame(cbind(gene,log_FC,log_pval,Significant))
          volcano_plot<-plot_ly(type='scatter', data = volcano_data, x = log_FC, y = log_pval, text = gene, mode = "markers", color = Significant) %>% layout(title=names(listDEGs)[i],xaxis=list(title="Fold Change",range =c(-5,5),tickvals=c(-5,-4,-3,-2,-1,0,1,2,3,4,5),ticktext=c('-32','-16','-8','-4','-2','1','2','4','8','16','32')),yaxis=list(title="-Log10 pvalue",range =c(0,15)))
          htmlwidgets::saveWidget(volcano_plot, paste0(workspace,"/volcano.html"))
        }
        return(list(listDEGs=listDEGs, norm_annotated=norm_annotated, pheno=pData(norm)))
      },
      error=function(cond) {
        return(outMsg)
      }
    ) 
    return(out)
  }
  return(catchAndReturn(outMsg))
  
  print("+++deg+++")
  sink(type='message')
}