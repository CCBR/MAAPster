library(shiny)
library(shinyjs)
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
#library(GSEA)
library(limma)
library(oligo)
library(gplots)
library(geneplotter)
library(multtest)
library(rgl)
library(rglwidget)
library(DT)
library(getopt)
library(annotate)
library(knitr)
library(reshape)
library(RColorBrewer)
library(mixOmics)
library(calibrate)
library(rmarkdown)
library(ggplot2)
library(ggfortify)


#setwd("/Users/valdezkm/Documents/MicroarrayPipeline/CodeInProgress/MicroArrayPipeline")
# 500 MB max upload size

options(shiny.maxRequestSize = 500*1024^2)

shinyServer(function(input, output) {
  
  url <- a("Manual", href="https://bioinformatics.cancer.gov/sites/default/files/course_material/microarray-pipeline-Btep-10032016.pptx")
  output$manu <- renderUI({
    tagList("Click to download:", url)
  })
  
  observeEvent(input$go, {
    # raw data
    raw=reactive(
      { 
        # system("rm *.txt")
        withProgress(message = 'Reading Raw data', detail = 'starting ...', value = 1/2, {
        # folder=path.expand(input$ProjectID)
        # dir.create(folder)
        # setwd(folder)
        myfiles=input$Indir
        # sort list by name
        # write.table(myfiles,"mycels.txt",sep="\t",row.names = F)
        myfiles=myfiles[order(myfiles$name),]
        write.table(myfiles,"mycels2.txt",sep="\t",row.names = F)
        cels = myfiles$datapath
        #print( cels )
        file1=input$pheno
        pd<-read.AnnotatedDataFrame(file1$datapath,header=TRUE,row.name="SampleName" ,sep="\t")
        celfiles <- read.celfiles(cels, phenoData=pd)
        # write.table(pData(celfiles),"celfiles2.txt",sep="\t",col.names=NA)
        # write.table(myfiles,"celfrominput.txt",sep="\t")
        cat(celfiles@annotation,file="annotation.txt")
        if (length(which(myfiles$name != rownames(pData(celfiles)))) > 0 ) {
          #cat("Please sort your phenotype on sample name and upload it again. \n")
          info("Please sort your phenotype on sample name and upload it again. Leaving...")
          stopApp(-1)
        }
        if (celfiles@annotation!="pd.hg.u133.plus.2" & celfiles@annotation!="pd.mogene.2.0.st" & celfiles@annotation!="pd.hugene.2.0.st" & celfiles@annotation!="pd.clariom.s.human.ht" & celfiles@annotation!="pd.clariom.s.human" & celfiles@annotation!="pd.clariom.s.mouse.ht" & celfiles@annotation!="pd.clariom.s.mouse") {
        #if (celfiles@annotation!="pd.hg.u133.plus.2" & celfiles@annotation!="pd.mogene.2.0.st") {
        #  if (celfiles@annotation!="pd.hg.u133.plus.2") {
               #cat("Please sort your phenotype on sample name and upload it again. \n")
          info(paste0("Affymetrix platform: ",celfiles@annotation," NOT supported. Leaving..."))
          stopApp(-1)
        }
        
        #validate(
        #  need(length(which(myfiles$name != rownames(pData(celfiles)))) == 0, "Please sort your phenotype on sample name and upload it again")
        #)
        celfiles
        })
      }
    )
    # norm data
    norm=reactive(
      {
        withProgress(message = 'Normalization', detail = 'starting ...', value = 1, {
        # if (input$Platform=="h133p2") {
        #if (raw()@annotation=="pd.hg.u133.plus.2") {
        if (raw()@annotation=="pd.hg.u133.plus.2" | raw()@annotation=="pd.clariom.s.human.ht" | raw()@annotation=="pd.clariom.s.human" | raw()@annotation=="pd.clariom.s.mouse.ht" | raw()@annotation=="pd.clariom.s.mouse") {
        celfiles.rma =rma(raw(), background=TRUE, normalize=TRUE, subset=NULL)
        } else {
          celfiles.rma =rma(raw(), background=TRUE, normalize=TRUE, subset=NULL, target="core")  
        }
        })
      }
    )
    # raw qc data
    qc=reactive(
      {
        withProgress(message = 'Fitting probe level model', detail = 'starting ...', value = 1, {
          celfiles.qc =fitProbeLevelModel(raw())
        
        })
         
      }
    )
    # list of DEG
    deg=reactive(
      {
         ##-------------
        withProgress(message = 'Computing differentially expressed genes', value = 0, {
          facs <- factor(pData(raw())$SampleGroup)
          labfacs=levels(facs)
          nbfacs=length(labfacs)
          file1=input$const
          contra=read.delim(file1$datapath)
          nb=dim(contra)[1]
          cons=c()
          for (k in 1:nb) {
            if ((contra[k,1] %in% labfacs) & (contra[k,2] %in% labfacs) )
            { 
              cons=c(cons,paste(contra[k,1],"-",contra[k,2],sep="")) 
            } else {
              cat("One of the groups in contrasts file at line :",k+1,"does not match a group in phenotype file..Quitting!!!\n")
              print( contra )
              stopApp(-1)
            }
          }
          
          myfactor <- factor(pData(norm())$SampleGroup)
          design1 <- model.matrix(~0+myfactor)
          colnames(design1) <- levels(myfactor)
          
          fit1 <- lmFit(norm(),design1)
          contrast.matrix <- makeContrasts(contrasts=cons,levels=design1)
          
          fit2 <- contrasts.fit(fit1, contrast.matrix)
          ebayes.fit2=eBayes(fit2) # smooths the std error
          incProgress(0.25, detail = 'Limma model fitted')
          # #EXTRACTING ALL GENES FOR EACH CONTRAST
          
          ##ANNOTATE PROBESET IDS FROM ANNOTATION PACKAGE FROM BIOCONDUCTOR
          ## load libraries as sources of annotation
          
          #library(mogene20sttranscriptcluster.db)
          #if (input$Platform=="mst2") {
          #if (raw()@annotation=="pd.mogene.2.0.st") {  
          #  Annot <- data.frame(ACCNUM=sapply(contents(mogene20sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mogene20sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(mogene20sttranscriptclusterGENENAME), paste, collapse=", "))
          #} else {
           # if (input$Platform=="h133p2") {
          #   if (raw()@annotation=="pd.hg.u133.plus.2") {
          #     Annot <- data.frame(ACCNUM=sapply(contents(hgu133plus2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133plus2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133plus2GENENAME), paste, collapse=", "))
          #  } 
          #} 
          
          if (raw()@annotation=="pd.mogene.2.0.st") {  
            Annot <- data.frame(ACCNUM=sapply(contents(mogene20sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mogene20sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(mogene20sttranscriptclusterGENENAME), paste, collapse=", "))
          } else {
           # if (input$Platform=="h133p2") {
             if (raw()@annotation=="pd.hg.u133.plus.2") {
               Annot <- data.frame(ACCNUM=sapply(contents(hgu133plus2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133plus2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133plus2GENENAME), paste, collapse=", "))
            } else {
              if (raw()@annotation=="pd.hugene.2.0.st") {
                Annot <- data.frame(ACCNUM=sapply(contents(hugene20sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hugene20sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hugene20sttranscriptclusterGENENAME), paste, collapse=", "))
              } else {
                if (raw()@annotation=="pd.clariom.s.human.ht") {
                  Annot <- data.frame(ACCNUM=sapply(contents(clariomshumanhttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(clariomshumanhttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(clariomshumanhttranscriptclusterGENENAME), paste, collapse=", "))
                } else {
                  if (raw()@annotation=="pd.clariom.s.mouse.ht") {
                    Annot <- data.frame(ACCNUM=sapply(contents(clariomsmousehttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(clariomsmousehttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(clariomsmousehttranscriptclusterGENENAME), paste, collapse=", "))
                  } else {
                    if (raw()@annotation=="pd.clariom.s.mouse") {
                      Annot <- data.frame(ACCNUM=sapply(contents(clariomsmousetranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(clariomsmousetranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(clariomsmousetranscriptclusterGENENAME), paste, collapse=", "))
                    } else {
                      if (raw()@annotation=="pd.clariom.s.human") {
                        Annot <- data.frame(ACCNUM=sapply(contents(clariomshumantranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(clariomshumantranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(clariomshumantranscriptclusterGENENAME), paste, collapse=", "))
                      }
                    }
                  }
                }
              }
            }
          }
          
          incProgress(0.25, detail = 'Preparing for pathway analysis')
          mylist=vector("list",nb)
             
          for (i in 1:nb)
          {
            
            all.genes.con = topTable(ebayes.fit2, coef = i, number=nrow(ebayes.fit2))
            
            # Merge data frames together (like a database table join)
            
            all <- merge(all.genes.con, Annot,by.x=0, by.y=0, all.x=T)
            all=all[order(all$P.Value),]
            colnames(all)[1]="probsetID"
          
            #L2P pathway starts here KV
            
            iup=which(all$P.Value<0.05 & all$logFC>=0)
            idw=which(all$P.Value<0.05 & all$logFC<0)
            fin.up=all[iup,]
  
            if (length(iup) > 500)
            {
              fin.up=fin.up[order(fin.up$P.Value),]
              fin.up=fin.up[1:500,]
            }
            #x2=rownames(fin.up)
            #gup=apply(array(as.character(x2)),1,function(z) unlist(strsplit(z, "\\|"))[2])
            
            fin.dw=all[idw,]
            if (length(idw) > 500)
            {
              fin.dw=fin.dw[order(fin.dw$P.Value),]
              fin.dw=fin.dw[1:500,]
            }
            #x2=rownames(fin.dw)
            #gdw=apply(array(as.character(x2)),1,function(z) unlist(strsplit(z, "\\|"))[2])
            
                 
            if (raw()@annotation=="pd.hg.u133.plus.2" | raw()@annotation=="pd.hugene.2.0.st" | raw()@annotation=="pd.clariom.s.human.ht" | raw()@annotation=="pd.clariom.s.human") 
            {
              cat(fin.up$SYMBOL,file=(paste0(input$ProjectID,'_',cons[i],'_Top500_Up.txt')), sep='\n')
              cat(fin.dw$SYMBOL,file=(paste0(input$ProjectID,'_',cons[i],'_Top500_Down.txt')),sep='\n')
            }
            else
            {
              cat(fin.up$SYMBOL,file=paste0(cons[i],"_Top500temp_Up.txt"),sep='\n')
              cat(fin.dw$SYMBOL,file=paste0(cons[i],"_Top500temp_Dw.txt"),sep='\n')
            
            system(paste0("cat ",cons[i],"_Top500temp_Up.txt | grep -v \"^NA\" | ./m2h | grep -v XXXX | cut -f2 -d\" \" >",input$ProjectID,'_',cons[i],"_Top500_Up.txt"))
            system(paste0("cat ",cons[i],"_Top500temp_Dw.txt | grep -v \"^NA\" | ./m2h | grep -v XXXX | cut -f2 -d\" \" >",input$ProjectID,'_',cons[i],"_Top500_Down.txt"))
            }
            system(paste0("cat ",input$ProjectID,'_',cons[i],"_Top500_Up.txt |sort | uniq | ./l2p >",input$ProjectID,'_',cons[i],"_Pathways_Up.txt"))
            system(paste0("cat ",input$ProjectID,'_',cons[i],"_Top500_Down.txt |sort | uniq | ./l2p >",input$ProjectID,'_',cons[i],"_Pathways_Down.txt"))
            
            addUpCol = read.delim(paste0(input$ProjectID,'_',cons[i],"_Pathways_Up.txt"), sep = '\t')
            addDwCol = read.delim(paste0(input$ProjectID,'_',cons[i],"_Pathways_Down.txt"), sep = '\t')
            
            colnames(addUpCol)=c("pval","fdr","ratio","nb.hits","nb.genes.path","nb.user.genes","tot.back.genes","path_id","source","description","type","gene.list")
            colnames(addDwCol)=c("pval","fdr","ratio","nb.hits","nb.genes.path","nb.user.genes","tot.back.genes","path_id","source","description","type","gene.list")
            addUpCol = addUpCol[order(addUpCol$pval),]
            addDwCol = addDwCol[order(addDwCol$pval),]
            addUpCol = addUpCol[,c(8,9,10,11,1,2,3,12,4,5,6,7)]
            addDwCol = addDwCol[,c(8,9,10,11,1,2,3,12,4,5,6,7)]
            write.table(addUpCol, file = paste0(input$ProjectID,'_',cons[i],"_Pathways_Up.txt"), sep = '\t', row.names = F)
            write.table(addDwCol, file = paste0(input$ProjectID,'_',cons[i],"_Pathways_Down.txt"), sep = '\t', row.names = F)
            
            # Write out to a file:
            all$FC = ifelse(all$logFC<0, -1/(2^all$logFC), 2^all$logFC)
           
            all = all[,c(9,1,8,10,11,2,5,6,3,4,7)]
            write.table(all,file=paste(input$ProjectID,"_",cons[i],"_all_genes.txt",sep=""),sep="\t",row.names=F)
            # cat("Contrast: ",i," done \n")
            
            mylist[[i]]=all
            ## end for
            
          }
          all <- merge(exprs(norm()), Annot,by.x=0, by.y=0, all.x=T)
          write.table(all,file=paste(input$ProjectID,"_normalized_data.txt",sep=""),sep="\t",row.names=F)
          #  
          names(mylist)=cons
          incProgress(0.5, detail = 'DEG done')
          mylist
          
        
          
        })
        ##-------------
      }
    )
  
    
    observeEvent(input$rep, {
       
      withProgress(message = 'Generating HTML report', detail = 'starting ...', value = 1, {
       # out <- render('../report_ver4.Rmd','html_document',paste(input$ProjectID,"_","report.html",sep=""),getwd(),getwd())
        out <- render('report_ver4.Rmd','html_document',paste(input$ProjectID,"_","report.html",sep=""))
      })
      
    })
    # processing all outputs
     output$projectid=renderText({paste("Project ID: ",input$ProjectID)})
     output$rawhist=renderPlot(
       {
         hist(raw(),which="all", main =" Raw Samples distribution")
       }
     )
     
    
     
     ## pca 2
     output$pca2d=renderPlot(
       {
         withProgress(message = 'Generating PCA', detail = 'starting ...', value = 1, {
           # myfactor <- factor(pData(norm())$SampleGroup)
           tedf= t(exprs(norm()))
       
           #removes zero  variances (issue with small sample sizes)
           if (length(which(apply(tedf, 2, var)==0)) >= 0){
             tedf = tedf[ , apply(tedf, 2, var) != 0]
           }
           
           rownames(tedf)=pData(norm())$SampleID
           # tedf1 = data.frame(tedf)
           pr1=prcomp(tedf,scale.=T)
           ff <- factor(pData(norm())$SampleGroup)
           dd=cbind(tedf,group=as.character(ff))
           
           pc1.var=100*round(((pr1$sdev)**2)[1]/sum((pr1$sdev)**2),digits=2) # %var pc1 
           pc2.var=100*round(((pr1$sdev)**2)[2]/sum((pr1$sdev)**2),digits=2) # % var pc2
           
           xLab=paste("PC1 - ",pc1.var," % of variation",sep="")
           yLab=paste("PC2 - ",pc2.var," % of variation",sep="")
           
           autoplot(pr1,data=dd, colour = 'group', label = T, xlab=xLab, ylab=yLab)  
         })
       }
       
     )
     ##
     ## pca 2
#     output$pca2d=renderPlot(
#       {
#         # myfactor <- factor(pData(norm())$SampleGroup)
#         tedf= t(exprs(norm()))
         # tedf1 = data.frame(tedf)
#         pr1=prcomp(tedf,scale.=T)
#         igp1=which(pData(norm())$SampleGroup==labfacs2[1])
#         xrange=range(pr1$x[,1])
#         yrange=range(pr1$x[,2])
         
#         pc1.var=100*round(((pr1$sdev)**2)[1]/sum((pr1$sdev)**2),digits=2) # %var pc1 
#         pc2.var=100*round(((pr1$sdev)**2)[2]/sum((pr1$sdev)**2),digits=2) # % var pc2
         
#         xlab=paste("PC1 - ",pc1.var," % of variation",sep="")
#         ylab=paste("PC2 - ",pc2.var," % of variation",sep="")
#         plot(pr1$x[igp1,1],pr1$x[igp1,2],col=1,main="2D PCA",pch=19,ylab=ylab,xlab=xlab,xlim=xrange,ylim=yrange)
#         for (i in 2:nbfacs2)
#         {
#           igp=which(pData(norm())$SampleGroup==labfacs2[i])
#           points(pr1$x[igp,1],pr1$x[igp,2],col=i,pch=19)
#           
#         }
#         legend('bottom', labfacs2,lty=1, col=1:nbfacs2, bty='n', cex=.8)
#       }
#     )
     ##
     output$rawbox=renderPlot(
       {
         boxplot(raw(), col=c(3,2),which="all", main="Boxplots before normalization",las=2,names=pData(raw())$SampleID)
        }
     )
     output$rle=renderPlot(
       {
         RLE(qc(), main="RLE plot",names=pData(raw())$SampleID)
       }
     )
     output$nuse=renderPlot(
       {
         NUSE(qc(), main="NUSE plot",names=pData(raw())$SampleID)
       }
     )
     output$rawmaplot=renderUI({
      #withProgress(message = 'Generating Raw Maplot', detail = 'starting ...', value = 0, {
       facs <- factor(pData(raw())$SampleGroup)
       labfacs=levels(facs)
       nbfacs=length(labfacs)
       plot_output_list <- lapply(1:nbfacs, function(i) {
         plotname <- paste("plot", i, sep="")
         plotOutput(plotname, height = 800, width = 800)
     })
       # Convert the list to a tagList - this is necessary for the list of items
       # to display properly.
       do.call(tagList, plot_output_list)
       #})
     })
     
     
         facs <- factor(pData(raw())$SampleGroup)
         labfacs=levels(facs)
         nbfacs=length(labfacs)
         #par(mfrow=c(1,nbfacs))
         
         for (i in 1:nbfacs) {
           local({
             my_i <- i
             plotname <- paste("plot", my_i, sep="")
           # MA plots are then used to visualize intensity-dependent ratio for each group
           igp=which(pData(raw())$SampleGroup==labfacs[i])
           output[[plotname]] <- renderPlot({
              withProgress(message = 'Generating Raw Maplot', detail = paste0('Plot ', my_i, ' starting ...'), value = (my_i/nbfacs), {
                MAplot(raw()[,igp],pairs=TRUE,plotFun=smoothScatter,main="MVA plot before normalization", labels=raw()[,igp]$SampleID)  
           })
           })
         })
  }
     
         
       ## MVAplot after normalization
         output$normaplot=renderUI({
           # withProgress(message = 'Generating Normalized Maplot', detail = 'starting ...', value = 0, {
           facs2 <- factor(pData(norm())$SampleGroup)
           labfacs2=levels(facs2)
           nbfacs2=length(labfacs2)
           plot_output_list2 <- lapply(1:nbfacs2, function(i) {
             plotname2 <- paste("plota", i, sep="")
             plotOutput(plotname2, height = 800, width = 800)
           })
           # Convert the list to a tagList - this is necessary for the list of items
           # to display properly.
           do.call(tagList, plot_output_list2)
          # })
         })
         
         
         facs2 <- factor(pData(norm())$SampleGroup)
         labfacs2=levels(facs2)
         nbfacs2=length(labfacs2)
         #par(mfrow=c(1,nbfacs))
         for (i in 1:nbfacs2) {
           local({
             my_i <- i
             plotname2 <- paste("plota", my_i, sep="")
             # MA plots are then used to visualize intensity-dependent ratio for each group
             igp=which(pData(norm())$SampleGroup==labfacs2[i])
             output[[plotname2]] <- renderPlot({
               withProgress(message = 'Generating Raw Maplot', detail = paste0('Plot ', my_i, ' starting ...'), value = (my_i/nbfacs2), {
               MAplot(norm()[,igp],pairs=TRUE,plotFun=smoothScatter,main="MVA plot after RMA Normalization", labels=norm()[,igp]$SampleID) # 
             })
             })
           })
         }
          
         
       ## end mvaplat after normalization
     output$rmabox=renderPlot(
       {
         boxplot(norm(),col=c(3,2), main="Boxplots after RMA normalization",las=2,names=pData(norm())$SampleID)
       }
     )
     output$heatmap=renderPlot(
       {
         # no std
         mat=as.matrix(dist(t(exprs(norm()))))
         rownames(mat)=pData(norm())$SampleID
         colnames(mat)=rownames(mat)
         heatmap.2(mat, trace="none", margin=c(10,10))
       }
     )
     output$rmahist=renderPlot(
       {
         hist(norm(), main ="Distribution after Normalization")
       }
     )
     output$deg=DT::renderDataTable(DT::datatable(
       {
        dat = deg()[[input$NumContrasts]]
        dat = dat[,-6]
        
        if (is.na(input$pval) & is.na(input$fc)) {   
          dat
        } else if (is.na(input$pval))  {
          dat = dat[(abs(as.numeric(dat[,5])) >= input$fc),]
        } else if (is.na(input$fc)) {
          dat = dat[(as.numeric(dat[,6]) <= input$pval),]
        } else {
          dat = dat[(as.numeric(dat[,6]) <= input$pval & abs(as.numeric(dat[,5])) >= input$fc),]
          dat
        }
        # deg()[[1]]
       }, caption =paste0("contrast: ",names(deg())[input$NumContrasts])
     )
     )
     
     output$topUp=DT::renderDataTable(DT::datatable(
       {
        callDEG = deg()[[input$NumContrasts]]
        topUp = read.delim(paste0(input$ProjectID,'_',names(deg())[input$NumContrasts],"_Pathways_Up.txt"), sep = '\t', header = T)
        
        if (!is.na(input$pathPval)) {
        topUp = topUp[(as.numeric(topUp[,5]) <= input$pathPval),]
        } else {
        topUp
        }
       } , caption=paste0("Pathways for the top 500 Upregulated Genes: ", names(deg())[input$NumContrasts]),
            options = list(columnDefs = list(list(targets = 8, 
               render = JS("function(data, type, row, meta) {",
                   "return type === 'display' && data.length > 30 ?",
                   "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
                   "}")))), callback = JS('table.page(3).draw(false);')
     )
     )
     output$topDown=DT::renderDataTable(DT::datatable(
       {
        callDEG = deg()[[input$NumContrasts]]
        topDw = read.delim(paste0(input$ProjectID,'_',names(deg())[input$NumContrasts],"_Pathways_Down.txt"), sep = '\t', header = T)
        
        if (!is.na(input$pathPval)) {
        topDw = topDw[(as.numeric(topDw[,5]) <= input$pathPval),]
        } else {
          topDw
        }
       } , caption=paste0("Pathways for the top 500 Downregulated Genes: ", names(deg())[input$NumContrasts]),
            options = list(columnDefs = list(list(targets = 8, 
              render = JS("function(data, type, row, meta) {",
                    "return type === 'display' && data.length > 30 ?",
                    "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
                    "}")))), callback = JS('table.page(3).draw(false);')
     )
     )

       
     
    
     
     
#     output$kegg=DT::renderDataTable(DT::datatable(
#       {
         ##
         
#         dat=deg()[[input$NumContrasts]]
         ## Hypergeometric Tests
         ## P value cutoff
#         dat.i <- which(umeric(dat[,5]) < input$pval & abs(as.numeric(dat[,2])) >= input$fc)
#         dat.s <- dat[dat.i,]
         ## get gene list
#         gene  <- as.character(dat.s[,9])
         # if (input$Platform=="mst2") {
#         if (raw()@annotation=="pd.mogene.2.0.st" | raw()@annotation=="pd.clariom.s.mouse.ht" | raw()@annotation=="pd.clariom.s.mouse") {  
#           data(KEGG_mm); mykegg=unique(unlist(GSEA.db$pathway.list)); ngene=intersect(gene,mykegg)
#           xx=doGSEA(db="KEGG_mm", gene=ngene, filter.num=0, fdr=T)
           
#         } else {
           # if (input$Platform=="h133p2") {
#           if (raw()@annotation=="pd.hg.u133.plus.2" | raw()@annotation=="pd.hugene.2.0.st" | raw()@annotation=="pd.clariom.s.human.ht" | raw()@annotation=="pd.clariom.s.human") {  
#             data(KEGG); mykegg=unique(unlist(GSEA.db$pathway.list)); ngene=intersect(gene,mykegg)
#             xx=doGSEA(db="KEGG", gene=ngene, filter.num=0, fdr=T)
#           }
#         }
#         write.table(xx,file=paste(input$ProjectID,"_kegg_results_contrast_",input$NumContrasts,"_",input$pval,"_pvalue_",input$fc,"_logFC.txt",sep=""),sep="\t",row.names=F)
#         xx 
         ##
     # } , caption =paste0("KEGG for contrast: ",names(deg())[input$NumContrasts])
     #)
     #)

#     output$go=DT::renderDataTable(DT::datatable(
#       {
         ##
         
#         dat=deg()[[input$NumContrasts]]
         ## Hypergeometric Tests
         ## P value cutoff
#         dat.i <- which(as.numeric(dat[,5]) < input$pval & abs(as.numeric(dat[,2])) >= input$fc)
#         dat.s <- dat[dat.i,]
         ## get gene list
#         gene  <- as.character(dat.s[,9])
         # doGSEA(db="GO_mm", gene=gene, filter.num=0, fdr=T);
         # if (input$Platform=="mst2") {
#         if (raw()@annotation=="pd.mogene.2.0.st" | raw()@annotation=="pd.clariom.s.mouse.ht" | raw()@annotation=="pd.clariom.s.mouse") {
#           data(GO_mm); mygo=unique(unlist(GSEA.db$pathway.list)); ngene=intersect(gene,mygo)
#           yy=doGSEA(db="GO_mm", gene=ngene, filter.num=0, fdr=T)
#         } else {
           # if (input$Platform=="h133p2") {
#           if (raw()@annotation=="pd.hg.u133.plus.2" | raw()@annotation=="pd.hugene.2.0.st" | raw()@annotation=="pd.clariom.s.human.ht" | raw()@annotation=="pd.clariom.s.human") {
#             data(GO); mygo=unique(unlist(GSEA.db$pathway.list)); ngene=intersect(gene,mygo) 
#             yy=doGSEA(db="GO", gene=ngene, filter.num=0, fdr=T)
#           }
#         };
         ##
#         write.table(yy,file=paste(input$ProjectID,"_go_results_contrast_",input$NumContrasts,"_",input$pval,"_pvalue_",input$fc,"_logFC.txt",sep=""),sep="\t",row.names=F)
#         yy
#       }  , caption =paste0("GO for contrast: ",names(deg())[input$NumContrasts])
#     )
#     )
     
     

     
     ##
     output$downloadTables <- downloadHandler(
       filename = function() {paste(input$ProjectID,"_tables.zip",sep="")},
       content = function(file) {
          # mytables=list.files(pattern="\\.txt$")
          mytables=list.files(pattern=paste(input$ProjectID,".*.txt",sep=""))
          zip(file, mytables, flags = "-r9X", extras = "", zip = Sys.getenv("R_ZIPCMD", "zip"))
         
       }
     )
     ##
    # output$downloadReport <- downloadHandler(
     #  filename = 'my-report.html',
    #   content = function(file) {
         # src <- normalizePath('MouseTestData_report_ver2.rmd')
         # src <- normalizePath('report_ver4.Rmd')
      #   out <- normalizePath('report_ver4.html')
         # temporarily switch to the temp dir, in case you do not have write
         # permission to the current working directory
         #
         # owd <- setwd(tempdir())
         # on.exit(setwd(owd))
         # file.copy(src, 'report_ver4.Rmd',overwrite = T)
         
         # library(rmarkdown)
         # out <- render('report_ver4.Rmd','html_document') #, output_file = "Report.html", output_dir = normalizePath(input$Indir))
         # out <- render('MouseTestData_report_ver4.Rmd',html_document())
         # out <- render('report.rmd') #, output_file = "ReportSimple.html", output_dir = normalizePath(input$Indir))
     #    file.rename(out, file)
    #   }
    # )
     ##
     output$downloadReport <- downloadHandler(
       filename = function() {paste(input$ProjectID,"_Report.zip",sep="")},
       content = function(file) {
         #mytables=list.files(pattern="\\.html$")
         mytables=list.files(pattern=paste(input$ProjectID,".*.html",sep=""))
         
         
         zip(file, mytables, flags = "-r9X", extras = "", zip = Sys.getenv("R_ZIPCMD", "zip"))
         
       }
     )
     
  })
  
})



