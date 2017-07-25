library(shiny)
library(shinyjs)
library(pd.mogene.2.0.st)
library(mogene20sttranscriptcluster.db)
library(pd.hg.u133.plus.2)
library(hgu133plus2.db)
library(GSEA)
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

# setwd("/Users/elloumif/Documents/shiny")
# 500 MB max upload size
options(shiny.maxRequestSize = 500*1024^2)

shinyServer(function(input, output) {
  observeEvent(input$go, {
    # raw data
    raw=reactive(
      { 
        # system("rm *.txt")
        withProgress(message = 'Reading Raw data', detail = 'starting ...', value = 0, {
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
        if (celfiles@annotation!="pd.hg.u133.plus.2" & celfiles@annotation!="pd.mogene.2.0.st") {
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
        withProgress(message = 'Normalization', detail = 'starting ...', value = 0, {
        # if (input$Platform=="h133p2") {
        if (raw()@annotation=="pd.hg.u133.plus.2") {
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
        withProgress(message = 'Fitting probe level model', detail = 'starting ...', value = 0, {
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
            }
            else 
            {
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
          if (raw()@annotation=="pd.mogene.2.0.st") {  
            Annot <- data.frame(ACCNUM=sapply(contents(mogene20sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mogene20sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(mogene20sttranscriptclusterGENENAME), paste, collapse=", "))
          } else {
           # if (input$Platform=="h133p2") {
             if (raw()@annotation=="pd.hg.u133.plus.2") {
               Annot <- data.frame(ACCNUM=sapply(contents(hgu133plus2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133plus2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133plus2GENENAME), paste, collapse=", "))
            } 
          } 
          incProgress(0.25, detail = 'Annotation done')
          mylist=vector("list",nb)
          
          for (i in 1:nb)
          {
            
            all.genes.con = topTable(ebayes.fit2, coef = i, number=nrow(ebayes.fit2))
            
            
            
            # Merge data frames together (like a database table join)
            
            all <- merge(all.genes.con, Annot,by.x=0, by.y=0, all.x=T)
            all=all[order(all$P.Value),]
            colnames(all)[1]="probsetID"
            
            # Write out to a file:
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
       
      withProgress(message = 'Generating HTML report', detail = 'starting ...', value = 0, {
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
     output$manu=renderUI(
       {
         "Manual coming soon"
         
         
       }
     )
     ## pca 2
     output$pca2d=renderPlot(
       {
         withProgress(message = 'Generating PCA', detail = 'starting ...', value = 0, {
         # myfactor <- factor(pData(norm())$SampleGroup)
         tedf= t(exprs(norm()))
         rownames(tedf)=pData(norm())$SampleID
         # tedf1 = data.frame(tedf)
         pr1=prcomp(tedf,scale.=T)
         ff <- factor(pData(norm())$SampleGroup)
         dd=cbind(tedf,group=as.character(ff))
         autoplot(pr1,data=dd, colour = 'group', label = T)  
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
           MAplot(raw()[,igp],pairs=TRUE,plotFun=smoothScatter,main="MVA plot before normalization", labels=raw()[,igp]$SampleID) # 
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
               MAplot(norm()[,igp],pairs=TRUE,plotFun=smoothScatter,main="MVA plot after RMA Normalization", labels=norm()[,igp]$SampleID) # 
            
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
       #{
         ##
        deg()[[input$NumContrasts]] , caption =paste0("contrast: ",names(deg())[input$NumContrasts])
         
        # deg()[[1]]
         ##
       #}
     )
     )
     output$kegg=DT::renderDataTable(DT::datatable(
       {
         ##
         
         dat=deg()[[input$NumContrasts]]
         ## Hypergeometric Tests
         ## P value cutoff
         dat.i <- which(as.numeric(dat[,5]) < input$pval & abs(as.numeric(dat[,2])) >= input$fc)
         dat.s <- dat[dat.i,]
         ## get gene list
         gene  <- as.character(dat.s[,9])
         # if (input$Platform=="mst2") {
         if (raw()@annotation=="pd.mogene.2.0.st") {  
         data(KEGG_mm); mykegg=unique(unlist(GSEA.db$pathway.list)); ngene=intersect(gene,mykegg)
         xx=doGSEA(db="KEGG_mm", gene=ngene, filter.num=0, fdr=T)
         
         } else {
           # if (input$Platform=="h133p2") {
           if (raw()@annotation=="pd.hg.u133.plus.2") {  
             data(KEGG); mykegg=unique(unlist(GSEA.db$pathway.list)); ngene=intersect(gene,mykegg)
             xx=doGSEA(db="KEGG", gene=ngene, filter.num=0, fdr=T)
           }
         }
         write.table(xx,file=paste(input$ProjectID,"_kegg_results_contrast_",input$NumContrasts,"_",input$pval,"_pvalue_",input$fc,"_logFC.txt",sep=""),sep="\t",row.names=F)
         xx 
         ##
       } , caption =paste0("KEGG for contrast: ",names(deg())[input$NumContrasts])
     )
     )
     
     output$go=DT::renderDataTable(DT::datatable(
       {
         ##
         
         dat=deg()[[input$NumContrasts]]
         ## Hypergeometric Tests
         ## P value cutoff
         dat.i <- which(as.numeric(dat[,5]) < input$pval & abs(as.numeric(dat[,2])) >= input$fc)
         dat.s <- dat[dat.i,]
         ## get gene list
         gene  <- as.character(dat.s[,9])
         # doGSEA(db="GO_mm", gene=gene, filter.num=0, fdr=T);
         # if (input$Platform=="mst2") {
         if (raw()@annotation=="pd.mogene.2.0.st") {
           data(GO_mm); mygo=unique(unlist(GSEA.db$pathway.list)); ngene=intersect(gene,mygo)
           yy=doGSEA(db="GO_mm", gene=ngene, filter.num=0, fdr=T)
         } else {
           # if (input$Platform=="h133p2") {
             if (raw()@annotation=="pd.hg.u133.plus.2") {
               data(GO); mygo=unique(unlist(GSEA.db$pathway.list)); ngene=intersect(gene,mygo) 
             yy=doGSEA(db="GO", gene=ngene, filter.num=0, fdr=T)
           }
         };
         ##
         write.table(yy,file=paste(input$ProjectID,"_go_results_contrast_",input$NumContrasts,"_",input$pval,"_pvalue_",input$fc,"_logFC.txt",sep=""),sep="\t",row.names=F)
         yy
       }  , caption =paste0("GO for contrast: ",names(deg())[input$NumContrasts])
     )
     )
     
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



