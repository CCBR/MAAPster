library(shiny)
library(shinyjs)
library(GEOquery)
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
library(shinyRGL)
library(plotly)
library(htmltools)
library(heatmaply)
library(Biobase)
library(GSVA)
library(GSEABase)
library(pheatmap)
library(viridis)
library(dendsort)

#setwd('/Users/valdezkm/Documents/MicroarrayPipeline/CodeInProgress/MicroArrayPipeline')

shinyServer(function(input, output) {
  
  v <- reactiveValues(data = NULL, platform=NULL)
  k <- reactiveValues()

  observeEvent(
    input$button, 
    isolate({
      shinyjs::show("hide1")
    })
  )
  observeEvent(
    input$button,
    isolate({
      shinyjs::show("hide2")
    })
  )
  observeEvent(
    input$CELbutton,
    isolate({
      shinyjs::show("hide1")
    })
  )
  observeEvent(
    input$CELbutton,
    isolate({
      shinyjs::show("hideCEL")
    })
  )
  observeEvent(
    input$button2,
    isolate({
      shinyjs:: show("hide3")
    })
  )
  observeEvent(
    input$addrow,
    isolate({
      shinyjs:: show("hideAnalysis")
    })
  )
  observeEvent(
    input$analyze,
    isolate({
      shinyjs:: show("hideResults")
    })
  )
  observeEvent(
    input$test,
    isolate({
      shinyjs:: show("hideContrast")
    })
  )
  observeEvent(
    input$analyze,
    isolate({
      shinyjs:: show("hideDownloads")
    })
  )
  
  
  observeEvent(
    input$CELbutton,
    isolate({
      #CEL file input
      myfiles=input$Indir
      cels = myfiles$name
    
      mytableCEL=matrix("",length(cels),1)
      colnames(mytableCEL)=c("file name")

      for (k in 1:length(cels))
      {
        mytableCEL[k,]<-c(cels[k])
      }
      mytableCEL <- data.frame(mytableCEL)
      
      mytableCEL$group <- ""
      v$data <- mytableCEL
      
      output$mytableCEL = DT::renderDataTable({
        if (is.null(v$data)) return()
        if (is.null(v$platform)) warning()
        DT::datatable(v$data, options = list(lengthMenu = c(2,4,6,8), pageLength = 8))
    })
    })
  )
 
  observeEvent(
    input$button, 
    isolate({
      withProgress(message = 'Loading files...', value = 0.25, {
        if (input$gseid=="8 digit GSE code") return()
        id=input$gseid
        id = gsub(" ","",id,fixed=TRUE)       
        
        incProgress(0.25)
        

        #gds <- getGEO(input$gseid, GSEMatrix = F,getGPL=T,AnnotGPL=T)
        gds <- getGEO(id, GSEMatrix = F,getGPL=T,AnnotGPL=T)
        
        gds
        
        mytable=matrix("",length(GSMList(gds)),3)
        colnames(mytable)=c("gsm","title","description")
        
        
        for (k in 1:length(GSMList(gds)))
        {
          if (is.null(Meta(GSMList(gds)[[k]])$description)) {     #error handling, some records don't have descriptions
            mytable[k,] <-c(Meta(GSMList(gds)[[k]])$geo_accession[1], Meta(GSMList(gds)[[k]])$title[1], 'No data available')
          } else {
            mytable[k,] <-c(Meta(GSMList(gds)[[k]])$geo_accession[1], Meta(GSMList(gds)[[k]])$title[1], Meta(GSMList(gds)[[k]])$description[1])
          }
        }
        # }
        
        mytable <- data.frame(mytable)
        
        mytable$group <- ""
        v$data <- mytable
        
        incProgress(0.25)
        output$mytable = DT::renderDataTable({
          if (is.null(v$data)) return()
          if (is.null(v$platform)) warning()
          DT::datatable(v$data, options = list(lengthMenu = c(2,4,6,8), pageLength = 8))
      })
      })
    })
  )
  
  observeEvent(
    input$button2,
    isolate({
      
      if (is.null(v$data)) return()
      if (is.null(input$group1)) return()
      
      #print(v$data)
      #print( "Error??2")

      #error handling names that conflict with DEG script
      tempGroup = input$group1
      tempGroup = make.names(tempGroup)
      
      v$data[input$mytable_rows_selected, "group" ] <- tempGroup
      v$data[input$mytableCEL_rows_selected, "group"] <- tempGroup
      
      #print( "Error??3")
    
      output$mytableCEL = DT::renderDataTable({
        
        if (is.null(v$data)) return()
        if (is.null(v$platform)) warning()
        DT::datatable(v$data, options = list(lengthMenu = c(2,4,6,8,10), pageLength = 8))
      })
      
      output$mytable = DT::renderDataTable({
        
        if (is.null(v$data)) return()
        if (is.null(v$platform)) warning()
        DT::datatable(v$data, options = list(lengthMenu = c(2,4,6,8,10), pageLength = 8))
      })
    })
  )
  
  observe({
    isolate({
      output$ui <- renderUI({
        # if (is.null(input$number))
        #   return()
            textInput("group1", "Please enter a group name (with no spaces)")
      })
    })
  })
  
  observe({
    isolate({
      output$choice1 <- renderUI({
        g = unique(v$data$group)
        selectizeInput("selectIn1","Select:",choices = g)
      })
    })
  })
  
  observe({
    isolate({
      output$choice2 <- renderUI({
        g = unique(v$data$group)
        selectizeInput("selectIn2","Versus Reference Group:",choices = g)
      })
    })
  })
  
  # observeEvent(input$number, {
  #              isolate({
  #                output$ui <- renderUI({
  #                  if (is.null(input$number))
  #                    return()
  # 
  #                  switch(input$number,
  #                         "1" = selectInput("group1", "Please select a group",
  #                                           choices = c("Group_1" = "Group_1"),
  #                                           selected = "Group_1"
  #                         ),
  #                         "2" = selectInput("group1", "Please select a group",
  #                                           choices = c("Group_1" = "Group_1",
  #                                                       "Group_2" = "Group_2"),
  #                                           selected = "Group_1"
  #                         ),
  #                         "3" = selectInput("group1", "Please select a group",
  #                                           choices = c("Group_1" = "Group_1",
  #                                                       "Group_2" = "Group_2",
  #                                                       "Group_3" = "Group_3"),
  #                                           selected = "Group_1"
  #                         ),
  #                         "4" = selectInput("group1", "Please select a group",
  #                                           choices = c("Group_1" = "Group_1",
  #                                                       "Group_2" = "Group_2",
  #                                                       "Group_3" = "Group_3",
  #                                                       "Group_4" = "Group_4"),
  #                                           selected = "Group_1"
  #                         ),
  #                         "5" = selectInput("group1", "Please select a group",
  #                                           choices = c("Group_1" = "Group_1",
  #                                                       "Group_2" = "Group_2",
  #                                                       "Group_3" = "Group_3",
  #                                                       "Group_4" = "Group_4",
  #                                                       "Group_5" = "Group_5"),
  #                                           selected = "Group_1"
  #                         ),
  #                         "6" = selectInput("group1", "Please select a group",
  #                                           choices = c("Group_1" = "Group_1",
  #                                                       "Group_2" = "Group_2",
  #                                                       "Group_3" = "Group_3",
  #                                                       "Group_4" = "Group_4",
  #                                                       "Group_5" = "Group_5",
  #                                                       "Group_6" = "Group_6"),
  #                                           selected = "Group_1"
  #                         ),
  #                         "7" = selectInput("group1", "Please select a group",
  #                                           choices = c("Group_1" = "Group_1",
  #                                                       "Group_2" = "Group_2",
  #                                                       "Group_3" = "Group_3",
  #                                                       "Group_4" = "Group_4",
  #                                                       "Group_5" = "Group_5",
  #                                                       "Group_6" = "Group_6",
  #                                                       "Group_7" = "Group_7"),
  #                                           selected = "Group_1"
  #                         ),
  #                         "8" = selectInput("group1", "Please select a group",
  #                                           choices = c("Group_1" = "Group_1",
  #                                                       "Group_2" = "Group_2",
  #                                                       "Group_3" = "Group_3",
  #                                                       "Group_4" = "Group_4",
  #                                                       "Group_5" = "Group_5",
  #                                                       "Group_6" = "Group_6",
  #                                                       "Group_7" = "Group_7",
  #                                                       "Group_8" = "Group_8"),
  #                                           selected = "Group_1"
  #                         ),
  #                         "9" = selectInput("group1", "Please select a group",
  #                                           choices = c("Group_1" = "Group_1",
  #                                                       "Group_2" = "Group_2",
  #                                                       "Group_3" = "Group_3",
  #                                                       "Group_4" = "Group_4",
  #                                                       "Group_5" = "Group_5",
  #                                                       "Group_6" = "Group_6",
  #                                                       "Group_7" = "Group_7",
  #                                                       "Group_8" = "Group_8",
  #                                                       "Group_9" = "Group_9"),
  #                                           selected = "Group_1"
  #                         ),
  #                         "10" = selectInput("group1", "Please select a group",
  #                                           choices = c("Group_1" = "Group_1",
  #                                                       "Group_2" = "Group_2",
  #                                                       "Group_3" = "Group_3",
  #                                                       "Group_4" = "Group_4",
  #                                                       "Group_5" = "Group_5",
  #                                                       "Group_6" = "Group_6",
  #                                                       "Group_7" = "Group_7",
  #                                                       "Group_8" = "Group_8",
  #                                                       "Group_9" = "Group_9",
  #                                                       "Group_10" = "Group_10"),
  #                                           selected = "Group_1"
  #                         )
  #                  )
  #                })
  #              })
  #        })

  observeEvent(input$addrow,{
    isolate({
      newLine1 = c(input$selectIn1)
      newLine2 = c(input$selectIn2)
      
      k$k1 <- rbind(k$k1, newLine1)
      k$k2 <- rbind(k$k2, newLine2) 
      })
    })
  output$contrastTable <- renderTable({
    data.frame(paste0(k$k1,'  vs  ', k$k2))
  }, colnames = FALSE)
  
  observe({
    isolate({
      output$displayContrast <- renderUI({
        selectizeInput("NumContrasts", label=h6("Choose contrast to show"),choices=paste0(k$k1, ' vs ', k$k2), width="150px")
      })
    })
  })
  
  
  
  observeEvent(
    input$analyze, {
      isolate({
      raw <- reactive(
        {
      withProgress(message = 'Loading files...', value = 0.25, {
        
      id=input$gseid
        
      if (id=='8 digit GSE code'){
        myfiles = input$Indir
        cels = myfiles$datapath
        Pheno = v$data
        row.names(Pheno) = Pheno$file.name
        SampleName = myfiles$name
        pd = AnnotatedDataFrame(Pheno)
        celfiles = read.celfiles(cels, phenoData = pd)
        colnames(pData(celfiles))[1] = 'SampleID'  
        
      } else {
        
      id = gsub(" ","",id,fixed=TRUE) 
      system(paste0('rm *.CEL.gz'))        #removes previous CEL files
      getGEOSuppFiles(id, makeDirectory = T, baseDir = getwd())
      fileID = paste0(id, '_RAW.tar')
      #system(paste0('tar -xvf', fileID))
      untar(paste0(getwd(),'/',id,'/',fileID))
      incProgress(0.25)
      
      #cels = paste0(Pheno$gsm,'_',Pheno$title,'.CEL.gz')   #adds filename
      Pheno = v$data
      SampleName = list.files(pattern = '/*CEL.gz')    #list contents of new directory with zipped CEL files
      
      if (length(grep('*CEL*',SampleName,ignore.case = T)) == 0) {
        info("Raw files must be CEL files")
      }
      rownames(Pheno) = Pheno$title
      cels = SampleName
      
      incProgress(0.25)
      
      pd = AnnotatedDataFrame(Pheno)
      celfiles = read.celfiles(cels, phenoData = pd)
      colnames(pData(celfiles))[2] = 'SampleID'    
      }
      
      cat(celfiles@annotation,file="annotation.txt")
      
      if (celfiles@annotation!="pd.hg.u133.plus.2" & celfiles@annotation!="pd.mogene.2.0.st" & celfiles@annotation!="pd.hugene.2.0.st" & celfiles@annotation!="pd.clariom.s.human.ht" & celfiles@annotation!="pd.clariom.s.human" & celfiles@annotation!="pd.clariom.s.mouse.ht" & celfiles@annotation!="pd.clariom.s.mouse" & celfiles@annotation!='pd.mouse430.2' & celfiles@annotation!='pd.hg.u133a' & celfiles@annotation!='pd.hugene.1.0.st.v1' & celfiles@annotation!='pd.mogene.1.0.st.v1' & celfiles@annotation!='pd.hg.u133a.2' & celfiles@annotation!='pd.huex.1.0.st.v2' & celfiles@annotation!='pd.hg.u219' & celfiles@annotation!='pd.mg.u74av2' & celfiles@annotation!='pd.mouse430a.2' & celfiles@annotation!='pd.moe430a' & celfiles@annotation!='pd.hg.u95av2' & celfiles@annotation!='pd.hta.2.0' & celfiles@annotation!='pd.moex.1.0.st.v1' & celfiles@annotation!='pd.hg.u133b' & celfiles@annotation!='pd.hugene.1.1.st.v1') {
        #cat("Please sort your phenotype on sample name and upload it again. \n")
        info(paste0("Affymetrix platform: ",celfiles@annotation," NOT supported. Leaving..."))
        stopApp(-1)
      }
      incProgress(0.25)
      
      celfiles
      })
})
      
      norm=reactive(
        {
          withProgress(message = 'Normalization', detail = 'starting ...', value = 0, {
            if (raw()@annotation=="pd.hg.u133.plus.2" | raw()@annotation=="pd.clariom.s.human.ht" | raw()@annotation=="pd.clariom.s.human" | raw()@annotation=="pd.clariom.s.mouse.ht" | raw()@annotation=="pd.clariom.s.mouse" | raw()@annotation=='pd.mouse430.2' | raw()@annotation=='pd.hg.u133a' | raw()@annotation=='pd.hg.u133a.2' | raw()@annotation=='pd.hg.u219' | raw()@annotation=='pd.mg.u74av2' | raw()@annotation=='pd.mouse430a.2' | raw()@annotation=='pd.moe430a' | raw()@annotation=='pd.hg.u95av2' | raw()@annotation=='pd.hg.u133b') {
              incProgress(0.5)
              celfiles.rma =rma(raw(), background=TRUE, normalize=TRUE, subset=NULL)
            } else {
              incProgress(0.5)
              celfiles.rma =rma(raw(), background=TRUE, normalize=TRUE, subset=NULL, target="core")
            }
          })
        })
      qc=reactive(
        {
          withProgress(message = 'Fitting probe level model', detail = 'starting ...', value = 0.75, {
            validate(
              need(raw()@annotation!= "pd.mogene.1.0.st.v1", 'NUSE and RLE plots unavailable for this platform.')
            )
            celfiles.qc=fitProbeLevelModel(raw())
          })
        }
      )
      # list of DEG
      deg=reactive(
        {
          ##-------------
          withProgress(message = 'Analysis starting...', value = 0, {
            facs <- factor(pData(raw())$group)
            labfacs=levels(facs)
            nbfacs=length(labfacs)

            contra=data.frame(k$k1,k$k2)
            nb=dim(contra)[1]
            cons=c()
            #validate(
            #  need((contra[k,1] %in% labfacs) & (contra[k,2] %in% labfacs), "One of the groups in contrast file does not match a group in phenotype file. Make sure names match and upload again. 
            #Once correct file is entered, 'Computing differentially expressed genes' message will display.")
            #)
            
            for (k in 1:nb) {
              #if ((contra[k,1] %in% labfacs) & (contra[k,2] %in% labfacs) )
              #{
                cons=c(cons,paste(contra[k,1],"-",contra[k,2],sep=""))
              #} else {
              #  cat("One of the groups in contrasts file at line :",k+1,"does not match a group in phenotype file..Quitting!!!\n")
              #  info('One of the groups in contrast file does not match a group in phenotype file. Make sure names match and upload again.')
              #  print( contra )
              #  stopApp(-1)
              #}
            }
            
            myfactor <- factor(pData(norm())$group)
            design1 <- model.matrix(~0+myfactor)
            colnames(design1) <- levels(myfactor)
            
            fit1 <- lmFit(norm(),design1)
            contrast.matrix <- makeContrasts(contrasts=cons,levels=design1)
            
            fit2 <- contrasts.fit(fit1, contrast.matrix)
            ebayes.fit2=eBayes(fit2) # smooths the std error
            incProgress(0.25)
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
                        } else {
                          if (raw()@annotation=="pd.mouse430.2") {
                            Annot <- data.frame(ACCNUM=sapply(contents(mouse4302ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mouse4302SYMBOL), paste, collapse=", "), DESC=sapply(contents(mouse4302GENENAME), paste, collapse=", "))
                          } else {
                            if (raw()@annotation=='pd.hg.u133a') {
                              Annot <- data.frame(ACCNUM=sapply(contents(hgu133aACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133aSYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133aGENENAME), paste, collapse=", "))
                            } else {
                              if (raw()@annotation=='pd.hugene.1.0.st.v1') {
                                Annot <- data.frame(ACCNUM=sapply(contents(hugene10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hugene10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hugene10sttranscriptclusterGENENAME), paste, collapse=", "))
                              } else {
                                if (raw()@annotation=='pd.mogene.1.0.st.v1') {
                                  Annot <- data.frame(ACCNUM=sapply(contents(mogene10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mogene10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(mogene10sttranscriptclusterGENENAME), paste, collapse=", "))
                                } else {
                                  if (raw()@annotation=='pd.hg.u133a.2') {
                                    Annot <- data.frame(ACCNUM=sapply(contents(hgu133a2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133a2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133a2GENENAME), paste, collapse=", "))
                                  } else {
                                    if (raw()@annotation=='pd.huex.1.0.st.v2') {
                                      Annot <- data.frame(ACCNUM=sapply(contents(huex10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(huex10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(huex10sttranscriptclusterGENENAME), paste, collapse=", "))
                                    } else {
                                      if (raw()@annotation=='pd.hg.u219') {
                                        Annot <- data.frame(ACCNUM=sapply(contents(hgu219ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu219SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu219GENENAME), paste, collapse=", "))
                                      } else {
                                        if (raw()@annotation=='pd.ht.hg.u133.plus.pm') {
                                          Annot <- data.frame(ACCNUM=sapply(contents(hgu133plus2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133plus2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133plus2GENENAME), paste, collapse=", "))
                                        } else {
                                          if (raw()@annotation=='pd.mg.u74av2') {
                                            Annot <- data.frame(ACCNUM=sapply(contents(mgu74av2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mgu74av2SYMBOL), paste, collapse=", "), DESC=sapply(contents(mgu74av2GENENAME), paste, collapse=", "))
                                          } else {
                                            if (raw()@annotation=='pd.mouse430a.2') {
                                              Annot <- data.frame(ACCNUM=sapply(contents(mouse430a2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mouse430a2SYMBOL), paste, collapse=", "), DESC=sapply(contents(mouse430a2GENENAME), paste, collapse=", "))
                                            } else {
                                              if (raw()@annotation=='pd.moe430a') {
                                                Annot <- data.frame(ACCNUM=sapply(contents(moe430aACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(moe430aSYMBOL), paste, collapse=", "), DESC=sapply(contents(moe430aGENENAME), paste, collapse=", "))
                                              } else {
                                                if (raw()@annotation=='pd.hg.u95av2') {
                                                  Annot <- data.frame(ACCNUM=sapply(contents(hgu95av2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu95av2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu95av2GENENAME), paste, collapse=", "))
                                                } else {
                                                  if (raw()@annotation=='pd.hta.2.0') {
                                                    Annot <- data.frame(ACCNUM=sapply(contents(hta20transcriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hta20transcriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hta20transcriptclusterGENENAME), paste, collapse=", "))
                                                  } else {
                                                    if (raw()@annotation=='pd.moex.1.0.st.v1') {
                                                      Annot <- data.frame(ACCNUM=sapply(contents(moex10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(moex10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(moex10sttranscriptclusterGENENAME), paste, collapse=", "))
                                                    } else {
                                                      if (raw()@annotation=='pd.hg.u133b') {
                                                        Annot <- data.frame(ACCNUM=sapply(contents(hgu133bACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133bSYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133bGENENAME), paste, collapse=", "))
                                                      } else {
                                                        if (raw()@annotation=='pd.hugene.1.1.st.v1') {
                                                          Annot <- data.frame(ACCNUM=sapply(contents(hugene11sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hugene11sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hugene11sttranscriptclusterGENENAME), paste, collapse=", "))
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
            
            incProgress(0.25)
            mylist=vector("list",nb)
            
            for (i in 1:nb)
            {
              
              all.genes.con = topTable(ebayes.fit2, coef = i, number=nrow(ebayes.fit2))
              
              # Merge data frames together (like a database table join)
              
              all <- merge(all.genes.con, Annot,by.x=0, by.y=0, all.x=T)
              all=all[order(all$P.Value),]
              colnames(all)[1]="probsetID"
              
              #add fold change and rearrange columns
              all$FC = ifelse(all$logFC<0, -1/(2^all$logFC), 2^all$logFC)
              all = all[,c(9,1,8,10,11,2,5,6,3,4,7)]
              
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
            
            #mylist
            list(mylist=mylist, Annot=Annot, cons=cons, design1=design1, nb=nb)
          })
          ##-------------
        }
      ) #end DEG
      
      pathways=reactive(
        {
          #L2P pathway starts here KV
          k$all = cbind(paste0(k$k1, ' vs ', k$k2))
          num = which(input$NumContrasts==k$all[,1])
          
          all = deg()$mylist[[num]]
          #up=vector("list",nb)
          #dw=vector("list",nb)
          
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
          
          #running volcano plot prior to pathways inexplicably turns SYMBOL into a factor so:
          fin.up$SYMBOL = as.character(fin.up$SYMBOL)
          fin.dw$SYMBOL = as.character(fin.dw$SYMBOL)
          
          #x2=rownames(fin.dw)
          #gdw=apply(array(as.character(x2)),1,function(z) unlist(strsplit(z, "\\|"))[2])
          
          if (raw()@annotation=="pd.hg.u133.plus.2" | raw()@annotation=="pd.hugene.2.0.st" | raw()@annotation=="pd.clariom.s.human.ht" | raw()@annotation=="pd.clariom.s.human" | raw()@annotation=='pd.hg.u133a' | raw()@annotation=='pd.hugene.1.0.st.v1' | raw()@annotation=='pd.hg.u133a.2' | raw()@annotation=='pd.huex.1.0.st.v2' | raw()@annotation=='pd.hg.u219' | raw()@annotation=='pd.ht.hg.u133.plus.pm' | raw()@annotation=='pd.hg.u95av2' | raw()@annotation=='pd.hta.2.0' | raw()@annotation=='pd.hg.u133b' | raw()@annotation=='pd.hugene.1.1.st.v1') 
          {
            cat(fin.up$SYMBOL,file=(paste0(input$ProjectID,'_',names(deg()$mylist[num]),'_Top500_Up.txt')), sep='\n')
            cat(fin.dw$SYMBOL,file=(paste0(input$ProjectID,'_',names(deg()$mylist[num]),'_Top500_Down.txt')),sep='\n')
          }
          else
          {
            cat(fin.up$SYMBOL,file=paste0(names(deg()$mylist[num]),"_Top500temp_Up.txt"),sep='\n')
            cat(fin.dw$SYMBOL,file=paste0(names(deg()$mylist[num]),"_Top500temp_Dw.txt"),sep='\n')
            
            system(paste0("cat ",names(deg()$mylist[num]),"_Top500temp_Up.txt | grep -v \"^NA\" | ./m2h | grep -v XXXX | cut -f2 -d\" \" >",input$ProjectID,'_',names(deg()$mylist[num]),"_Top500_Up.txt"))
            system(paste0("cat ",names(deg()$mylist[num]),"_Top500temp_Dw.txt | grep -v \"^NA\" | ./m2h | grep -v XXXX | cut -f2 -d\" \" >",input$ProjectID,'_',names(deg()$mylist[num]),"_Top500_Down.txt"))
          }
          system(paste0("cat ",input$ProjectID,'_',names(deg()$mylist[num]),"_Top500_Up.txt |sort | uniq | ./l2p >",input$ProjectID,'_',names(deg()$mylist[num]),"_Pathways_Up.txt"))
          system(paste0("cat ",input$ProjectID,'_',names(deg()$mylist[num]),"_Top500_Down.txt |sort | uniq | ./l2p >",input$ProjectID,'_',names(deg()$mylist[num]),"_Pathways_Down.txt"))
          
          addUpCol = read.delim(paste0(input$ProjectID,'_',names(deg()$mylist[num]),"_Pathways_Up.txt"), sep = '\t')
          addDwCol = read.delim(paste0(input$ProjectID,'_',names(deg()$mylist[num]),"_Pathways_Down.txt"), sep = '\t')
          
          colnames(addUpCol)=c("pval","fdr","ratio","nb.hits","nb.genes.path","nb.user.genes","tot.back.genes","path_id","source","description","type","gene.list")
          colnames(addDwCol)=c("pval","fdr","ratio","nb.hits","nb.genes.path","nb.user.genes","tot.back.genes","path_id","source","description","type","gene.list")
          addUpCol = addUpCol[order(addUpCol$pval),]
          addDwCol = addDwCol[order(addDwCol$pval),]
          addUpCol = addUpCol[,c(8,9,10,11,1,2,3,12,4,5,6,7)]
          addDwCol = addDwCol[,c(8,9,10,11,1,2,3,12,4,5,6,7)]
          write.table(addUpCol, file = paste0(input$ProjectID,'_',names(deg()$mylist[num]),"_Pathways_Up.txt"), sep = '\t', row.names = F)
          write.table(addDwCol, file = paste0(input$ProjectID,'_',names(deg()$mylist[num]),"_Pathways_Down.txt"), sep = '\t', row.names = F)
          
          list(up=addUpCol,dw=addDwCol)
        } 
      ) #end pathways
      
      ssGSEA=reactive(
        {
          withProgress(message = 'Performing ssGSEA', detail = 'may take a couple of minutes ...', value = 0, {
          ssgs =  merge(exprs(norm()),deg()$Annot,by.x=0, by.y=0, all.x=T)
          ssgs = ssgs[ssgs$SYMBOL!='NA',]
          ssgs = subset(ssgs, select=-c(ACCNUM,DESC,Row.names))
          ssgs = aggregate(.~SYMBOL,data=ssgs,mean)                               #aggregate duplicate probes by mean
          rownames(ssgs) = ssgs$SYMBOL
          ssgs = subset(ssgs, select=-c(SYMBOL))
          ssgs = as.matrix(ssgs)
    
          gset = getGmt(input$geneSet)                                            #input gene set 
          ssgsResults = gsva(ssgs, gset, method='ssgsea')                         #run GSVA
          
          incProgress(amount = 0.50, detail = 'Performing differential expression analysis of pathways...')
          
          fit1 = lmFit(ssgsResults,deg()$design1)
          contrast.matrix = makeContrasts(contrasts=deg()$cons,levels=deg()$design1)
          fit2 = contrasts.fit(fit1,contrast.matrix)
          ebayes.fit2 = eBayes(fit2)
          
          myPathways=vector("list",deg()$nb)
          for (i in 1:deg()$nb)
          {
            all.pathways = topTable(ebayes.fit2, coef=i, number=nrow(ebayes.fit2))
            all.pathways = all.pathways[order(abs(all.pathways$logFC),decreasing=T),]
            colnames(all.pathways)[2] = 'EnrichmentScore'
            myPathways[[i]] = all.pathways
          }
          names(myPathways)=deg()$cons
          list(myPathways=myPathways,ssgsResults=ssgsResults)
          })
        }
      )
      
      
      
  #Processing all outputs
  
  ####creates a list of colors specific to each group
  fs = factor(pData(raw())$group)
  #fs = factor(pData(raw())$SampleGroup)
  lFs=levels(fs)
  numFs=length(lFs)
  colors = list()
  for (i in 1:numFs){
    colors[which(fs==lFs[i])] = i*5
  }
  colors = unlist(colors)
  ####end
  
  output$projectid=renderText({paste("Project ID: ",input$ProjectID)})
  output$rawhist=renderPlot(
    {
      hist(raw(),which="all", main =" Raw Samples distribution")
    }
  )
  
  ###Beginning raw maplot###
  output$rawmaplot=renderUI({
    #facs <- pData(raw())$SampleID
    facs <- pData(raw())$SampleID
    nbfacs=length(facs)
    plot_output_list <- lapply(1:nbfacs, function(i) {
      plotname <- paste("plot", i, sep="")
      plotOutput(plotname, height = 400, width = 600)
    })
    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, plot_output_list)
  })
  
  facs <- pData(raw())$SampleID
  nbfacs=length(facs)

  for (i in 1:nbfacs) {
    local({
      my_i <- i
      plotname <- paste("plot", my_i, sep="")
      output[[plotname]] <- renderPlot({
        withProgress(message = 'Generating Raw Maplot', detail = paste0('Plot ', my_i, ' ...'), value = (my_i/nbfacs), {
          MAplot(raw(),which=my_i,plotFun=smoothScatter,refSamples=c(1:nbfacs), main='', cex=2)
        })
      })
    })
  }
  ###end raw maplot###
  
  output$rawbox=renderPlot(
    {
      par(mar=c(8,4,4,2))
      boxplot(raw(), col=colors, which="all", main="Boxplots before normalization",las=2,names=pData(raw())$SampleID)
    }
  )
  output$rle=renderPlot(
    {
      par(mar=c(8,4,4,2))
      RLE(qc(), main="RLE plot",names=pData(raw())$SampleID, col=colors, las=2)
    }
  )
  output$nuse=renderPlot(
    {
      par(mar=c(8,4,4,2))
      NUSE(qc(), main="NUSE plot",names=pData(raw())$SampleID, col=colors, las=2)
    }
  )
  output$rmahist=renderPlot(
    {
      hist(norm(), main ="Distribution after Normalization")
    }
  )
  
  ## MVAplot after normalization
  output$normaplot=renderUI({
    facs2 <- pData(norm())$SampleID
    nbfacs2=length(facs2)
    plot_output_list2 <- lapply(1:nbfacs2, function(i) {
      plotname2 <- paste("plota", i, sep="")
      plotOutput(plotname2, height = 400, width = 600)
    })
    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, plot_output_list2)
  })
  
  facs2 <- pData(norm())$SampleID
  nbfacs2=length(facs2)
  
  for (i in 1:nbfacs2) {
    local({
      my_i <- i
      plotname2 <- paste("plota", my_i, sep="")
      # MA plots are then used to visualize intensity-dependent ratio for each group
      output[[plotname2]] <- renderPlot({
        withProgress(message = 'Generating Normalized Maplot', detail = paste0('Plot ', my_i, ' ...'), value = my_i/nbfacs2, {
          MAplot(norm(),which=my_i,plotFun=smoothScatter,refSamples=c(1:nbfacs2),main='', cex=2)
        })
      })
    })
  }
  output$rmabox=renderPlot(
    {
      par(mar=c(8,4,4,2))
      boxplot(norm(),col=colors, main="Boxplots after RMA normalization",las=2,names=pData(norm())$SampleID)
    }
  )
  ## end mvaplat after normalization
  
  ## pca 3D
  output$pca3d=renderRglwidget(
    {
      withProgress(message = 'Generating PCA', detail = 'starting ...', value = 0.5, {
        tedf= t(exprs(norm()))
        
        #removes zero  variances (issue with small sample sizes)
        if (length(which(apply(tedf, 2, var)==0)) >= 0){
          tedf = tedf[ , apply(tedf, 2, var) != 0]
        }
        
        pca=prcomp(tedf, scale. = T)
        incProgress(amount = 0.25, detail = 'determining variance ...')
        rgl.open(useNULL=T)
        bg3d('white')
        plot3d(pca$x[,1:3],col=colors, type='s',size=2)
        group.v=as.vector(pData(norm())$SampleID)
        text3d(pca$x, pca$y, pca$z, group.v, cex=0.6, adj=1.5)
        par3d(mouseMode = "trackball")
        rglwidget()
      })
    }
  )
  output$heatmap=renderPlotly(
    {
      mat=as.matrix(dist(t(exprs(norm()))))
      rownames(mat)=pData(norm())$SampleID
      colnames(mat)=rownames(mat)
      heatmaply(mat,margins = c(80,120,60,40),colorRampPalette(colors = c("red", "yellow")))
    }
  )
  output$deg=DT::renderDataTable(DT::datatable(
    {
      #turn NumContrasts back into a number
      k$all = cbind(paste0(k$k1, ' vs ', k$k2))
      num = which(input$NumContrasts==k$all[,1])
   
      dat = deg()$mylist[[num]]
      dat = dat[,-6]
      dat[,6:7] = format(dat[,6:7], scientific = TRUE)
      
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
    }, caption =paste0("contrast: ",names(deg()$mylist)[num])
  )
  )
  output$topUp=DT::renderDataTable(DT::datatable(
    {
      k$all = cbind(paste0(k$k1, ' vs ', k$k2))
      num = which(input$NumContrasts==k$all[,1])
      
      #deg()[[input$NumContrasts]]
      #topUp = read.delim(paste0(input$ProjectID,'_',names(deg())[input$NumContrasts],"_Pathways_Up.txt"), sep = '\t', header = T)
      #dev.off(which = plotly)
      topUp = pathways()$up
      
      if(is.na(input$pathPval)) {
        topUp
      } else {
        topUp = topUp[(as.numeric(topUp[,5]) <= input$pathPval),]
      }
      topUp
    } , caption=paste0("Pathways for the top 500 Upregulated Genes: ", names(deg()$mylist)[num]),
    options = list(columnDefs = list(list(targets = 8,
                                          render = JS("function(data, type, row, meta) {",
                                                      "return type === 'display' && data.length > 30 ?",
                                                      "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
                                                      "}")))), callback = JS('table.page(3).draw(false);')
  )
  )
  output$topDown=DT::renderDataTable(DT::datatable(
    {
      #callDEG = deg()[[input$NumContrasts]]
      #topDw = read.delim(paste0(input$ProjectID,'_',names(deg())[input$NumContrasts],"_Pathways_Down.txt"), sep = '\t', header = T)
      
      k$all = cbind(paste0(k$k1, ' vs ', k$k2))
      num = which(input$NumContrasts==k$all[,1])
      
      topDw = pathways()$dw
      
      if (is.na(input$pathPval)) {
        topDw
      } else {
        topDw = topDw[(as.numeric(topDw[,5]) <= input$pathPval),]
      }
      topDw
      
    } , caption=paste0("Pathways for the top 500 Downregulated Genes: ", names(deg()$mylist)[num]),
    options = list(columnDefs = list(list(targets = 8, 
                                          render = JS("function(data, type, row, meta) {",
                                                      "return type === 'display' && data.length > 30 ?",
                                                      "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
                                                      "}")))), callback = JS('table.page(3).draw(false);')
  )
  )
  output$volcano=renderPlotly(
    {
      withProgress(message = 'Generating Volcano Plot', detail = 'starting ...', value = 0, {
        k$all = cbind(paste0(k$k1, ' vs ', k$k2))
        num = which(input$NumContrasts==k$all[,1])
        
        dat=deg()$mylist[[num]]
        #dat=deg()[[input$NumContrasts]]
        
        dat = dat[dat$SYMBOL!='NA',]
        log_FC=dat$logFC
        log_pval=-log10(dat$P.Value)
        Significant=rep("NotSignificant",length(log_FC))
        Significant[which(dat$P.Value<0.05 & abs(dat$logFC)>=1)]="AbsLogFoldChange>1 & PValue<0.05"
        Significant[which(dat$P.Value<0.05 & abs(dat$logFC)<1)]="PValue<0.05"
        Significant[which(dat$P.Value>=0.05 & abs(dat$logFC)>=1)]="AbsLogFoldChange>1"
        gene=dat$SYMBOL
        volcano_data=as.data.frame(cbind(gene,log_FC,log_pval,Significant))
        incProgress(0.50)
        plot_ly(type='scatter', data = volcano_data, x = log_FC, y = log_pval, text = gene, mode = "markers", color = Significant) %>% layout(title=paste0('Volcano plot for: ',names(deg()$mylist)[num]),xaxis=list(title="Fold Change",range =c(-5,5),tickvals=c(-5,-4,-3,-2,-1,0,1,2,3,4,5),ticktext=c('-32','-16','-8','-4','-2','1','2','4','8','16','32')),yaxis=list(title="-Log10 pvalue",range =c(0,15)))
      })
    })
  output$ssgsea=DT::renderDataTable(DT::datatable(
    {
      k$all = cbind(paste0(k$k1, ' vs ', k$k2))
      num = which(input$NumContrasts==k$all[,1])
      
      ssGSEA()$myPathways[[num]]
    }
  ))
  output$ssHeatmap=renderPlot(
    {
      k$all = cbind(paste0(k$k1, ' vs ', k$k2))
      num = which(input$NumContrasts==k$all[,1])
      each = ssGSEA()$myPathways[[num]]

      sampleColumns = c(which(v$data$group==k$k2[num]),which(v$data$group==k$k1[num]))          #subset columns (samples) for user input contrast
      paths = ssGSEA()$ssgsResults[rownames(ssGSEA()$ssgsResults) %in% rownames(each)[1:50],]   #subset diff exprs pathways for user input contrast
      paths = paths[,sampleColumns]
      
      matCol = data.frame(group=v$data[sampleColumns,2])
      rownames(matCol) = v$data[sampleColumns,1]
      matColors = list(group = unique(colors[sampleColumns]))
      names(matColors$group) = unique(v$data[sampleColumns,2])
      
      pheatmap(paths,color=inferno(10),annotation_col=matCol,annotation_colors=matColors,drop_levels=TRUE,fontsize=7, main='Enrichment Scores for Top 50 Differentially Expressed ssGSEA Pathways')
     
    }
  )
  #observeEvent(input$rep, {
    # withProgress(message = 'Generating HTML report', detail = 'starting ...', value = 0.5, {
    #   # out <- render('../report_ver4.Rmd','html_document',paste(input$ProjectID,"_","report.html",sep=""),getwd(),getwd())
    #   out <- rmarkdown::render('report_ver4.Rmd','html_document',paste(input$ProjectID,"_","report.html",sep=""))
    # })
  #})
  output$projectid=renderText({paste("Project ID: ",input$ProjectID)})
  output$downloadReport <- downloadHandler(
    filename = function() {paste(input$ProjectID,"_Report.zip",sep="")},
    content = function(file) {
      withProgress(message = 'Generating HTML report', detail = 'starting ...', value = 0.5, {
      out <- rmarkdown::render('report_ver4.Rmd','html_document',paste(input$ProjectID,"_","report.html",sep=""))
      #mytables=list.files(pattern="\\.html$")
      mytables=list.files(pattern=paste(input$ProjectID,".*.html",sep=""))
      zip(file, mytables, flags = "-r9X", extras = "", zip = Sys.getenv("R_ZIPCMD", "zip"))
      })
    }
    )
  output$downloadTables <- downloadHandler(
    filename = function() {paste(input$ProjectID,"_tables.zip",sep="")},
    content = function(file) {
      withProgress(message = 'Preparing download', detail = 'starting ...', value = 0.5, {
      # mytables=list.files(pattern="\\.txt$")
      mytables=list.files(pattern=paste(input$ProjectID,".*.txt",sep=""))
      zip(file, mytables, flags = "-r9X", extras = "", zip = Sys.getenv("R_ZIPCMD", "zip"))
      })
    }
  )
  })
})
})


