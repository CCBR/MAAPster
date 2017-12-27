library(shiny)
library(DT)
library(Biobase)
library(shinyjs)
library(rglwidget)
library(plotly)

options(shiny.maxRequestSize = 500*1024^2)

shinyUI(
  fluidPage(
    tags$head(
      tags$style(
        HTML(".shiny-notification {
             height: 100px;
             width: 800px;
             position:fixed;
             font-weight: 500;
             font-size: 20px;
             background-color: #c1cdcd;
             top: calc(40% - 50px);;
             left: calc(50% - 400px);;
             }"
        )
      )
    ),
    useShinyjs(),
    
    titlePanel("CCBR Microarray analysis workflow", windowTitle="CCBR Microarray analysis workflow"),
    h5("(For Affymetrix human and mouse data)"),
    
    sidebarPanel(
      selectInput('analysisType', "Choose type of analysis",
                  c('Upload CEL files' = "CEL", 'Analyze GEO data' = 'GEO', ""),
                  selected = ''
      )
    ),
      
    conditionalPanel(
      condition = "input.analysisType == 'CEL'",
      fluidRow(),
      column(2,
        textInput("ProjectID", label=h6("Project ID:"), value="CCBR", width="150px")
      ),
      column(3,
        fileInput("Indir", label=h6("Select CEL files"),multiple =T),
        h5("Choose short descriptive file names with no spaces (eg Ctl1.CEL, KO1.CEL)")
      ),
      br(),
      br(),
      actionButton(inputId="CELbutton", label="Display"),
      br(),
      br()
    ),
    
    br(),
    br(),
    
    conditionalPanel(
      condition = "input.analysisType == 'GEO'",
      fluidRow(),
      column(2,
        textInput("ProjectID", label=h6("Project ID:"), value="CCBR", width="150px")
      ),
      column(2,
        textInput("gseid", label= h6("Accession Code"), value="8 digit GSE code", width="150px")
        ),
      br(),
      br(),
      actionButton(inputId="button", label="Display")
    ),
    
    
    br(),
    br(),
    br(),
    br(),
    
    


    shinyjs::hidden(
      div(id= "hide1",
          # fluidRow(
          #   column(3,
          #          selectInput("number", "Number of groups:",
          #                      c("1"="1", "2"="2", "3"="3", "4"="4", "5"="5", "6"="6", "7"="7", "8"="8", "9"="9", "10"="10", ""),
          #                      selected=""
          #          )
          #          #actionButton("button3", "Enter")
          #   )),
          column(3, 
                 uiOutput("ui"),
                 h6("(Name may be changed for compatibility with R)"),
                 br(),
                 h5("Select samples below and click \"Define\" to add to group"),
                 actionButton("button2", "Define")
      ))
    ),
    br(),
    
    
    #shinyjs::hidden(
    #  div(id='hide1',
      # conditionalPanel(
      # condition = 'input.number == 1',
      #     fluidRow(),
      #       column(2,
      #         #lapply(1:input$number, function(i) {     
      #         selectInput("group1", "Please select a group",
      #                     #choices = paste0('Group', i),
      #                     choices = c("Group_1" = "Group_1"),
      #                     selected = "Group_1"
      #         #)}
      #       ),
      #       actionButton("button2", "Define")
      # )
      # ),
      # 
    #   conditionalPanel(
    #   condition = 'input.number == 2',
    #   fluidRow(),
    #   column(2,
    #          #lapply(1:input$number, function(i) {     
    #          selectInput("group1", "Please select a group",
    #                      #choices = paste0('Group', i),
    #                      choices = c("Group_1" = "Group_1", "Group_2" = "Group_2"),
    #                      selected = "Group_1"
    #                      #)}
    #          ),
    #          actionButton("button2", "Define")
    #   )
    # ),
    
      
    br(),
    shinyjs::hidden(
      div(id= "hide2",
          br(),
          br(),
          fluidRow(column(10, wellPanel(DT:: dataTableOutput("mytable"))))
        )),
    
    shinyjs::hidden(
      div(id="hideCEL",
          br(),
          br(),
          fluidRow(column(10, wellPanel(DT:: dataTableOutput("mytableCEL"))))
      )),
    
    # shinyjs:: hidden(
    #   div(id= "hide3",
    #       mainPanel(
    #         actionButton("test","Contrast"))
    #   )),
    
   shinyjs:: hidden(
     div(id= "hide3",
         mainPanel(
           h4("Contrasts: "),
           uiOutput("choice1"),
           uiOutput("choice2"),
         mainPanel(
           fluidRow(align='Top',
                    column(2,
                           actionButton("addrow","Add Contrast")
                    )),
           mainPanel(
             tableOutput("contrastTable"))
         ))),
     
 
     # div(id= "hide3",
     #      # conditionalPanel(
     #      #   condition = "input.number == '1'",
     #        mainPanel(
     #          h4("Contrasts: "),
     #          selectizeInput("selectIn1", "Select", 
     #                         choices = c("Group_1" = "Group_1",
     #                                     "Group_2" = "Group_2",
     #                                     "Group_3" = "Group_3",
     #                                     "Group_4" = "Group_4",
     #                                     "Group_5" = "Group_5",
     #                                     "Group_6" = "Group_6",
     #                                     "Group_7" = "Group_7",
     #                                     "Group_8" = "Group_8",
     #                                     "Group_9" = "Group_9",
     #                                     "Group_10" = "Group_10")),
     #          selectizeInput("selectIn2", "Versus",
     #                         choices = c("Group_1" = "Group_1",
     #                                     "Group_2" = "Group_2",
     #                                     "Group_3" = "Group_3",
     #                                     "Group_4" = "Group_4",
     #                                     "Group_5" = "Group_5",
     #                                     "Group_6" = "Group_6",
     #                                     "Group_7" = "Group_7",
     #                                     "Group_8" = "Group_8",
     #                                     "Group_9" = "Group_9",
     #                                     "Group_10" = "Group_10"))),
     #     mainPanel(
     #         fluidRow(align='Top',
     #                  column(2,
     #                         actionButton("addrow","Add Contrast")
     #                  )),
     #          
     #    mainPanel(
     #        tableOutput("contrastTable"))
     # )),
    
    shinyjs::hidden(
      div(id='hideAnalysis',
        sidebarLayout(position = 'left',
          sidebarPanel(
            width=12,
            fluidRow(align='Top',
              # column(3,
              #        numericInput("NumContrasts", label=h6("Choose contrast to show"),value="1", width="150px")),
              column(2,
                     uiOutput("displayContrast")),
              column(2,
                     numericInput("pval", label=h6("P-value threshold for DEGs"),value=0.05, width="200px")),
              column(2,
                     numericInput("fc", label=h6("FC threshold for DEGs"),value=1.5, width="150px")),
              column(2, 
                     numericInput("pathPval", label=h6("P-value threshold for pathways"),value=0.05, width="200px")),
              column(4,
                     selectInput('geneSet', label=h6("Choose Gene Set for ssGSEA"),
                                 c("H: Hallmark Gene Sets"="h.all.v6.1.symbols.gmt", "C1: Positional Gene Sets"="c1.all.v6.1.symbols.gmt", "C2: Curated Gene Sets"="c2.all.v6.1.symbols.gmt", 
                                   "C3: Motif Gene Sets"="c3.all.v6.1.symbols.gmt", "C4: Computational Gene Sets"="c4.all.v6.1.symbols.gmt","C5: GO gene sets"="c5.all.v6.1.symbols.gmt", 
                                   "C6: Oncogenic Signatures"="c6.all.v6.1.symbols.gmt", "C7: Immunologic Signatures"="c7.all.v6.1.symbols.gmt"), selected="h.all.v6.1.symbols.gmt"))
            ),
            actionButton(inputId="analyze",label="Start")),
        mainPanel(
          br(),br(),br(),
          navbarPage(title = "Results",
                     navbarMenu (title="Pre-normalization QC plots",
                                 tabPanel("Histogram",plotOutput("rawhist")),
                                 tabPanel("Maplots", uiOutput("rawmaplot")),
                                 tabPanel("Boxplots", plotOutput("rawbox")),
                                 tabPanel("RLE",plotOutput("rle")),
                                 tabPanel("NUSE",plotOutput("nuse"))
                     ),
                    navbarMenu (title="Post-normalization plots",
                                tabPanel("Histogram",plotOutput("rmahist")),
                                tabPanel("Maplots",uiOutput("normaplot")),
                                tabPanel("Boxplots",plotOutput("rmabox")),
                                tabPanel("3D-PCA",rglwidgetOutput("pca3d")),
                                tabPanel("Interactive Heatmap",plotlyOutput("heatmap"))
                    ),
                    navbarMenu (title="DEG-Enrichments-tables",
                                tabPanel("Differentially Expressed Genes",DT::dataTableOutput("deg")),
                                tabPanel("Pathways for Upregulated Genes",DT::dataTableOutput("topUp")),
                                tabPanel("Pathways for Downregulated Genes",DT::dataTableOutput("topDown")),
                                tabPanel("Interactive Volcano Plot",plotly::plotlyOutput('volcano'))
                    ),
                    navbarMenu (title='Single Sample GSEA',
                                tabPanel("Enriched Pathways",DT::dataTableOutput('ssgsea')),
                                tabPanel("Pathway Heatmap",plotOutput("ssHeatmap", width='100%', height='800px'))
                                )
          ),br(),br(),br(),br()
        )
        )
      )
    ),
    shinyjs::hidden(
      div(id='hideDownloads',
          mainPanel(
            #actionButton(inputId="rep",label="Generate Report"),
            downloadButton('downloadReport', label = 'Download Report'),
            downloadButton('downloadTables', label = 'Download Tables')
          )
      )
    )
   )
  )
)
 
 

  
  
