library(shiny)
library(shinyjs)
library(shinyRGL)
library(rgl)
library(plotly)
 
shinyUI(fluidPage(
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
      ),
    tags$style(
      HTML(".shiny-output-error-validation {
            color: red;
            font-size: 20px;
      }"
      )
    )
    ),
  useShinyjs(),
  titlePanel("CCBR Microarray analysis workflow", windowTitle="CCBR Microarray analysis workflow"),
  h5("(For Affymetrix human and mouse data)"),
  sidebarLayout(
    sidebarPanel(
      # tags$link(rel = 'stylesheet', type = 'text/css', href = 'styles.css'),
                width=12,
                 fluidRow(align="Top",
                          
                          column(2,
                                 textInput("ProjectID", label=h6("Project ID:"), value="ccbr", width="150px")
                                 # radioButtons("Platform", label=h6("Select platform"), choices=c("hgu133plus2" = "h133p2", "Mouse.gene.2.0.st" = "mst2"),selected="mst2")
                          ),
                          column(3,
                                 # textInput("Indir", label=h6("Path to input directory"), value="/Users/elloumif/Documents/cels6/cels",width="300px")
                                 fileInput("Indir", label=h6("Select CEL files"),multiple =T)
                          ),
                          
                          column(3,
                                 fileInput("pheno", label=h6("Choose phenotype file"))
                          ),
                          column(2,
                                 fileInput("const", label=h6("Choose contrast file")),
                                 numericInput("NumContrasts", label=h6("Choose contrast to show"),value="1", width="150px")
                                 #numericInput("pval", label=h6("KEGG/GO Enrichment Pvalue threshold"),value="0.05", width="150px"),
                                 #sliderInput("fc",label=h6("KEGG/GO Enrichment FC threshold"),min=0.5,max=3.5,value=1.5,step=0.5)
                          ),
                          column(2,
                                 numericInput("pval", label=h6("P-value threshold for DEGs"),value=0.05, width="150px"),
                                 numericInput("fc", label=h6("FC threshold for DEGs"),value=1.5, width="150px"),
                                 numericInput("pathPval", label=h6("P-value threshold for pathways"),value=0.05, width="150px")
                                 # sliderInput("fc",label=h6("KEGG/GO Enrichment FC threshold"),min=0.5,max=3.5,value=1.5,step=0.5)
                          )
                          #       numericInput("NumContrasts", label=h6("Which contrast to show"),value="1", width="150px"),
                           #      numericInput("fdr", label=h6("FDR for Pathway enrichment"),value="0.05", width="150px")
                          #)
                 ),
                
            
                # submitButton(text="Click to assign samples to groups and create contrasts")
                actionButton(inputId="go",label="Start"),
                actionButton(inputId="rep",label="Generate Report"),
                downloadButton('downloadReport', label = 'Download Report'),
                downloadButton('downloadTables', label = 'Download Tables')
    ),
    mainPanel(
      navbarPage(title = "Results",
                  #title = "Microarray",
                  #tabPanel("Results"," "),
                  navbarMenu (title="Pre-normalization QC plots",
                    tabPanel("Histogram",plotOutput("rawhist")),
                    tabPanel("Maplots",uiOutput("rawmaplot")),
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
                  navbarMenu (title="Help",
                              tabPanel("Manual",uiOutput("manual"))
                  )
                  
      )
      # end solution2
    )
  )
  ## div(style="display:inline-block",submitButton("Generate PDF report"), style="float:center")
))
