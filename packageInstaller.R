InstalledPackage <- function(package) 
{
  available <- suppressMessages(suppressWarnings(sapply(package, require, quietly = TRUE, character.only = TRUE, warn.conflicts = FALSE)))
  missing <- package[!available]
  if (length(missing) > 0) return(FALSE)
  return(TRUE)
}

CRANChoosen <- function()
{
  return(getOption("repos")["CRAN"] != "@CRAN@")
}

UsePackage <- function(package, defaultCRANmirror = "http://cran.at.r-project.org") 
{
  if(!InstalledPackage(package))
  {
    if(!CRANChoosen())
    {       
      chooseCRANmirror()
      if(!CRANChoosen())
      {
        options(repos = c(CRAN = defaultCRANmirror))
      }
    }
    
    suppressMessages(suppressWarnings(install.packages(package)))
    if(!InstalledPackage(package)) return(FALSE)
  }
  return(TRUE)
}

libraries = c("rgl","rglwidget","DT","getopt",
              "calibrate","shinyRGL","htmltools",
              "viridis","dendsort","shiny","shinyjs","gplots","knitr",
              "reshape","RColorBrewer","mixOmics","rmarkdown","ggplot2","ggfortify",
              "heatmaply","plotly","pheatmap")
for(library in libraries) 
{ 
  if(!UsePackage(library))
  {
    stop("Error!", library)
  }
}

listOfBiocPackages = c("GEOquery","pd.clariom.s.rat","clariomsrattranscriptcluster.db",
                       "pd.mogene.2.0.st","mogene20sttranscriptcluster.db","pd.hg.u133.plus.2",
                       "hgu133plus2.db","pd.hugene.2.0.st","hugene20sttranscriptcluster.db",
                       "pd.clariom.s.human.ht","clariomshumanhttranscriptcluster.db",
                       "pd.clariom.s.human","clariomshumantranscriptcluster.db","pd.clariom.s.mouse.ht",
                       "clariomsmousehttranscriptcluster.db","pd.clariom.s.mouse","clariomsmousetranscriptcluster.db",
                       "hugene10sttranscriptcluster.db","pd.mogene.1.0.st.v1","mogene10sttranscriptcluster.db",
                       "pd.hg.u133a.2","hgu133a2.db","pd.huex.1.0.st.v2","huex10sttranscriptcluster.db",
                       "pd.hg.u219","hgu219.db","pd.mg.u74av2","mgu74av2.db","pd.mouse430a.2",
                       "mouse430a2.db","pd.moe430a","moe430a.db","pd.hg.u95av2","hgu95av2.db","pd.hta.2.0",
                       "hta20transcriptcluster.db","pd.moex.1.0.st.v1","moex10sttranscriptcluster.db","pd.hg.u133b",
                       "hgu133b.db","pd.hugene.1.1.st.v1","hugene11sttranscriptcluster.db","pd.mogene.1.1.st.v1",
                       "mogene11sttranscriptcluster.db","pd.hugene.2.1.st","hugene21sttranscriptcluster.db","oligo","geneplotter","multtest","Biobase",
                       "GSVA","GSEABase","annotate","limma")

notInstalled = setdiff(listOfBiocPackages, rownames(installed.packages()))
if( length(notInstalled) ) {
  BiocManager::install(notInstalled, update = FALSE, version = "3.8")
}
