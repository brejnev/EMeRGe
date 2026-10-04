## ------------------------------------------------------------ ##
##                      *  EMeRGe v1.0 *                         ##   
## Installation of required R libraries                         ## 
## Execute this script '' in R Studio or use Rscript command    ##
## ------------------------------------------------------------ ##

get_os <- function(){
        sysinf <- Sys.info()
        if (!is.null(sysinf)){
                os <- sysinf['sysname']
                if (os == 'Darwin')
                        os <- "osx"
        } else { ## mystery machine
                os <- .Platform$OS.type
                if (grepl("^darwin", R.version$os))
                        os <- "osx"
                if (grepl("linux-gnu", R.version$os))
                        os <- "linux"
        }
        tolower(os)
}

os=get_os()
if (os=="windows"){
        ## Installation for all user run the program (Rstudio) as administrator
        libpaths=.libPaths()[2]
        libpaths=c(libpaths,"E:/RHome/R-4.0.5/lib")
        .libPaths(libpaths)
}

## Cran packages:
list.of.packages=c("shiny",
                   "shinyjs",
                   "shinyWidgets",
                   "shinyalert",
                   "shinyFiles",
                   "fs",
                   "qs",
                   "ellipse",
                   "pdftools", 
                   "gridExtra",
                   "grid",
                   "kableExtra", 
                   "FactoMineR", 
                   "factoextra", 
                   "pheatmap",
                   "Cairo",
                   "magick",
                   "igraph",
                   "reshape",      
                   "Matrix",       
                   "ggplot2",      
                   "ggpubr",       
                   "crayon",
                   "png",
                   "plotly",
                   "scatterplot3d",
                   "rhandsontable", 
                   "RColorBrewer",
                   "colorRamps",
                   "readxl",
                   "readr",
                   "RJSONIO",
                   "rlang",
                   "Rserve",
                   "dplyr",
                   "tidyverse",
                   "htmltools",
                   "webshot",
                   "flashClust")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)!=0){
        cat("Installing new packages from CRAN: please wait ...\n")
        cat(new.packages,"\n")
        install.packages(new.packages)
        
}else(cat("All required R libraries are installed."))



## Biocoductor packages


list.of.packages=c("pathview","org.Hs.eg.db","org.Mm.eg.db","KEGGgraph","Rgraphviz")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)!=0){
        cat("Installing new packages from Bioconductor: please wait ...\n")
        cat(new.packages,"\n")
        
        if (!requireNamespace("BiocManager", quietly = TRUE))
                install.packages("BiocManager")
        
        BiocManager::install(new.packages )
        
}else(cat("All required Bioconductor libraries are installed."))


### Recent versison of dplyr (might be required)

if (packageVersion("devtools") < 1.6) {
        install.packages("devtools")
}
library(devtools)
devtools::install_github("hadley/lazyeval")
devtools::install_github("hadley/dplyr")

