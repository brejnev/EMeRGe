#### app.R : The code below launches the application ##### ----  

## Detects the operating system which the program is running on.
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

## This is specific for the MNP PC (to locate the libraries)
if (get_os()=="windows"){
        libpaths=.libPaths()[2]
        libpaths=c("E:/RHome/R-4.0.5/lib", libpaths)
        .libPaths(libpaths)
}

## Loading the user interface and the server functions
suppressPackageStartupMessages({
        source("ui.R")
        source("server.R")
})

## Global settings and variables
if (get_os()=="windows"){
        options(stringsAsFactors = F, browser= "C:/Program Files (x86)/Google/Chrome/Application/chrome.exe")
}else{
        options(stringsAsFactors = F, browser="/Applications/Google Chrome.app/Contents/Frameworks/Google Chrome Framework.framework/Versions/86.0.4240.75/Helpers/Google Chrome Helper (GPU).app (0xd2c)")
}

shinyApp(ui=ui, server=server)
