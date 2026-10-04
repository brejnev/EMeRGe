#### server.R : The code below controls the execution / functionality of the application ##### ----

## Import R libraries
library(shinyjs)
library(shiny)
library(kableExtra)
library(rhandsontable)
library(readxl)

## Importing all program modules
source("idxUtils.R")
source("molarUtils.R")
source("midaUtils.R")

options(scipen = 999)


## Server function
server <- function(input, output,session) {
        
        
        # Some global variables ----
        .on.public.web<<-FALSE
        api.base <<- "http://api.xialab.ca"
        app.Path<<-getwd()
       
        # idxMOLAR global variables ----
        methods_regex <<- "HILIC_[AB]_|LIPID_|T3pos_|PHpos_|PHneg_|PHTCA_"
        methods_regex2<<-gsub("_\\|","\\|",methods_regex); methods_regex2=gsub("_$","",methods_regex2)
        
        idx_PB_Blanks_pattern<<-"Blank|pb|qc|pos.id|neg.id|sample_"
        
        # Other common global variables
        currModule<<-""
        
        MID_FE_global<<-NULL
        TotalLabel_FE_global<<-NULL
        PoolSize_global<<-NULL
        MID_PoolSize_global<<-NULL 
        TotalLabel_PoolSize_global<<-NULL
        
       
        #Initial setting
        if (file.exists("www/Plot.pdf")){file.remove("www/Plot.pdf")}
        if (file.exists("www/Kable.csv")){file.remove("www/Kable.csv")}

        output$docview <- renderUI({tags$iframe(style="height:700px; width:100%; scrolling=yes", src="EMeRGe.pdf")})

        volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
        
        shinyFileChoose(input, "file_inputData", roots = volumes, session = session,filetypes = c('',"tsv","csv","xls","xlsx"))
        shinyFileChoose(input, "file_PDF", roots = volumes, session = session,filetypes = c('', "pdf"))
        shinyFileChoose(input, "file_Kable", roots = volumes, session = session,filetypes = c('', "csv","tsv", "xls","xlsx") )

     
        # Save settings button  ----
        ### Saves the settings for group pairwise comparisons and statistical tests to be performed
        ### Save pairwise comparisons
        observeEvent(input$btn_Save, {
                if (is.null(input$hot_Groups)){
                        shinyalert("Settings","Please load the input file first.",type="error")
                }else{
                        studySettings[["samplesDF"]]=hot_to_r(input$hot_Samples)
                        studySettings[["groupsDF"]]=hot_to_r(input$hot_Groups)
                        studySettings[["orderDF"]]=hot_to_r(input$hot_Order)
                        
                        oT=studySettings[["orderDF"]] 
                        
                        ## check the order
                        if (length(which(duplicated(oT[,"Order"])))>0){shinyalert(paste0("       Invalid group order : ", oT[which(duplicated(oT[,"Order"])),"Order"], " is dublicated"), type="warning"); return()}
                        if (length(which(is.na(oT[,"Order"])))>0){shinyalert("       Invalid input : please correct the value next to ", rownames(oT)[which(is.na(oT[,"Order"]))], type="warning");return()}
                        ## check indicated group for exclusion 
                        if (length(which(!tolower(unique(oT[,"Included"]))%in% c("y","yes","n","no")))>0){
                                shinyalert("           Error: please use Yes or No, or Y or N to indicate what group to include/exclude")
                                return()
                        }
                        ## check if which have at least two group checked 
                        if (length(grep("y|yes",tolower(oT[,"Included"])))<2){
                                shinyalert("           Error: atleast two groups need to be checked")
                                return()
                        }
                        ## check whether there are at least two groups
                      
                        ## check the sample normalization factors 
                        sT=studySettings$samplesDF
                        if (length(which(is.na(as.numeric(sT[,"Normalization_factors"]))))>0){
                                shinyalert("           Error: non-numeric data was entered in the normalization factors table (Study settings)")
                                return()
                        }
                        
                        ## check if which have at least two group checked 
                        if (length(grep("y|yes",tolower(sT[,"Included"])))< 4){
                            shinyalert("           Error: atleast two groups each with two samples need to be checked")
                            return()
                        }
                        
                        ## check indicated sample for exclusion 
                        if (length(which(!tolower(unique(sT[,"Included"]))%in% c("y","yes","n","no")))>0){
                                shinyalert("           Error: please use Yes or No, or Y or N to indicate what sample to include/exclude")
                                return()
                        }
                        
                        
                        studySettings[["factorsDF"]]=hot_to_r(input$hot_Factors)
                        
                        save(studySettings, file=paste0(inPath,"/",prefix,"_studySettings.RData"))
                        #View(studySettings$samplesDF)
                        
                        shinyalert("Settings", "Study settings saved succesfully",type="success")
                }
        })        
        
        
        ### Get the input path | read data ----
        output$dirPath <- renderText({
                 ret=""
                 if (is.integer(input$file_inputData)) {
                         ret="   No .xlsx file is selected ..."
                 } else {
                        
                         ## Some global variables
                         inPathFile=parseFilePaths(volumes, input$file_inputData)
                         inputDataFile<<-inPathFile[[length(inPathFile)]]
                         
                         inPath<<-dirname(inputDataFile)
                         
                         all.ext="_HILIC_A\\.xlsx$|_HILIC_B\\.xlsx$|_LIPID\\.xlsx$|_T3pos\\.xlsx$|_PHpos\\.xlsx$|_PHneg\\.xlsx$|_PHTCA\\.xlsx$|_MOLAR_RawData\\.xlsx$|_MOLAR_Metadata\\.xlsx$|_MOLAR_CmpdNotes\\.xlsx$|_MIDA_RawData\\.xlsx$|_MIDA_Metadata\\.xlsx$|_MIDA_CmpdNotes\\.xlsx$"
                         
                         if (length(grep(pattern="\\.xlsx$", basename(inputDataFile)))==0){
                                 shinyalert("Input error",paste0("The input file must an excel file with extension '.xlsx' or 'xls'."),type="error")
                                 ret="Please choose the proper input file for your study."
                         }else if (length(grep(pattern=all.ext, basename(inputDataFile)))==0){
                                 shinyalert(paste0("The only allowed input file suffixes:\n", gsub("\\|","\n",gsub("\\\\|\\$","",all.ext))))
                         }else{
                                currModule<<-""
                                
                                if  ( length(grep(pattern="MIDA_RawData\\.xlsx$|MIDA_Metadata\\.xlsx$|MIDA_CmpdNotes\\.xlsx$", basename(inputDataFile)))> 0){currModule<<-"MIDA"}
                                if  ( length(grep(pattern="MOLAR_RawData\\.xlsx$|MOLAR_Metadata\\.xlsx$|MOLAR_CmpdNotes\\.xlsx$", basename(inputDataFile)))> 0){currModule<<-"MOLAR"}
                                if  ( length(grep(pattern=methods_regex2, basename(inputDataFile)))> 0){currModule<<-"idxMOLAR"} 
                                if  ( currModule==""){ ret=paste0(" Error : wrong input name: For MOLAR/MIDA file ends with MIDA_RawData.xlsx  or MOLAR_RawData.xlsx"); shinyalert(ret,type="error"); return(ret)}
                                
                                Dis.EnableWidgets_all(currModule_=currModule)
                                
                                
                                #idxMOLAR file parsing ---- 
                                if (currModule=="idxMOLAR"){
                                    ## Reading and checking the input files'
                                    prefix<<-strsplit(basename(inputDataFile),".xlsx")[[1]]
                                    
                                          
                                    ## Reads input files and reformat sample names removing ( "_S01", "_plasma_R","_CSF_C" etc) from sample names
                                    ## Place data in compData.A compData.B
                                    cD<<-read_INPUT_files_Idx(inPath_=inPath,  prefix_=prefix, studyType_=input$idxMolarRun.Type, studyMethod_=input$idxMolarRun.Method )
                                    if (length(cD)==1){ shinyalert(cD,type="error",size= "m" );return()}
                                    
                                    ## Reformat sample names if the replicate number if missing (_1.raw, _2.raw, _3.raw etc)
                                    ## generate "Replicate column" [ combined name based on study factors: Factor1,Factor2,Factor3 ...[eg: Fact1-Fact2-Fact3_1 ], also pick up QC and Blank sample names as they are, and save table
                                    
                                    
                                    ## Remove blank and QC and other unwanted samples names
                                    sampleNames_<-as.character(get_sampleNameList_idx(cD$compData_List[[1]]))         
                                    ixx=grep(pattern=idx_PB_Blanks_pattern,sampleNames_,ignore.case=T)  
                                    if (length(ixx>0)){sampleNames_=sampleNames_[-ixx]} 
                                    
                                    ## Get sample groups and study factors
                                    sgroups_sfactors<<-get_splitSampleNames_df_idx(sampleNames_)
                                }
                                
                                
                                # MOLAR file parsing ----        
                                if (currModule=="MOLAR"){
                                    
                                    prefix<<-gsub(pattern="_MOLAR_RawData\\.xlsx$|_MOLAR_CmpdNotes\\.xlsx$|_MOLAR_Metadata\\.xlsx$" ,"", basename(inputDataFile))
                                    
                                    cD<<-read_INPUT_files_Molar(inPath_=inPath,  prefix_=prefix)
                                    if (length(cD)==1){ shinyalert(cD,type="error",size= "m" );return()} ## if any problem with input return the error messages
                                    
                                    
                                    ret<-get_sampleNameList_Mida_Molar(InputData_=cD$InputData, metaData_=cD$metaData,currModule_=currModule)
                                    if (length(ret)==1){shinyalert(ret);return(ret)}else{cS=ret} # if any problem with sample name between data and metadata file return a error message
                                    
                                    
                                    ## Remove blank and QC and unlabeled samples names
                                    ixx1=grep("blank|qc|mc-pb|mc-db",cS$snames,ignore.case = T)  #, grep("qc",tolower(sampleNames)), grep("pos.id",tolower(sampleNames)), grep("neg.id",tolower(sampleNames)))
                                    ixx2=grep("blank|qc|mc-pb|mc-db",cS$new_snames,ignore.case = T)
                                    ixx=unique(c(ixx1,ixx2))
                                    if (length(ixx>0)){sampleNames_=cS$new_snames[-ixx]}else{sampleNames_=cS$new_snames}
                                    
                                    
                                    ## Get sample groups and study factors
                                    sgroups_sfactors<<-get_splitSampleNames_df_Mida_Molar(sampleNames_)
                                    
                                }
                                
                                
                                # MIDA file parsing ----
                                if (currModule=="MIDA"){
                                        
                                        prefix<<-gsub(pattern="_MIDA_RawData\\.xlsx$|_MIDA_CmpdNotes\\.xlsx$|_MIDA_Metadata\\.xlsx$" ,"", basename(inputDataFile))
                                        
                                        cD<<-read_INPUT_files_Mida(inPath_=inPath,  prefix_=prefix)
                                        if (length(cD)==1){ shinyalert(cD,type="error",size= "m" );return()}## if any problem with input return the error messages
                                        
                             
                                        ret<-get_sampleNameList_Mida_Molar(InputData_=cD$InputData, metaData_=cD$metaData, currModule_=currModule)
                                        if (length(ret)==1){shinyalert(ret);return(ret)}else{cS=ret} # if any problem with sample name between data and metadata file return a error message
                                                
                                        
                                        ## Remove blank and QC and unlabled samples names
                                        ixx1=grep("blank|qc|mc-pb|mc-db|unlabeled|unlab",cS$snames,ignore.case = T)  #, grep("qc",tolower(sampleNames)), grep("pos.id",tolower(sampleNames)), grep("neg.id",tolower(sampleNames)))
                                        ixx2=grep("blank|qc|mc-pb|mc-db|unlabeled|unlab",cS$slabeling ,ignore.case = T )
                                        ixx3=grep("blank|qc|mc-pb|mc-db|unlabeled|unlab",cS$new_snames,ignore.case = T)
                                        ixx=unique(c(ixx1,ixx2,ixx3))
                                        if (length(ixx>0)){sampleNames_=cS$new_snames[-ixx]}else{sampleNames_=cS$new_snames} 
                                        
                 
                                        ## Get sample groups and study factors
                                        sgroups_sfactors<<-get_splitSampleNames_df_Mida_Molar(sampleNames_)
                                        
                                }
                                
                                
                                        
                                nGroups<<-length(unique(sgroups_sfactors$group))
                                
                                
                                message("   * Number of groups detected : ",nGroups  )
                                
                                output$hot_Samples=NULL
                                output$hot_Groups=NULL
                                output$hot_Factors=NULL
                                
                                         
                                         
                                # If more than two groups detected; save settings/ view setting matrix : pairwise comparison matrix
                                if (nGroups>1){
                                                 ## get the metadata file and generate the groups and factor DFs
                                                 
                                                 studySettings<<-get_Group_Samples_Factors_DF_all(sgroups_sfactors) 
                                                 
                                                 
                                                 cat(paste0("Current working directory: ",inPath,"\n"))
                                                 save(studySettings, file=paste0(inPath,"/",prefix,"_studySettings.RData"))
                                                 
                                               
                                                 ## Samples
                                                 
                                                 rHWs=30 ## Adjust the row width
                                                
                                                 output$hot_Samples<-renderRHandsontable({
                                                         rhandsontable(studySettings$samplesDF, rowHeaderWidth=rHWs, stretchH="all")%>%
                                                                 hot_cols(renderer = "function (instance, td, row, col, prop, value, cellProperties) {
                                                                 
                                                                   Handsontable.renderers.TextRenderer.apply(this, arguments);
                                                                   
                                                                            if (col == 0) {
                                                                                      td.style.background = 'lightgrey';
                                                                                      td.style.color = 'black'
                                                                                      cellProperties.readOnly = true;
                                                                             }else{
                                                                                      td.style.background = 'white';
                                                                                      cellProperties.readOnly = false;
                                                                             }
                                                                          
                                                                   }")
                                                 })
                                                 
                                        
                                                 ## Groups
                                                 
                                                 rHWg=30+max(sapply(rownames(studySettings$groupsDF),nchar))*8 ## Adjust the row width
                                                 
                                                 output$hot_Groups<-renderRHandsontable({
                                                         rhandsontable(studySettings$groupsDF, rowHeaderWidth=rHWg,stretchH="all")%>%
                                                                 hot_cols(renderer = "function (instance, td, row, col, prop, value, cellProperties) {
                                                                Handsontable.renderers.CheckboxRenderer.apply(this, arguments);
                                                                
                                                                             if (row == col) {
                                                                                      td.style.background = 'darkgrey';
                                                                                      cellProperties.readOnly = true;
                                                                             }else{
                                                                                      td.style.background = 'white';
                                                                                      cellProperties.readOnly = false;
                                                                             }
                                                                             
                                                                   }")
                                                 })
                                                 
                                                 output$hot_Order<-renderRHandsontable({
                                                         rhandsontable(studySettings$orderDF, rowHeaderWidth=rHWg, colHeaderWidth=8)%>%
                                                                 hot_cols(renderer = "function (instance, td, row, col, prop, value, cellProperties) {
                                                                 
                                                                        Handsontable.renderers.NumericRenderer.apply(this, arguments);
                                                                        
                                                                 }")
                                                 })
                                                 
                                                 
                                                 ## Group order
                                                 
                                                 rHWf=max(sapply(rownames(studySettings$factorsDF),nchar))*9 ## Adjust the row width
                                                 #cat(rHWf,"\n")
                                                 
                                                
                                                 output$hot_Factors<-renderRHandsontable({
                                                         rhandsontable(studySettings$factorsDF, rowHeaderWidth=rHWf,stretchH="all")%>%
                                                                 hot_cols(renderer = "function (instance, td, row, col, prop, value, cellProperties) {
                                                                Handsontable.renderers.CheckboxRenderer.apply(this, arguments);
                                                                   }")
                                                         
                                                 })
                                    
                                    }else{
                                                #output$dirPath <- renderText(paste0("    * Input error: There is only ",nGroups," sample group, the analysis can not be beformed"))
                                                cat(paste0("   * Input error: There is only ",nGroups," sample group, the analysis can not be beformed\n"))         
                                                return()
                                 }
                                         
                                 ## return the full path of metadata file
                                if (currModule=="idxMOLAR"){files_=unlist(cD$prefixM_List)}else{files_=unlist(cD$prefix_List)}
                                
                                ret<-paste0("Module: ",currModule," | Files (",length(files_),"): ",paste(files_,collapse=", ")," | Path: ",inPath)
                                
                                # Remove table global variable, table   
                                remove_tb_globals()
                                output$out_Table<-DT::renderDataTable(DT::datatable({data.frame(0)}))
                                output$table_Name<-renderText("")

                                  
                         }
                         
                 }
                 
                 ret
        })


        ## Action Button for running idxMOLAR pipeline ----
        observeEvent(input$btn_idxMolar, {
             
                if (currModule==""){
                        shinyalert(paste0("Error: no input files selected!"), type = "warning")
                        
                }else if (currModule!="idxMOLAR"){
                        shinyalert(paste0("Error: wrong button, please use appropriate button to run ",currModule,"."), type = "warning")
                
                }else if (!exists("inputDataFile")) {
                        shinyalert("Input error: no file for analysis was selected, please select your raw data file to run the program.", type = "warning")
                        
                }else if (!file.exists(paste0(inPath,"/",prefix,"_studySettings.RData")) & (nGroups>2)) {
                        shinyalert("Settings error : no study settings were selected, please select you settings and save.", type = "warning")
                }else{
                       
                     
                        shinyjs::disable("file_inputData")
                        shinyjs::disable("btn_idxMolar")
                        shinyjs::disable("btn_Mida")
                        shinyjs::disable("btn_Molar")
                        shinyjs::disable("file_inputData")
                        shinyjs::disable("btn_Save")
                        shinyjs::disable("file_PDF")                    #disable some controls when the program is running
                        shinyjs::disable("file_Kable")
                        shinyjs::disable("btn_Documentation")
                        
                        shinyjs::disable("btn_MID_FE")
                        shinyjs::disable("btn_TotalLabel_FE")
                        shinyjs::disable("btn_PoolSize")
                        shinyjs::disable("btn_MID_PoolSize")
                        shinyjs::disable("btn_TotalLabel_PoolSize")
                        
                        ## Removes table global variables, and the table displayed on interface
                        remove_tb_globals()
                        output$out_Table<-DT::renderDataTable(DT::datatable({data.frame(0)}))
                        output$table_Name<-renderText("")
                        
                        
                        ## Removes log and warning files if they exist
                        logF=paste0(inPath,"/",prefix,"_log.txt")
                        
                        if (file.exists(logF)){unlink(logF)}
                        
                        ## Calls idxMOLAR pipeline 's main function ----
                        withCallingHandlers(
                                

                                mainFunIdxMOLAR(

                                        inFile_=inputDataFile,
                                        prefix_=prefix,
                                        inPath_=inPath,
                                        
                                        topMetaboFDR_=input$idxMolarTopMetaboliteFDR,
                                        normMethod_=input$idxMolarNormalizationMethod,
                                        transMethod_=input$idxMolarTransformationMethod,
                                        scalingMethod_=input$idxMolarScalingMethod,
                                        
                                        KeggPathwayAnalysisCheck_=input$idxMolarKeggPathwayAnalysisCheck,
                                        KeggPathwayEnrichMethod_= input$idxMolarKeggPathwayEnrichMethod,
                                        KeggPathwayOrganism_=input$idxMolarKeggPathwayOrganism,
                                        KeggPathwayFDR_=input$idxMolarKeggPathwayFDR,
                                       
                                        run.Type_=input$idxMolarRun.Type,
                                        run.Method_=input$idxMolarRun.Method

                                        ),

                                
                                ## Print messages to the console    ## write messages into a log file   ## write warning into a log file
                                message = function(m) { shinyjs::html("log", m$message, add=TRUE)},
                                message=function(m){write_lines(m$message, file =logF, sep="",na = "NA", append = TRUE)}
                                
                                
                        )
                        shinyjs::enable("file_inputData")
                        shinyjs::enable("btn_idxMolar")
                        shinyjs::enable("btn_Mida")
                        shinyjs::enable("btn_Molar")
                        shinyjs::enable("btn_Save")
                        shinyjs::enable("file_PDF")                    #enable some controls when the program is running
                        shinyjs::enable("file_Kable")
                        shinyjs::enable("btn_Documentation")
                        
                        if (input$idxMolarRun.Type=="profiling"){
                                shinyjs::enable("btn_PoolSize")
                        }else if (input$idxMolarRun.Type=="tracer"){
                                shinyjs::enable("btn_MID_FE")
                                shinyjs::enable("btn_TotalLabel_FE")
                                shinyjs::enable("btn_PoolSize")
                                shinyjs::enable("btn_MID_PoolSize")
                                shinyjs::enable("btn_TotalLabel_PoolSize")
                        }

                }
        })


        
        # Action Button for running MOLAR pipeline----
        observeEvent(input$btn_Molar, {
                
                if (currModule==""){
                        shinyalert(paste0("Error: no input files selected!"), type = "warning")
                        
                }else if (currModule!="MOLAR"){
                        shinyalert(paste0("Error: wrong button, please use appropriate button to run ",currModule,"."), type = "warning")        
                }else{
                        ## Main Execution
                       
                        shinyjs::disable("file_inputData")
                        shinyjs::disable("btn_idxMolar")
                        shinyjs::disable("btn_Mida")
                        shinyjs::disable("btn_Molar")
                        shinyjs::disable("file_inputData")
                        shinyjs::disable("btn_Save")
                        shinyjs::disable("file_PDF")                    
                        shinyjs::disable("file_Kable")
                        shinyjs::disable("btn_Documentation")
                        
                        shinyjs::disable("btn_MID_FE")
                        shinyjs::disable("btn_TotalLabel_FE")
                        shinyjs::disable("btn_PoolSize")
                        shinyjs::disable("btn_MID_PoolSize")
                        shinyjs::disable("btn_TotalLabel_PoolSize")
                        
                        
                        ## Removes table global variables, and the table displayed on interface
                        remove_tb_globals()
                        output$out_Table<-DT::renderDataTable(DT::datatable({data.frame(0)}))
                        output$table_Name<-renderText("")
                        
                        
                        ## Removes log and warning files if they exist
                        logF=paste0(inPath,"/",prefix,"_log.txt")
                        wrnF=paste0(inPath,"/",prefix,"_warnings.txt")
                        if (file.exists(logF)){unlink(logF)}
                      
                        ## Calls MOLAR pipeline 's main function ----
                        withCallingHandlers(
                            
                            
                            mainFunMolar(
                                
                                    prefix_=prefix,
                                    inPath_=inPath,
                                    
                                    QC.Blank.thr_=input$molarQC.Blank.thresh,
                                    NA.thr_=input$molarNA.thresh,
                                    imput.method_=input$molarImputationMethod,
                                    corr.method_=input$molarCorrectionMethod,
                                    QC.span_=input$molarQC.span,
                                    P.degree_=input$molarP.degree,
                                    
                                    preproSignalCorr.Flg_=input$molarPreprocSignalCorrection,
                                    signalCorr.Flg_=input$molarSignalCorrCheck,
                                    qc.blank.Flg_=input$molarQC.Blank.threshCheck,
                                    NA.thr.Flg_=input$molarQC.Blank.threshCheck,
                                    imputation.Flg_=input$molarImputationMethodCheck,
                                    
                                    topMetaboFDR_= input$idxMolarTopMetaboliteFDR,
                                    normMethod_= input$idxMolarNormalizationMethod,
                                    transMethod_= input$idxMolarTransformationMethod,
                                    scalingMethod_= input$idxMolarScalingMethod,
                                    
                                    KeggPathwayAnalysisCheck_=input$idxMolarKeggPathwayAnalysisCheck,
                                    KeggPathwayEnrichMethod_=input$idxMolarKeggPathwayEnrichMethod,
                                    KeggPathwayOrganism_=input$idxMolarKeggPathwayOrganism,
                                    KeggPathwayFDR_=input$idxMolarKeggPathwayFDR,
                                    run.Type_="profiling"
                                
                            ),
                            
                            
                            ## Print messages to the console    ## write messages into a log file   ## write warning into a log file
                            message = function(m) { shinyjs::html("log", m$message, add=TRUE)},
                            message=function(m){write_lines(m$message, file=logF, sep="",na = "NA", append = TRUE)}
                           
                        )
                        shinyjs::enable("file_inputData")
                        shinyjs::enable("btn_idxMolar")
                        shinyjs::enable("btn_Mida")
                        shinyjs::enable("btn_Molar")
                        shinyjs::enable("btn_Save")
                        shinyjs::enable("file_PDF")                    #enable some controls when the program is running
                        shinyjs::enable("file_Kable")
                        shinyjs::enable("btn_Documentation")
                        
                        shinyjs::enable("btn_PoolSize")
                        
                        
                                           
                }
                
                
                
        })
        
        
        # Action Button for running MIDA pipeline ----
        observeEvent(input$btn_Mida, {
            
                if (currModule==""){
                        shinyalert(paste0("Error: no input files selected!"), type = "warning")
                    
                }else if (currModule!="MIDA"){
                         shinyalert(paste0("Error: wrong button, please use appropriate button to run ",currModule,"."), type = "warning")        
                }else{
                    
                    
                        ## Disables controls during execution
                        shinyjs::disable("file_inputData")
                        shinyjs::disable("btn_idxMolar")
                        shinyjs::disable("btn_Mida")
                        shinyjs::disable("btn_Molar")
                        shinyjs::disable("file_inputData")
                        shinyjs::disable("btn_Save")
                        shinyjs::disable("file_PDF")                    
                        shinyjs::disable("file_Kable")
                        shinyjs::disable("btn_Documentation")
              
                        shinyjs::disable("btn_MID_FE")
                        shinyjs::disable("btn_TotalLabel_FE")
                        shinyjs::disable("btn_PoolSize")
                        shinyjs::disable("btn_MID_PoolSize")
                        shinyjs::disable("btn_TotalLabel_PoolSize")
                        
                        
                        ## Removes table global variable, and the table displayed on interface
                        remove_tb_globals()
                        output$out_Table<-DT::renderDataTable(DT::datatable({data.frame(0)}))
                        output$table_Name<-renderText("")
                        
                        
                        ## Removes log and warning files if they exist
                        logF=paste0(inPath,"/",prefix,"_log.txt")
                        if (file.exists(logF)){unlink(logF)}
                        
                        
                        ## Calls MIDA pipeline 's main function ----
                        withCallingHandlers(
                            
                            mainFunMida(
                                inputDataFile_=inputDataFile,
                                prefix_=prefix,
                                inPath_=inPath,
                                
                                baseLineDet_= input$midaBaseLDet,
                                FEpercThresh_ = input$midaFEpercThresh,
                                MinLabelThresh_ = input$midaMinLabelThresh,
                                
                                topMetaboFDR_=input$idxMolarTopMetaboliteFDR,
                                normMethod_=input$idxMolarNormalizationMethod,
                                transMethod_=input$idxMolarTransformationMethod,
                                scalingMethod_=input$idxMolarScalingMethod,
                                
                                KeggPathwayAnalysisCheck_=input$idxMolarKeggPathwayAnalysisCheck,
                                KeggPathwayEnrichMethod_=input$idxMolarKeggPathwayEnrichMethod,
                                KeggPathwayOrganism_=input$idxMolarKeggPathwayOrganism,
                                KeggPathwayFDR_=input$idxMolarKeggPathwayFDR,
                                
                                run.Type_="tracer"
                            ),
                            
                            
                            ## Print messages to the console    ## write messages into a log file   ## write warning into a log file
                            message = function(m) { shinyjs::html("log", m$message, add=TRUE)},
                            message=function(m){write_lines(m$message, file=logF, sep="",na = "NA", append = TRUE)}
                        )
                        
                        
                        # Enabling controls during exectution
                        shinyjs::enable("file_inputData")
                        shinyjs::enable("btn_idxMolar")
                        shinyjs::enable("btn_Mida")
                        shinyjs::enable("btn_Molar")
                        shinyjs::enable("file_inputData")
                        shinyjs::enable("btn_Save")
                        shinyjs::enable("file_PDF")                    
                        shinyjs::enable("file_Kable")
                        shinyjs::enable("btn_Documentation")
                        
                        shinyjs::enable("btn_MID_FE")
                        shinyjs::enable("btn_TotalLabel_FE")
                        shinyjs::enable("btn_PoolSize")
                        shinyjs::enable("btn_MID_PoolSize")
                        shinyjs::enable("btn_TotalLabel_PoolSize")
                        
                    
                }
            
        })
        
        
        
        # Action Button for running Seahorse Bio-Energetics pipeline ----
        observeEvent(input$btn_SeahorseBE, {
                shinyalert("MSG","The Seahorse Bio-Energetics pipeline is still under implementation.",type="info")
        })
        
        
        
        # Action Button for running Metabolomics Cage pipeline ----
        observeEvent(input$btn_MetaboCage, {
                shinyalert("MSG","The Metabolomic Cage pipeline is still under implementation.",type="info")
        })
        
        
        # Disabling and enabling widgets ----
        
        ## Disable/Enable the whole QC-correction procedure
        observeEvent(input$molarPreprocSignalCorrection,{
                if (input$molarPreprocSignalCorrection==FALSE){
                        shinyjs::disable("molarQC.Blank.thresh") 
                        shinyjs::disable("molarNA.thresh") 
                        shinyjs::disable("molarImputationMethod") 
                        shinyjs::disable("molarCorrectionMethod")  
                        shinyjs::disable("molarQC.span")  
                        shinyjs::disable("molarP.degree") 
                        
                        shinyjs::disable("molarQC.Blank.threshCheck")
                        shinyjs::disable("molarImputationMethodCheck")
                        shinyjs::disable("molarNA.threshCheck")
                        shinyjs::disable("molarSignalCorrCheck")
                        
                }else if (input$molarPreprocSignalCorrection==TRUE){
                        shinyjs::enable("molarQC.Blank.thresh")
                        shinyjs::enable("molarNA.thresh") 
                        shinyjs::enable("molarImputationMethod") 
                        shinyjs::enable("molarCorrectionMethod")  
                        shinyjs::enable("molarQC.span")  
                        shinyjs::enable("molarP.degree") 
                        
                        shinyjs::enable("molarQC.Blank.threshCheck")
                        shinyjs::enable("molarImputationMethodCheck")
                        shinyjs::enable("molarNA.threshCheck")
                        shinyjs::enable("molarSignalCorrCheck")
                }
        },ignoreInit = FALSE)
        
        ### Disable/Enable QC-blank  based filtration controls
        observeEvent(input$molarQC.Blank.threshCheck,{
                if (input$molarQC.Blank.threshCheck==FALSE){
                        shinyjs::disable("molarQC.Blank.thresh")  
                }else if (input$molarQC.Blank.threshCheck==TRUE){
                        shinyjs::enable("molarQC.Blank.thresh")  
                }
        },ignoreInit = FALSE)
        
        
        ### Disable/Enable NA  based filtration controls
        observeEvent(input$molarNA.threshCheck,{
                if (input$molarNA.threshCheck==FALSE){
                        shinyjs::disable("molarNA.thresh")  
                }else if (input$molarNA.threshCheck==TRUE){
                        shinyjs::enable("molarNA.thresh")  
                }
        },ignoreInit = FALSE)
        
        ### Disable/Enable imputation controls
        observeEvent(input$molarImputationMethodCheck,{
                if (input$molarImputationMethodCheck==FALSE){
                        shinyjs::disable("molarImputationMethod")  
                }else if (input$molarImputationMethodCheck==TRUE){
                        shinyjs::enable("molarImputationMethod")  
                }
        },ignoreInit = FALSE)
        
        
        ### Disable/Enable signal correction controls
        observeEvent(input$molarSignalCorrCheck,{
                if (input$molarSignalCorrCheck==FALSE){
                        shinyjs::disable("molarCorrectionMethod")  
                        shinyjs::disable("molarQC.span")  
                        shinyjs::disable("molarP.degree")  
                }else if (input$molarSignalCorrCheck==TRUE){
                        shinyjs::enable("molarCorrectionMethod")  
                        shinyjs::enable("molarQC.span")  
                        shinyjs::enable("molarP.degree")  
                }
        },ignoreInit = FALSE)
        
        
        
        ## Interactive plot visualization
        ## NOTE: This functionality is NOT yet implemented and be removed, if there is interest in adding interactive plot visualization 
        output$distPlot <- renderPlot({
                
                x    <- faithful$waiting
                bins <- seq(min(x), max(x), length.out = input$bins + 1)
                
                hist(x, breaks = bins, col = "#75AADB", border = "white",
                     xlab = "Waiting time to next eruption (in mins)",
                     main = "Histogram of waiting times")
                
        })

        
        ## Table visualization ----
        ## MID_FE
        observeEvent(input$btn_MID_FE, {
                if (is.null(MID_FE_global)){shinyalert("The MID_FE table not yet produced")
                    output$out_Table<-DT::renderDataTable(DT::datatable({data.frame(0)}))
                    output$table_Name<-renderText("");return()}
            
                output$out_Table <- DT::renderDataTable(DT::datatable({  
                    tt=format_DT(MID_FE_global)
                    tt
                }, 
                filter = 'top', extensions = c('Buttons', 'Scroller'),
                options = list(scrollY = 650, scrollX = 500, deferRender = TRUE, scroller = TRUE, buttons = list('excel', list(extend = 'colvis', targets = 0, visible = FALSE)),dom = 'lBfrtip',fixedColumns = TRUE)), )     
                output$table_Name<-renderText(paste0("MID_FE - ",prefix))
                })
        
        ## PoolSize
        observeEvent(input$btn_PoolSize, {
                if (is.null(PoolSize_global)){shinyalert("The PoolSize table not yet produced")
                    output$out_Table<-DT::renderDataTable(DT::datatable({data.frame(0)}))
                    output$table_Name<-renderText("");return()}
            
                output$out_Table <- DT::renderDataTable(DT::datatable({
                    tt=format_DT(PoolSize_global)
                    tt
                }, 
                filter = 'top', extensions = c('Buttons', 'Scroller'),
                options = list(scrollY = 650, scrollX = 500, deferRender = TRUE, scroller = TRUE, buttons = list('excel', list(extend = 'colvis', targets = 0, visible = FALSE)),dom = 'lBfrtip',fixedColumns = TRUE)), )    
                output$table_Name<-renderText(paste0("PoolSize - ",prefix))
                })
        
        ## TotalLabel_FE
        observeEvent(input$btn_TotalLabel_FE, {
                if (is.null(TotalLabel_FE_global)){shinyalert("The TotalLabel_FE table not yet produced")
                    output$out_Table<-DT::renderDataTable(DT::datatable({data.frame(0)}))
                    output$table_Name<-renderText("");return()}
            
            
                output$out_Table <- DT::renderDataTable(DT::datatable({
                    tt=format_DT(TotalLabel_FE_global)
                    tt
                }, 
                filter = 'top', extensions = c('Buttons', 'Scroller'),
                options = list(scrollY = 650, scrollX = 500, deferRender = TRUE, scroller = TRUE, buttons = list('excel', list(extend = 'colvis', targets = 0, visible = FALSE)),dom = 'lBfrtip',fixedColumns = TRUE)), )     
                output$table_Name<-renderText(paste0("TotalLabel_FE - ",prefix))
                })
        
        ## MID_PoolSize
        observeEvent(input$btn_MID_PoolSize, {
                if (is.null(MID_PoolSize_global)){shinyalert("The MID_PoolSize table not yet produced");return()}
                output$out_Table <- DT::renderDataTable(DT::datatable({
                    tt=format_DT(MID_PoolSize_global)
                    tt
                }, 
                filter = 'top', extensions = c('Buttons', 'Scroller'),
                options = list(scrollY = 650, scrollX = 500, deferRender = TRUE, scroller = TRUE, buttons = list('excel', list(extend = 'colvis', targets = 0, visible = FALSE)),dom = 'lBfrtip',fixedColumns = TRUE)), )     
                output$table_Name<-renderText(paste0("MID_PoolSize - ",prefix))
                })
        
        ## TotalLabel_PoolSize
        observeEvent(input$btn_TotalLabel_PoolSize, {
                if (is.null(TotalLabel_PoolSize_global)){shinyalert("The TotalLabel_PoolSize table not yet produced");return()}
                output$out_Table <- DT::renderDataTable(DT::datatable({
                    tt=format_DT(TotalLabel_PoolSize_global)
                    tt
                }, 
                filter = 'top', extensions = c('Buttons', 'Scroller'),
                options = list(scrollY = 650, scrollX = 500, deferRender = TRUE, scroller = TRUE, buttons = list('excel', list(extend = 'colvis', targets = 0, visible = FALSE)),dom = 'lBfrtip',fixedColumns = TRUE)), )     
                output$table_Name<-renderText(paste0("TotalLabel_PoolSize - ",prefix))
                })
            
        
        
        
        ## Report visualization  (It opens the pdf report)
        observe({
                req(input$file_PDF)
                output$pdfview <- renderUI({
                        #tags$iframe(style="height:700px; width:100%; scrolling=yes", src="")
                        inPathPlot=parseFilePaths(volumes, input$file_PDF)
                        
                        if(file.exists("www/Plot.pdf")) {file.remove("www/Plot.pdf")}
                        
                        file.copy(inPathPlot[[length(inPathPlot)]],"www/Plot.pdf", overwrite = T)
                        
                        if(file.exists("www/Plot.pdf")){
                                tags$iframe(style="height:700px; width:100%; scrolling=yes", src="Plot.pdf")
                        }
                        
                })
                
        })
        
       
        ## View the  the documentation
        observeEvent(input$btn_Documentation, {
                output$docview <-  renderUI({tags$iframe(style="height:700px; width:100%; scrolling=yes", src="EMeRGe.pdf")})
        })
        
        
        
        
        ## Cancel button to interrupt the analysis
        observeEvent(input$btn_Cancel, {
                shinyalert("MSG","This functionality still under implementation.",type="info")
        }
        
        )


        ## Automatically stop a Shiny app when closing the browser tab
        session$onSessionEnded(stopApp)

}

