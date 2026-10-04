
#### ui.R : The code build the graphical user interface of the application ##### ----

## Import R libraries
library(shinyWidgets)
library(shinyalert)
library(shinyFiles)
library(rhandsontable)
library(fs)

## Java Script function
js <- "
$(document).ready(function(){
          var objDiv = document.getElementById('log');
          // create an observer instance
          var observer = new MutationObserver(function(mutations) {
                        objDiv.scrollTop = objDiv.scrollHeight - objDiv.clientHeight;
          });
          // configuration of the observer
          var config = {childList: true};
          // observe objDiv
          observer.observe(objDiv, config);
})
"

ui<-fluidPage(
        #includeCSS('www/mycss.css'),
        useShinyalert(),
        shinyjs::useShinyjs(),
        tags$head(tags$script(HTML(js))),

        ## Header ----
        h5("VAI Metabolism -  Energetics and MEtabolomics Report GEnerator (EMeRGe v1.0)" ,style='padding:5px; color:  #0046E7',img(src="VAI_Logo.jpg" ,align = "right",height=30,width=105)),

        tabsetPanel(

                # Console & Command buttons ----
                tabPanel("Console",
                        
                         br(),

                         ## Command buttons
                         fluidRow(
                                 column(12,align="center",
                                        shinyFilesButton("file_inputData", "Open file", "Please an input raw dat file", multiple = FALSE,style='padding:5px; font-size:80%; color: #0046E7;width: 90px'),
                                        actionButton("btn_idxMolar", "idxMOLAR",style='padding:5px; font-size:80%; color:  #0046E7; width: 90px'),
                                        actionButton("btn_Molar", "MOLAR",style='padding:5px; font-size:80%; color:  #0046E7; width: 90px'),
                                        actionButton("btn_Mida", "MIDA",style='padding:5px; font-size:80%; color:  #0046E7; width: 90px'),
                                        actionButton("btn_SeahorseBE", "SeahorseBE",style='padding:5px; font-size:80%; color:  #0046E7; width: 90px'),
                                        actionButton("btn_MetaboCage", "MetaboCage",style='padding:5px; font-size:80%; color:  #0046E7; width: 90px'),
                                        actionButton("btn_Cancel", "Cancel",style='padding:5px; font-size:80%; color:  #0046E7; width: 90px')
                                 )
                         ),
                        
                        ## Console
                        tags$style(HTML("#log{height:700px; font-size:85%; background-color: #0046E7; color:#ffffff}")), # VAI Blue #003882 #00265d
                        
                        verbatimTextOutput("log",placeholder=T),
                        
                        ## Working directory textbox 
                        tags$style(HTML("#dirPath{font-size:75%; height:100%; text-align:left; color:#0046E7}")),
                        verbatimTextOutput("dirPath", placeholder = TRUE),
                        icon = icon("laptop")
                ),

                
                ## General Settings tab ----
                tabPanel("General settings",

                         sidebarPanel(

                               tags$style(".well {background-color:#0046E7; color:#ffffff;}"), ##00265d
                               h4(tags$b("Description:")),
                               
                               h4(tags$b("MIDA")),
                               br(),
                               tags$b("Correction method"),
                               p("This allows to use labeled samples or theoretical mass isotopologue distribution values"),
                               
                               tags$b("Baseline detection level (default=1000)"),
                               p("Base line detection threshold at which compounds are filtered out before analysis. By default the method excludes compounds 
                      if 2/3 of unlabeled samples [M+0] responses are below 1000, or if all labeled sample entries are < 1000"),
                               
                               tags$b("Unlabeled FE % threshold (2/3)"),
                               p("FE percentage threshold (default 5), it used to filter compounds with unlabeled samples with  FE <95 or >105, after correction for natural abundance"),
                                    
                               tags$b("Min % labeling threshold (2/3)"),
                               p("Minimum threshold for total labeling (sum 13C or sum of M+1,M+2, ..., M+n) for a compound to be reported"),          
                               br(),
                         
                               h4(tags$b("MOLAR")),
                               
                               tags$b("QC/Blank filtration threshold"),
                               p("Removes metabolites whose average QC response is less than “X” times the average blank response. Uncheck the box to skip QC/Blank filtration."),
                               
                               tags$b("Per group missing values filtration threshold"),
                               p("Removes metabolites which contain “X”% or more missing values in any group. Uncheck the box to turn off missing value filtration."),
                               
                               tags$b("Imputation method"),
                               p("Sets the method to impute response for missing values. Uncheck the box to skip imputation."),
                               
                               tags$b("Signal drift correction method"),
                               p("Select signal correction method for within batch signal drift correction. LOESS signal drift performs best with at least 8-10 QC runs, if many QC runs are available, Random Forest may provide improved correction. Uncheck the box to skip signal drift correction."),
                               
                               
                               tags$b("QC Span"),
                               p("If using LOESS Signal Correction, select a value between 0.5-0.75 which indicates the amount of smoothing to be applied to the correction function. A higher value increases smoothing, but requires more QC data points."),
                               
                               tags$b("Polynomial degree"),
                               p("If using LOESS signal correction method, it determines the type of the model to use (ie. Linear = 1, Quadratic = 2 etc.)."),
                               
                               br(),
                               
                               h4(tags$b("Statistical Analyses setting")),
                               
                               tags$b("FDR threshold (ANOVA, T-test)"),
                               p("Sets the threshold for metabolites significantly different between groups (ANOVA or T-test)"),
                                 
                               tags$b("Row normalisation"),
                               p("Sets the data normalization method"),
                                 
                               tags$b("Data transformation"),
                               p("Sets the data transformation method "),
                                 
                               tags$b("Data scaling"),
                               p("Sets the data scaling method ")
        
                               
                        ),

                        mainPanel(
                           
                           column(6,
                                  
                                  ## MIDA specific settings ----
                                  br(),
                                  fluidRow(
                                    
                                    column(8,
                                      h4("MIDA - Isotope natural abundance correction"),
                                      radioButtons("midaAnalysis", "Correction method:",
                                                   c("Using Unlabeled sample" = "Unlabeled samples","Theoretical MIDs [Not yet implemented]" = "Theoretical MIDs")
                                      ), 
                                      
                                      numericInput("midaBaseLDet", label="Baseline detection:", value = 1000,  min = 500, max = 3000, step = 100),
                                      numericInput("midaFEpercThresh", label="Unlabeled FE % threshold (2/3):", value = 5,  min = 5, max = 20, step = .5),
                                      numericInput("midaMinLabelThresh", label="Min % labeling threshold  (2/3):", value = 5,  min = 0, max = 20, step = .5)
                                    
                                    )
                                  ),
                                  
                                  
                                  ## MOLAR specific settings ----
                                  fluidRow(
                                    column(8,h4("MOLAR - QC-based signal correction (StatTarget)")),
                                    column(1,checkboxInput("molarPreprocSignalCorrection",NULL, value = TRUE, width = NULL))
                                    
                                  ),
                                  
                                  fluidRow(
                                    column(8,numericInput("molarQC.Blank.thresh", label="QC/Blank filtration threshold:", value = 5,  min = 0, max = 20, step = 1)),
                                    column(1,checkboxInput("molarQC.Blank.threshCheck",NULL, value = TRUE, width = NULL))
                                  ),
                                  
                                  fluidRow(
                                    
                                    column(8, numericInput("molarNA.thresh", label="Per group missing value filtration threshold:", value = .35,  min = .1, max = .50, step = .05)),
                                    column(1,checkboxInput("molarNA.threshCheck",NULL, value = TRUE, width = NULL))
                                  ),
                                  
                                  fluidRow(
                                    column(8,radioButtons("molarImputationMethod", "Imputation method:",
                                                          c("Half of the minimum" = "minHalf","Minimum" = "min","Median" = "median","K-nearest neigbors" = "knn"))),
                                    column(1,checkboxInput("molarImputationMethodCheck",NULL, value = TRUE, width = NULL))
                                  ),
                                  
                            
                                  fluidRow(
                                    column(8,
                                           radioButtons("molarCorrectionMethod", "Signal drift correction method:",
                                                        c("QC-based Robust Loess Regression" = "QCRLSC", "QC-based Random Forest Classifier" = "QCRFSC")
                                           ),
                                           
                                           numericInput("molarQC.span", label="QC span", value = .75,  min = .50, max = .80, step = .05),
                                           numericInput("molarP.degree", label="Polynomial degree", value = 1,  min = 1, max = 2, step = 1)
                                    ),
                                    column(1,checkboxInput("molarSignalCorrCheck",NULL, value = TRUE, width = NULL)),
                                  )
                                  
                                  
                                  
                           ),  
                           
                           
                           ## Statistical analysis with MetaboAnalyst ----
                           column(6,
                                  br(),
                                  h4("Statistical analyses (MetaboAnalystR)"),
                                  selectInput("idxMolarTopMetaboliteFDR", label="FDR threshold (ANOVA, T-test)", choices=list(0.1, 0.05,0.01,0.005,0.001,0.0005,0.0001,0.00005,.00001,0.000005,.000001), selected =0.1, multiple = FALSE, selectize = TRUE, width = NULL, size = NULL),
                        
                                  radioButtons("idxMolarNormalizationMethod", "Row normalization:",
                                               c("None"="None","Quantile Normalization" = "QuantileNorm","Probabilistic Quotient Normalization without ref sample" = "ProbNormT","Probabilistic Quotient Normalization with ref sample" = "ProbNormF",
                                                 "Normalization by a reference feature" = "CompNorm", "Normalization to constant sum"="SumNorm", "Normalization to sample median" = "MedianNorm", "Normalization by a sample-specific factor"="SpecNorm")
                                  ),

                                  radioButtons("idxMolarTransformationMethod", "Data transformation:",
                                               c("Log transformation"="LogNorm","Cubic Root Transformation" = "CrNorm")
                                  ),

                                  radioButtons("idxMolarScalingMethod", "Data scaling:",
                                               c("Mean-centering"="MeanCenter", "Autoscaling"="AutoNorm","Pareto scaling" = "ParetoNorm","Range scaling"="RangeNorm")
                                  )

                           )

                         ),
                         icon = icon("sliders-h")
                ),

                
                
                ## Study Settings tab ----
                tabPanel("Study settings",
                                
                         sidebarPanel(
                                   tags$style(".well {background-color:#0046E7; color:#ffffff;}"), ##00265d
        
                                   h4(tags$b("Description:")),
                                   
                                   
                                   tags$b("IDX – Method:"),
                                   
                                   p("Indicate if input file contains data from multiple analytical methods or only one."),
                                   
                                   tags$b("IDX – Study type:"), 
                                   
                                   p("Indicate if the input data is from a metabolite profiling experiment or tracer analysis."), 
                                   
                                   
                                   tags$b("Group pairwise  comparisons"),
                                   p("Selects the group pairwise comparisons to be perfomed on the input dataset."),
                          
                                   tags$b("Group order"), 
                                   p("- This table allows to exclude any specific group by changing Y to N (case insensitive)"),
                                   p("- The “order” column allows changing the ordering of group on final plots."), 
                                   
                                   tags$b("Sample names and normalization factors"), 
                                   p("-	Enter “N” if one or more samples should be excluded from analysis (for samples which may be an outlier). Each sample group must have a minimum of two biological replicates. Values for removed samples are imputation using average of the other replicates in the group."),
                                   p("- Normalization factor for one or more samples (e.g. used to account for number of cells."),
                                   
                                   tags$b("Study factor analysis"), 
                                   p("Select combinations of study factors which should be compared to each other or analyzed separately."), 
                                   
                                   tags$b("Metabolite Set Enrichment Analysis"),
                                   p(tags$u("(QEA) Quantitative Enrichment Analysis)"),"- metabolite set enrichment analysis is based on the globalTest algorithm."),
                                   p(tags$u("ORA (Over Representation Analysis)"),"- metabolite set over representation test is based on cumulative hypergeometic distribution."),
                                   
                                   tags$b("Top enriched Kegg pathway (FDR)"),
                                   p("Sets the FDR threshold for enriched pathways between two groups."),
                             
                                   tags$b("Organism"),
                                   p("Sets the organism under study to get the correct Kegg pathway libary."),


                         ),

                         mainPanel(
                                  br(),
                                  radioButtons("idxMolarRun.Method", "IDX - Method:",
                                                c("Multiple methods (HILIC_A, HILIC_B ...)"="MultiMethods","One method" = "OneMethod")
                                  ),
                              
                                  radioButtons("idxMolarRun.Type", "IDX - Study type:",
                                                c("IDX Metabolite profiling"="profiling","IDX Tracer analysis" = "tracer")
                                  ),
                                   
                                  br(),
                                  actionButton("btn_Save", "Save"),
                                   
                                   
                                  # Group pairwise comparison
                                  fluidRow(
                                        column(8,h4("Group pairwise comparison"))
                                  ),
                                  rHandsontableOutput("hot_Groups"),
                                   
                                  
                                  # Group order
                                  fluidRow(
                                        column(8,h4("Group order"))
                                  ),
                                  rHandsontableOutput("hot_Order"),
                                  
                                  
                                  #  Sample names and normalization factor
                                  fluidRow(
                                        column(8,h4("Sample names and normalization factors (Only for PoolSize)"))
                                  ),
                                  rHandsontableOutput("hot_Samples"),
                                  
                                  
                                  fluidRow(
                                        column(8,h4("Study factor analysis"))
                                  ),
                                  rHandsontableOutput("hot_Factors"),
                                              
                                  br(),
                                  fluidRow(
                                    column(1,   checkboxInput("idxMolarKeggPathwayAnalysisCheck", NULL,value = TRUE, width = NULL), align="left"),
                                    column(8,   h4("Pathway Analysis"))
                                  ), 
                                  selectInput("idxMolarKeggPathwayEnrichMethod", "Metabolites Set Enrichment Analysis", choices=list("QEA","ORA"), selected = "QEA", multiple = FALSE, selectize = TRUE, width = NULL, size = NULL),
                                  
                                  selectInput("idxMolarKeggPathwayFDR", label="Top enriched Kegg pathways (FDR)", choices=list(0.1, 0.05, 0.01, 0.005,0.001,0.0005,0.0001,0.00005,.00001,0.000005,.000001), selected = 0.1, multiple = FALSE, selectize = TRUE, width = NULL, size = NULL),
                                   
                                  selectInput("idxMolarKeggPathwayOrganism", "Organism", choices=list("mmu", "hsa",""), selected = "mmu", multiple = FALSE, selectize = TRUE, width = NULL, size = NULL),
                                   
                                  br()
                                   
                         ),
                         icon = icon("sliders-h")
                ),

                ## View tables  ----
                tabPanel("Tables",
                         
                         br(),

                         column(12,align="center",
                                actionButton("btn_MID_FE", "MID FE",style='padding:5px; font-size:80%; color: #0046E7; width: 90px'),
                                actionButton("btn_TotalLabel_FE", "T.Label FE",style='padding:5px; font-size:80%; color: #0046E7; width: 90px'),
                                actionButton("btn_PoolSize", "PoolSize",style='padding:5px; font-size:80%; color: #0046E7; width: 90px'),
                                actionButton("btn_MID_PoolSize", "MID PoolSize",style='padding:5px; font-size:80%; color: #0046E7; width: 90px'),
                                actionButton("btn_TotalLabel_PoolSize", "T.Label PoolSize",style='padding:5px; font-size:80%; color: #0046E7; width: 90px'),
                                br(),
                                br(),
                                tags$style(HTML("#table_Name{height:100%, font-weight: bold , background-color:#FFFFFF}")), 
                                verbatimTextOutput("table_Name", placeholder = FALSE)
                         ),
                         
                         DT::dataTableOutput("out_Table"),
                         
                         icon = icon("fal fa-table")
      
                ),

                
                # View pdf report  ----
                tabPanel("Report",
                         
                         column(12,align="center",
                                shinyFilesButton("file_PDF", "Select PDF", "Please select a file", multiple = FALSE,style='padding:5px; font-size:80%; color:  #0046E7;width: 90px')
                         ),
                         
                         uiOutput("pdfview"),
                         
                         icon = icon("far fa-file-pdf")
                ),
                
                
                ## View manual  ----
                tabPanel("Manual",
                         
                         column(12,align="center",
                                actionButton("btn_Documentation", "Open",style='padding:5px; font-size:80%; color: #0046E7; width: 90px')
                         ),
                         uiOutput("docview"),
                         
                         #tags$iframe(style="height:700px; width:100%; scrolling=yes", src="EMeRGe.pdf"),
                         
                         icon = icon("book")
                )
        )
)