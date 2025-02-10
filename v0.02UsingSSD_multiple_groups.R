require("fastqcr") # ok
require("DESeq2") # ok
require("tximport") # ok
require("DT") # ok
require("shiny") # ok 
require("rhandsontable") # ok
require("shiny") # ok 
require("shinyWidgets") #ok
require("shinyjs") #ok 
require("shinyFiles") #ok 
require("ggplot2") # ok
require("tibble") #ok
require("dplyr") #ok
require("patchwork") #ok
require("shinymanager") #ok
require("shinyBS") #ok
require("shinyscreenshot") #ok
require("shinyalert") #Ok
require("decoupleR") #ok
require("tibble") #ok
require("tidyr") #ok
require("pheatmap") #ok
require("ggrepel") #ok
require("OmnipathR") #ok 
require("clusterProfiler") #ok 
require("openxlsx") #ok
require(enrichR)
require(crosstalk)
#require('org.Hs.eg.db')
#require('org.Mm.eg.db')
require(plotly)
options(shiny.maxRequestSize=10000*1024^2)

salmon <- 'salmon '
bbduk.sh <- 'bbduk.sh '
#fastqc_conda <- '/home/zakiev/miniconda3/bin/conda run -n bioinfo fastqc '
fastqc_powershell <- "fastqc -t 12 "

credentials <- data.frame(
  user = c("shiny", "shinymanager", "emma", 'fabrice', 'aurelia','emile','maeva','aneta','elodie','manon','cedric','yann','theo'), # mandatory
  password = c("flaviallab", "flavial", "emma", 'fabrice','aurelia','emile','maeva','aneta','elodie','manon','cedric','yann','theo'), # mandatory
  #start = c("2019-04-15"), # optinal (all others)
  expire = c(rep(NA, 13)),
  admin = c(rep(NA, 13)),
  comment = "Simple and secure entification mechanism 
  for single ‘Shiny’ applications.",
  stringsAsFactors = FALSE
)

extractNamesFromFastqcrReports <- function(lst){
  toRet <- NULL;
  for (i in 1:length(lst))
    toRet <- c(toRet, as.character(lst[[i]]$basic_statistics[1,2]))
  return(toRet)
}

fcListFromDF <- function(x, decreasing =T){
  toRet <- x$log2FoldChange
  names(toRet) <- x$SYMBOL
  return(sort(toRet, decreasing))
}

ConvertRownamesEnsembl2Symbol <- function(matr, reference_df=gene2tx2name_GRCm39_u_systematic_u_onlyGeneNames){
  possibleGenes <- reference_df[rownames(matr),]
  droppedNAs <- possibleGenes[which(!is.na(possibleGenes[,2])),]
  matr <- matr[which(!is.na(possibleGenes[,2])),]
  uniqs_genes <- droppedNAs[which(!duplicated(droppedNAs[,2])),]
  matr <- matr[which(!duplicated(droppedNAs[,2])),]
  rownames(matr) <- uniqs_genes[,2]
  return(matr)
}

ConvertEnsemblStringToSymbols <- function(vec, reference_df=gene2tx2name_GRCm39_u_systematic_u_onlyGeneNames){
  vv <- strsplit(vec, split = "/")[[1]]
  if(length(vv)>1){
    vv2 <- unique(reference_df[vv,2])
    return(paste(vv2, collapse = "/"))
  } else {
    return(reference_df[vec,2])
  }
}

Convert8thColumnBackToSymbols <- function(mat, reference_df=gene2tx2name_GRCm39_u_systematic_u_onlyGeneNames){
  mat[,8] <- sapply(mat[,8], ConvertEnsemblStringToSymbols, reference_df)
  return(mat)
}


Convert11thColumnBackToSymbols <- function(mat, reference_df=gene2tx2name_GRCm39_u_systematic_u_onlyGeneNames){
  mat[,11] <- sapply(mat[,11], ConvertEnsemblStringToSymbols, reference_df)
  return(mat)
}

strsplits <- function(x, splits, ...)
{
  for (split in splits)
  {
    x <- unlist(strsplit(x, split, ...))
  }
  return(x[!x == ""]) # Remove empty values
}

shiny_home <- '/windows/Users/flavial/Documents/Shiny_0.9.0/ShinyData/'
house <- '/media/minicluster/Data/FASTQ/'

ui <- fluidPage(
  #verbatimTextOutput("auth_output"),
  fluidRow(
    h1(strong("File Uploading (optional)"), style = "font-size:30px;"),
    column(width = 3,
           textInput(inputId = "whereToSavefastq", label = "Provide a name for the folder where to save your fastq files", value = "")),
    column(width = 3,
           fileInput(inputId = "uploadFASTQ", label ="Upload your own FASTQ files", multiple = T, placeholder = '.gz or .fastq.gz files', accept = c('.gz', '.fastq'))),
    column(width = 1, tags$img(src = 'Possible_logo_30pct.png')),
    style='margin-bottom:0px;border:1px solid; padding: 10px;'
  ),
  fluidRow(
    column(width = 3,
           h1(strong("Analysis"), style = "font-size:30px;"),
           textInput(inputId = "ProjectName", label = "A new project name?", value = "..."),
           pickerInput(
             inputId = "referenceGenomeChoice", 
             #label = "Select your organism",
             choices = c("Human GRCh38"=1,"Mouse GRCm39"=2,"Mouse C57B6"=3,"Mouse MM129"=4), selected = 2,
             options = list(
               style = "btn-danger")
           ),
           shinyFilesButton(id = 'infiles', label = 'Select .fastq files', title = 'Please select fastq or fastq.gz
                       files', multiple = T)
    ),
    #column(width = 1, br(), "or"),
    
    column(width = 2, h1(strong("Load previous analysis"), style = "font-size:30px;"), br(), shinyDirButton(id = 'infolder', label = 'Load previous project?', title = 'Please select a 
                       folder')),
    column(width = 6, br(), verbatimTextOutput('previousprojname')),
    style='margin-bottom:30px;border:1px solid; padding: 10px;'
  ),
  useShinyjs(),
  textOutput("selectedfiles_text"),
  verbatimTextOutput("selectedfiles"),
  textOutput("filenameDiagnostics"),
  textOutput("filenameDiagnosticsWarning"),
  verbatimTextOutput("selectedfiles_multiplexedReads_and_salmon_pairedEndReads"),
  br(),
  textOutput("GroupsPromptText"),
  textOutput("GroupsPromptText2"),
  textOutput("GroupsPromptText3"),
  br(),
  fluidRow(column(width = 6, rHandsontableOutput('GroupsPrompt')),
           column(width = 6,  textOutput("GroupsPromptTextFork"),
                  textInput(inputId = 'folder2fork2', label = 'Please specify new folder name'),
                  uiOutput('ForkingPrompt'),
                  actionButton(inputId = 'validateforking', label = 'Validate choice'),
                  verbatimTextOutput('forkingfeedback')
           )),
  br(),
  uiOutput("button_save_colDatt"),
  uiOutput("button_rebase_colDatt"),
  ### the stuff below this point is after the button "confirm"
  
  tableOutput('dynamic_log'),
  
  textOutput("fastqc_finished_text"),
  
  br(),
  textOutput("Fastqc_Stats_Name"),
  tabsetPanel(id='fastqcstatspanel', type = "tabs", 
              tabPanel(title = "General stats", DT::DTOutput("fastqc_stats")),
              tabPanel(title = "Per base sequence quality", plotOutput("fastqc_sequence_q", height = '2000px')),
              tabPanel(title = "Overrepresented sequences", DT::DTOutput("fastqc_overrepseq"))#,
              #tabPanel(title = "Adapter Content", plotOutput("fastqc_adapters"))
  ),
  br(),
  textOutput("Aftertrimming"),
  br(),
  tabsetPanel(id='fastqcstatspanel_AfterTrimming', type = "tabs",
              tabPanel(title = "General stats", DT::DTOutput("fastqc_stats_AfterTrimming")),
              tabPanel(title = "Per base sequence quality", plotOutput("fastqc_sequence_q_AfterTrimming", height = '2000px')) ,
              tabPanel(title = "Overrepresented sequences", DT::DTOutput("fastqc_overrepseq_AfterTrimming"))#,
              #tabPanel(title = "Adapter Content", plotOutput("fastqc_adapters_AfterTrimming"))
  ),
  br(),
  br(),
  textOutput("mapping_rates_header"),
  textOutput("mapping_rates_subheader"),
  #DT::DTOutput("salmon_log"),
  rHandsontableOutput("salmon_log"),
  br(),
  #textOutput("Analyzing"),
  br(),
  textOutput("text_TPM"),
  br(),
  tags$head(tags$style(HTML(".search {float: right;}"))),
  br(),
  tags$input(type = "text", id = "mySearch", placeholder = "Search"),
  DT::DTOutput("TPM"),
  br(),
  br(),
  br(),
  textOutput("text_DEGs"),
  textOutput("text_which_is_Control_in_DEG"),
  #uiOutput('multitab_DEG'),
  DT::DTOutput("DESeq_DEGs"),
  uiOutput('DESeq_DEGsMultitab'),
  textOutput("PCA_title"),
  #plotOutput("PCA", height = "900px"),
  uiOutput('PCA_2tab'),
  #imageOutput("PCA_image", height = "100%", width = "100%"),
  textOutput("Volcano_title"),
  plotlyOutput("Volcano", height = '1000px'),
  uiOutput('VolcanoMultitab', height='1000px`'),
  
  br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),
  
  textOutput("Enrichments_title"),
  br(),
  #tabsetPanel(id='Enrichments_panel', type = "tabs",
  #            tabPanel("Dotplots", 
  #                     tabsetPanel(id = 'tabset1', tabPanel("by overrepresentation test", uiOutput('GO_dotplot_multitab'), plotOutput("GO_dotplot", height = '2000px')), 
  #                                 tabPanel("by GSEA",plotOutput("GSEA_dotplot", height = '2000px')))),
  #            tabPanel("Tables of terms", 
  #                     tabsetPanel(id = 'tabset2', tabPanel("by overrepresentation test", uiOutput('GO_multitab'),  DT::DTOutput("GO")), 
  #                                 tabPanel("by GSEA",DT::DTOutput("GSEA"))))
  #),
  tabsetPanel(id='Enrichments_panel', type = "tabs",
              tabPanel("Dotplots",
                       textInput('showCategory', label = 'Show how many top enriched terms?', value = "50", width='15%'),
                       uiOutput('GO_dotplot_multitab'), plotOutput("GO_dotplot", height = '2000px')),
              tabPanel("Tables of terms", uiOutput('GO_multitab'),  DT::DTOutput("GO"))
  ),
  br(),
  textOutput("Dorothea_title"),
  uiOutput("DorotheaMultitab"),
  plotOutput("Dorothea", height = "750px", width = "750px"),
  br(),
  
  textOutput("GSVA_title"),
  fluidRow(
    column(8, textInput('genes_for_gsva', label = '', value = "", placeholder = 'Enter/paste a list of genes in SYMBOL and/or ENSEMBL formats', width='100%')),
    column(1, div( style = "margin-top: 20px;", actionButton(inputId = "runGSVA", "Calculate")))
  ),
  verbatimTextOutput("GSVA_genes_matching"),
  plotOutput("GSVAplot", height = '1200px'),
  DT::DTOutput("GSVAtable"),
  
  
  textOutput("enrichR_title"),
  fluidRow(
    column(8, textInput('genes_for_enrichR', label = '', value = "", placeholder = 'Enter/paste a list of genes in SYMBOL and/or ENSEMBL formats', width='100%')),
    column(1, div( style = "margin-top: 20px;", actionButton(inputId = "runenrichR", "Calculate")))
  ),
  verbatimTextOutput("enrichR_genes_matching"),
  shinycssloaders::withSpinner(DT::DTOutput("enrichRtable")),
  
  #downloadButton("report", "Generate report")
  #screenshotButton(download = TRUE, id = "plot", filename = "poopity.scoop")
  #actionButton("toSnap", "Take a snapshot of the current state"),
  downloadButton("downloadData", "Download"),
  downloadButton("downloadDotPlot", "Download"),
  downloadButton("downloadDotPlots", "Download"),
  
  tags$head(tags$style("#salmon_path{
                           border: dotted grey; border-radius: 15px; 
  }
                        #bbduk_path{
       border: dotted grey; border-radius: 15px; 
                        }
                      #fastqc_finished_text{color: black;
                                 font-size: 22px;
                                 /*font-style: italic;*/
                                 font-weight: bold;
                                 }
                       #Aftertrimming{color: black;
                                 font-size: 22px;
                                 /*font-style: italic;*/
                                 font-weight: bold;
                                 }
                       #Analyzing{color: black;
                                 font-size: 22px;
                                 /*font-style: italic;*/
                                 font-weight: bold;
                                 }
                       #GroupsPromptText{color: red;
                                 font-size: 22px;
                                 /*font-style: italic;*/
                                 font-weight: bold;
                       }
                      #GroupsPromptText2{color: black;
                                 font-size: 14px;
                                 /*font-style: italic;*/
                                 font-weight: bold;
                      }
                      #GroupsPromptText3{color: grey;
                                 font-size: 14px;
                                 /*font-style: italic;*/
                      }  
                                 #GroupsPromptTextFork{color: black;
                                 font-size: 22px;
                                 /*font-style: italic;*/
                                 /*font-weight: bold;*/
                       }
                                 
                       #text_TPM{color: black;
                                 font-size: 22px;
                                 /*font-style: italic;*/
                                 font-weight: bold;
                       }
                       #PCA_title{color: black;
                                 font-size: 22px;
                                 font-weight: bold;
                      }
                      #or{color: black;
                                 font-size: 16px;
                      }
                      #or2{color: black;
                                 font-size: 16px;
                      }
                      #orSalmon{color: black;
                                 font-size: 16px;
                      }
                                 #button_salmon{color: black;
                                 font-size: 22px;
                                 /*font-style: italic;*/
                                 font-weight: bold;}
                                 #mapping_rates_header{color: black;
                                 font-size: 22px;
                                 /*font-style: italic;*/
                                 font-weight: bold;}
                                 #mapping_rates_subheader{color: grey;
                                 font-size: 16px;
                                 /*font-style: italic;*/
                                 font-weight: bold;}
                                 #text_DEGs{color: black;
                                 font-size: 22px;
                                 /*font-style: italic;*/
                                 font-weight: bold;}
                                 #Enrichments_title{color: black;
                                 font-size: 22px;
                                 /*font-style: italic;*/
                                 font-weight: bold;}
                                 #Dorothea_title{color: black;
                                 font-size: 22px;
                                 /*font-style: italic;*/
                                 font-weight: bold;}
                                 #GSVA_title{color: black;
                                 font-size: 22px;
                                 /*font-style: italic;*/
                                 font-weight: bold;}
                                 #enrichR_title{color: black;
                                 font-size: 22px;
                                 /*font-style: italic;*/
                                 font-weight: bold;}
                                 #Volcano_title{color: black;
                                 font-size: 22px;
                                 /*font-style: italic;*/
                                 font-weight: bold;}
                                 #textDetectedFastQC{color: black;
                                 font-size: 16px;}
                                 #fastqc_sequence_q_caption{
                                 color: grey;
                                 font-size: 12px;
                                 }
                                 #fastqc_after_trimming_sequence_q_caption{
                                 color: grey;
                                 font-size: 12px;
                                 }
                                 #salmon_log_where_saved{
                                 color: grey;
                                 font-size: 12px;
                                 }
                                 #text_which_is_Control_in_DEG{
                                 color: grey;
                                 font-size: 16px;
                                 }
                                 #filenameDiagnostics{color: black;
                                 font-size: 16px;
                                 /*font-style: italic;*/
                                 font-weight: bold;}
                                 #filenameDiagnosticsWarning{color: red;
                                 font-size: 14px;
                                 }
                       ")
  )
  
  
  
  
)


callback <- '
$("div.search").append($("#mySearch"));
$("#mySearch").on("keyup redraw", function(){
  var splits = $("#mySearch").val().split(" ").filter(function(x){return x !=="";})
  var searchString = "(" + splits.join("|") + ")";
  table.search(searchString, true).draw(true);
});
'

server <- function(input, output, session){
  tempStorage <- tempdir()
  
  observe({
    if (is.null(input$uploadFASTQ)) return()
    usr <- res_auth$user
    print("printin input$uploadFASTQ")
    print(input$uploadFASTQ)
    print("printin input$uploadFASTQ$datapath")
    print(input$uploadFASTQ$datapath)
    print("printin destination files")
    fold <- input$whereToSavefastq
    dest_dir <- paste0(house, usr, '/', fold)
    if(!dir.exists(dest_dir)) dir.create(dest_dir, recursive = T)
    #dest <- paste0(house, usr, '/',basename(input$uploadFASTQ$datapath))
    dest <- paste0(dest_dir, '/',input$uploadFASTQ$name)
    print(dest)
    file.copy(input$uploadFASTQ$datapath, dest, overwrite = T)
    file.remove(input$uploadFASTQ$datapath)
  })
  print("printing input")
  print(input)
  now <- gsub(pattern = ":", replacement = "-", x = gsub(pattern = " ", replacement = "_", Sys.time()))
  
  observe({
    #usr <- session$user
    #usr <- Sys.getenv("LOGNAME")
    usr <- res_auth$user
    toProjName <- paste0(usr,"_analysis_",now)
    updateTextInput(session = session, inputId = "ProjectName", value = toProjName)
  })
  
  proceedWithLoad <- reactiveVal(0)
  referenceGenomeChoice <- reactiveVal(2)
  
  observeEvent(input$referenceGenomeChoice, handlerExpr = {
    referenceGenomeChoice(input$referenceGenomeChoice)
  })
  
  output$previousprojname <- renderText({
    req(input$infolder)
    if(is.null(input$infolder)) return({})
    inFolder <- paste0(house,paste0(input$infolder[[1]], collapse = "/"))
    print('line 321 printing inFolder')
    print(inFolder)
    objfile <-paste0(inFolder, "/txi.RDS")
    
    print('looking for the txi.RDS at destination:')
    
    if(file.exists(objfile)){
      toOutput <- as.character(paste0('Loading project ', paste0(input$infolder[[1]], collapse = "/")))
      proceedWithLoad(1)
      fastqc_finished(1)
      bbduk_finished(1)
      val <- readRDS(paste0(house, paste0(input$infolder[[1]], collapse = "/"), '/referenceGenomeChoice.RDS'))
      print('printing the stored reference genome index')
      print(val)
      referenceGenomeChoice(val)
    } else 
    {
      toOutput <- "no project results detected at selected destination"
    }
    return(toOutput)
  })
  
  ProjFolder <- reactive({
    var <- as.character(input$ProjectName)
    var_corrected <- gsub(pattern = ' ', replacement = '_', x=var)
    var_corrected <- gsub(pattern = '.', replacement = '_', x=var_corrected, fixed=T)
    return(var_corrected)
  })
  
  
  detected_filenames_anomaly <- reactiveVal(0)
  fastqc_finished <- reactiveVal(0)
  bbduk_finished <- reactiveVal(0)
  salmon_finished <- reactiveVal(0)
  bbduk_and_fastqc_reused <- reactiveVal(0)
  started_listening <- reactiveVal(0)
  multiplegroups <- reactiveVal(0)
  
  res_auth <- secure_server(
    check_credentials = check_credentials(credentials),
    timeout = 0
  )
  
  output$auth_output <- renderPrint({
    reactiveValuesToList(res_auth)
  })
  
  output$Aftertrimming <- renderText({
    if(length(input$infiles) <= 1) return({})
    if(fastqc_finished()==0) 
      return({})
    return("Trimmed counts QC:")
  }
  )
  files_selected <- observe({
    shinyFileChoose(
      input, 
      "infiles", 
      roots = c(home = house),
      filetypes = c('', 'gz','fastq', 'fq')
    )
    req(input$infiles)
    if (is.null(input$infiles))
      return(NULL) 
    return(parseFilePaths(c(home = house), input$infiles)$datapath)
  })
  previous_project_selected <- observe({
    shinyDirChoose(
      input, 
      "infolder", 
      roots=c(home=house)
    )
    req(input$infolder)
    if (is.null(input$infolder))
      return(NULL)
    return(parseFilePaths(c(home = house), input$infolder)$datapath)
  })
  files <- reactive({input$infiles})  
  fls <- reactive({
    fl <- as.character(parseFilePaths(c(home = house),files())$datapath)
    print('line 393 ok, printing fls')
    print(fl)
    return(fl)
  })
  
  ProjFolderFull <- reactive({
    print('entering ProjFolderFull calculation')
    if(proceedWithLoad()==0){
      req(input$infiles)
      print('printing fls()')
      print(fls())
      print('line 388 ok, printing paste0(dirname(fls()[1]),ProjFolder())')
      toRet <- paste0(dirname(fls()[1]),'/',ProjFolder())
      print(toRet)
      return(toRet)
    } else {
      message('line 391!!!! ok')
      unl <- unlist(input$infolder)
      toPrint <- paste0(unl[-length(unl)], collapse = '/')
      print('printing toPrint')
      print(toPrint)
      toRet <- paste0(house, paste0(unl[-length(unl)], collapse = '/'))
      print('printing input$infolder')
      print(toRet)
      return(toRet)
    }
  })
  
  trimmed_folder <- reactive({
    req(input$infiles)
    paste0(ProjFolderFull(), '/trimmed/')
  })
  
  observe(
    {
      if(length(input$infiles) <= 1) return({})
      logfile <- file.path(ProjFolderFull(), 'dynamic_log.txt')
      print('printing logfile')
      print(logfile)
      system(command=paste0('mkdir -p ', ProjFolderFull(), '/trimmed'))
      system(paste0('touch ', logfile))
      system(paste0('echo "Progress.." > ', logfile))
      dynamic_log <- reactiveFileReader(1000, session, logfile, read.csv)
      output$dynamic_log <- renderTable(
        as.data.frame(tail(dynamic_log(), 20)),
        options=list(scrollY=TRUE)
      )
      started_listening(1)
    })
  
  observe({
    
  })
  
  indexPath <- reactive(
    {
      #req(input$referenceGenomeChoice)
      #Human GRCh38"=1,"Mouse GRCm39"=2,"Mouse C57B6"=3,"Mouse MM129"=4
      #print(input$referenceGenomeChoice)
      #if(input$referenceGenomeChoice==1) indexPath <- paste0(shiny_home,'/Homo_sapiens_GRCh38')
      #if(input$referenceGenomeChoice==2) indexPath <- paste0(shiny_home,'/GRCm39_mm_index')
      if(referenceGenomeChoice()==1) indexPath <- '/home/minicluster/BulkRNAsuite/GRCh38_index_SAF'
      if(referenceGenomeChoice()==2) indexPath <- '/home/minicluster/BulkRNAsuite/GRCm39_index_SAF'
      if(referenceGenomeChoice()==3) indexPath <- paste0(shiny_home,'/c57_b6_index')
      if(referenceGenomeChoice()==4) indexPath <- paste0(shiny_home,'/MM129_mm_index')
      return(indexPath)
    })
  txdbPath <- reactive(
    {req(input$referenceGenomeChoice)
      #Human GRCh38"=1,"Mouse GRCm39"=2,"Mouse C57B6"=3,"Mouse MM129"=4
      #print(input$referenceGenomeChoice)
      if(referenceGenomeChoice()==1) txdbPath <- paste0(house, '/txdb_hsapiens.RDS')
      if(referenceGenomeChoice()==2) txdbPath <- paste0(house, '/txdb_GRCm39.RDS')
      if(referenceGenomeChoice()==3) txdbPath <- paste0(house, '/txdb_c57bl6.RDS')
      if(referenceGenomeChoice()==4) txdbPath <- paste0(house, '/txdb_MM129.RDS')
      return(txdbPath)
    })
  
  output$selectedfiles_text <- renderText({
    req(input$infiles)
    if(length(input$infiles) <= 1) return({})
    message <-  "The following files will be analyzed"
    return(message)
  })
  
  output$selectedfiles <- renderText({
    req(input$infiles)
    if(length(input$infiles) <= 1) return({})
    return(paste(fls(), collapse = "\n"))
  })
  
  detected_multiplexedReads <- reactiveVal(0)
  detected_pairedEndReads <- reactiveVal(0)
  observe({
    req(input$infiles)
    if(length(input$infiles) <= 1) return({})
    if(length(grep(pattern = "_L002", x = fls()))>0){
      print("line 435 detected multiplexed reads")
      detected_multiplexedReads(1)
    }
    if(length(grep(pattern = "_R1", x = fls()))>0){
      print("line 436 detected pairedEndReads")
      detected_pairedEndReads(1)
    }
  })
  
  demultiplexed_fls <- reactive({
    req(input$infiles)
    message("line 445 ok")
    if(length(input$infiles) <= 1) return({})
    if(detected_multiplexedReads()==0) return(fls())
    else{
      reference_tracks <- grep(pattern = "L001", x = fls(), value = T)
      files2Concat <- list()
      resulting_tracks <- list()
      for (i in 1:length(reference_tracks)){
        reference_track <- reference_tracks[i]
        patt <- gsub(pattern = "L001", replacement = "L00*", reference_track)
        resulting_track <- gsub(pattern = "L001", replacement = "", reference_track)
        resulting_track <- gsub(pattern = "__", replacement = "_", resulting_track)
        resulting_tracks[[i]] <- resulting_track
      }
    }
    resulting_tracks <- unlist(resulting_tracks)
    resulting_tracks <- paste0(dirname(resulting_tracks[1]),'/',ProjFolder(),'/',basename(resulting_tracks))
    print("line 406 ok")
    print(resulting_tracks)
    return(resulting_tracks)
  })
  
  observe({
    message("line 469 ok")
    if(length(input$infiles) <= 1) return({})
    print("line 470 ok")
    toRet <- NULL
    
    resulting_tracks <- demultiplexed_fls()
    
    if(detected_multiplexedReads()==1){
      message <- "Detected multiplexed reads which will be concatenated first:\n\n"
      reference_tracks <- grep(pattern = "L001", x = fls(), value = T)
      files2Concat <- list()
      for (i in 1:length(reference_tracks)){
        reference_track <- reference_tracks[i]
        patt <- gsub(pattern = "L001", replacement = "L00*", reference_track)
        files2Concat[[i]] <- grep(pattern = glob2rx(patt), x = fls(), value = T)
      }
      
      print("line 418 ok")
      toRet <- paste0(message,
                      paste0(files2Concat[[1]], collapse = "\n"),
                      "\ninto\n",resulting_tracks[[1]], "\n===")
      if(length(reference_tracks)<2){
        detected_filenames_anomaly(1)
        #return("Cannot have only one sample!")
      } 
      print("line 427 ok")
      for (i in 2:length(reference_tracks)){
        toRet <- paste0(toRet, "\nas well as: \n",
                        paste0(files2Concat[[i]], collapse = "\n"),
                        "\ninto:\n",resulting_tracks[[i]], "\n===")
      }
      
      resulting_tracks_unlisted <- unlist(resulting_tracks)
      print("printing 489 line:")
      #print(fls)
      #return(fls)
    }
    if(detected_pairedEndReads()==1){
      
      print('line 499 ok')
      if(detected_multiplexedReads()==1) message <- "Paired-end reads in the data after concatenation of multiplexed lines which will be treated in pairs:\n" 
      if(detected_multiplexedReads()==0) message <- "Detected paired-end reads in the data which will be treated in pairs:\n" 
      
      R1_tracks <- grep(pattern = "_R1", x = demultiplexed_fls(), value = T)
      R1_tracks_swapped <- gsub(pattern = "_R1", replacement = "_R2", R1_tracks)
      R2_tracks <- grep(pattern = "_R2", x = demultiplexed_fls(), value = T)
      R2_tracks_swapped <- gsub(pattern = "_R2", replacement = "_R1", R2_tracks)
      
      R1_matched_tracks <- R1_tracks[which(R1_tracks_swapped %in% R2_tracks)]
      R1_tracks[which(!R1_tracks_swapped %in% R2_tracks)] -> R1_unbalanced_tracks
      R2_tracks[which(!R2_tracks_swapped %in% R1_tracks)] -> R2_unbalanced_tracks
      
      if(length(R1_unbalanced_tracks)>0|length(R2_unbalanced_tracks)>0) 
        detected_filenames_anomaly(1)
      
      R1_R2_matched_tracks_list <- list()
      for (i in 1:length(R1_matched_tracks)){
        R1_R2_matched_tracks_list[[i]] <- c(R1_matched_tracks[i], gsub(pattern = "_R1", replacement = "_R2", R1_matched_tracks[i]))
      }
      if(length(R1_unbalanced_tracks)>0){
        for (i in 1:length(R1_unbalanced_tracks)){
          R1_R2_matched_tracks_list <- append(R1_R2_matched_tracks_list, list(c(R1_unbalanced_tracks[i], paste0('[warning!]: file ',  gsub(pattern = "_R1", replacement = "_R2", R1_unbalanced_tracks[i]), ' was expected, but not provided, are you sure you selected all paired-end reads?'))))
        }
      }
      if(length(R2_unbalanced_tracks)>0){
        for (i in 1:length(R2_unbalanced_tracks)){
          R1_R2_matched_tracks_list <- append(R1_R2_matched_tracks_list, list(c(paste0('[warning!]: file ',  gsub(pattern = "_R2", replacement = "_R1", R2_unbalanced_tracks[i]), ' was expected, but not provided, are you sure you selected all paired-end reads?'), R2_unbalanced_tracks[i])))
        }
      }
      if((length(R1_tracks)<2) | (length(R2_tracks)<2)){
        detected_filenames_anomaly(1)
        print("R1_tracks:")
        print(R1_tracks)
        print("R2_tracks:")
        print(R2_tracks)
        toReturn <- "You have only one (complete paired-end) sample, need at least 4 samples (a duplicate from control and a duplicate from experiment)!"
      }
      if(length(R1_tracks)!=length(R2_tracks)){
        detected_filenames_anomaly(1)
        toReturn <- "unequal number of R1 and R2 tracks selected, are you sure you have selected both R1 and R2 reads for each sample?"
      }
      
      R1_R2_matched_tracks_list -> files2Pair
      
      toRet2 <- paste0(message,
                       paste0(files2Pair[[1]], collapse = "\n"), "\n===")
      for (i in 2:length(R1_R2_matched_tracks_list)){
        toRet2 <- paste0(toRet2, "\nas well as: \n",
                         paste0(files2Pair[[i]], collapse = "\n"),
                         "\n===")
      }
      if((length(R1_tracks)>=2) & (length(R2_tracks)>=2) & (length(R1_tracks)==length(R2_tracks))) #& (length(reference_tracks)>=2))
        toReturn <- paste0(toRet, "\n******************************\n\n", toRet2)
      print("wahoo")
      output$selectedfiles_multiplexedReads_and_salmon_pairedEndReads<- renderText({toReturn})
    }
  })
  
  output$filenameDiagnostics <- renderText({
    req(input$infiles)
    if(length(input$infiles) <= 1) return({})
    if((detected_multiplexedReads()==1)|(detected_pairedEndReads()==1)){
      return("Automatic detection of multiplexed and/or paired-end reads:")
    }
  })
  
  output$filenameDiagnosticsWarning <- renderText({
    req(input$infiles)
    if(length(input$infiles) <= 1) return({})
    if((detected_multiplexedReads()==1)|(detected_pairedEndReads()==1)){
      if(detected_filenames_anomaly()==1)
        return("Detected anomalies in the provided files, see details below!")
    }
  })
  
  
  observe({
    req(input$groups_specified)
    print("line 516 ok")
    #if(length(input$infiles) <= 1) return({})
    system(command=paste0('rm ', dirname(fls()[1]), "/*_Finished"))
    #system(command=paste0('mkdir -p ', ProjFolderFull(), '/trimmed'))
    if(detected_multiplexedReads()==1){
      reference_tracks <- grep(pattern = "L001", x = fls(), value = T)
      for (i in 1:length(reference_tracks)){
        reference_track <- reference_tracks[i]
        resulting_track <- gsub(pattern = "L001", replacement = "", x = reference_track)
        resulting_track <- gsub(pattern = "__", replacement = "_", x = resulting_track)
        reference_track <- gsub(pattern = "L001", replacement = "L00*", x = reference_track)
        resulting_track_full <- paste0(ProjFolderFull(),"/", basename(resulting_track))
        cmmnd <- paste0("cat ", reference_track, " > ", resulting_track_full)
        system(command=cmmnd)
      }
    }
  })
  
  observe({
    print("line 527 was ok")
    if(length(input$infiles) <= 1) return({})
    print("line 607 was ok")
    if(fastqc_finished()==0) return({})
    print("line 608 was ok")
    if(detected_multiplexedReads()==1){
      files2bbduk <- demultiplexed_fls()
    } else {
      files2bbduk <- fls()
    }
    print("kool!:\n")
    print(files2bbduk)
    
    
    if(detected_pairedEndReads()==1){
      cmmnd <- paste0("; filename_dir=$(dirname $filename1); filename_small1=$(basename $filename1); filename_small2=$(basename $filename2); extension=${filename_small1##*.}; filename_small_sans_ext1=${filename_small1%.*}; filename_small_sans_ext2=${filename_small2%.*}; ", bbduk.sh, " in1=${filename1} in2=${filename2} out1=", tempStorage, "/${filename_small_sans_ext1}.trimmed.fastq.gz out2=", tempStorage, "/${filename_small_sans_ext2}.trimmed.fastq.gz ref=/home/minicluster/miniconda3/opt/bbmap-39.06-0/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo ziplevel=3 -Xmx1g")
      ref<- grep(pattern = "_R1", x = files2bbduk, value = T)
      match_ref <- gsub("_R1", replacement = "_R2", x = ref)
      
      if(length(ref)>6){
        
        for (i in 1:6){
          flname <- c(ref[i], match_ref[i])
          print("printing the exact command that is passed to bbduk")
          arggg <- paste0("(filename1=", flname[1], "; filename2=",flname[2],cmmnd,"; touch ${filename1}_Finished)")
          print(arggg)
          system(command=arggg, wait = F)
          Sys.sleep(2)
        }
        Sys.sleep(240)
        
        if(length(ref)>12){
          for (i in 7:12){
            flname <- c(ref[i], match_ref[i])
            print("printing the exact command that is passed to bbduk")
            arggg <- paste0("(filename1=", flname[1], "; filename2=",flname[2],cmmnd,"; touch ${filename1}_Finished)")
            print(arggg)
            system(command=arggg, wait = F)
            Sys.sleep(2)
          }
          Sys.sleep(240)
          for (i in 13:length(ref)){
            flname <- c(ref[i], match_ref[i])
            print("printing the exact command that is passed to bbduk")
            arggg <- paste0("(filename1=", flname[1], "; filename2=",flname[2],cmmnd,"; touch ${filename1}_Finished)")
            print(arggg)
            system(command=arggg, wait = F)
            Sys.sleep(2)
          }
        } else {
          for (i in 7:length(ref)){
            flname <- c(ref[i], match_ref[i])
            print("printing the exact command that is passed to bbduk")
            arggg <- paste0("(filename1=", flname[1], "; filename2=",flname[2],cmmnd,"; touch ${filename1}_Finished)")
            print(arggg)
            system(command=arggg, wait = F)
            Sys.sleep(2)
          }
        }
        
      } else {
        for (i in 1:length(ref)){
          flname <- c(ref[i], match_ref[i])
          print("printing the exact command that is passed to bbduk")
          arggg <- paste0("(filename1=", flname[1], "; filename2=",flname[2],cmmnd,"; touch ${filename1}_Finished)")
          print(arggg)
          system(command=arggg, wait = F)
          Sys.sleep(2)
        }
      }
      
      
      filenames2monitor <- paste0(ref, "_Finished")
      print("printing filenames2monitor")
      print(filenames2monitor)
      while (!all(file.exists(filenames2monitor))) {
        Sys.sleep(1)
      }
      Sys.sleep(20)
      print("printing the exact command that is passed to fastqc")
      
      arggg <- paste0(fastqc_powershell, paste(unlist(final_fls()), collapse = " "), ' -o ',  trimmed_folder(),' 2>> ', ProjFolderFull(), '/dynamic_log.txt')
      print(arggg)
      system(arggg)
    } else {
      cmmnd <- paste0("; filename_dir=$(dirname $filename); filename_small=$(basename $filename); extension=${filename_small##*.}; filename_small_sans_ext=${filename_small%.*}; ", bbduk.sh," in=${filename} out=", tempStorage, "/${filename_small_sans_ext}.trimmed.fastq.gz ref=/home/minicluster/miniconda3/opt/bbmap-39.06-0/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo ziplevel=3 -Xmx1g")
      for (i in 1:length(files2bbduk)){
        flname <- files2bbduk[[i]]
        print("printing the exact command that is argggg")
        argggg <- paste0("(filename=", flname, cmmnd,"; touch ${filename}_Finished)")
        system(command=argggg, wait=F)
      }
      filenames2monitor <- paste0(unlist(files2bbduk), "_Finished")
      print("printing filenames2monitor")
      print(filenames2monitor)
      while (!all(file.exists(filenames2monitor))) {
        Sys.sleep(1)
      }
      Sys.sleep(10)
      print("printing the exact command that is argggg")
      argggg <- paste0(fastqc_powershell, paste(unlist(final_fls()), collapse = " "), ' -o ',  trimmed_folder(),' 2>> ', ProjFolderFull(), '/dynamic_log.txt')
      print(argggg)
      system(argggg)
    }
    bbduk_finished(1)
    print("bbduk_finished?")
    print(bbduk_finished())
  })
  
  final_fls <- reactive({
    print("line 559 ok")
    req(input$infiles)
    if(length(input$infiles) <= 1) return({})
    toReturn <- list()
    
    print('printing ProjFolderFull()')
    print(ProjFolderFull())
    print('printing trimmed_folder()')
    print(trimmed_folder())
    print('printing demulti')
    print(demultiplexed_fls())
    
    inputfiles_names <- paste0(trimmed_folder(),tools::file_path_sans_ext(basename(demultiplexed_fls())),'.trimmed.fastq.gz' )
    inputfiles <- paste0(tempStorage,'/',tools::file_path_sans_ext(basename(demultiplexed_fls())),'.trimmed.fastq.gz' )
    print("line and 559 ok")
    
    if(detected_pairedEndReads()==1){
      print("line after 559 ok")
      reference_initial_files <- grep(pattern = "_R1", x = inputfiles, value = T)
      reference_initial_files_names <- grep(pattern = "_R1", x = inputfiles_names, value = T)
      matched_reference_initial_files <- gsub(pattern = "_R1", replacement = "_R2", x = reference_initial_files)
      for (i in 1:length(reference_initial_files)){
        toReturn[[i]] <- c(reference_initial_files[i], matched_reference_initial_files[i])
      }
      names(toReturn) <- gsub("_R1", replacement = "", x = reference_initial_files_names)
    } else {
      print("else line after 559 ok")
      reference_initial_files <- inputfiles
      for (i in 1:length(reference_initial_files)){
        toReturn[[i]] <- inputfiles[i]
      }
      names(toReturn) <- inputfiles_names
    }
    print("line 630 ok, printing final_fls() below")
    print(toReturn)
    return(toReturn)
  })
  
  fastqc_initial_already_exist <- reactiveVal(0)
  trimmed_fastq_already_exist <- reactiveVal(0)
  fastqc_after_trimming_already_exist <- reactiveVal(0)
  salmon_alignments_already_exist <- reactiveVal(0)
  
  #observeEvent(input$infiles,{
  # TOFIX checking 4 salmon results
  #  if(length(input$infiles) <= 1) return({})
  
  #  flz <- final_fls()
  #  files2check <- c(paste0(tools::file_path_sans_ext(flz), "_fastqc.zip"), paste0(tools::file_path_sans_ext(flz), "_fastqc.html"))
  
  #  if(all(file.exists(files2check)))
  #    salmon_alignments_already_exist(1)
  #}
  #)
  
  output$button_save_colDatt<-renderUI({
    if(proceedWithLoad()==1) return({})
    req(input$infiles)
    if(length(input$infiles) <= 1) return({})
    actionButton(inputId = 'groups_specified', label = "Validate groups and start analysis", icon("paper-plane")) 
  })
  
  observeEvent( eventExpr = input$groups_specified, {
    #Sys.sleep(300)
    print('saving referenceGenomeChoice value to referenceGenomeChoice.RDS')
    print('the value is ')
    print(referenceGenomeChoice())
    saveRDS(referenceGenomeChoice(), paste0(ProjFolderFull(), '/referenceGenomeChoice.RDS'))
    print("line 650 ok, printing fls()")
    print(fls())
    print('printing demultiplexed_fls()')
    print(demultiplexed_fls())
    print('printing demultiplexed_fls()')
    print(demultiplexed_fls())
    print("printing final_fls()")
    print(final_fls())
    files2fastqc <- unlist(demultiplexed_fls())
    print('printing colDatt$data')
    print(colDatt$data)
    print('startingfastqc')
    startingfastqc_command <- paste0(fastqc_powershell, paste(files2fastqc, collapse = " "), ' -o ',  ProjFolderFull(),' 2>> ', ProjFolderFull(), '/dynamic_log.txt; touch ', ProjFolderFull(), '/1stStageFinished')
    print(startingfastqc_command)
    system(startingfastqc_command, wait = F)
    while(!file.exists(paste0(ProjFolderFull(), '/1stStageFinished'))){
      Sys.sleep(1)
      print('waiting for the 1st stage to Finish')
    }
    #cutext <- function(x){ y = gsub(pattern = '.fastq', replacement = '', x); return(gsub(pattern = '.gz', replacement = '', y))}
    #movingcommand <- function(x) paste0("file='", cutext(x), "'; dir2move2='",ProjFolderFull(),"';fl2mv=${file}_fastqc.zip; mv $fl2mv $dir2move2")
    #print("printing one moving command:")
    #i=1; print(movingcommand(files2fastqc[i]))
    #for (i in 1:length(files2fastqc)) system(command = movingcommand(files2fastqc[i]))
    fastqc_finished(1)
  })
  
  
  output$fastqc_finished_text <- renderText({
    if(fastqc_finished()==0) 
      return({})
    "Finished FastQC, here are the stats:"
  })
  
  fastqc_report_before_trimming <- reactive({
    if(proceedWithLoad()==0){
      fls <- demultiplexed_fls()
      fls <- unlist(fls)
      fls <- gsub('.fastq', replacement = '', x = fls)
      fls <- gsub('.gz', replacement = '', x = fls)
      fls <- gsub('.fq', replacement = '', x = fls)
      fastqc.zip_files <- paste0(ProjFolderFull(),'/', basename(fls),"_fastqc.zip")
      print('printing fastqc.zip_files')
      print(fastqc.zip_files)
      fastqc.zip_files_read <-list()
      for (i in 1:length(fastqc.zip_files)){
        fastqc.zip_files_read[[i]] <- suppressMessages(qc_read(file = fastqc.zip_files[[i]]))
      }
      saveRDS(fastqc.zip_files_read, file = paste0(ProjFolderFull(),'/fastqc_report_before_trimming.RDS'))
      return(fastqc.zip_files_read)
    } else {
      print("reading fastqc_report_before_trimming.RDS file")
      return(readRDS(paste0(ProjFolderFull(),'/fastqc_report_before_trimming.RDS')))
    }
  })
  
  output$fastqc_stats <- DT::renderDT({
    if(fastqc_finished()==0) 
      return({})
    basicstats <- data.frame(NULL)
    frbt <- fastqc_report_before_trimming()
    for (i in 1:length(frbt)){
      dedup <- round(100-frbt[[i]]$total_deduplicated_percentage,2)
      basicstats <- rbind(basicstats, c(as.data.frame(frbt[[i]]$basic_statistics)[c(1,4,5,6,7),2], dedup))
    }
    colnames(basicstats) <- c(as.data.frame(frbt[[1]]$basic_statistics)[c(1,4,5,6,7),1], "%Duplicated")
    
    return(basicstats)
  })
  output$fastqc_overrepseq <- DT::renderDT({
    if(fastqc_finished()==0) 
      return({})
    
    frbt <- fastqc_report_before_trimming()
    samplenames <- tools::file_path_sans_ext(extractNamesFromFastqcrReports(frbt))
    
    overseq<-list(); for (i in 1:length(frbt)){
      overseq[[i]] <- qc_plot(frbt[[i]], "Overrepresented sequences")
      if(is.null(dim(overseq[[i]]))) overseq[[i]] <- tibble(Sequence="none", Count=NULL,Percentage=NULL,"Possible Source"="none", Sample=samplenames[i]) else {
        overseq[[i]] <- overseq[[i]] %>% add_column(Sample=samplenames[i])
      } 
    }
    
    overseq <- overseq %>% bind_rows()
    
    d1 <- DT::datatable(overseq,
                        rownames = FALSE, 
                        extensions = 'RowGroup', 
                        options = list(rowGroup = list(dataSrc=c(2)),
                                       pageLength = 100,
                                       columnDefs = list(list(visible=FALSE, targets=c(2)))
                        ))
    return(d1)
  })
  
  output$fastqc_sequence_q <- renderPlot({
    if(fastqc_finished()==0) 
      return({})
    p_pbsq <- list(); 
    samplenames <- tools::file_path_sans_ext(extractNamesFromFastqcrReports(fastqc_report_before_trimming()))
    p_pbsq[[1]] <- qc_plot(fastqc_report_before_trimming()[[1]], 
                           modules = 'Per base sequence quality') + ggtitle(label="Per base sequence quality", 
                                                                            subtitle =samplenames[1] ); 
    for (i in 2:length(fastqc_report_before_trimming()))
      p_pbsq[[i]] <- qc_plot(fastqc_report_before_trimming()[[i]], modules = 'Per base sequence quality') + ggtitle(label="",subtitle = samplenames[i]); 
    return(Reduce(f = '+', p_pbsq))
  })
  
  fastqc_report_after_trimming <- reactive({
    if(salmon_finished()==0) return({})
    #if(proceedWithLoad()==0){
    print("line 895 was ok")
    fastqc.zip_files <- list.files(path = trimmed_folder(), pattern = "fastqc.zip", full.names = T)
    print("line 897 was ok")
    print(fastqc.zip_files)
    print("line 898 was ok")
    fastqc.zip_files_read <-list()
    for (i in 1:length(fastqc.zip_files)){
      fastqc.zip_files_read[[i]] <- suppressMessages(qc_read(file = fastqc.zip_files[[i]]))
    }
    saveRDS(fastqc.zip_files_read, file = paste0(ProjFolderFull(),'/fastqc_report_after_trimming.RDS'))
    return(fastqc.zip_files_read)
    #} else {
    #  print("reading fastqc_report_after_trimming.RDS file")
    #  return(readRDS(paste0(ProjFolderFull(),'/fastqc_report_after_trimming.RDS')))
    #}
  })
  
  observeEvent(input$validateforking,{
    
    if(length(input$folder2fork2)==0) output$forkingfeedback <- renderText("Can't have empty folder name")
    print('printing folder2fork2')
    print(input$folder2fork2)
    print('printing ForkingPrompt')
    print(input$samples2fork)
    
    foldr <- gsub(pattern = ' ', replacement = '_', x = input$folder2fork2)
    foldr <- gsub(pattern = '.', replacement = '_', x = foldr, fixed = T)
    foldr <- gsub(pattern = '!', replacement = '_', x = foldr, fixed = T)
    
    subst <- which(colDatt$data$Sample %in% input$samples2fork)
    colData_subset <- colDatt$data[subst,]
    txi_subset <- txi()
    txi_subset$abundance <- txi_subset$abundance[,subst]
    txi_subset$counts <- txi_subset$counts[,subst]
    txi_subset$length <- txi_subset$length[,subst]
    
    rownames(colData_subset) <- 1:nrow(colData_subset)
    
    dir.create(path = paste0(ProjFolderFull(),'/../',foldr))
    
    saveRDS(txi_subset, file = paste0(ProjFolderFull(),'/../',foldr,'/txi.RDS'))
    saveRDS(colData_subset, file = paste0(ProjFolderFull(),'/../',foldr,'/colData.RDS'))
    saveRDS(referenceGenomeChoice(), file = paste0(ProjFolderFull(),'/../',foldr,'/referenceGenomeChoice.RDS'))
    
    output$forkingfeedback <- renderText(paste0("~ Successfully copied files ", paste(input$samples2fork, collapse = ', '), ' into a folder called ', foldr))
  }
  
  )
  
  output$fastqc_stats_AfterTrimming <- DT::renderDT({
    if(fastqc_finished()==0) 
      return({})
    if(bbduk_finished()==0) 
      return({})
    basicstats <- data.frame(NULL)
    
    for (i in 1:length(fastqc_report_after_trimming())){
      dedup <- round(100-fastqc_report_after_trimming()[[i]]$total_deduplicated_percentage,2)
      basicstats <- rbind(basicstats, c(as.data.frame(fastqc_report_after_trimming()[[i]]$basic_statistics)[c(1,4,5,6,7),2], dedup))
    }
    colnames(basicstats) <- c(as.data.frame(fastqc_report_after_trimming()[[1]]$basic_statistics)[c(1,4,5,6,7),1], "%Duplicated")
    
    return(basicstats)
  })
  output$fastqc_overrepseq_AfterTrimming <- DT::renderDT({
    if(fastqc_finished()==0) 
      return({})
    if(bbduk_finished()==0) 
      return({})
    samplenames <- tools::file_path_sans_ext(extractNamesFromFastqcrReports(fastqc_report_after_trimming()))
    
    overseq<-list(); for (i in 1:length(fastqc_report_after_trimming())){
      overseq[[i]] <- qc_plot(fastqc_report_after_trimming()[[i]], "Overrepresented sequences")
      if(is.null(dim(overseq[[i]]))) overseq[[i]] <- tibble(Sequence="none", Count=NULL,Percentage=NULL,"Possible Source"="none", Sample=samplenames[i]) else {
        overseq[[i]] <- overseq[[i]] %>% add_column(Sample=samplenames[i])
      } 
    }
    
    overseq <- overseq %>% bind_rows()
    
    d1 <- DT::datatable(overseq,
                        rownames = FALSE, 
                        extensions = 'RowGroup', 
                        options = list(rowGroup = list(dataSrc=c(2)),
                                       pageLength = 100,
                                       columnDefs = list(list(visible=FALSE, targets=c(2)))
                        )
    )
    return(d1)
  })
  output$fastqc_sequence_q_AfterTrimming <- renderPlot({
    if(fastqc_finished()==0)
      return({})
    if(bbduk_finished()==0) 
      return({})
    p_pbsq <- list(); 
    samplenames <- tools::file_path_sans_ext(extractNamesFromFastqcrReports(fastqc_report_after_trimming()))
    p_pbsq[[1]] <- qc_plot(fastqc_report_after_trimming()[[1]], 
                           modules = 'Per base sequence quality') + ggtitle(label="Per base sequence quality", 
                                                                            subtitle =samplenames[1] ); 
    for (i in 2:length(fastqc_report_after_trimming()))
      p_pbsq[[i]] <- qc_plot(fastqc_report_after_trimming()[[i]], modules = 'Per base sequence quality') + ggtitle(label="",subtitle = samplenames[i]); 
    return(Reduce(f = '+', p_pbsq))
  })
  
  
  observe({
    if(salmon_finished()==0){
      shinyjs::hide(id='GSVA_title')
      shinyjs::hide(id='genes_for_gsva')
      shinyjs::hide(id='runGSVA')
      shinyjs::hide(id='enrichR_title')
      shinyjs::hide(id='genes_for_enrichR')
      shinyjs::hide(id='runenrichR')
      shinyjs::hide(id='ForkingPrompt')
      shinyjs::hide(id='GroupsPromptTextFork')
      shinyjs::hide(id='validateforking')
      shinyjs::hide(id='folder2fork2')
      shinyjs::hide(id='forkingfeedback')
    } else {
      shinyjs::show(id='GSVA_title')
      shinyjs::show(id='genes_for_gsva')
      shinyjs::show(id='runGSVA')
      shinyjs::show(id='enrichR_title')
      shinyjs::show(id='genes_for_enrichR')
      shinyjs::show(id='runenrichR')
      shinyjs::show(id = "fastqcstatspanel")
      shinyjs::show(id = "fastqcstatspanel_AfterTrimming")
      shinyjs::show(id='ForkingPrompt')
      shinyjs::show(id='GroupsPromptTextFork')
      shinyjs::show(id='validateforking')
      shinyjs::show(id='folder2fork2')
      shinyjs::show(id='forkingfeedback')
    }
  })
  
  observe({
    
    if(multiplegroups()==0){
      shinyjs::hide(id = "DESeq_DEGsMultitab")
      shinyjs::hide(id = "VolcanoMultitab")
      shinyjs::hide(id = "GO_multitab")
      shinyjs::hide(id = "GO_dotplot_multitab")
      shinyjs::hide(id = "DorotheaMultitab")
      
      shinyjs::show(id = "DESeq_DEGs")
      shinyjs::show(id = "Volcano")
      shinyjs::show(id = "GO")
      shinyjs::show(id = "GO_dotplot")
      shinyjs::show(id = "Dorothea")
    } else {
      
      shinyjs::show(id = "DESeq_DEGsMultitab")
      shinyjs::show(id = "VolcanoMultitab")
      shinyjs::show(id = "GO_multitab")
      shinyjs::show(id = "GO_dotplot_multitab")
      shinyjs::show(id = "DorotheaMultitab")
      
      shinyjs::hide(id = "DESeq_DEGs")
      shinyjs::hide(id = "Volcano")
      shinyjs::hide(id = "GO")
      shinyjs::hide(id = "GO_dotplot")
      shinyjs::hide(id = "Dorothea")
      
      
    }
  })
  
  observeEvent(input$ProjectName, {
    shinyjs::hide(id = "fastqcstatspanel")
    shinyjs::hide(id = "fastqcstatspanel_AfterTrimming")
    shinyjs::hide(id = "mySearch")
    shinyjs::hide(id = "downloadData")
    shinyjs::hide(id = "downloadDotPlot")
    shinyjs::hide(id = "downloadDotPlots")
  })
  observeEvent(input$groups_specified, {
    print("line 909 ok")
    shinyjs::hide(id = "button_save_colDatt")
    shinyjs::show(id = "fastqcstatspanel")
    shinyjs::show(id = "mySearch")
  })
  
  output$button_rebase_colDatt<-renderUI({
    while(salmon_finished()==0) return({})
    actionButton(inputId = 'groups_rebase', label = "Modify the background factor and/or sample groupings?", icon("paper-plane"))
  })
  
  output$button_fork_project <-renderUI({
    while(salmon_finished()==0) return({})
    actionButton(inputId = 'fork', label = "Fork some samples into a different folder?")
  })
  
  
  
  observe({
    req(input$groups_specified)
    print("line 913 ok")
    if(bbduk_finished()==1){ 
      print("bbduk's finished yepcock")
      shinyjs::show(id = "fastqcstatspanel_AfterTrimming")
    }
  })
  
  
  observeEvent(input$ProjectName, {
    shinyjs::hide(id = "Enrichments_panel")
    shinyjs::hide(id = "tabset1")
    shinyjs::hide(id = "tabset2")
    shinyjs::hide(id = "multitab_DEG")
    shinyjs::hide(id = "multitab_Volcano")
    shinyjs::hide(id = "multitab_Enrichments")
    shinyjs::hide(id = "multitab_Dorothea")
  })
  observeEvent(c(input$groups_specified,input$infolder), {
    shinyjs::show(id = "Enrichments_panel")
    shinyjs::show(id = "tabset1")
    shinyjs::show(id = "tabset2")
  })
  
  output$mapping_rates_header<- renderText({
    if(salmon_finished()==0) return({})
    #if(proceedWithLoad()==1) return({})
    return("Finished alignment, here are the mapping rates:")
  })
  
  output$mapping_rates_subheader<- renderText({
    if(salmon_finished()==0) return({})
    #if(proceedWithLoad()==1) return({})
    return("A respectable mapping rate would rarely go below 80%")
  })
  
  observe({
    if(bbduk_finished()==0) return({})
    initial_files <- final_fls()
    if(detected_pairedEndReads()==1){
      ##
      for (i in 1:length(initial_files)){
        flname1 <- initial_files[[i]][1]
        flname2 <- initial_files[[i]][2]
        flname <- gsub(pattern = "_R1", replacement = "", x = flname1)
        flname_base <- basename(flname)
        arggggg <- paste0("filename=",flname,"; ", salmon, " quant -i ", indexPath()," -l A -1 ", flname1," -2 ", flname2," -p 16 --validateMappings -o ", trimmed_folder(),'/', flname_base, "_quant")
        print('printing the arggggg passed to salmon')
        print(arggggg)
        system(command=arggggg)
      }
    } else {
      for (i in 1:length(initial_files)){
        flname <- initial_files[[i]]
        flname_base <- basename(flname)
        arggggg <- paste0("filename=",flname,"; ", salmon, " quant -i ", indexPath()," -l A -r ", flname," -p 16 --validateMappings -o ", trimmed_folder(),'/', flname_base, "_quant")
        print('printin the arggggg passed to salmon')
        print(arggggg)
        system(command=arggggg)
      }
    }
    print("salmon is officially finished!!!!!!")
    salmon_finished(1)
    ### removing temp files
    unlink(tempStorage, recursive=TRUE)
    
  })
  
  txi <- reactive({
    if(proceedWithLoad()==0){
      print("line 888 ok")
      if(salmon_finished()==0) return ({})
      filelist <- names(final_fls())
      TPM_filelist <- NULL;
      for (i in 1:length(filelist))
        TPM_filelist <- c(TPM_filelist, paste0(filelist[i],"_quant/quant.sf"))
      print("line 899 ok")
      print("printing TPM_filelist")
      names(TPM_filelist) <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(filelist)))
      print(TPM_filelist)
      txdb <- readRDS(txdbPath())
      txi_ <- tximport::tximport(TPM_filelist, type = "salmon", tx2gene = txdb, ignoreTxVersion = T)
      saveRDS(txi_, file = paste0(ProjFolderFull(),'/txi.RDS'))
      return(txi_)
    } else {
      salmon_finished(1)
      print("reading txi.RDS file at location:")
      loc <- paste0(ProjFolderFull(),'/txi.RDS')
      print(loc)
      return(readRDS(loc))
    }
    
  })
  
  txi_tpms <- reactive({
    if(salmon_finished()==0) return ({})
    print("line 987 ok")
    res<-txi()$abundance
    res <- data.frame(ENSEMBL=rownames(res), res)
    annots <- AnnotationDbi::select(OrgDeeBee(), keys=rownames(res), 
                                    columns="SYMBOL", keytype="ENSEMBL")
    result <- merge(annots, res, by.x="ENSEMBL", by.y="ENSEMBL")
    result <- rbind(c('-','-',colDatt$data$Group),result)
    
    openxlsx::write.xlsx(result, file = paste0(ProjFolderFull(), '/TPMs.xlsx'))
    saveRDS(result, file = paste0(ProjFolderFull(), '/txi_tpms.RDS'))
    
    res_counts<-txi()$counts
    res_counts <- data.frame(ENSEMBL=rownames(res_counts), res_counts)
    annots2 <- AnnotationDbi::select(OrgDeeBee(), keys=rownames(res_counts), 
                                     columns="SYMBOL", keytype="ENSEMBL")
    result_counts <- merge(annots2, res_counts, by.x="ENSEMBL", by.y="ENSEMBL")
    result_counts <- rbind(c('-','-',colDatt$data$Group),result_counts)
    openxlsx::write.xlsx(result_counts, file = paste0(ProjFolderFull(), '/Counts.xlsx'))
    
    return(result)
  })
  
  output$salmon_log <- renderRHandsontable({
    #if(proceedWithLoad()==0){
    if(salmon_finished()==0) return({})
    print("line 999 ok")
    toPrint <- NULL
    reactiveFileLog <- list()
    filelist <- names(final_fls())
    
    flpth <- paste0(filelist,"_quant/logs/salmon_quant.log")
    print('Printing flpth:')
    print(flpth)
    texts <- NULL
    rates <- NULL
    for (i in 1:length(flpth)){
      reactiveFileLog[[i]] <- reactiveFileReader(1000, session, filePath = flpth[i], read.table, sep = "\t")
      texts <- c(texts, grep(pattern = "Mapping rate", x = as.vector(reactiveFileLog[[i]]())[[1]], value = T))
      rates <- c(rates, sapply(strsplit(x = texts[i], split = "=", fixed = T), function(x) x[2]))
    }
    saveRDS(filelist, file = paste0(ProjFolderFull(),"/filelist.RDS"))
    saveRDS(rates, file = paste0(ProjFolderFull(),"/rates.RDS"))
    print("line 1061 was ok")
    #return(data.frame(AlignmentLog=grep(pattern = "Index contained|Mapping rate|warning", toPrint, value = T )))
    openxlsx::write.xlsx(x = data.frame(File=filelist, Rate=rates), file = 'alignment_rates.xlsx')
    rhandsontable(data.frame(File=filelist, Rate=rates), width = 2000) %>%
      hot_context_menu(
        customOpts = list(
          csv = list(name = "Download to CSV",
                     callback = htmlwidgets::JS(
                       "function (key, options) {
                         var csv = csvString(this, sep=',', dec='.');

                         var link = document.createElement('a');
                         link.setAttribute('href', 'data:text/plain;charset=utf-8,' +
                           encodeURIComponent(csv));
                         link.setAttribute('download', 'data.csv');

                         document.body.appendChild(link);
                         link.click();
                         document.body.removeChild(link);
                       }"))))
    
    #} else {
    #  filelist <- readRDS(file = paste0(ProjFolderFull(),"/filelist.RDS"))
    #  rates <- readRDS(file = paste0(ProjFolderFull(),"/rates.RDS"))
    #  openxlsx::write.xlsx(x = data.frame(File=filelist, Rate=rates), file = 'alignment_rates.xlsx')
    #  print("line 1061 was ok")
    #  #return(data.frame(AlignmentLog=grep(pattern = "Index contained|Mapping rate|warning", toPrint, value = T )))
    #  rhandsontable(data.frame(File=filelist, Rate=rates), width = 2000) %>%
    #    hot_context_menu(
    #      customOpts = list(
    #        csv = list(name = "Download to CSV",
    #                   callback = htmlwidgets::JS(
    #                     "function (key, options) {
    #                     var csv = csvString(this, sep=',', dec='.');
    #
    #                       var link = document.createElement('a');
    #                       link.setAttribute('href', 'data:text/plain;charset=utf-8,' +
    #                         encodeURIComponent(csv));
    #                       link.setAttribute('download', 'data.csv');
    #
    #                      document.body.appendChild(link);
    #                       link.click();
    #                       document.body.removeChild(link);
    #                    }"))))
    #}
  }
  )
  output$text_TPM<- renderText({
    if(salmon_finished()==0) return({})
    return("Transcripts-per-million (TPM) values:")
  })
  output$TPM <- DT::renderDT(datatable(txi_tpms(), extensions = 'Buttons', 
                                       options = list(
                                         dom = "Bl<'search'>rtip",
                                         buttons = c('copy', 'csv', 'excel', 'pdf', 'print')), callback = JS(callback)),
                             server=FALSE, rownames=FALSE)
  
  output$GroupsPromptText <- renderText({
    if(length(input$infiles) <= 1) return({}) else {
      "Please allocate samples to groups and select the control group:"
    }
  })
  output$GroupsPromptText2 <- renderText({
    if(length(input$infiles) <= 1) return({}) else {
      "the middle and right columns are editable"
    }
  })
  output$GroupsPromptText3 <- renderText({
    if(length(input$infiles) <= 1) return({}) else {
      "Note: levels of factors in the design should contain only letters, numbers, '_' and '.', as these are safe characters
  for column names in R"
    }
  })
  
  output$GroupsPromptTextFork <- renderText({
    if(salmon_finished()==0) return({}) else {
      "Do you want to take a subset of this run's samples and put them into another folder?"
    }
  })
  
  dat <- data.frame(Sample=c("sample1",'sample2'), Group= c("group?",'group?'))
  colDatt <- reactiveValues(data = dat)
  
  observeEvent(input$infiles, {
    if(length(input$infiles) <= 1) return({})
    filelist <- names(final_fls())
    df <- data.frame(Sample=tools::file_path_sans_ext(tools::file_path_sans_ext(basename(filelist))), Group=rep("Group?", length(filelist)), Control = rep(F, length(filelist)))
    colDatt$data <- df
  })
  
  #observeEvent(input$button_fastqc, colDatt$data <- hot_to_r(colDatt$data))
  
  output$GroupsPrompt <- renderRHandsontable({
    if(proceedWithLoad()==0){
      if(length(input$infiles) <= 1) return({})
      colDatt$data -> datt
      print('printing and writing datt')
      if(groups_specified()==1){
        write.table(x = as.data.frame(datt), file = paste0(ProjFolderFull(),"/colData.tsv"), sep = "\t")
        saveRDS(datt, file = paste0(ProjFolderFull(),"/colData.RDS"))
      }
      if(length(unique(colDatt$data$Group)) > 2) multiplegroups(1)
      print(datt)
      rhandsontable(datt) %>%
        hot_col("Sample", readOnly = TRUE)
    } else {
      datt <- readRDS(file = paste0(ProjFolderFull(),"/colData.RDS"))
      rownames(datt) <- NULL
      print('loading datt for rHandsonTable')
      print(datt)
      colDatt$data <- datt
      if(length(unique(datt$Group)) > 2) multiplegroups(1)
      print('printing the current multiplegroups value line 1386')
      salmon_finished(1)
      shinyjs::show(id = "mySearch")
      print(multiplegroups())
      rhandsontable(datt)
    }
  })
  
  output$ForkingPrompt <- renderUI({
    if(salmon_finished()==0) return({})
    selectInput(inputId = "samples2fork", label = "Samples to fork into another directory", choices = colDatt$data$Sample, multiple = T)
  }
  )
  
  groups_specified <- reactiveVal(0)
  observeEvent(input$groups_specified,{
    groups_specified(1)
  }
  )
  observeEvent(input$button_save_colDatt,{
    groups_specified(1)
  }
  )
  
  detected_salmon_results <- reactiveVal(0)
  
  #observeEvent(c(input$start_bbduk, input$reuse_bbduk),{
  #  fl <- names(final_fls())
  #print(fl)
  #  print("line 1146 was ok")
  #  files2check <- c(paste0(fl, "_quant/quant.sf"))
  #print(files2check)
  #  print("line 1147 was ok")
  #  if(all(file.exists(files2check)))
  #    detected_salmon_results(1)
  #})
  
  observeEvent(input$groups_specified, colDatt$data <- hot_to_r(input$GroupsPrompt))
  observeEvent(input$groups_rebase, {colDatt$data <- hot_to_r(input$GroupsPrompt); saveRDS(colDatt$data, file = paste0(ProjFolderFull(),"/colData.RDS"))})
  
  txi_deseq <- reactive({
    #if(proceedWithLoad()==0){
    if(salmon_finished()==0) return({})
    print('diagnostic colDatt$data from txi_deseq')
    print(colDatt$data)
    toRet <- DESeq2::DESeqDataSetFromTximport(txi(), colData = colDatt$data, design = ~Group)
    saveRDS(toRet, file = paste0(ProjFolderFull(),'/txi_deseq.RDS'))
    return(toRet)
    #} else {
    #  print("reading txi_deseq.RDS file")
    #  return(readRDS(paste0(ProjFolderFull(),'/txi_deseq.RDS')))
    #}
    
  })
  
  txi_deseq_deseq <- reactive({
    #if(proceedWithLoad()==0){
    if(salmon_finished()==0) return({})
    txi_deseq() -> matr
    colDatt$data -> colDat
    reff <- colDat[which(colDat[,3]),2][1]
    #reff <- matr$Group[which(colDat[,3])]
    #print(reff)
    if(length(which(colDat[,3]))>0) matr$Group <- relevel(x = matr$Group, ref = reff)
    smallestGroupSize <- floor(nrow(colDat))
    keep <- rowSums(counts(matr) >= 10) >= smallestGroupSize
    matr <- matr[keep,]
    toRet <- DESeq2::DESeq(matr)
    saveRDS(toRet, file = paste0(ProjFolderFull(),'/txi_deseq_deseq.RDS'))
    return(toRet)
    #}
    #else {
    #  print("reading txi_deseq_deseq.RDS file")
    #  return(readRDS(paste0(ProjFolderFull(),'/txi_deseq_deseq.RDS')))
    #}
  })
  res_txi_deseq <- reactive({
    #if(proceedWithLoad()==0){
    if(salmon_finished()==0) return({})
    if(length(resultsNames(txi_deseq_deseq()))<=2){
      toRet <- results(txi_deseq_deseq())
      toRet$log10padj <- log10(toRet$padj)
      toRet$Significance <- 'Non-significant'
      toRet$Significance[toRet$padj<0.05] <- 'Significant (padj < 0.05)'
      
      toRet2 <- data.frame(ENSEMBL=rownames(toRet), toRet[,1:2], abs_log2FoldChange=abs(toRet[,2]),toRet[,3:ncol(toRet)])
      annots <- AnnotationDbi::select(OrgDeeBee(), keys=rownames(toRet2), 
                                      columns="SYMBOL", keytype="ENSEMBL")
      
      toRet3 <- merge(annots, toRet2, by.x="ENSEMBL", by.y="ENSEMBL")
      toRet3 <- toRet3[order(toRet3$padj),]
      
    }  else {
      toRet3 <- toRet <- list()
      for (i in 2:length(resultsNames(txi_deseq_deseq()))){
        toRet[[(i-1)]] <- results(txi_deseq_deseq(), name = resultsNames(txi_deseq_deseq())[i] )
        toRet[[(i-1)]]$log10padj <- log10( toRet[[(i-1)]]$padj)
        toRet[[(i-1)]]$Significance <- 'Non-significant'
        toRet[[(i-1)]]$Significance[ toRet[[(i-1)]]$padj<0.05] <- 'Significant (padj < 0.05)'
        
        toRet2 <- data.frame(ENSEMBL=rownames(toRet[[(i-1)]]), toRet[[(i-1)]][,1:2], abs_log2FoldChange=abs(toRet[[(i-1)]][,2]),toRet[[(i-1)]][,3:ncol(toRet[[(i-1)]])])
        annots <- AnnotationDbi::select(OrgDeeBee(), keys=rownames(toRet2), 
                                        columns="SYMBOL", keytype="ENSEMBL")
        
        toRet3[[(i-1)]] <- merge(annots, toRet2, by.x="ENSEMBL", by.y="ENSEMBL")
        toRet3[[(i-1)]] <- toRet3[[(i-1)]][order(toRet3[[(i-1)]]$padj),]
        
      }
      names(toRet3) <- resultsNames(txi_deseq_deseq())[2:length(resultsNames(txi_deseq_deseq()))]
    }
    
    saveRDS(toRet3, file = paste0(ProjFolderFull(),'/res_txi_deseq.RDS'))
    openxlsx::write.xlsx(toRet3, file = paste0(ProjFolderFull(),'/DEGs_full.xlsx'))
    return(toRet3)
    #}  else {
    #  print("reading txi_deseq.RDS file")
    #  return(readRDS(paste0(ProjFolderFull(),'/res_txi_deseq.RDS')))
    #}
  })
  res_DEGs_txi_deseq <- reactive({
    #if(proceedWithLoad()==0){
    if(salmon_finished()==0) return({})
    toRet <- res_txi_deseq()
    print('printing class of res_txi_deseq()')
    print(class(toRet))
    if(class(toRet)=="data.frame"){
      toRet <- toRet[which(toRet$padj<0.05),]
      if(nrow(toRet)==0) return({as.data.frame("oops, none significant")})
      toRet <- toRet[order(toRet$padj),]
      openxlsx::write.xlsx(as.data.frame(toRet), file = paste0(ProjFolderFull(),'/DEGs.xlsx'))
    } else {
      for (i in 1:length(toRet)){
        toRet[[i]] <- toRet[[i]][which(toRet[[i]]$padj<0.05),]
        if(nrow(toRet[[i]])==0){
          toRet[[i]] <- "oops, none significant" 
        }
        else {
          toRet[[i]] <- toRet[[i]][order(toRet[[i]]$padj),]
        } 
      }
      openxlsx::write.xlsx(toRet, file = paste0(ProjFolderFull(),'/DEGs.xlsx'))
    }
    
    saveRDS(toRet, file = paste0(ProjFolderFull(),'/res_DEGs_txi_deseq.RDS'))
    return(toRet)
    #}
    #else {
    #  print("reading res_DEGs_txi_deseq.RDS file")
    #  return(readRDS(paste0(ProjFolderFull(),'/res_DEGs_txi_deseq.RDS')))
    #}
  })
  
  
  output$text_DEGs <- renderText({
    if(salmon_finished()==0) return({})
    "List of differentially expressed genes:"
  })
  output$text_which_is_Control_in_DEG <- renderText({
    if(salmon_finished()==0) return({})
    paste0("The control (background) is ", levels(txi_deseq_deseq()$Group)[1])
  })
  
  n_contrasts <- reactive({
    if(salmon_finished()==0) return({})
    return(length(resultsNames(res_DEGs_txi_deseq())))
  })
  
  
  
  output$DESeq_DEGs <- DT::renderDT({
    #req(input$groups_specified)
    if(salmon_finished()==0) return({})
    if(multiplegroups()==1) return({})
    
    datatable(as.data.frame(res_txi_deseq()), filter = 'top', extensions = 'Buttons', 
              options = list(
                dom = "Bl<'search'>rtip",
                buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                searchCols = list(NULL, NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL, NULL, NULL, list(search = 'Significant (padj < 0.05)'))
              ),
              callback = JS(callback))
  },
  server=FALSE, rownames = FALSE
  
  )
  
  output$DESeq_DEGsMultitab <- renderUI({
    #req(input$groups_specified)
    if(salmon_finished()==0) return({})
    if(multiplegroups()==0) return({})
    nTabs = length(res_DEGs_txi_deseq())
    
    myTabs = lapply(1: nTabs, function(x){
      tabPanel(paste(strsplit(names(res_txi_deseq()[x]), split = "_")[[1]][2:4], collapse=' '), 
               renderDT({datatable(res_txi_deseq()[[x]], filter = 'top', extensions = 'Buttons', 
                                   options = list(dom = "Blrtip",
                                                  buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                                  searchCols = list(NULL, NULL, NULL, NULL,
                                                                    NULL, NULL, NULL, NULL, NULL, NULL, NULL, list(search = 'Significant (padj < 0.05)'))
                                   ), callback = JS(callback),
               )},
               server=FALSE, rownames = FALSE))});
    return(do.call(tabsetPanel, myTabs))
  })
  
  output$PCA_title<- renderText({
    if(salmon_finished()==0) return({})
    return("PCA plot:")
  })
  
  output$Enrichments_title<- renderText({
    if(salmon_finished()==0) return({})
    return("Enriched terms:")
  })
  
  #output$PCA <- renderPlot({
  #  if(salmon_finished()==0) return({})
  #  p <- plotPCA(DESeq2::vst(txi_deseq()), intgroup="Group", ntop = nrow(txi()$abundance))
  #  p <- p + ggrepel::geom_label_repel(label=colnames(txi_deseq()), max.overlaps = 30)+ theme_classic()
  #  png(paste0(ProjFolderFull(), '/PCA_plot_300dpi.png'), 
  #      width = 12, 
  #      height = 12,
  #      res = 300, units = 'in')
  #  print(p)
  #  dev.off()
  #  return(p)
  #})
  
  output$PCA_2tab <- renderUI({
    if(salmon_finished()==0) return({})
    p <- list()
    
    dsq <- DESeq2::vst(txi_deseq())
    
    p[[1]] <- plotPCA(dsq, intgroup="Group", ntop = nrow(txi()$abundance))
    p[[1]] <- p[[1]] + ggrepel::geom_label_repel(label=colnames(txi_deseq()), max.overlaps = 30)+ theme_classic()
    png(paste0(ProjFolderFull(), '/PCA_plot_allGenes_300dpi.png'), 
        width = 12, 
        height = 12,
        res = 300, units = 'in')
    print(p[[1]])
    dev.off()
    
    p[[2]] <- plotPCA(dsq, intgroup="Group")
    p[[2]] <- p[[2]] + ggrepel::geom_label_repel(label=colnames(txi_deseq()), max.overlaps = 30)+ theme_classic()
    png(paste0(ProjFolderFull(), '/PCA_plot_top500HVG_300dpi.png'), 
        width = 12, 
        height = 12,
        res = 300, units = 'in')
    print(p[[2]])
    dev.off()
    
    nTabs = 2; tytl <- c('Top 500 HVG','All genes')
    myTabs = lapply(1: nTabs, function(x){tabPanel(tytl[x], renderPlot(p[[x]], height = 1200))});
    
    return(do.call(tabsetPanel, myTabs))
  })
  
  
  output$Volcano_title<- renderText({
    if(salmon_finished()==0) return({})
    return("Volcano plot:")
  })
  output$Volcano <- renderPlotly({
    if(salmon_finished()==0) return({})
    print('printing the current multiplegroups value line 1583')
    print(multiplegroups())
    if(multiplegroups()==1) return({})
    res_txi_deseq() -> df
    print('printing head of df just to check')
    print(head(df))
    df[which(is.na(df$SYMBOL)),'SYMBOL'] <- df[which(is.na(df$SYMBOL)),'ENSEMBL']
    df$SignificanceLevel <- 'NS'
    df[which(df$padj < 0.05 & abs(df$log2FoldChange) > 0.5 ),"SignificanceLevel"] <- "Significant&FoldChange"
    df[which(df$padj > 0.05 & abs(df$log2FoldChange) > 0.5 ),"SignificanceLevel"] <- "FoldChange"
    df[which(df$padj < 0.05 & abs(df$log2FoldChange) < 0.5 ),"SignificanceLevel"] <- "Significant"
    
    
    tytl <-''
    
    
    df.df <- as.data.frame(df)
    df.df$padj[which(df.df$padj==0)] <- min(df.df$padj[which(df.df$padj!=0)])/10000
    tx <- highlight_key(df.df, ~SYMBOL)
    
    # initiate a plotly object
    base <- plot_ly(tx, height=1000) %>%
      add_trace(x = ~log2FoldChange, 
                y = ~-log10(padj), 
                text = ~SYMBOL, mode = 'markers', 
                color = ~SignificanceLevel,  
                colors = c("#32a852","black", "blue", "red"),
                hovertemplate = paste('<b>%{text}</b><br>', '-log10(FDR): %{y:.2f}<br>','log2FC: %{x:.2f}'), 
                showlegend = T) %>%
      add_trace(data = df.df %>% 
                  filter(SignificanceLevel=='Significant&FoldChange') %>% 
                  top_n(-20, wt=padj), 
                x = ~log2FoldChange, 
                y = ~-log10(padj), 
                text = ~SYMBOL, mode = 'text',  textposition = "topright",  
                showlegend = T, name = 'Annotations') %>% 
      layout(xaxis=list(showgrid=F), yaxis=list(showgrid=F), title=tytl)
    
    
    # create a time series of median house price
    hlght <- highlight(
      base, 
      on = "plotly_click", 
      selectize = TRUE, 
      dynamic = TRUE, 
      persistent = TRUE,
      opacityDim = 0.07
    )
    #setwd(ProjFolderFull())
    #orca(hlght, file = paste0('Volcano_plot.svg'), width = 900)
    
    return(hlght)
  })
  
  output$VolcanoMultitab<- renderUI({
    #req(input$groups_specified)
    
    if(salmon_finished()==0) return({})
    if(multiplegroups()==0) return({})
    p <- list()
    tytl <- list()
    for (i in 1:length(res_txi_deseq())){
      res_txi_deseq()[[i]] -> df
      df[which(is.na(df$SYMBOL)),'SYMBOL'] <- df[which(is.na(df$SYMBOL)),'ENSEMBL']
      df$SignificanceLevel <- 'NS'
      df[which(df$padj < 0.05 & abs(df$log2FoldChange) > 0.5 ),"SignificanceLevel"] <- "Significant&FoldChange"
      df[which(df$padj > 0.05 & abs(df$log2FoldChange) > 0.5 ),"SignificanceLevel"] <- "FoldChange"
      df[which(df$padj < 0.05 & abs(df$log2FoldChange) < 0.5 ),"SignificanceLevel"] <- "Significant"
      
      library(crosstalk)
      library(plotly)
      tytl[[i]] <- paste(strsplit(names(res_txi_deseq()[i]), split = "_")[[1]][2:4], collapse=' ')
      
      
      df.df <- as.data.frame(df)
      df.df$padj[which(df.df$padj==0)] <- min(df.df$padj[which(df.df$padj!=0)])/10000
      tx <- highlight_key(df.df, ~SYMBOL)
      
      # initiate a plotly object
      base <- plot_ly(tx, height=1000) %>%
        add_trace(x = ~log2FoldChange, 
                  y = ~-log10(padj), 
                  text = ~SYMBOL, mode = 'markers', 
                  color = ~SignificanceLevel,  
                  colors = c("#32a852","black", "blue", "red"),
                  hovertemplate = paste('<b>%{text}</b><br>', '-log10(FDR): %{y:.2f}<br>','log2FC: %{x:.2f}'), 
                  showlegend = T) %>%
        add_trace(data = df.df %>% 
                    filter(SignificanceLevel=='Significant&FoldChange') %>% 
                    top_n(-20, wt=padj), 
                  x = ~log2FoldChange, 
                  y = ~-log10(padj), 
                  text = ~SYMBOL, mode = 'text',  textposition = "topright",  
                  showlegend = T, name = 'Annotations') %>% 
        layout(xaxis=list(showgrid=F), yaxis=list(showgrid=F), title=tytl)
      
      
      # create a time series of median house price
      p[[i]] <- highlight(
        base, 
        on = "plotly_click", 
        selectize = TRUE, 
        dynamic = TRUE, 
        persistent = TRUE,
        opacityDim = 0.07
      )
      setwd(ProjFolderFull())
      #orca(p[[i]], file = paste0('/Volcano_plot_',gsub(pattern = ' ', replacement = '_', tytl[[i]]),'.svg'), width = 900)
    }
    
    
    
    nTabs = length(res_txi_deseq())
    myTabs = lapply(1: nTabs, function(x){tabPanel(tytl[[x]], renderPlotly(p[[x]]))});
    
    return(do.call(tabsetPanel, myTabs))
  })
  
  m_t2g_most_symbol <- reactive({
    if(referenceGenomeChoice()!=1) 
      return(readRDS(paste0(house, '/m_t2g_most_symbol.RDS'))) else {
        return(readRDS(paste0(house, '/m_t2g_most_symbol_human.RDS')))
      }
  }
  )
  
  net <- reactive({
    if(referenceGenomeChoice()!=1){
      return(readRDS(paste0(house, '/dorothea_net_mouse.RDS')))} else {
        return(readRDS(paste0(house, '/dorothea_net_human.RDS')))
      }
  })
  
  
  
  OrgDeeBee <- reactive({
    print('printing input$referenceGenomeChoice')
    print(input$referenceGenomeChoice)
    
    print('printing referenceGenomeChoice()')
    print(referenceGenomeChoice())
    
    if(referenceGenomeChoice()!=1){
      return(org.Mm.eg.db::org.Mm.eg.db)} else {
        return(org.Hs.eg.db::org.Hs.eg.db)
      }
  }
  )
  
  GSEA_result <- reactive({
    if(proceedWithLoad()==0)
    {
      if(salmon_finished()==0) return({})
      print("printing fcListFromDF(res_DEGs_txi_deseq())")
      w1 <- fcListFromDF(res_DEGs_txi_deseq())
      toKeep <- which(!is.na(names(w1)) & !duplicated(names(w1)))
      w2 <- w1[toKeep,]
      print(w2)
      toRet <- GSEA(geneList = w2, TERM2GENE = m_t2g_most_symbol())
      openxlsx::write.xlsx(as.data.frame(toRet), file = paste0(ProjFolderFull(),'/GSEAs.xlsx'))
      saveRDS(toRet, file = paste0(ProjFolderFull(),'/GSEA_result.RDS'))
      return(
        toRet
      )} else {
        print("reading GSEA_result.RDS file")
        return(readRDS(paste0(ProjFolderFull(),'/GSEA_result.RDS')))
      }
  })
  GO_result <- reactive({
    #req(input$groups_specified)
    if(salmon_finished()==0) return({})
    #if(proceedWithLoad()==0){
    toRet <- list()
    toWrite <- list()
    tytl <- list()
    if(class(res_DEGs_txi_deseq())=='list'){
      # multiple groups
      for (i in 1:length(res_DEGs_txi_deseq())){
        tytl[[i]] <- paste(strsplit(names(res_txi_deseq()[i]), split = "_")[[1]][2:4], collapse=' ')
        if(is.null(dim(res_DEGs_txi_deseq()[[i]]))){
          toRet[[i]] <- NULL
          toWrite[[i]] <- NULL
        } else {
          toRet[[i]] <- enrichGO(gene = res_DEGs_txi_deseq()[[i]][,2], 
                                 keyType = "SYMBOL", 
                                 OrgDb = OrgDeeBee(), 
                                 ont = "BP", 
                                 pAdjustMethod = "BH", 
                                 qvalueCutoff = 0.05, 
                                 readable = TRUE)
          
          toWrite[[i]] <- as.data.frame(toRet[[i]])
          toWrite[[i]]$geneID <- gsub(pattern='/', replacement=' ', toWrite[[i]]$geneID)
          
        }
        
      }
      names(toRet) <- unlist(tytl)
      saveRDS(toRet, file = paste0(ProjFolderFull(),'/GO_result.RDS'))
      names(toRet) <- gsub(pattern = 'Group ', replacement = '', x = names(toRet))
      names(toRet) <- substr(names(toRet), start = 1, stop = 30)
      openxlsx::write.xlsx(toWrite, file = paste0(ProjFolderFull(),'/GOs.xlsx'))
      return(toRet)
      
    } else {
      # 1 vs 1 
      toRet <- enrichGO(gene = res_DEGs_txi_deseq()[,2], 
                        keyType = "SYMBOL", 
                        OrgDb = OrgDeeBee(), 
                        ont = "BP", 
                        pAdjustMethod = "BH", 
                        qvalueCutoff = 0.05, 
                        readable = TRUE)
      saveRDS(toRet, file = paste0(ProjFolderFull(),'/GO_result.RDS'))
      toWrite <- as.data.frame(toRet)
      toWrite$geneID <- gsub(pattern='/', replacement=' ', toWrite$geneID)
      openxlsx::write.xlsx(toWrite, file = paste0(ProjFolderFull(),'/GOs.xlsx'))
      return(toRet)
    }
    #}  else {
    #  print("reading GO_result.RDS file")
    #  return(readRDS(paste0(ProjFolderFull(),'/GO_result.RDS')))
    #}
    
    
  })
  
  
  output$enrichR_title<- renderText({
    if(salmon_finished()==0) return({})
    return("Calculate deeper enrichment against the enrichR database")
  })
  
  enrichR_string <- reactive({
    req(input$runenrichR)
    if(length(input$genes_for_enrichR)==0) return(NULL)
    toRet <- strsplits(input$genes_for_enrichR, c(" ", ",", "/"))
    return(toRet)
  })
  
  enrichR_string_matched <- reactiveVal('')
  
  enrichR_result <- reactive({
    req(input$runenrichR)
    if(is.null(enrichR_string())) return({})
    
    abund <- txi_tpms()
    groups <- abund[1,3:ncol(abund)]
    abund <- abund[-1,]
    
    abund <- abund[which(!duplicated(abund$ENSEMBL)),]
    rownames(abund) <- abund$ENSEMBL
    #matching_genes <- unique(enrichR_string()[which(enrichR_string() %in% unlist(abund[,1:2]))])
    #non_matching_genes <- unique(enrichR_string()[which(!enrichR_string() %in% unlist(abund[,1:2]))])
    
    matching_genes <- unique(enrichR_string()[which(toupper(enrichR_string()) %in% toupper(unlist(abund[,1:2])))])
    non_matching_genes <- unique(enrichR_string()[which(!toupper(enrichR_string()) %in% toupper(unlist(abund[,1:2])))])
    
    if(length(matching_genes) < 2){
      return({})
    } else {
      enrichR_string_matched(matching_genes)
      
      #which(unlist(abund[,1:2]) %in% enrichR_string()) -> v1
      #which(toupper(unlist(abund[,1:2])) %in% toupper(enrichR_string())) -> v1
      #ifelse(test=v1>nrow(abund), yes=v1-nrow(abund), no=v1) -> v2
      
      # signature <- rownames(abund)[unique(v2)]
      #abund <- abund[,-2]
      #abund <- abund[,-1]
      #print(signature)
      #print(head(abund))
      signature <- matching_genes
      
      dbs <- listEnrichrDbs()
      dbs_selected <- dbs$libraryName[dbs$categoryId %in% c(2,3,4,5)]
      #saveRDS(signature, file = paste0(ProjFolderFull(),'/signature.RDS'))
      enriched <- enrichr(signature, dbs_selected)
      
      enriched_clean <- Filter(function(x) !(nrow(x)==0), enriched);
      enriched_clean_named <- Map(cbind, group = names(enriched_clean), enriched_clean)
      enrichR_res <- dplyr::bind_rows(enriched_clean_named)
      saveRDS(enrichR_res, file = paste0(ProjFolderFull(),'/enrichR_res.RDS'))
      
      enrichR_res <- enrichR_res %>% dplyr::filter(Adjusted.P.value < 0.05)
      enrichR_res[,1] <- as.factor(enrichR_res[,1])
      
      openxlsx::write.xlsx(x = enrichR_res, file = paste0(ProjFolderFull(),'/enrichR_table.xlsx'))
      saveRDS(enrichR_res, file = paste0(ProjFolderFull(),'/enrichR_result.RDS'))
      return(enrichR_res)
    }
  })
  
  
  output$enrichRtable<- DT::renderDT({
    req(input$runenrichR)
    if(is.null(enrichR_result())) return(NULL)
    as.data.frame(enrichR_result())
  }, extensions = 'Buttons', 
  options = list(
    dom = 'Bflrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
  server=FALSE, rownames=T, filter='top')
  
  output$enrichR_genes_matching<- renderText({
    if(salmon_finished()==0) return({})
    req(input$runenrichR)
    if(is.null(enrichR_string())) return("List of genes shouldn't be empty")
    
    
    if(length(enrichR_string_matched())<2){
      return(paste0('Only 1 gene was correctly mapped: ', enrichR_string_matched())) 
    } else {
      nnmpdgns <- setdiff(enrichR_string(),enrichR_string_matched())
      if(length(nnmpdgns)==0){
        return(paste0('All genes were mapped correctly! The genes were: ', paste(enrichR_string_matched(), collapse = ",")))
      } else {
        return(paste0('Correctly mapped genes: ', paste(enrichR_string_matched(), collapse = ","), '\n',
                      'Non-mapped genes: ', paste(nnmpdgns, collapse = ","), '\n',
                      "(If the list of unmapped genes is too long, try to replace them with their ENSEMBL notation)"
        ))
      }
      
    }
  })
  
  
  
  
  output$GSEA<- renderDT({
    if(salmon_finished()==0) return({})
    as.data.frame(GSEA_result())
  }, extensions = 'Buttons', 
  options = list(
    dom = 'Bfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
  server=FALSE, rownames=FALSE)
  
  output$GO<- renderDT({
    if(salmon_finished()==0) return({})
    if(multiplegroups()==1) return({})
    varbl <- as.data.frame(GO_result())
    varbl$geneID <- gsub(pattern='/', replacement=' ', varbl$geneID)
    varbl
  }, extensions = 'Buttons', 
  options = list(
    dom = 'Bfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
  server=FALSE, rownames=FALSE)
  
  output$GSEA_dotplot<- renderPlot({
    if(salmon_finished()==0) return({})
    p <- dotplot(GSEA_result(), showCategory=input$showCategory)
    png(paste0(ProjFolderFull(),'/GSEA_dotplot_300dpi.png'), 
        width = 10, 
        height = 20,
        res = 300, units = 'in')
    print(p)
    dev.off()
    return(p)
  })
  
  
  
  output$GO_multitab <- renderUI({
    #req(input$groups_specified)
    if(salmon_finished()==0) return({})
    if(multiplegroups()==0) return({})
    
    nTabs = length(GO_result())
    
    toOutput <- list()
    for (o in 1:nTabs){ 
      toOutput[[o]] <- as.data.frame(GO_result()[[o]])
      toOutput[[o]]$geneID <- gsub(pattern='/', replacement=' ', toOutput[[o]]$geneID)
    }
    
    myTabs = lapply(1: nTabs, function(x){tabPanel(paste(strsplit(names(res_txi_deseq()[x]), split = "_")[[1]][2:4], collapse=' ')
                                                   , renderDT(datatable(toOutput[[x]], filter = 'top', extensions = 'Buttons', 
                                                                        options = list(dom = "Blrtip",
                                                                                       buttons = c('copy', 'csv', 'excel', 'pdf', 'print'))),
                                                              server=FALSE, rownames = FALSE))});
    return(do.call(tabsetPanel, myTabs))
  })
  
  output$GO_dotplot<- renderPlot({
    if(salmon_finished()==0) return({})
    if(multiplegroups()==1) return({})
    ncat <- as.numeric(input$showCategory)
    
    p <- dotplot(GO_result(), showCategory=50, x='p.adjust', decreasing=F)
    png(paste0(ProjFolderFull(),'/GO_dotplot_300dpi.png'), 
        width = 10, 
        height = 20,
        res = 300, units = 'in')
    print(p)
    dev.off()
    
    if(ncat!=50){
      p_custom <- dotplot(GO_result(), showCategory=ncat, x='p.adjust', decreasing=F)
      png(paste0(ProjFolderFull(),'/GO_dotplot_custom_number_of_Terms_300dpi.png'), 
          width = 10, 
          height = round(20*ncat/50, 1),
          res = 300, units = 'in')
      print(p_custom)
      dev.off()
      click("downloadDotPlot")
      return(p_custom)
    }
    return(p)
  })
  
  output$GO_dotplot_multitab <- renderUI({
    #req(input$groups_specified)
    if(salmon_finished()==0) return({})
    if(multiplegroups()==0) return({})
    p <- list()
    tytl <- list()
    ncat <- as.numeric(input$showCategory)
    for (i in 1:length(res_txi_deseq())){
      tytl[[i]] <- paste(strsplit(names(res_txi_deseq()[i]), split = "_")[[1]][2:4], collapse=' ')
      if(is.null(GO_result()[[i]])){ p[[i]] <- ggplot()+theme_void()} else {
        if(ncat!=50){
          p[[i]] <- dotplot(GO_result()[[i]], showCategory=ncat, x='p.adjust', decreasing=F)
          flnm <- paste0(ProjFolderFull(),'/GO_dotplot_custom_number_of_Terms_300dpi_', tytl[[i]],'.png')
          cat('printing flnm')
          cat(flnm)
          png(filename = flnm, 
              width = 10, 
              height = round(20*ncat/50, 1),
              res = 300, units = 'in')
          print(p[[i]])
          dev.off()
        } else {
          p[[i]] <- dotplot(GO_result()[[i]], showCategory=50, x='p.adjust', decreasing=F)
          flnm <- paste0(ProjFolderFull(),'/GO_dotplot_300dpi_', tytl[[i]],'.png')
          cat('printing flnm')
          cat(flnm)
          png(filename = flnm, 
              width = 10, 
              height = 20,
              res = 300, units = 'in')
          print(p[[i]])
          dev.off()
        }
        
      }
    }
    
    if(ncat!=50){
      files2zip <- list.files(path = ProjFolderFull(), pattern = 'GO_dotplot_custom_number_of_Terms', full.names = F, recursive = F)
      print('line 2117 printing files2zip')
      print(files2zip)
      zip::zip(root = ProjFolderFull(), zipfile = 'Plots_custom.zip', 
               files=files2zip, recurse = F, include_directories = F)
      click('downloadDotPlots')
    }
    
    nTabs = length(res_txi_deseq())
    myTabs = lapply(1: nTabs, function(x){tabPanel(tytl[[x]], 
                                                   renderPlot(p[[x]], height = 2000))});
    
    return(do.call(tabsetPanel, myTabs))
  })
  
  
  
  output$GSVA_title<- renderText({
    if(salmon_finished()==0) return({})
    return("GSVA")
  })
  
  gsva_string <- reactive({
    req(input$runGSVA)
    if(length(input$genes_for_gsva)==0) return(NULL)
    toRet <- strsplits(input$genes_for_gsva, c(" ", ",", "/"))
    return(toRet)
  })
  
  gsva_string_matched <- reactiveVal('')
  
  gsva_result <- reactive({
    req(input$runGSVA)
    if(is.null(gsva_string())) return({})
    
    abund <- txi_tpms()
    groups <- abund[1,3:ncol(abund)]
    abund <- abund[-1,]
    
    abund <- abund[which(!duplicated(abund$ENSEMBL)),]
    rownames(abund) <- abund$ENSEMBL
    #matching_genes <- unique(gsva_string()[which(gsva_string() %in% unlist(abund[,1:2]))])
    #non_matching_genes <- unique(gsva_string()[which(!gsva_string() %in% unlist(abund[,1:2]))])
    
    matching_genes <- unique(gsva_string()[which(toupper(gsva_string()) %in% toupper(unlist(abund[,1:2])))])
    non_matching_genes <- unique(gsva_string()[which(!toupper(gsva_string()) %in% toupper(unlist(abund[,1:2])))])
    
    if(length(matching_genes) < 2){
      return({})
    } else {
      gsva_string_matched(matching_genes)
      
      #which(unlist(abund[,1:2]) %in% gsva_string()) -> v1
      which(toupper(unlist(abund[,1:2])) %in% toupper(gsva_string())) -> v1
      ifelse(test=v1>nrow(abund), yes=v1-nrow(abund), no=v1) -> v2
      
      signature <- rownames(abund)[unique(v2)]
      abund <- abund[,-2]
      abund <- abund[,-1]
      #print(signature)
      #print(head(abund))
      
      gsva <- GSVA::gsva(as.matrix(dplyr::mutate_all(abund, function(x) as.numeric(as.character(x)))),
                         gset.idx.list = list(signature))
      gsva <- rbind(groups, gsva)
      gsva <- t(gsva)
      colnames(gsva) <- c('Group','Value')
      gsva <- as.data.frame(gsva)
      gsva[,2] <- as.numeric(gsva[,2])
      
      openxlsx::write.xlsx(x = gsva, file = paste0(ProjFolderFull(),'/gsva_table.xlsx'))
      saveRDS(gsva, file = 'gsva_result.RDS')
      return(gsva)
    }
  })
  
  output$GSVAplot<- renderPlot({
    if(salmon_finished()==0) return({})
    req(input$runGSVA)
    if(is.null(gsva_result())) return({})
    
    fntsize=12
    p1 <- ggplot(gsva_result(), aes(x=Group, y=Value, fill=Group)) + 
      geom_boxplot(width=0.4)+
      #geom_point(position = position_jitterdodge(jitter.width = 0.0, dodge.width = 0.8))+
      ggtitle(paste0("GSVA analysis"), subtitle = "on TPM counts") +  
      xlab("") +
      ylab("GSVA value")+
      theme_classic() +
      theme(legend.position = "none") + theme(axis.text = element_text(size = fntsize))+scale_x_discrete(guide = guide_axis(n.dodge = 2))
    
    #flnm <- paste0("GSVAsignature.png",paste(gsva_string_matched(), collapse = ","), ".png")
    flnm <- paste0("GSVAsignature.png")
    #print('printing flnm')
    #print(flnm)
    png(filename = flnm, width = 8, height = 6, units = "in", res =300)
    print(p1)
    dev.off()
    
    
    
    
    return(p1)
  })
  
  output$GSVAtable<- DT::renderDT({
    req(input$runGSVA)
    if(is.null(gsva_result())) return(NULL)
    as.data.frame(gsva_result())
  }, extensions = 'Buttons', 
  options = list(
    dom = 'Bfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
  server=FALSE, rownames=T)
  
  output$GSVA_genes_matching<- renderText({
    if(salmon_finished()==0) return({})
    req(input$runGSVA)
    if(is.null(gsva_string())) return("List of genes shouldn't be empty")
    
    
    if(length(gsva_string_matched())<2){
      return(paste0('Only 1 gene was correctly mapped: ', gsva_string_matched())) 
    } else {
      nnmpdgns <- setdiff(gsva_string(),gsva_string_matched())
      if(length(nnmpdgns)==0){
        return(paste0('All genes were mapped correctly! The genes were: ', paste(gsva_string_matched(), collapse = ",")))
      } else {
        return(paste0('Correctly mapped genes: ', paste(gsva_string_matched(), collapse = ","), '\n',
                      'Non-mapped genes: ', paste(nnmpdgns, collapse = ","), '\n',
                      "(If the list of unmapped genes is too long, try to replace them with their ENSEMBL notation)"
        ))
      }
      
    }
  })
  
  output$Dorothea_title<- renderText({
    if(salmon_finished()==0) return({})
    return("Dorothea regulons:")
  })
  
  output$Dorothea <- renderPlot({
    #req(input$groups_specified)
    if(salmon_finished()==0) return({})
    if(multiplegroups()==1) return({})
    #if(proceedWithLoad()==0){
    deg <- res_DEGs_txi_deseq()
    #counts <- txi()$abundance
    
    #TOFIX add dorothea per-sample composition
    
    net <- net()
    
    
    #counts[which(rowSums(counts)!=0),] -> meaningful_counts
    #meaningful_counts[rownames(deg),] -> meaningful_counts
    deg -> meaningful_deg
    
    print('if rownames(meaningful_counts) are %in% gene2tx')
    #print(all(rownames(meaningful_counts) %in% gene2tx2name_GRCm39_u_systematic_u_onlyGeneNames[,1]))
    
    print('rownames(meaningful_counts):')  
    #print(head(rownames(meaningful_counts)))
    
    print('rownames(meaningful_deg):')  
    print(head(rownames(meaningful_deg)))
    
    #res1 <- data.frame(ENSEMBL=rownames(meaningful_counts), meaningful_counts)
    #annots1 <- AnnotationDbi::select(OrgDeeBee, keys=res1$ENSEMBL, 
    #                                 columns="SYMBOL", keytype="ENSEMBL")
    #meaningful_counts_symbol <- merge(annots1, res1, by.x="ENSEMBL", by.y="ENSEMBL")
    
    #rownames(meaningful_deg) <- meaningful_deg$ENSEMBL 
    #res2 <- data.frame(ENSEMBL=rownames(meaningful_deg), meaningful_deg)
    #annots2 <- AnnotationDbi::select(OrgDeeBee, keys=res2$ENSEMBL, 
    #                                 columns="SYMBOL", keytype="ENSEMBL")
    meaningful_deg_symbol <- meaningful_deg
    
    #meaningful_counts_symbol2 <- meaningful_counts_symbol[which(!is.na(meaningful_deg_symbol$stat)),]
    meaningful_deg_symbol2 <- meaningful_deg_symbol[which(!is.na(meaningful_deg_symbol$stat)),]
    meaningful_deg_symbol2 <- meaningful_deg_symbol2[!is.na(meaningful_deg_symbol2$SYMBOL),]
    
    #toKeep <- which(!is.na(meaningful_counts_symbol2$SYMBOL));
    #meaningful_counts_symbol2 <- meaningful_counts_symbol2[toKeep,];
    #toKeep2 <- which(!duplicated(meaningful_counts_symbol2$SYMBOL));
    #meaningful_counts_symbol2 <- meaningful_counts_symbol2[toKeep2,-1]
    #rownames(meaningful_counts_symbol2) <- meaningful_counts_symbol2$SYMBOL
    #meaningful_counts_symbol2 <- meaningful_counts_symbol2[,-1];
    
    #meaningful_deg_symbol2 <- meaningful_deg_symbol2[toKeep,]
    #meaningful_deg_symbol2 <- meaningful_deg_symbol2[toKeep2,]
    #rownames(meaningful_deg_symbol2) <- rownames(meaningful_counts_symbol2)
    
    #saveRDS(meaningful_counts_symbol2, file = 'meaningful_counts_symbol2.RDS')
    #saveRDS(meaningful_deg_symbol2, file = 'meaningful_deg_symbol2.RDS')
    
    #rownames(meaningful_counts_symbol2) <- toupper(rownames(meaningful_counts_symbol2))
    rownames(meaningful_deg_symbol2) <- make.unique(toupper(meaningful_deg_symbol2$SYMBOL))
    
    
    #sample_acts <- decoupleR::run_ulm(mat=meaningful_counts_symbol2, net=net, .source='source', .target='target',
    #                                  .mor='mor', minsize = 5)
    contrast_acts <- decoupleR::run_ulm(mat=meaningful_deg_symbol2[,'stat', drop=F], net=net, .source='source', .target='target',
                                        .mor='mor', minsize = 5)
    
    
    
    
    n_tfs <- 25
    f_contrast_acts <- contrast_acts %>%
      mutate(rnk = NA)
    msk <- f_contrast_acts$score > 0
    f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
    f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))
    tfs <- f_contrast_acts %>%
      arrange(rnk) %>%
      head(n_tfs) %>%
      pull(source)
    f_contrast_acts <- f_contrast_acts %>%
      filter(source %in% tfs)
    
    n_tfs <- 100
    f_contrast_acts_100 <- contrast_acts %>%
      mutate(rnk = NA)
    msk_100 <- f_contrast_acts_100$score > 0
    f_contrast_acts_100[msk_100, 'rnk'] <- rank(-f_contrast_acts_100[msk_100, 'score'])
    f_contrast_acts_100[!msk_100, 'rnk'] <- rank(-abs(f_contrast_acts_100[!msk_100, 'score']))
    tfs_100 <- f_contrast_acts_100 %>%
      arrange(rnk) %>%
      head(n_tfs) %>%
      pull(source)
    f_contrast_acts_100 <- f_contrast_acts_100 %>%
      filter(source %in% tfs)
    
    saveRDS(f_contrast_acts, file = paste0(ProjFolderFull(),'/dorothea_f_contrast_acts.RDS'))
    saveRDS(f_contrast_acts_100, file = paste0(ProjFolderFull(),'/dorothea_f_contrast_acts_100.RDS'))
    as.data.frame(f_contrast_acts_100)->f_contrast_acts_100_df
    #names(f_contrast_acts_100_df) <- NULL
    openxlsx::write.xlsx(x = f_contrast_acts_100_df, file = paste0(ProjFolderFull(),'/dorothea_enrichments.xlsx'))
    # Plot
    library(ggplot2)
    g <- ggplot(f_contrast_acts, aes(x = reorder(source, score), y = score)) + 
      geom_bar(aes(fill = score), stat = "identity") +
      scale_fill_gradient2(low = "darkblue", high = "indianred", 
                           mid = "whitesmoke", midpoint = 0) + 
      theme_minimal() +
      theme(axis.title = element_text(face = "bold", size = 12),
            axis.text.x = 
              element_text(angle = 45, hjust = 1, size =10, face= "bold"),
            axis.text.y = element_text(size =10, face= "bold"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) +
      xlab("TF")
    png(paste0(ProjFolderFull(),'/Dorothea_regulons_300dpi.png'), 
        width = 8, 
        height = 8,
        res = 300, units = 'in')
    print(g)
    dev.off()
    
    
    files2zip <- list.files(path = ProjFolderFull(), pattern = 'png|xlsx|svg|html|quants.sf', full.names = F, recursive = F)
    print('printing files2zip')
    print(files2zip)
    
    fls2 <- list.files(path = paste0(ProjFolderFull(), '/trimmed'), pattern = 'html|sf', full.names = F, recursive = T)
    files2zip_2 <- paste0('trimmed/', fls2)
    
    print('printing files2zip_2')
    print(files2zip_2)
    
    print('printing length(fls2)')
    print(length(fls2))
    
    if(length(fls2)==0){
      filescombined <-files2zip
    } else {
      filescombined <- c(files2zip, files2zip_2)
    }
    
    
    
    zip::zip(root = ProjFolderFull(), zipfile = 'results.zip', files=filescombined, recurse = F, include_directories = F)
    click("downloadData")
    return(g)
    #}
    #else {
    #  print("loading dorothea_f_contrast_acts.RDS file")
    #  f_contrast_acts <- readRDS(paste0(ProjFolderFull(),'/dorothea_f_contrast_acts.RDS'))
    #  library(ggplot2)
    #  g <- ggplot(f_contrast_acts, aes(x = reorder(source, score), y = score)) + 
    #    geom_bar(aes(fill = score), stat = "identity") +
    #    scale_fill_gradient2(low = "darkblue", high = "indianred", 
    #                         mid = "whitesmoke", midpoint = 0) + 
    #    theme_minimal() +
    #    theme(axis.title = element_text(face = "bold", size = 12),
    #          axis.text.x = 
    #            element_text(angle = 45, hjust = 1, size =10, face= "bold"),
    #          axis.text.y = element_text(size =10, face= "bold"),
    #          panel.grid.major = element_blank(), 
    #          panel.grid.minor = element_blank()) +
    #    xlab("TF")
    #  
    #  files2zip <- list.files(path = ProjFolderFull(), pattern = 'png|xlsx', full.names = F, recursive = F)
    #  #files2zip <- dir(ProjFolderFull(), full.names = TRUE, recursive = F)
    #  #files2zip <- grep(pattern = 'trimmed|results.zip', x = files2zip, value = T, invert = T)
    #  print('zipping files2zip')
    #  print(files2zip)
    #  zip::zip(root = ProjFolderFull(), zipfile = 'results.zip', files = files2zip, recurse = F, include_directories = F)
    #  click("downloadData")
    #  return(g)
    #}
  })
  
  output$DorotheaMultitab<- renderUI({
    #req(input$groups_specified)
    message('line 2279 ok')
    if(salmon_finished()==0) return({})
    message('line 2281 ok')
    if(multiplegroups()==0) return({})
    message('line 2283 ok')
    g <- list()
    tytl <- list()
    net <- net() 
    f_contrast_acts_100_list <- list()
    f_contrast_acts_list <- list()
    
    message('line 2287 ok')
    
    #if(proceedWithLoad()==0){
    for (i in 1:length(res_txi_deseq())){
      tytl[[i]] <- paste(strsplit(names(res_txi_deseq()[i]), split = "_")[[1]][2:4], collapse=' ')
      deg <- res_DEGs_txi_deseq()[[i]]
      if(is.null(dim(deg))){
        g[[i]] <- ggplot()+theme_void()
      } else {
        
        #TOFIX only select certain samples 
        #TOFIX add dorothea per-sample composition
        counts <- txi()$abundance
        
        #counts[which(rowSums(counts)!=0),] -> meaningful_counts
        #meaningful_counts[rownames(deg),] -> meaningful_counts
        deg -> meaningful_deg
        
        #print('if rownames(meaningful_counts) are %in% gene2tx')
        #print(all(rownames(meaningful_counts) %in% gene2tx2name_GRCm39_u_systematic_u_onlyGeneNames[,1]))
        
        print('rownames(meaningful_counts):')  
        #print(head(rownames(meaningful_counts)))
        
        #print('rownames(meaningful_deg):')  
        #print(head(rownames(meaningful_deg)))
        
        #res1 <- data.frame(ENSEMBL=rownames(meaningful_counts), meaningful_counts)
        #annots1 <- AnnotationDbi::select(OrgDeeBee, keys=res1$ENSEMBL, 
        #                                 columns="SYMBOL", keytype="ENSEMBL")
        #meaningful_counts_symbol <- merge(annots1, res1, by.x="ENSEMBL", by.y="ENSEMBL")
        
        #rownames(meaningful_deg) <- meaningful_deg$ENSEMBL 
        #
        #res2 <- data.frame(ENSEMBL=rownames(meaningful_deg), meaningful_deg)
        #annots2 <- AnnotationDbi::select(OrgDeeBee, keys=res2$ENSEMBL, 
        #                                 columns="SYMBOL", keytype="ENSEMBL")
        #meaningful_deg_symbol <- merge(annots2, res2, by.x="ENSEMBL", by.y="ENSEMBL")
        meaningful_deg_symbol <- meaningful_deg
        
        #meaningful_counts_symbol2 <- meaningful_counts_symbol[which(!is.na(meaningful_deg_symbol$stat)),]
        meaningful_deg_symbol2 <- meaningful_deg_symbol[which(!is.na(meaningful_deg_symbol$stat)),]
        meaningful_deg_symbol2 <- meaningful_deg_symbol2[!is.na(meaningful_deg_symbol2$SYMBOL),]
        
        #toKeep <- which(!is.na(meaningful_counts_symbol2$SYMBOL));
        #meaningful_counts_symbol2 <- meaningful_counts_symbol2[toKeep,];
        #toKeep2 <- which(!duplicated(meaningful_counts_symbol2$SYMBOL));
        #meaningful_counts_symbol2 <- meaningful_counts_symbol2[toKeep2,-1]
        #rownames(meaningful_counts_symbol2) <- meaningful_counts_symbol2$SYMBOL
        #meaningful_counts_symbol2 <- meaningful_counts_symbol2[,-1];
        
        #meaningful_deg_symbol2 <- meaningful_deg_symbol2[toKeep,]
        #meaningful_deg_symbol2 <- meaningful_deg_symbol2[toKeep2,]
        #rownames(meaningful_deg_symbol2) <- rownames(meaningful_counts_symbol2)
        
        #saveRDS(meaningful_counts_symbol2, file = 'meaningful_counts_symbol2.RDS')
        #saveRDS(meaningful_deg_symbol2, file = 'meaningful_deg_symbol2.RDS')
        
        # rownames(meaningful_counts_symbol2) <- toupper(rownames(meaningful_counts_symbol2))
        meaningful_deg_symbol2$SYMBOL -> vec
        vec[is.na(vec)] <- 'none'
        rownames(meaningful_deg_symbol2) <- make.unique(toupper(vec))
        
        
        #  sample_acts <- decoupleR::run_ulm(mat=meaningful_counts_symbol2, net=net, .source='source', .target='target',
        #                                    .mor='mor', minsize = 5)
        contrast_acts <- decoupleR::run_ulm(mat=meaningful_deg_symbol2[,'stat', drop=F], net=net, .source='source', .target='target',
                                            .mor='mor', minsize = 5)
        
        n_tfs <- 25
        f_contrast_acts <- contrast_acts %>%
          mutate(rnk = NA)
        msk <- f_contrast_acts$score > 0
        f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
        f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))
        tfs <- f_contrast_acts %>%
          arrange(rnk) %>%
          head(n_tfs) %>%
          pull(source)
        f_contrast_acts <- f_contrast_acts %>%
          filter(source %in% tfs)
        
        n_tfs <- 100
        f_contrast_acts_100 <- contrast_acts %>%
          mutate(rnk = NA)
        msk_100 <- f_contrast_acts_100$score > 0
        f_contrast_acts_100[msk_100, 'rnk'] <- rank(-f_contrast_acts_100[msk_100, 'score'])
        f_contrast_acts_100[!msk_100, 'rnk'] <- rank(-abs(f_contrast_acts_100[!msk_100, 'score']))
        tfs_100 <- f_contrast_acts_100 %>%
          arrange(rnk) %>%
          head(n_tfs) %>%
          pull(source)
        f_contrast_acts_100 <- f_contrast_acts_100 %>%
          filter(source %in% tfs)
        
        f_contrast_acts_100_list[[i]] <- f_contrast_acts_100
        f_contrast_acts_list[[i]] <- f_contrast_acts
        
        message('line 2386 ok')
        
        # Plot
        library(ggplot2)
        g[[i]] <- ggplot(f_contrast_acts, aes(x = reorder(source, score), y = score)) + 
          geom_bar(aes(fill = score), stat = "identity") +
          scale_fill_gradient2(low = "darkblue", high = "indianred", 
                               mid = "whitesmoke", midpoint = 0) + 
          theme_minimal() +
          theme(axis.title = element_text(face = "bold", size = 12),
                axis.text.x = 
                  element_text(angle = 45, hjust = 1, size =10, face= "bold"),
                axis.text.y = element_text(size =10, face= "bold"),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank()) +
          xlab("TF")
        png(paste0(ProjFolderFull(),'/Dorothea_regulons_',tytl[[i]],'_300dpi.png'), 
            width = 8, 
            height = 8,
            res = 300, units = 'in')
        print(g[[i]])
        dev.off()
      }
    }
    
    message('line 2411 ok')
    names(f_contrast_acts_100_list) <- tytl
    names(f_contrast_acts_list) <- tytl
    
    saveRDS(f_contrast_acts_list, file = paste0(ProjFolderFull(),'/dorothea_f_contrast_acts.RDS'))
    saveRDS(f_contrast_acts_100_list, file = paste0(ProjFolderFull(),'/dorothea_f_contrast_acts_100.RDS'))
    
    f_contrast_acts_100_list_toPrint <- f_contrast_acts_100_list
    names(f_contrast_acts_100_list_toPrint) <- substr(x = names(f_contrast_acts_100_list_toPrint), start = 1, stop = 30)
    openxlsx::write.xlsx(x = f_contrast_acts_100_list_toPrint, file = paste0(ProjFolderFull(),'/dorothea_enrichments.xlsx'))
    saveRDS(g, file = paste0(ProjFolderFull(),'/dorothea_plots.RDS'))
    
    files2zip <- list.files(path = ProjFolderFull(), pattern = 'png|xlsx|svg|html|quants.sf', full.names = F, recursive = F)
    print('printing files2zip')
    print(files2zip)
    
    fls2 <- list.files(path = paste0(ProjFolderFull(), '/trimmed'), pattern = 'html|sf', full.names = F, recursive = T)
    files2zip_2 <- paste0('trimmed/', fls2)
    
    print('printing files2zip_2')
    print(files2zip_2)
    
    print('printing length(fls2)')
    print(length(fls2))
    
    if(length(fls2)==0){
      filescombined <-files2zip
    } else {
      filescombined <- c(files2zip, files2zip_2)
    }
    
    
    
    zip::zip(root = ProjFolderFull(), zipfile = 'results.zip', files=filescombined, recurse = F, include_directories = F)
    click("downloadData")
    
    
    #} else {
    #  print("loading dorothea_plots.RDS file")
    #  g <- readRDS(paste0(ProjFolderFull(), '/dorothea_plots.RDS'))
    #  
    #  files2zip <- list.files(path = ProjFolderFull(), pattern = 'png|xlsx', full.names = F, recursive = F)
    #  #files2zip <- dir(ProjFolderFull(), full.names = TRUE, recursive = F)
    #  #files2zip <- grep(pattern = 'trimmed|results.zip', x = files2zip, value = T, invert = T)
    #  print('zipping files2zip')
    #  print(files2zip)
    #  zip::zip(root = ProjFolderFull(), zipfile = 'results.zip', files = files2zip, recurse = F, include_directories = F)
    #  click("downloadData")
    #  
    #}
    #  
    #  for (j in 1:length(res_txi_deseq())){
    #    tytl[[j]] <- strsplit(res_txi_deseq()[[j]]@elementMetadata[5,2], split = ": ")[[1]][2]
    #  }
    
    
    nTabs = length(g)
    myTabs = lapply(1: nTabs, function(x){tabPanel(tytl[[x]], renderPlot(g[[x]], height = 1500))});
    
    return(do.call(tabsetPanel, myTabs))
  })
  
  
  reportFinished <-reactiveVal(0)
  observeEvent(input$button_save_colDatt,{
    if(salmon_finished()==0) return({})
    reportFinished(1)
    delay(10000, click("toSnap"))
  }
  )
  
  output$downloadData <- downloadHandler(
    filename = 'results.zip', 
    content = function(file) {
      file.copy(paste0(ProjFolderFull(),'/results.zip'), file)
    },
    contentType = "application/zip"
  )
  
  output$downloadDotPlot <- downloadHandler(
    filename = 'GO_dotplot_custom_number_of_Terms_300dpi.png', 
    content = function(file) {
      file.copy(paste0(ProjFolderFull(),'/GO_dotplot_custom_number_of_Terms_300dpi.png'), file)
    },
    contentType = 'image/png'
  )
  output$downloadDotPlots <- downloadHandler(
    filename = 'Plots_custom.zip', 
    content = function(file) {
      file.copy(paste0(ProjFolderFull(),'/Plots_custom.zip'), file)
    },
    contentType = "application/zip"
  )
  
}

ui <- secure_app(ui)
shinyApp(ui, server)
