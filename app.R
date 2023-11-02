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
#require('org.Hs.eg.db')
#require('org.Mm.eg.db')
options(shiny.maxRequestSize=10000*1024^2)

salmon <- '/home/zakiev/miniconda3/bin/conda run -n bioinfo salmon '
bbduk.sh <- 'source ~/.bashrc; bbduk.sh '
#fastqc_conda <- '/home/zakiev/miniconda3/bin/conda run -n bioinfo fastqc '
fastqc_powershell <- "pwsh C:/Users/flavial/fastqc/testing.ps1 "
wsl <- "wsl -e bash -c '"
gs <- function(x){ 
  y <- gsub(pattern = 'D:/', replacement = '/mnt/d/', x)
  return(gsub(pattern = 'C:/', replacement = '/mnt/c/', y))
}


credentials <- data.frame(
  user = c("shiny", "shinymanager", "emma", 'fabrice', 'aurelia','emile','maeva','aneta','elodie','manon','debora','baptiste'), # mandatory
  password = c("flaviallab", "flavial", "emma", 'fabrice','aurelia','emile','maeva','aneta','elodie','manon','debora','baptiste'), # mandatory
  #start = c("2019-04-15"), # optinal (all others)
  expire = c(NA, NA, NA),
  admin = c(FALSE, TRUE, FALSE),
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



ui <- fluidPage(
  #verbatimTextOutput("auth_output"),
  fluidRow(
    h1(strong("File Uploading (optional)"), style = "font-size:30px;"),
    column(width = 3,
           textInput(inputId = "whereToSavefastq", label = "Provide a name for the folder where to save your fastq files", value = "")),
    column(width = 3,
           fileInput(inputId = "uploadFASTQ", label ="Upload your own FASTQ files", multiple = T, placeholder = '.gz or .fastq.gz files', accept = c('.gz', '.fastq'))),
    column(width = 1, tags$img(src = "Possible_logo_30pct.png")),
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
  br(),
  rHandsontableOutput('GroupsPrompt'),
  br(),
  uiOutput("button_save_colDatt"),
  ### the stuff below this point is after the button "confirm"
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
  DT::DTOutput("TPM"),
  br(),
  br(),
  br(),
  textOutput("text_DEGs"),
  textOutput("text_which_is_Control_in_DEG"),
  #uiOutput('multitab_DEG'),
  DT::DTOutput("DESeq_DEGs"),
  textOutput("PCA_title"),
  plotOutput("PCA", height = "900px"),
  #imageOutput("PCA_image", height = "100%", width = "100%"),
  textOutput("Volcano_title"),
  plotOutput("Volcano", height = '900px'),
  #uiOutput('multitab_Volcano'),
  textOutput("Enrichments_title"),
  br(),
  tabsetPanel(id='Enrichments_panel', type = "tabs",
              tabPanel("Dotplots", 
                       tabsetPanel(id = 'tabset1', tabPanel("by overrepresentation test", plotOutput("GO_dotplot", height = '2000px')), 
                                   tabPanel("by GSEA",plotOutput("GSEA_dotplot", height = '2000px')))),
              tabPanel("Tables of terms", 
                       tabsetPanel(id = 'tabset2', tabPanel("by overrepresentation test", DT::DTOutput("GO")), 
                                   tabPanel("by GSEA",DT::DTOutput("GSEA"))))
  ),
  br(),
  textOutput("Dorothea_title"),
  plotOutput("Dorothea", height = "750px", width = "750px"),
  #downloadButton("report", "Generate report")
  #screenshotButton(download = TRUE, id = "plot", filename = "poopity.scoop")
  #actionButton("toSnap", "Take a snapshot of the current state"),
  downloadButton("downloadData", "Download"),
  
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

server <- function(input, output, session){
  shiny_home <- c('C:/Users/flavial/Documents/Shiny_0.9.0/ShinyData/')
  house <- 'D:/FASTQ/'
  
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
      salmon_finished(1)
      proceedWithLoad(1)
      fastqc_finished(1)
      bbduk_finished(1)
    } else 
    {
      toOutput <- "no project results detected at selected destination"
    }
    return(toOutput)
  })
  
  ProjFolder <- reactive({
    as.character(input$ProjectName)
  })
  
  #output$confirm<-renderUI({
  #  req(isTruthy(input$infiles) | isTruthy(input$infolder))
  #  if(length(input$infiles) <= 1){
  #    if(!is.null(input$infolder)){ return({})} else {
  #      labl <- "Load the existing results"
  #    }
  #  } else {labl <- "Start analysis"}
  #  if(fastqc_finished()==1) return({})
  #  actionButton(inputId = 'confirm', label = labl, icon("paper-plane")) 
  #})
  
  detected_filenames_anomaly <- reactiveVal(0)
  fastqc_finished <- reactiveVal(0)
  bbduk_finished <- reactiveVal(0)
  salmon_finished <- reactiveVal(0)
  bbduk_and_fastqc_reused <- reactiveVal(0)
  
  
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
      toRet <- paste0('D:/FASTQ/', paste0(unl[-length(unl)], collapse = '/'))
      print('printing input$infolder')
      print(toRet)
      return(toRet)
    }
  })
  
  trimmed_folder <- reactive({
    req(input$infiles)
    paste0(ProjFolderFull(), '/trimmed/')
  })
  
  indexPath <- reactive(
    {req(input$referenceGenomeChoice)
      #Human GRCh38"=1,"Mouse GRCm39"=2,"Mouse C57B6"=3,"Mouse MM129"=4
      #print(input$referenceGenomeChoice)
      if(input$referenceGenomeChoice==1) indexPath <- paste0(shiny_home,'/Homo_sapiens_GRCh38')
      if(input$referenceGenomeChoice==2) indexPath <- paste0(shiny_home,'/GRCm39_mm_index')
      if(input$referenceGenomeChoice==3) indexPath <- paste0(shiny_home,'/c57_b6_index')
      if(input$referenceGenomeChoice==4) indexPath <- paste0(shiny_home,'/MM129_mm_index')
      return(indexPath)
    })
  txdbPath <- reactive(
    {req(input$referenceGenomeChoice)
      #Human GRCh38"=1,"Mouse GRCm39"=2,"Mouse C57B6"=3,"Mouse MM129"=4
      #print(input$referenceGenomeChoice)
      if(input$referenceGenomeChoice==1) txdbPath <- paste0(shiny_home, '/txdb_hsapiens.RDS')
      if(input$referenceGenomeChoice==2) txdbPath <- paste0(shiny_home, '/txdb_GRCm39.RDS')
      if(input$referenceGenomeChoice==3) txdbPath <- paste0(shiny_home, '/txdb_c57bl6.RDS')
      if(input$referenceGenomeChoice==4) txdbPath <- paste0(shiny_home, '/txdb_MM129.RDS')
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
    if(length(grep(pattern = "L002", x = fls()))>0){
      print("line 435 detected multiplexed reads")
      detected_multiplexedReads(1)
    }
    if(length(grep(pattern = "R1", x = fls()))>0){
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
      
      R1_tracks <- grep(pattern = "R1", x = demultiplexed_fls(), value = T)
      R1_tracks_swapped <- gsub(pattern = "R1", replacement = "R2", R1_tracks)
      R2_tracks <- grep(pattern = "R2", x = demultiplexed_fls(), value = T)
      R2_tracks_swapped <- gsub(pattern = "R2", replacement = "R1", R2_tracks)
      
      R1_matched_tracks <- R1_tracks[which(R1_tracks_swapped %in% R2_tracks)]
      R1_tracks[which(!R1_tracks_swapped %in% R2_tracks)] -> R1_unbalanced_tracks
      R2_tracks[which(!R2_tracks_swapped %in% R1_tracks)] -> R2_unbalanced_tracks
      
      if(length(R1_unbalanced_tracks)>0|length(R2_unbalanced_tracks)>0) 
        detected_filenames_anomaly(1)
      
      R1_R2_matched_tracks_list <- list()
      for (i in 1:length(R1_matched_tracks)){
        R1_R2_matched_tracks_list[[i]] <- c(R1_matched_tracks[i], gsub(pattern = "R1", replacement = "R2", R1_matched_tracks[i]))
      }
      if(length(R1_unbalanced_tracks)>0){
        for (i in 1:length(R1_unbalanced_tracks)){
          R1_R2_matched_tracks_list <- append(R1_R2_matched_tracks_list, list(c(R1_unbalanced_tracks[i], paste0('[warning!]: file ',  gsub(pattern = "R1", replacement = "R2", R1_unbalanced_tracks[i]), ' was expected, but not provided, are you sure you selected all paired-end reads?'))))
        }
      }
      if(length(R2_unbalanced_tracks)>0){
        for (i in 1:length(R2_unbalanced_tracks)){
          R1_R2_matched_tracks_list <- append(R1_R2_matched_tracks_list, list(c(paste0('[warning!]: file ',  gsub(pattern = "R2", replacement = "R1", R2_unbalanced_tracks[i]), ' was expected, but not provided, are you sure you selected all paired-end reads?'), R2_unbalanced_tracks[i])))
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
    system2("powershell", args=paste0(wsl, 'rm ', gs(dirname(fls()[1])), "/*_Finished'"))
    system2("powershell", args=paste0('mkdir ', ProjFolderFull(), '/trimmed -ea 0'))
    if(detected_multiplexedReads()==1){
      reference_tracks <- grep(pattern = "L001", x = fls(), value = T)
      for (i in 1:length(reference_tracks)){
        reference_track <- reference_tracks[i]
        resulting_track <- gsub(pattern = "L001", replacement = "", x = reference_track)
        resulting_track <- gsub(pattern = "__", replacement = "_", x = resulting_track)
        reference_track <- gsub(pattern = "L001", replacement = "L00*", x = reference_track)
        resulting_track_full <- paste0(ProjFolderFull(),"/", basename(resulting_track))
        cmmnd <- paste0("type ", reference_track, " > ", resulting_track_full)
        system2("powershell", args=cmmnd)
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
      cmmnd <- paste0("; filename_dir=$(dirname $filename1); filename_small1=$(basename $filename1); filename_small2=$(basename $filename2); extension=${filename_small1##*.};  ", bbduk.sh, " in1=${filename1} in2=${filename2} out1=", gs(trimmed_folder()), "${filename_small1}.trimmed.${extension} out2=", gs(trimmed_folder()), "${filename_small2}.trimmed.${extension} ref=/home/zakiev/miniconda3/envs/bioinfo/opt/bbmap-39.01-1/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo")
      ref<- grep(pattern = "R1", x = files2bbduk, value = T)
      match_ref <- gsub("R1", replacement = "R2", x = ref)
      
      for (i in 1:length(ref)){
        flname <- c(ref[i], match_ref[i])
        print("printing the exact command that is passed to bbduk")
        arggg <- paste0(wsl, "filename1=", gs(flname[1]), "; filename2=",gs(flname[2]),cmmnd,"; touch ${filename1}_Finished'")
        print(arggg)
        system2("powershell", args=arggg, wait = F)
      }
      filenames2monitor <- paste0(ref, "_Finished")
      print("printing filenames2monitor")
      print(filenames2monitor)
      while (!all(file.exists(filenames2monitor))) {
        Sys.sleep(1)
      }
      Sys.sleep(20)
      print("printing the exact command that is passed to fastqc")
      #arggg <- paste0("pwsh C:/Users/flavial/fastqc/testing.ps1 (ls '", trimmed_folder(),"*.fastq*').fullname")
      #arggg <- c(,"C:/Users/flavial/fastqc/testing.ps1", paste0("(ls '", trimmed_folder(),"*.fastq*').fullname"))
      arggg <- paste0(fastqc_powershell, paste(unlist(final_fls()), collapse = " "))
      print(arggg)
      system(arggg)
    } else {
      cmmnd <- paste0("; filename_dir=$(dirname $filename); filename_small=$(basename $filename); extension=${filename_small##*.}; ", bbduk.sh," in=${filename} out=", gs(trimmed_folder()), "${filename_small}.trimmed.${extension} ref=/home/zakiev/miniconda3/envs/bioinfo/opt/bbmap-39.01-1/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo")
      for (i in 1:length(files2bbduk)){
        flname <- files2bbduk[[i]]
        print("printing the exact command that is argggg")
        argggg <- paste0(wsl, "filename=", gs(flname), cmmnd,"; touch ${filename}_Finished'")
        system2("powershell", args=argggg, wait=F)
      }
      filenames2monitor <- paste0(unlist(files2bbduk), "_Finished")
      print("printing filenames2monitor")
      print(filenames2monitor)
      while (!all(file.exists(filenames2monitor))) {
        Sys.sleep(1)
      }
      Sys.sleep(10)
      print("printing the exact command that is argggg")
      argggg <- paste0(fastqc_powershell, paste(unlist(final_fls()), collapse = " "))
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
    
    inputfiles <- paste0(trimmed_folder(),basename(demultiplexed_fls()),'.trimmed.', tools::file_ext(demultiplexed_fls()[1]))
    print("line and 559 ok")
    
    if(detected_pairedEndReads()==1){
      print("line after 559 ok")
      reference_initial_files <- grep(pattern = "R1", x = inputfiles, value = T)
      matched_reference_initial_files <- gsub(pattern = "R1", replacement = "R2", x = reference_initial_files)
      for (i in 1:length(reference_initial_files)){
        toReturn[[i]] <- c(reference_initial_files[i], matched_reference_initial_files[i])
      }
      names(toReturn) <- gsub("_R1", replacement = "", x = reference_initial_files)
    } else {
      print("else line after 559 ok")
      reference_initial_files <- inputfiles
      for (i in 1:length(reference_initial_files)){
        toReturn[[i]] <- inputfiles[i]
      }
      names(toReturn) <- toReturn
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
    req(input$infiles)
    if(length(input$infiles) <= 1) return({})
    #validate(
    #  need(length(which(colDatt$data[,3]==T))<=1, "Please check only one checkbox or don't select at all if you want the first group (alphabetically ordered) to be the control")
    #)
    actionButton(inputId = 'groups_specified', label = "Validate groups and start analysis", icon("paper-plane")) 
  })
  
  #observeEvent(input$start_fastqc,{
  #  fastqc_finished(1)
  #}
  #)
  
  #observeEvent(input$start_salmon,{
  #  salmon_finished(1)
  #})
  #observeEvent(input$reuse_salmon,{
  #  salmon_finished(1)
  #})
  
  observeEvent(input$groups_specified,{
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
    startingfastqc_command <- paste0(fastqc_powershell, paste(files2fastqc, collapse = " "))
    print(startingfastqc_command)
    system(startingfastqc_command)
    cutext <- function(x){ y = gsub(pattern = '.fastq', replacement = '', x); return(gsub(pattern = '.gz', replacement = '', y))}
    movingcommand <- function(x) paste0("$file='", cutext(x), "'; $dir2move2='",ProjFolderFull(),"';$fl2mv=($file+'_fastqc.zip'); mv $fl2mv $dir2move2")
    print("printing one moving command:")
    i=1; print(movingcommand(files2fastqc[i]))
    for (i in 1:length(files2fastqc)) system2(command = 'powershell', args=movingcommand(files2fastqc[i]))
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
    if(proceedWithLoad()==0){
      print("line 895 was ok")
      fls <- final_fls()
      fls <- unlist(fls)
      fls <- gsub('.fastq', replacement = '', x = fls)
      fls <- gsub('.gz', replacement = '', x = fls)
      fls <- gsub('.fq', replacement = '', x = fls)
      
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
    } else {
      print("reading fastqc_report_after_trimming.RDS file")
      return(readRDS(paste0(ProjFolderFull(),'/fastqc_report_after_trimming.RDS')))
    }
  })
  
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
  
  
  observeEvent(input$ProjectName, {
    shinyjs::hide(id = "fastqcstatspanel")
    shinyjs::hide(id = "fastqcstatspanel_AfterTrimming")
  })
  observeEvent(input$groups_specified, {
    print("line 909 ok")
    if(fastqc_finished()==1)
      shinyjs::show(id = "fastqcstatspanel")
  })
  observe({
    req(input$groups_specified)
    print("line 913 ok")
    if(bbduk_finished()==1){ 
      print("bbduk's finished yepcock")
      shinyjs::show(id = "fastqcstatspanel_AfterTrimming")
    }
  })
  observe({
    req(input$infolder)
    shinyjs::show(id = "fastqcstatspanel")
    shinyjs::show(id = "fastqcstatspanel_AfterTrimming")
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
  #observe({
  #  req(input$groups_specified)
  #  print("line 913 ok")
  #  if(nrow()){ 
  #    print("bbduk's finished yepcock")
  #    shinyjs::show(id = "")
  #  }
  #})
  
  
  
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
        arggggg <- paste0(wsl, "filename=",gs(flname),"; ", salmon, " quant -i ", gs(indexPath())," -l A -1 ", gs(flname1)," -2 ", gs(flname2)," -p 16 --validateMappings -o ${filename}_quant", "'")
        print('printing the arggggg passed to salmon')
        print(arggggg)
        system2("powershell", args=arggggg)
      }
    } else {
      for (i in 1:length(initial_files)){
        flname <- initial_files[[i]]
        arggggg <- paste0(wsl, "filename=",gs(flname),"; ", salmon, " quant -i ", gs(indexPath())," -l A -r ", gs(flname)," -p 16 --validateMappings -o ${filename}_quant", "'")
        print('printin the arggggg passed to salmon')
        print(arggggg)
        system2("powershell", args=arggggg)
      }
    }
    print("salmon is officially finished!!!!!!")
    salmon_finished(1)
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
      names(TPM_filelist) <- tools::file_path_sans_ext(basename(filelist))
      txdb <- readRDS(txdbPath())
      txi_ <- tximport::tximport(TPM_filelist, type = "salmon", tx2gene = txdb, ignoreTxVersion = T)
      saveRDS(txi_, file = paste0(ProjFolderFull(),'/txi.RDS'))
      return(txi_)
    } else {
      salmon_finished(1)
      print("reading txi.RDS file")
      return(readRDS(paste0(ProjFolderFull(),'/txi.RDS')))
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
    openxlsx::write.xlsx(result, file = paste0(ProjFolderFull(), '/TPMs.xlsx'))
    return(result)
  })
  output$salmon_log <- renderRHandsontable({
    if(proceedWithLoad()==0){
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
      
    } else {
      filelist <- readRDS(file = paste0(ProjFolderFull(),"/filelist.RDS"))
      rates <- readRDS(file = paste0(ProjFolderFull(),"/rates.RDS"))
      openxlsx::write.xlsx(x = data.frame(File=filelist, Rate=rates), file = 'alignment_rates.xlsx')
      print("line 1061 was ok")
      #return(data.frame(AlignmentLog=grep(pattern = "Index contained|Mapping rate|warning", toPrint, value = T )))
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
    }
  }
  )
  output$text_TPM<- renderText({
    if(salmon_finished()==0) return({})
    return("Transcripts-per-million (TPM) values:")
  })
  output$TPM <- DT::renderDT(txi_tpms(), extensions = 'Buttons', 
                             options = list(
                               dom = 'Bfrtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
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
      print(datt)
      rhandsontable(datt) %>%
        hot_col("Sample", readOnly = TRUE)
    } else {
      datt <- readRDS(file = paste0(ProjFolderFull(),"/colData.RDS"))
      colDatt$data <- datt
      rhandsontable(datt)
    }
  })
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
  
  txi_deseq <- reactive({
    if(proceedWithLoad()==0){
      if(salmon_finished()==0) return({})
      toRet <- DESeq2::DESeqDataSetFromTximport(txi(), colData = colDatt$data, design = ~Group)
      saveRDS(toRet, file = paste0(ProjFolderFull(),'/txi_deseq.RDS'))
      return(toRet)
    } else {
      print("reading txi_deseq.RDS file")
      return(readRDS(paste0(ProjFolderFull(),'/txi_deseq.RDS')))
    }
    
  })
  txi_deseq_deseq <- reactive({
    if(proceedWithLoad()==0){
      if(salmon_finished()==0) return({})
      txi_deseq() -> matr
      colDatt$data -> colDat
      reff <- colDat[which(colDat[,3]),2]
      #reff <- matr$Group[which(colDat[,3])]
      #print(reff)
      if(length(which(colDat[,3]))>0) matr$Group <- relevel(x = matr$Group, ref = reff)
      smallestGroupSize <- floor(nrow(colDat))
      keep <- rowSums(counts(matr) >= 10) >= smallestGroupSize
      matr <- matr[keep,]
      toRet <- DESeq2::DESeq(matr)
      saveRDS(toRet, file = paste0(ProjFolderFull(),'/txi_deseq_deseq.RDS'))
      return(toRet)
    }
    else {
      print("reading txi_deseq.RDS file")
      return(readRDS(paste0(ProjFolderFull(),'/txi_deseq.RDS')))
    }
  })
  res_txi_deseq <- reactive({
    if(proceedWithLoad()==0){
      if(salmon_finished()==0) return({})
      if(length(resultsNames(txi_deseq_deseq()))<=2){
        toRet <- results(txi_deseq_deseq())
      }  else {
        toRet <- list()
        for (i in 2:length(resultsNames(txi_deseq_deseq()))){
          contr <- strsplit(resultsNames(dds_Deb_relevel_NTN1ko)[3], split = '_')[[1]]; 
          toRet[[(i-1)]] <- results(txi_deseq_deseq(), contrast = c('Group', contr[2], contr[4]))
        }
        
      }
      saveRDS(toRet, file = paste0(ProjFolderFull(),'/res_txi_deseq.RDS'))
      return(toRet)
    }
    else {
      print("reading txi_deseq.RDS file")
      return(readRDS(paste0(ProjFolderFull(),'/res_txi_deseq.RDS')))
    }
  })
  res_DEGs_txi_deseq <- reactive({
    if(proceedWithLoad()==0){
      if(salmon_finished()==0) return({})
      toRet <- res_txi_deseq()
      if(length(toRet)<500){
        toRet <- toRet[which(toRet$padj<0.05),]
        if(nrow(toRet)==0) return({as.data.frame("oops, none significant")})
        toRet <- toRet[order(toRet$padj),]
        toRet2 <- data.frame(ENSEMBL=rownames(toRet), toRet[,1:2], abs_log2FoldChange=abs(toRet[,2]),toRet[,3:ncol(toRet)])
        annots <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=rownames(toRet2), 
                                        columns="SYMBOL", keytype="ENSEMBL")
        
        toRet3 <- merge(annots, toRet2, by.x="ENSEMBL", by.y="ENSEMBL")
        openxlsx::write.xlsx(as.data.frame(toRet3), file = paste0(ProjFolderFull(),'/DEGs.xlsx'))
        saveRDS(toRet3, file = paste0(ProjFolderFull(),'/res_DEGs_txi_deseq.RDS'))
        return(toRet3)
      } else {
        toRet3 <- list()
        for (i in 1:length(toRet)){
          toRet[[i]] <- toRet[[i]][which(toRet[[i]]$padj<0.05),]
          if(nrow(toRet[[i]])==0) toRet[[i]] <- "oops, none significant"
          toRet[[i]] <- toRet[[i]][order(toRet[[i]]$padj),]
          toRet2 <- data.frame(ENSEMBL=rownames(toRet[[i]]), toRet[[i]][,1:2], abs_log2FoldChange=abs(toRet[[i]][,2]),toRet[[i]][,3:ncol(toRet[[i]])])
          annots <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=rownames(toRet2), 
                                          columns="SYMBOL", keytype="ENSEMBL")
          
          toRet3[[i]] <- merge(annots, toRet2, by.x="ENSEMBL", by.y="ENSEMBL")
        }
        names(toRet3) <- names(toRet)
        openxlsx::write.xlsx(toRet3, file = paste0(ProjFolderFull(),'/DEGs.xlsx'))
        saveRDS(toRet3, file = paste0(ProjFolderFull(),'/res_DEGs_txi_deseq.RDS'))
      }
    }
    else {
      print("reading res_DEGs_txi_deseq.RDS file")
      return(readRDS(paste0(ProjFolderFull(),'/res_DEGs_txi_deseq.RDS')))
    }
  })
  
  
  output$text_DEGs <- renderText({
    if(salmon_finished()==0) return({})
    "List of differentially expressed genes:"
  })
  output$text_which_is_Control_in_DEG <- renderText({
    if(salmon_finished()==0) return({})
    paste0("The control (background) is ", tail(strsplit(res_txi_deseq()@elementMetadata[5,2], split = " ")[[1]], 1))
  })
  
  n_contrasts <- reactive({
    if(salmon_finished()==0) return({})
    return(length(resultsNames(res_DEGs_txi_deseq())))
  })
  
  output$multitab_DEG <- renderUI(do.call(tabsetPanel, c(id='t',lapply(1:(n_contrasts()-1), function(i) {
    tabPanel(
      title=paste0("Contrast ", i), 
      DTOutput(paste0('a',i))
    )
  }))))
  
  
  #lapply(1:n, function(j) {
  #  if(n_contrasts() <=2) output[[paste0('a',j)]] <- NULL else {
  #    output[[paste0('a',j)]] <- DT::renderDT(res_DEGs_txi_deseq()[[j]], filter = 'top', extensions = 'Buttons', 
  #                                            options = list(
  #                                              dom = 'Bfrtip',
  #                                              buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
  #                                            server=FALSE, rownames = FALSE,
  #                                            caption = htmltools::tags$caption(
  #                                              style = 'caption-side: bottom; text-align: left;',
  #                                              'can also be found at  ', htmltools::em(paste0(ProjFolder(),'/DEGs.xlsx'))
  #                                            )
  #    )
  #  }
  #  return(output)
  #})
  
  
  output$multitab_Volcano <- renderUI(do.call(tabsetPanel, c(id='t2',lapply(1:(n_contrasts()-1), function(i) {
    tabPanel(
      title=paste0("Contrast ", i), 
      DTOutput(paste0('b',i))
    )
  }))))
  
  
  
  output$DESeq_DEGs <- DT::renderDT(res_DEGs_txi_deseq(), filter = 'top', extensions = 'Buttons', 
                                    options = list(
                                      dom = 'Bfrtip',
                                      buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
                                    server=FALSE, rownames = FALSE
  )
  
  output$PCA_title<- renderText({
    if(salmon_finished()==0) return({})
    return("PCA plot:")
  })
  
  output$Enrichments_title<- renderText({
    if(salmon_finished()==0) return({})
    return("Enriched terms:")
  })
  
  output$PCA<- renderPlot({
    if(salmon_finished()==0) return({})
    p <- plotPCA(DESeq2::vst(txi_deseq()), intgroup="Group")
    p <- p + ggrepel::geom_label_repel(label=colnames(txi_deseq()), max.overlaps = 30)+ theme_classic()
    png(paste0(ProjFolderFull(), '/PCA_plot_300dpi.png'), 
        width = 7, 
        height = 7,
        res = 300, units = 'in')
    print(p)
    dev.off()
    return(p)
  })
  
  
  output$Volcano_title<- renderText({
    if(salmon_finished()==0) return({})
    return("Volcano plot:")
  })
  output$Volcano<- renderPlot({
    if(salmon_finished()==0) return({})
    keyvals <- ifelse(
      ((res_txi_deseq()$log2FoldChange < 0) & (res_txi_deseq()$padj < 0.05)), 'red',
      ifelse(((res_txi_deseq()$log2FoldChange > 0) & (res_txi_deseq()$padj < 0.05)), 'green',
             'black'))
    
    keyvals[which(is.na(keyvals))]  <- 'black'
    
    toKeyvals <- data.frame(ENSEMBL=rownames(res_txi_deseq()))
    annot<-AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=rownames(res_txi_deseq()), 
                                 columns="SYMBOL", keytype="ENSEMBL")
    common_row_names <- intersect(toKeyvals[[1]], annot$ENSEMBL)
    toKeyvals <- merge(annot, toKeyvals, by.x="ENSEMBL", by.y="ENSEMBL")
    toKeyvals <- toKeyvals[!duplicated(toKeyvals[,1]),]
    names(keyvals) <- toKeyvals$SYMBOL
    names(keyvals)[keyvals == 'red'] <- 'downreg'
    names(keyvals)[keyvals == 'green'] <- 'upreg'
    names(keyvals)[keyvals == 'black'] <- 'NS'
    
    percentiles <- quantile(res_txi_deseq()$log2FoldChange, probs=c(0.005,0.995), na.rm =T);
    xlim <- c(-max(abs(percentiles)), max(abs(percentiles)))
    
    p <- EnhancedVolcano::EnhancedVolcano(res_txi_deseq(),
                                          lab = keyvals,
                                          x = 'log2FoldChange',
                                          y = 'padj',
                                          title = strsplit(res_txi_deseq()@elementMetadata[5,2], split = ": ")[[1]][2],
                                          pCutoff = 5e-2,
                                          FCcutoff = 100,
                                          pointSize = 2.0,
                                          labFace = 'bold',
                                          colCustom = keyvals,
                                          labSize = 4.0,
                                          gridlines.major = FALSE,
                                          gridlines.minor = FALSE, 
                                          subtitle = "",
                                          boxedLabels = TRUE,
                                          drawConnectors = T#,
                                          #ylim = ylim_NEO,
                                          #xlim = xlim)
    )
    png(paste0(ProjFolderFull(), '/Volcano_plot_300dpi.png'), 
        width = 15, 
        height = 15,
        res = 300, units = 'in')
    print(p)
    dev.off()
    return(p)
  })
  m_t2g_most_symbol <- reactive({
    if(input$referenceGenomeChoice!=1) 
      return(readRDS('~/Analyses/RNAseq/m_t2g_most_symbol.RDS')) else {
        return(readRDS('~/Analyses/RNAseq/m_t2g_most_symbol_human.RDS'))
      }
  }
  )
  
  OrgDeeBee <- reactive({
    if(input$referenceGenomeChoice!=1) 
      return(org.Mm.eg.db::org.Mm.eg.db) else {
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
    #print("printing res_DEGs_txi_deseq()[,2]")
    #print(res_DEGs_txi_deseq()[,2])
    if(proceedWithLoad()==0){
      if(salmon_finished()==0) return({})
      toRet <- enrichGO(gene = res_DEGs_txi_deseq()[,2], 
                        keyType = "SYMBOL", 
                        OrgDb = OrgDeeBee(), 
                        ont = "BP", 
                        pAdjustMethod = "BH", 
                        qvalueCutoff = 0.05, 
                        readable = TRUE)
      saveRDS(toRet, file = paste0(ProjFolderFull(),'/GO_result.RDS'))
      openxlsx::write.xlsx(as.data.frame(toRet), file = paste0(ProjFolderFull(),'/GOs.xlsx'))
      return(toRet)
    }  else {
      print("reading GO_result.RDS file")
      return(readRDS(paste0(ProjFolderFull(),'/GO_result.RDS')))
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
    as.data.frame(GO_result())
  }, extensions = 'Buttons', 
  options = list(
    dom = 'Bfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
  server=FALSE, rownames=FALSE)
  output$GSEA_dotplot<- renderPlot({
    if(salmon_finished()==0) return({})
    p <- dotplot(GSEA_result(), showCategory=50)
    png(paste0(ProjFolderFull(),'/GSEA_dotplot_300dpi.png'), 
        width = 10, 
        height = 20,
        res = 300, units = 'in')
    print(p)
    dev.off()
    return(p)
  })
  output$GO_dotplot<- renderPlot({
    if(salmon_finished()==0) return({})
    p <- dotplot(GO_result(), showCategory=50)
    png(paste0(ProjFolderFull(),'/GO_dotplot_300dpi.png'), 
        width = 10, 
        height = 20,
        res = 300, units = 'in')
    print(p)
    dev.off()
    return(p)
  })
  
  output$Dorothea_title<- renderText({
    if(salmon_finished()==0) return({})
    return("Dorothea regulons:")
  })
  
  output$Dorothea <- renderPlot({
    if(proceedWithLoad()==0){
      if(salmon_finished()==0) return({})
      deg <- res_txi_deseq()
      counts <- txi()$abundance
      
      #print('printing counts')
      #print(counts)
      
      #print("printing deg")
      #print(deg)
      
      #print('printing nrow(counts)')
      #print(nrow(counts))
      
      #print("printing nrow(deg)")
      #print(nrow(deg))
      
      #print("printing if rownames(counts)==rownames(deg)")
      #print(all(rownames(counts)==rownames(deg)))
      
      
      if(input$referenceGenomeChoice==1){
        net <- decoupleR::get_collectri(organism='human', split_complexes=FALSE)
      } else {
        net <- decoupleR::get_collectri(organism='mouse', split_complexes=FALSE)
      }
      
      counts[which(rowSums(counts)!=0),] -> meaningful_counts
      meaningful_counts[rownames(deg),] -> meaningful_counts
      deg -> meaningful_deg
      
      print('if rownames(meaningful_counts) are %in% gene2tx')
      #print(all(rownames(meaningful_counts) %in% gene2tx2name_GRCm39_u_systematic_u_onlyGeneNames[,1]))
      
      print('rownames(meaningful_counts):')  
      print(head(rownames(meaningful_counts)))
      
      print('rownames(meaningful_deg):')  
      print(head(rownames(meaningful_deg)))
      
      res1 <- data.frame(ENSEMBL=rownames(meaningful_counts), meaningful_counts)
      annots1 <- AnnotationDbi::select(OrgDeeBee(), keys=res1$ENSEMBL, 
                                       columns="SYMBOL", keytype="ENSEMBL")
      meaningful_counts_symbol <- merge(annots1, res1, by.x="ENSEMBL", by.y="ENSEMBL")
      
      res2 <- data.frame(ENSEMBL=rownames(meaningful_deg), meaningful_deg)
      annots2 <- AnnotationDbi::select(OrgDeeBee(), keys=res2$ENSEMBL, 
                                       columns="SYMBOL", keytype="ENSEMBL")
      meaningful_deg_symbol <- merge(annots2, res2, by.x="ENSEMBL", by.y="ENSEMBL")
      
      
      meaningful_counts_symbol2 <- meaningful_counts_symbol[which(!is.na(meaningful_deg_symbol$stat)),]
      meaningful_deg_symbol2 <- meaningful_deg_symbol[which(!is.na(meaningful_deg_symbol$stat)),]
      
      toKeep <- which(!is.na(meaningful_counts_symbol2$SYMBOL));
      meaningful_counts_symbol2 <- meaningful_counts_symbol2[toKeep,];
      toKeep2 <- which(!duplicated(meaningful_counts_symbol2$SYMBOL));
      meaningful_counts_symbol2 <- meaningful_counts_symbol2[toKeep2,-1]
      rownames(meaningful_counts_symbol2) <- meaningful_counts_symbol2$SYMBOL
      meaningful_counts_symbol2 <- meaningful_counts_symbol2[,-1];
      
      meaningful_deg_symbol2 <- meaningful_deg_symbol2[toKeep,]
      meaningful_deg_symbol2 <- meaningful_deg_symbol2[toKeep2,]
      rownames(meaningful_deg_symbol2) <- rownames(meaningful_counts_symbol2)
      
      #saveRDS(meaningful_counts_symbol2, file = 'meaningful_counts_symbol2.RDS')
      #saveRDS(meaningful_deg_symbol2, file = 'meaningful_deg_symbol2.RDS')
      
      rownames(meaningful_counts_symbol2) <- toupper(rownames(meaningful_counts_symbol2))
      rownames(meaningful_deg_symbol2) <- toupper(rownames(meaningful_deg_symbol2))
      
      
      sample_acts <- decoupleR::run_ulm(mat=meaningful_counts_symbol2, net=net, .source='source', .target='target',
                                        .mor='mor', minsize = 5)
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
      openxlsx::write.xlsx(x = as.data.frame(f_contrast_acts_100), file = 'dorothea_enrichments.xlsx')
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
      
      files2zip <- list.files(path = ProjFolderFull(), pattern = 'png|xlsx', full.names = F, recursive = F)
      #files2zip <- grep(pattern = 'trimmed|results.zip', x = files2zip, value = T, invert = T)
      print('zipping files2zip')
      print(files2zip)
      zip::zip(root = ProjFolderFull(), zipfile = 'results.zip', files = files2zip, recurse = F, include_directories = F)
      click("downloadData")
      return(g)
    }
    else {
      print("loading dorothea_f_contrast_acts.RDS file")
      f_contrast_acts <- readRDS(paste0(ProjFolderFull(),'/dorothea_f_contrast_acts.RDS'))
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
      
      files2zip <- list.files(path = ProjFolderFull(), pattern = 'png|xlsx', full.names = F, recursive = F)
      #files2zip <- dir(ProjFolderFull(), full.names = TRUE, recursive = F)
      #files2zip <- grep(pattern = 'trimmed|results.zip', x = files2zip, value = T, invert = T)
      print('zipping files2zip')
      print(files2zip)
      zip::zip(root = ProjFolderFull(), zipfile = 'results.zip', files = files2zip, recurse = F, include_directories = F)
      click("downloadData")
      return(g)
    }
  })
  
  reportFinished <-reactiveVal(0)
  observeEvent(input$button_save_colDatt,{
    if(salmon_finished()==0) return({})
    reportFinished(1)
    delay(10000, click("toSnap"))
  }
  )
  
  output$report <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "~/poopity_report.html",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      #params <- list(n = input$slider)
      params <- list(NA)
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  
  output$downloadData <- downloadHandler(
    filename = 'results.zip', 
    content = function(file) {
      file.copy(paste0(ProjFolderFull(),'/results.zip'), file)
    },
    contentType = "application/zip"
  )
  
  
  observeEvent(input$toSnap, {
    screenshot(download = F, server_dir = ProjFolderFull(), filename = "logOfCommands")
    #delay(1000, system2("powershell", args=""))
  })
  #session$allowReconnect("force")
}

ui <- secure_app(ui)
shinyApp(ui, server)