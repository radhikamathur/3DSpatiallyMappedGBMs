## R script for Shiny app (https://3d-gbms.shinyapps.io/search/)

library(pacman)

p_load(shiny, ggplot2, tidyverse, ggpubr, rgl, misc3d)
ModelPath <- './Models/'
PatientColors = c('P455' = "#ff7500", 'P475' = "#ae7000", "P498" = "#44cef6","P500" = "#1bd1a5", "P503" = "#FFD92F", "P519" = "#8d4bbb","P521" = "#ff0097", "P524" = "#BEBEBE", "P529" = "#0B0B45", "P530" = "#FF0000")
RNAcolors <- c('Classical' = "#3e5063", 'Mesenchymal' = '#0433ff', 'Proneural' = '#ff40ff', 'Neural' = '#ff9300')

# setwd('~/Dropbox/Postdoc/Papers/ATAC Paper/Shiny App/')
# p_load(openxlsx, EnvStats, janitor)
# dataPath <- '~/Dropbox/Postdoc/Papers/ATAC Paper/Revision/Supplementary Tables/'
# RNA <- read.xlsx(paste0(dataPath, 'TableS5 - Transcriptome analysis.xlsx'), sheet = 2)%>% pivot_longer(cols = !Gene, names_to = "Sample", values_to = "Signal")
# InputFile <-  read.xlsx(paste0(dataPath, 'TableS1 - Patient and sample information.xlsx'),sheet =3, startRow = 2, )
# SampleInfo <- InputFile%>%dplyr::select(Patient, Sample, Purity, Dist, RNAsubtype, ID, CEL)%>%mutate(Sample = as.factor(Sample))
# RNAModules <- read.xlsx(paste0(dataPath, 'TableS5 - Transcriptome analysis.xlsx'), fillMergedCells = TRUE, sheet = 3) %>%row_to_names(1)%>%mutate(ModuleID = paste0("R_", `RNA module (R_)`), Genes = as.numeric(Genes), PurityR = as.numeric(PurityR), DistR = as.numeric(DistR))%>%relocate(ModuleID)
# RNAModules <- RNAModules%>%select(ModuleID|Genes|`Label/Top GSEA result (FDR)`|13:ncol(RNAModules))%>% pivot_longer(cols = !c(ModuleID, Genes,`Label/Top GSEA result (FDR)`) , names_to = "Sample", values_to = "Signal")
# RNAModuleGenes <- read.xlsx(paste0(dataPath, 'TableS5 - Transcriptome analysis.xlsx'),sheet = 4)%>%mutate(ModuleID = paste0("R_", Module)) %>%select( -Module)
# ATACModules <- read.xlsx(paste0(dataPath, 'TableS6 - Chromatin landscape analysis.xlsx'), fillMergedCells = TRUE, sheet = 3) %>%row_to_names(1)%>%mutate(Peaks = as.numeric(Peaks), PurityR = as.numeric(PurityR), DistR = as.numeric(DistR))
# ATACModules <- ATACModules%>%select(ModuleID|Peaks|`Top GSEA result (FDR)`|`Top Motif (p-value)`|12:ncol(ATACModules))%>% pivot_longer(cols = !c(ModuleID, Peaks, `Top GSEA result (FDR)`, `Top Motif (p-value)`), names_to = "Sample", values_to = "Signal")
# ATACModuleRankedPeaks <- read.xlsx(paste0(dataPath, 'TableS6 - Chromatin landscape analysis.xlsx'), fillMergedCells = TRUE, sheet = 4)%>%select(-chr, -start, -end, -Genes_n)%>%mutate(Rank = as.integer(Rank))%>%filter(Rank<=10)
# ATACModuleTopLinkedGene <- read.xlsx(paste0(dataPath, 'TableS6 - Chromatin landscape analysis.xlsx'),sheet = 5)%>%mutate(Linkages = as.integer(Linkages))%>%group_by(ModuleID) %>%top_n(n=10)%>%ungroup()
# scATAC_geneactivity <- read.xlsx(paste0(dataPath, 'TableS7 - Single-cell analysis.xlsx'), sheet = 2)%>%mutate(Sample = as.factor(as.integer(Sample)))
# scATAC_peaks <- read.xlsx(paste0(dataPath, 'TableS7 - Single-cell analysis.xlsx'), sheet = 3)%>%mutate(Sample = as.factor(as.integer(Sample)))
# scATAC_geneactivity_longer <- scATAC_geneactivity%>% pivot_longer(cols = 8:ncol(scATAC_geneactivity), names_to = "ModuleID", values_to = "Score") %>%arrange(Score)
# scATAC_peaks_longer <- scATAC_peaks%>% pivot_longer(cols = 8:ncol(scATAC_peaks), names_to = "ModuleID", values_to = "Score") %>%arrange(Score)
# scATAC_overview <- scATAC_geneactivity %>%select(1:7)
# scATAC_all <- bind_rows('Gene Activity Score' = scATAC_geneactivity_longer, 'Read-in-peak score'= scATAC_peaks_longer, .id = "Type")%>%select(Cell, Type, ModuleID, Score)
# RNA <- RNA %>% separate(Sample, c("Patient", "Sample"))%>% mutate(Signal = as.numeric(Signal), Sample = as.factor(as.numeric(Sample)))
# RNAModuleInfo <- RNAModules %>% select(-Sample, -Signal)%>%unique()
# RNAModules <- RNAModules %>%select(ModuleID, Sample, Signal) %>% separate(Sample, c("Patient", "Sample"))%>% mutate(Signal = as.numeric(Signal), Sample = as.factor(as.numeric(Sample)))
# ATACModuleInfo <- ATACModules %>%select(-Sample, -Signal)%>%unique()
# ATACModuleInfo$`Top GSEA result (FDR)`[ATACModuleInfo$ModuleID == "L_salmon4"] <- "EGFR amp"
# ATACModules <- ATACModules %>% select(ModuleID, Sample, Signal) %>% separate(Sample, c("Patient", "Sample"))%>% mutate(Signal = as.numeric(Signal), Sample = as.factor(as.numeric(Sample)))
# FACETSResults <- read.xlsx(paste0(dataPath, 'TableS3 - Genomic characterization.xlsx'), sheet = 5)
# FACETSResults$Gene[FACETSResults$chrom == 7 & FACETSResults$start<55324313 & FACETSResults$end>55086714] <- "EGFR"
# FACETSResults$Gene[FACETSResults$chrom == 4 & FACETSResults$start<55164414 & FACETSResults$end>55095264] <- "PDGFRA"
# FACETSResults$Gene[FACETSResults$chrom == 1 & FACETSResults$start<204542871 & FACETSResults$end>204485511] <- "MDM4"
# FACETSResults$Gene[FACETSResults$chrom == 2 & FACETSResults$start<16087129 & FACETSResults$end>16080683] <- "MYCN"
# FACETSResults$Gene[FACETSResults$chrom == 9 & FACETSResults$start<21975132 & FACETSResults$end>21967751] <- "CDKN2A"
# FACETSResults$Gene[FACETSResults$chrom == 10 & FACETSResults$start<89728532 & FACETSResults$end>89623195] <- "PTEN"
# FACETSResults$Gene[FACETSResults$chrom == 13 & FACETSResults$start<48878049 & FACETSResults$end>49054207] <- "RB1"
# FACETSResults$Gene[FACETSResults$chrom == 17 & FACETSResults$start<7572927 & FACETSResults$end>7579569] <- "TP53"
# Fusions <- read.xlsx(paste0(dataPath, 'TableS4 - Structural variants and fusion genes.xlsx'), sheet = 3)%>%filter(confidence != "low")%>%mutate(Fusion = paste0(gene1, "-", gene2)) %>%filter(!str_detect(Fusion, "RP11"))%>%filter(!str_detect(Fusion, ","))%>%filter(!str_detect(Fusion, "\\."))%>%select(Patient, Sample, Fusion) %>% unique()%>%mutate(Signal = 1, Sample = as.factor(Sample))
# CopyNumber <- FACETSResults %>% filter(!is.na(Gene))%>%select(Patient, Sample, Gene, Signal = tcn.em)%>%mutate(Sample = as.factor(Sample))
# DistanceMatrix <- read.xlsx(paste0(dataPath, 'TableS2 - Pairwise inter-sample distances.xlsx'), sheet = 3)
# 
# save(DistanceMatrix, Fusions, CopyNumber, RNA, SampleInfo, RNAModuleInfo, ATACModuleInfo, RNAModules, RNAModuleGenes, scATAC_all,scATAC_overview, ATACModules, ATACModuleRankedPeaks, ATACModuleTopLinkedGene, file = "~/Dropbox/Postdoc/Papers/ATAC Paper/Shiny App/ShinyAppData.RData")

load("ShinyAppData.RData")

GeneList <- RNA %>%distinct(Gene)%>%arrange(Gene) %>% pull(Gene)
PatientList <- SampleInfo %>% filter(Patient != "P565") %>% distinct(Patient)%>%pull(Patient)
CNList <- CopyNumber %>%arrange(Gene) %>% distinct(Gene)%>%pull(Gene)
RNAModuleList <- RNAModules %>%distinct(ModuleID)%>%pull(ModuleID)
ATACModuleList <- ATACModules %>%distinct(ModuleID)%>% pull(ModuleID)
scATAC_modulelist <- scATAC_all %>%distinct(ModuleID) %>% arrange(ModuleID) %>% pull(ModuleID)

colorByFeatureMain <- function(vector){
  rbPal <- colorRampPalette(c("blue","red"))
  if (length(unique(vector)) > 1){
    mappedColors <- rbPal(length(vector))[as.numeric(cut(vector,breaks = length(vector)))][1:length(vector)]
  } else {
    mappedColors <- rep('blue', length(vector))
  }
  names(mappedColors) <- names(vector)
  return(mappedColors)
}

ui <- fluidPage(
  navbarPage(title = "Search in 3D spatially mapped GBMs", id = "navigation",
             
  tabPanel("About",tags$hr(),
    sidebarLayout(sidebarPanel(width = 3,
    tags$h4("About"),
    tags$div("Online platform for 3D visualization and data exploration designed to accompany the manuscript"),
    tags$br(),
    tags$b("'Glioblastoma evolution and heterogeneity from a 3D whole-tumor perspective'"),
    tags$br(),tags$br(),
    tags$div("Created by Radhika Mathur, Ph.D."),
    tags$div("Costello Lab, UCSF Department of Neurosurgery"),
    tags$em("radhika.mathur@ucsf.edu")),
    mainPanel(
    imageOutput('welcomeimage', height = 600),
    tags$span("This platform enables exploration of genomic, epigenomic, and transcriptomic datasets across 3D spatially mapped tumors from 10 IDH-WT GBM patients. Search by individual patient for 3D visualization, by gene or module of interest for summary figures across all patients, or explore our single-nucleus (sn-) ATAC dataset.")),
  )),
             
             
  tabPanel("3D visualization",
    sidebarLayout(sidebarPanel(width = 4,
      tags$h5("Search by patient for interactive 3D visualization of datasets"), tags$hr(),
      selectizeInput(inputId = "SpecificPatient", label = "Patient:", choices = PatientList, selected = "P521"),
      selectizeInput(inputId = "DataSet", label = "Visualize:", choices = c("Purity", "Distance", "Copy number", "Gene expression", "Fusions", "RNA module", "ATAC or linkage module"), selected = "Purity"),
      uiOutput("ReactiveSelection"),
      textOutput(outputId = "RNAModuleInfo3D", container = tags$i),
      textOutput(outputId = "ATACModuleInfo3D", container = tags$i),
      checkboxInput(inputId = "BrainPlot", label = "Show patient brain", value = FALSE, width = NULL)),
      
    mainPanel(
      rglwidgetOutput("model3D", height = "600"),
      htmlOutput("units"),
      imageOutput('image', height = "30", width="100"),
      htmlOutput('colorbartext'),
      tableOutput(outputId = "VisualizationTable"),
  ))),
             
  tabPanel("Gene expression",
    sidebarLayout(sidebarPanel(width = 2,
    tags$h5("Search by gene to see tumor-wide expression patterns"),tags$hr(),
      selectizeInput(inputId = "SpecificGene", label = "Gene: ", choices =NULL),
      textOutput(outputId = "genemoduleText", container = tags$i)),
    mainPanel(
      textOutput(outputId = "geneText", container = tags$h4),
      splitLayout(cellWidths = c("20%", "80%"), plotOutput(outputId = "genepatientPlot"),plotOutput(outputId = "genesamplePlot")),
      splitLayout(cellWidths = c("33%", "33%", "34%"), plotOutput(outputId = "genepurityPlot"), plotOutput(outputId = "genedistPlot"), plotOutput(outputId = "genesubtypePlot")),
  ))),

  tabPanel("Transcriptomic programs",
    sidebarLayout(sidebarPanel(width = 2,
      tags$h5("Search by RNA module to see tumor-wide expression patterns"),tags$hr(),
      selectizeInput(inputId = "SpecificRNAModule", label = "RNA Module: ", choices =RNAModuleList, selected = "R_blue")),
    mainPanel(
      textOutput(outputId = "RNAmoduleText", container = tags$h4),
      tableOutput(outputId = "genetable"),
      splitLayout(cellWidths = c("20%", "80%"), plotOutput(outputId = "patientPlot"),plotOutput(outputId = "samplePlot")),
      splitLayout(cellWidths = c("33%", "33%", "34%"), plotOutput(outputId = "purityPlot"), plotOutput(outputId = "distPlot"), plotOutput(outputId = "subtypePlot")),
      splitLayout(cellWidths = c("60%", "40%", "34%"), plotOutput(outputId = "RNAmatrixPlot"), plotOutput(outputId = "celPlot"))
  ))),

  tabPanel("Chromatin programs",
    sidebarLayout(sidebarPanel(width = 2,
          tags$h5("Search by ATAC or Linkage module to see tumor-wide activation patterns"),tags$hr(),
          selectizeInput(inputId = "SpecificATACModule", label = "ATAC or Linkage Module: ", choices =ATACModuleList, selected = "A_magenta4")),
      mainPanel(
         textOutput(outputId = "ATACmoduleText", container = tags$h4),
          splitLayout(cellWidths = c("70%", "30%"), tableOutput(outputId = "ATACtable_rankedpeaks"), tableOutput(outputId = "ATACtable_toplinkedgenes")),
          splitLayout(cellWidths = c("16%", "84%"), plotOutput(outputId = "ATACpatientPlot"),plotOutput(outputId = "ATACsamplePlot")),
          splitLayout(cellWidths = c("33%", "33%", "34%"), plotOutput(outputId = "ATACpurityPlot"), plotOutput(outputId = "ATACdistPlot"), plotOutput(outputId = "ATACsubtypePlot")),
          splitLayout(cellWidths = c("60%", "40%", "34%"), plotOutput(outputId = "ATACmatrixPlot"), plotOutput(outputId = "ATACcelPlot"))
    ))),

  tabPanel("Explore snATAC data",
    sidebarLayout(sidebarPanel(width = 2,
          tags$h5("Explore datasets at the single-cell level"),tags$hr(),
          selectizeInput(inputId = "scATAC_select", label = "Select data type: ", choices =c('Sample', 'Neoplastic status', 'Module'), selected = 'Sample'),
          uiOutput("scATAC_reactive"),
          textOutput(outputId = "sc_ModuleInfo", container = tags$i)),
      mainPanel(
          plotOutput(outputId = "exploreATAC", height = 700, width = 700),
          plotOutput(outputId = "sc_barcharts", width = 700)
  ))),
  

))

server <- function(input, output, session) {

  ## search by patient

  session$onSessionEnded(function() {
    stopApp()
  })
  
  observeEvent(input$DataSet, {
    updateSelectizeInput(session, "SelectedGene", choices = GeneList,  selected = "EGFR", server = TRUE)
  })
  
  
  output$ReactiveSelection <- renderUI({
    if (input$DataSet == "Gene expression"){selectizeInput("SelectedGene", "Gene: ", choices = NULL)}
    else if (input$DataSet == "RNA module"){selectizeInput("SelectedRNAModule", "Module: ", choices = RNAModuleList, selected = "R_brown")}
    else if (input$DataSet == "ATAC or linkage module"){selectizeInput("SelectedATACModule", "Module: ", choices = ATACModuleList, selected = "L_salmon4")}
    else if(input$DataSet == "Fusions"){
      req(input$SpecificPatient)
      FusionList <- Fusions %>%filter(Patient == input$SpecificPatient) %>%arrange(Fusion)%>%distinct(Fusion) %>%pull(Fusion)
      SelectedFusion <- FusionList[1]
      selectizeInput("SelectedFusion", "Fusion: ", choices = FusionList, selected = SelectedFusion)
      }
    else if(input$DataSet == "Copy number"){selectizeInput("SelectedCN", "Gene:", choices = CNList, selected = "EGFR")}
  })


  output$RNAModuleInfo3D <-renderText({
    req(input$SelectedRNAModule)
    if(input$DataSet == "RNA module"){
    ModuleInfo <- RNAModuleInfo %>%filter(ModuleID == input$SelectedRNAModule) %>% pull(`Label/Top GSEA result (FDR)`)
    ModuleGeneN <- RNAModuleInfo %>%filter(ModuleID == input$SelectedRNAModule)  %>% pull(`Genes`)
    paste0(ModuleGeneN, " constituent genes; ",ModuleInfo)}
  })
  
  output$ATACModuleInfo3D <- renderText({
    req(input$SelectedATACModule)
    if(input$DataSet == "ATAC or linkage module"){
    SpecificModuleLabel <-  ATACModuleInfo %>%filter(ModuleID == input$SelectedATACModule)%>%pull(`Top GSEA result (FDR)`)
    SpecificModulePeaks <-  ATACModuleInfo %>%filter(ModuleID == input$SelectedATACModule)%>%pull(Peaks)
    SpecificModuleMotif <-  ATACModuleInfo %>%filter(ModuleID == input$SelectedATACModule)%>%pull(`Top Motif (p-value)`)
    paste0(SpecificModulePeaks, " constituent peaks; enriched for ", SpecificModuleLabel, " and ", SpecificModuleMotif, " motif")}
  })
  
  dataValues <- reactive({
    req(input$DataSet, input$SpecificPatient)
    if (input$DataSet == "Purity"){
      dataValues <- SampleInfo %>%filter(Patient == input$SpecificPatient) %>% pull(Purity)
      names(dataValues) = SampleInfo %>%filter(Patient == input$SpecificPatient) %>% pull(Sample)
    }
    else if (input$DataSet == "Distance"){
      dataValues <- SampleInfo %>%filter(Patient == input$SpecificPatient) %>% filter(!is.na(Dist)) %>% pull(Dist)
      names(dataValues) = SampleInfo %>%filter(Patient == input$SpecificPatient) %>% filter(!is.na(Dist)) %>% pull(Sample)
    }
    else if (input$DataSet == "Fusions"){
      req(input$SelectedFusion)
      dataTables <- Fusions %>%filter(Fusion == input$SelectedFusion) 
      dataValues <- SampleInfo %>%filter(!is.na(RNAsubtype)) %>%filter(Patient == input$SpecificPatient)%>%select(Patient, Sample) %>%left_join(dataTables)%>% mutate(Signal = replace_na(Signal, 0))%>%mutate(Signal = as.numeric(Signal)) %>%pull(Signal)
      names(dataValues) = SampleInfo %>%filter(!is.na(RNAsubtype))%>%filter(Patient == input$SpecificPatient)%>%select(Patient, Sample) %>%left_join(dataTables) %>% pull(Sample)
      if (var(dataValues) == 0){dataValues <- append(dataValues, 0)}
    }
    else if (input$DataSet == "Gene expression"){
      req(input$SelectedGene)
      dataTables <- RNA %>% filter(Gene == input$SelectedGene & Patient == input$SpecificPatient) %>%filter(!is.na(Signal))
      dataValues <- dataTables%>%pull(Signal)
      names(dataValues) = dataTables%>% pull(Sample)
    }
    else if (input$DataSet == "RNA module"){
      req(input$SelectedRNAModule)
      dataTables <- RNAModules %>%filter(ModuleID == input$SelectedRNAModule & Patient == input$SpecificPatient) %>%filter(!is.na(Signal))
      dataValues <- dataTables%>%pull(Signal)
      names(dataValues) = dataTables%>% pull(Sample)
    }
    else if (input$DataSet == "ATAC or linkage module"){
      req(input$SelectedATACModule)
      dataTables <- ATACModules %>%filter(ModuleID == input$SelectedATACModule & Patient == input$SpecificPatient) %>%filter(!is.na(Signal))
      dataValues <- dataTables%>%pull(Signal)
      names(dataValues) = dataTables%>% pull(Sample)
    }
    else if (input$DataSet == "Copy number"){
      req(input$SelectedCN)
      dataTables <- CopyNumber %>%filter(Gene == input$SelectedCN & Patient == input$SpecificPatient) %>%filter(!is.na(Signal))
      dataValues <- dataTables%>%pull(Signal)
      names(dataValues) = dataTables%>% pull(Sample)
    }
   
    return(dataValues)
  })
  
  
  output$image <- renderImage({list(src = 'www/colorbar.png', height=30, width=100, align="center")}, deleteFile = F)
  
  output$colorbartext <- renderUI({ 
    min <- as.character(round(min(dataValues()),1))
    max <- as.character(round(max(dataValues()),1))
    HTML(paste0('<p>', min, ' &nbsp&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp ',max, '</p>'))
  })
  
  output$units <- renderUI({
    req(input$DataSet)
    if (input$DataSet=="Purity"){unitsForData <- "Purity ψ (proportion of tumor cells in sample)"} 
    else if (input$DataSet == "Distance") {unitsForData <- "Relative distance (d) from tumor centroid (0) to periphery (1)"}
    else if (input$DataSet == "Gene expression") {unitsForData <- "Gene expression in log2tpm+1"}
    else if (input$DataSet == "Fusions") {unitsForData <- "Fusion gene (0 = absent, 1 = present)"}
    else if (input$DataSet == "Copy number") {unitsForData <- "Total copy number"}
    else if (input$DataSet == "RNA module") {unitsForData <- "Average expression of constituent genes in log2tpm+1"}
    else if (input$DataSet == "ATAC or linkage module") {unitsForData <- "Average ATAC signal of constituent peaks in log2cpm+1"}
    HTML(as.character(unitsForData))
  })
  
  output$model3D <- renderRglwidget({
    req(input$DataSet, input$SpecificPatient)
    try(close3d(), silent = TRUE)
    sampleCoordinates <- readRDS(paste0(ModelPath, input$SpecificPatient, '/coordinates_samples.rds'))
    tumorModel <- readRDS(paste0(ModelPath, input$SpecificPatient, '/tumor_t2.rds'))
    brainModel<- readRDS(paste0(ModelPath, input$SpecificPatient,  '/coordinates_brain_periphery.rds'))
    voxelToMM <- readRDS(paste0(ModelPath, input$SpecificPatient,'/adj.rds'))
    dtemp <- dim(tumorModel)
    maximumDimension <- max(dtemp)
    dtempScaled <- dtemp/maximumDimension # the aspect that fits the original voxel matrix, scaled to 1
    aspectInMM <- dtempScaled * c(voxelToMM$x, voxelToMM$y, voxelToMM$z) # and this then converts the voxel dimensions to mm
    
    if (input$BrainPlot == T){
      message('Creating brain contour and plotting brain')
      plot3d(brainModel, alpha=0.01, col='#726665', axes=F,  xlab = "", ylab = "", zlab= "", aspect = aspectInMM)}
    
    message('Creating tumor contour and plotting tumor')
    tumor <- contour3d(tumorModel, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3], level = 1, alpha = .2, add = TRUE, draw = TRUE, color = 'yellow')

    if (input$SpecificPatient == "P530"){
      sampleCoordinates_frontal <- readRDS(paste0(ModelPath, 'P530 frontal/coordinates_samples.rds'))
      sampleCoordinates <- rbind(sampleCoordinates, sampleCoordinates_frontal)
      rownames(sampleCoordinates) <- 1:19
      
      tumorModel <- readRDS(paste0(ModelPath, 'P530 frontal/tumor_t2.rds'))
      voxelToMM <- readRDS(paste0(ModelPath, 'P530 frontal/adj.rds'))
      dtemp <- dim(tumorModel)
      maximumDimension <- max(dtemp)
      dtempScaled <- dtemp/maximumDimension # the aspect that fits the original voxel matrix, scaled to 1
      aspectInMM <- dtempScaled * c(voxelToMM$x, voxelToMM$y, voxelToMM$z) # and this then converts the voxel dimensions to mm
      
      message('Creating frontal tumor contour and plotting tumor')
      tumor2 <- contour3d(tumorModel, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3], level = 1, alpha = .2, add = TRUE, draw = TRUE, color = 'yellow')
    }
    
    colors <- colorByFeatureMain(dataValues())
    samplesToPlot <- rownames(sampleCoordinates)[rownames(sampleCoordinates) %in% names(colors)]
    sampleCoordinates <- sampleCoordinates[samplesToPlot,]
    sampleColors <- colors[samplesToPlot]
    points3d(x=sampleCoordinates[,1], y=sampleCoordinates[,2], z=sampleCoordinates[,3], level = 1, size = 7, color=sampleColors)
    text3d(x=sampleCoordinates[,1], y=sampleCoordinates[,2], z=sampleCoordinates[,3], texts = rownames(sampleCoordinates), cex=1, adj=-.3)
    
    if (input$BrainPlot == T){view3d(zoom = 0.5)} else{view3d()}
    
    
    rglwidget()

  })

  output$VisualizationTable <- renderTable({
    dataValues() %>% bind_rows
  })
  
  
  ### search by gene
  updateSelectizeInput(session, "SpecificGene", choices = GeneList,  selected = "AIF1", server = TRUE)
   
    output$geneText <- renderText({
      paste0(input$SpecificGene, " gene expression (log2tpm+1)")
    })

    output$genemoduleText <- renderText({
      req(input$SpecificGene %in% RNAModuleGenes$Gene)
      SpecificModule <-  RNAModuleGenes %>%filter(Gene == input$SpecificGene)%>%pull(ModuleID)
      ModuleInfo <- RNAModuleInfo %>%filter(ModuleID == SpecificModule)%>% pull(`Label/Top GSEA result (FDR)`)
      paste0("Constituent gene of ", SpecificModule, "; ", ModuleInfo)
    })

    RNA_SpecificGene <- reactive({
      req(input$SpecificGene)
      RNA_SpecificGene <- RNA %>% filter(Gene == input$SpecificGene)  %>% inner_join(SampleInfo)%>%rename(RNA=Signal)
      return(RNA_SpecificGene)
    })
      
    output$genepatientPlot <- renderPlot({
      ggplot(RNA_SpecificGene(), aes(x=Patient, y = RNA, fill = Patient))+
      geom_boxplot()+
      geom_point()+
      scale_fill_manual(values =PatientColors)+
      labs(title = "across patients", x = NULL, y = "Gene expression")+
      theme_bw(base_size = 14)+
      theme(legend.position = "none")+
      theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    })

    output$genesamplePlot <- renderPlot({
        ggplot(RNA_SpecificGene(), aes(x=Sample, y = RNA, fill = Patient))+
        geom_bar(stat = "identity", color = "black")+
        scale_fill_manual(values =PatientColors)+
        facet_grid(.~Patient, scales = "free", space = "free")+
        labs(title = "across samples", y = "Gene expression")+
        theme_bw(base_size = 14)+
        theme(strip.background = element_rect(colour="white", fill="white"))+
        theme(strip.text.y = element_text( size = 14, angle = 0))+
        theme(legend.position = "none")+
        theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    })

     output$genepurityPlot <- renderPlot({
    ggplot(RNA_SpecificGene(), aes(x=Purity, y = RNA, color = Patient, group = "none"))+
        geom_point()+
        labs(title = "by purity", x= "Sample purity (ψ)", y = "Gene Expression")+
        theme_bw(base_size = 14)+
        theme(strip.background = element_rect(colour="white", fill="white"))+
        theme(strip.text.y = element_text( size = 14, angle = 0))+
        theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
        geom_smooth(method='lm')+
        stat_cor(method = "pearson", size =5,label.x.npc = "center")+
        theme(legend.position = "none")+
        scale_color_manual(values = PatientColors)
    })

    output$genedistPlot <- renderPlot({
        ggplot(RNA_SpecificGene(), aes(x=Dist, y = RNA, color = Patient, group = "none"))+
        geom_point()+
        labs(title = "by location", x= "Relative distance from centroid (d)", y = "Gene Expression")+
        theme_bw(base_size = 14)+
        theme(strip.background = element_rect(colour="white", fill="white"))+
        theme(strip.text.y = element_text( size = 14, angle = 0))+
        theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
        geom_smooth(method='lm')+
        stat_cor(method = "pearson", size =5, label.x.npc = "center")+
        theme(legend.position = "none")+
        scale_color_manual(values = PatientColors)
    })

    output$genesubtypePlot <- renderPlot({
      ggplot(RNA_SpecificGene(), aes(x=RNAsubtype, y = RNA, color = Patient, group = RNAsubtype))+
      geom_boxplot()+
      geom_jitter()+
      labs(title = "by subtype", x= "RNAsubtype", y = "Gene Expression")+
      theme_bw(base_size = 14)+
      theme(strip.background = element_rect(colour="white", fill="white"))+
      theme(strip.text.y = element_text( size = 14, angle = 0))+
      theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
      scale_color_manual(values = PatientColors)+
      theme(legend.position = "none")+
      stat_compare_means(method = "anova", size =5)
    })

    ### search by RNA module

    output$RNAmoduleText <- renderText({
      SpecificModuleLabel <-  RNAModuleInfo %>%filter(ModuleID == input$SpecificRNAModule)%>%pull(`Label/Top GSEA result (FDR)`)
      SpecificModuleGenes <-  RNAModuleInfo %>%filter(ModuleID == input$SpecificRNAModule)%>%pull(Genes)
      paste0(input$SpecificRNAModule, " contains ", SpecificModuleGenes, " constituent genes; ", SpecificModuleLabel)
    })

    output$genetable <- renderTable({
      RNAModuleGenes %>%filter(ModuleID == input$SpecificRNAModule) %>%select(-ModuleID)%>%head(n=10)%>%mutate(Rank = as.integer(Rank))
    }, bordered = T)
    
    RNA_SpecificModule <- reactive({
      req(input$SpecificRNAModule)
      RNA_SpecificModule <- RNAModules %>% filter(ModuleID == input$SpecificRNAModule)%>%rename(RNA=Signal)%>% inner_join(SampleInfo)
      return(RNA_SpecificModule)
    })
    
    RNA_matrix <- reactive({
      req(input$SpecificRNAModule)
      RNA_matrix <- DistanceMatrix %>%filter(ModuleID == input$SpecificRNAModule)
      return(RNA_matrix)
    })

    output$patientPlot <- renderPlot({
        ggplot(RNA_SpecificModule(), aes(x=Patient, y = RNA, fill = Patient))+
        geom_boxplot()+
        geom_point()+
        scale_fill_manual(values =PatientColors)+
        theme_bw(base_size = 14)+
        labs(title = "across patients", x = NULL, y = "Average expression")+
        theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
        theme(legend.position = "none")+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    })

    output$samplePlot <- renderPlot({
        ggplot(RNA_SpecificModule(), aes(x=Sample, y = RNA, fill = Patient))+
        geom_bar(stat = "identity", color = "black")+
        scale_fill_manual(values =PatientColors)+
        facet_grid(.~Patient, scales = "free", space = "free")+
        labs(title = "across samples", y = NULL)+
        theme_bw(base_size = 14)+
        theme(strip.background = element_rect(colour="white", fill="white"))+
        theme(strip.text.y = element_text( size = 14, angle = 0))+
        theme(legend.position = "none")+
        theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    })

    output$purityPlot <- renderPlot({
        ggplot(RNA_SpecificModule(), aes(x=Purity, y = RNA, color = Patient, group = "none"))+
        geom_point()+
        labs(title = "by purity", x= "Sample purity (ψ)", y = "Average expression")+
        theme_bw(base_size = 14)+
        theme(strip.background = element_rect(colour="white", fill="white"))+
        theme(strip.text.y = element_text( size = 14, angle = 0))+
        theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
        geom_smooth(method='lm')+
        stat_cor(method = "pearson", size =5)+
        theme(legend.position = "none")+
        scale_color_manual(values = PatientColors)
    })

    output$distPlot <- renderPlot({
        ggplot(RNA_SpecificModule(), aes(x=Dist, y = RNA, color = Patient, group = "none"))+
        geom_point()+
        labs(title = "by location", x= "Relative distance from centroid (d)", y = "Average expression")+
        theme_bw(base_size = 14)+
        theme(strip.background = element_rect(colour="white", fill="white"))+
        theme(strip.text.y = element_text( size = 14, angle = 0))+
        theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
        geom_smooth(method='lm')+
        stat_cor(method = "pearson", size =5)+
        theme(legend.position = "none")+
        scale_color_manual(values = PatientColors)
    })

    output$subtypePlot <- renderPlot({
        ggplot(RNA_SpecificModule(), aes(x=RNAsubtype, y = RNA, color = Patient, group = RNAsubtype))+
        geom_boxplot()+
        geom_jitter()+
        labs(title = "by subtype", x= "RNAsubtype", y = "Average expression")+
        theme_bw(base_size = 14)+
        theme(strip.background = element_rect(colour="white", fill="white"))+
        theme(strip.text.y = element_text( size = 14, angle = 0))+
        theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
        scale_color_manual(values = PatientColors)+
        theme(legend.position = "none")+
        stat_compare_means(method = "anova", size =5)
    })
    
    output$RNAmatrixPlot <- renderPlot({
      ggplot(RNA_matrix(), aes(x=Distance, y = ValueDifference, color = Patient, group = "none"))+
        geom_point()+
        geom_smooth(method='lm', se = F)+
        stat_cor(method = "pearson")+
        scale_color_manual(values = PatientColors)+
        theme_bw(base_size = 14)+
        theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
        labs(title = "distribution within tumors", x = "Inter-sample distance", y = "Difference in average expression")
    })
    
    output$celPlot <- renderPlot({
      CEL_plot <-  RNA_SpecificModule() %>%filter(CEL != "Border")
      ggplot(CEL_plot, aes(x = CEL, y = RNA, color = RNAsubtype, group = CEL))+
        geom_boxplot(outlier.shape = NA)+
        geom_jitter()+
        theme_bw(base_size = 14)+
        scale_color_manual(values = RNAcolors)+
        labs(title = "contrast-enhancing lesion", x = "Contrast-enhancing",y = "Average expression" )+
        theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
        stat_compare_means(method = "t.test", label = "p.format")
    })

    ### search by ATAC module

    output$ATACmoduleText <- renderText({
      SpecificModuleLabel <-  ATACModuleInfo %>%filter(ModuleID == input$SpecificATACModule)%>%pull(`Top GSEA result (FDR)`)
      SpecificModulePeaks <-  ATACModuleInfo %>%filter(ModuleID == input$SpecificATACModule)%>%pull(Peaks)
      SpecificModuleMotif <-  ATACModuleInfo %>%filter(ModuleID == input$SpecificATACModule)%>%pull(`Top Motif (p-value)`)
      paste0(input$SpecificATACModule, " contains ", SpecificModulePeaks, " constituent peaks; ", SpecificModuleLabel, " and ", SpecificModuleMotif, " motif")
    })

    output$ATACtable_rankedpeaks <- renderTable({
      ATACModuleRankedPeaks %>%filter(ModuleID == input$SpecificATACModule)%>%select(-ModuleID)%>%rename(`Linked gene(s)` = Genes) 
    }, bordered = T)

    output$ATACtable_toplinkedgenes <- renderTable({
      ATACModuleTopLinkedGene %>%filter(ModuleID == input$SpecificATACModule) %>%select(-ModuleID)%>%head(n=10)%>%rename(`Top Gene`=Gene)
    }, bordered = T)
    
    ATAC_SpecificModule <- reactive({
      req(input$SpecificATACModule)
      ATAC_SpecificModule <- ATACModules %>% filter(ModuleID == input$SpecificATACModule) %>%rename(ATAC=Signal) %>% inner_join(SampleInfo)
      return(ATAC_SpecificModule)
    })
    
    ATAC_matrix <- reactive({
      req(input$SpecificATACModule)
      ATAC_matrix <- DistanceMatrix %>%filter(ModuleID == input$SpecificATACModule)
      return(ATAC_matrix)
    })

    output$ATACpatientPlot <- renderPlot({
        ggplot(ATAC_SpecificModule(), aes(x=Patient, y = ATAC, fill = Patient))+
        geom_boxplot()+
        geom_point()+
        scale_fill_manual(values =PatientColors)+
        labs(title = "across patients", x = NULL, y = "Average ATAC signal")+
        theme_bw(base_size = 14)+
        theme(legend.position = "none")+
        theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    })

    output$ATACsamplePlot <- renderPlot({
        ggplot(ATAC_SpecificModule(), aes(x=Sample, y = ATAC, fill = Patient))+
        geom_bar(stat = "identity", color = "black")+
        scale_fill_manual(values =PatientColors)+
        facet_grid(.~Patient, scales = "free", space = "free")+
        labs(title = "across samples", y = NULL)+
        theme_bw(base_size = 14)+
        theme(strip.background = element_rect(colour="white", fill="white"))+
        theme(strip.text.y = element_text( size = 14, angle = 0))+
        theme(legend.position = "none")+
        theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    })

    output$ATACpurityPlot <- renderPlot({
        ggplot(ATAC_SpecificModule(), aes(x=Purity, y = ATAC, color = Patient, group = "none"))+
        geom_point()+
        labs(title = "by purity", x= "Sample purity (ψ)", y = "Average ATAC signal")+
        theme_bw(base_size = 14)+
        theme(strip.background = element_rect(colour="white", fill="white"))+
        theme(strip.text.y = element_text( size = 14, angle = 0))+
        theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
        geom_smooth(method='lm')+
        stat_cor(method = "pearson", size =5)+
        theme(legend.position = "none")+
        scale_color_manual(values = PatientColors)
    })

    output$ATACdistPlot <- renderPlot({
        ggplot(ATAC_SpecificModule(), aes(x=Dist, y = ATAC, color = Patient, group = "none"))+
        geom_point()+
        labs(title = "by location", x= "Relative distance from centroid (d)", y = "Average ATAC signal")+
        theme_bw(base_size = 14)+
        theme(strip.background = element_rect(colour="white", fill="white"))+
        theme(strip.text.y = element_text( size = 14, angle = 0))+
        theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
        geom_smooth(method='lm')+
        stat_cor(method = "pearson", size =5)+
        theme(legend.position = "none")+
        scale_color_manual(values = PatientColors)
    })

    output$ATACsubtypePlot <- renderPlot({
      ATAC_SpecificModule_subtype <- ATAC_SpecificModule() %>% filter(!is.na(RNAsubtype))
        ggplot(ATAC_SpecificModule_subtype, aes(x=RNAsubtype, y = ATAC, color = Patient, group = RNAsubtype))+
        geom_boxplot()+
        geom_jitter()+
        labs(title = "by subtype", x= "RNAsubtype", y = "Average ATAC signal")+
        theme_bw(base_size = 14)+
        theme(strip.background = element_rect(colour="white", fill="white"))+
        theme(strip.text.y = element_text( size = 14, angle = 0))+
        theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
        scale_color_manual(values = PatientColors)+
        theme(legend.position = "none")+
        stat_compare_means(method = "anova", size =5)
    })

    output$ATACmatrixPlot <- renderPlot({
      ggplot(ATAC_matrix(), aes(x=Distance, y = ValueDifference, color = Patient, group = "none"))+
        geom_point()+
        geom_smooth(method='lm', se = F)+
        stat_cor(method = "pearson")+
        scale_color_manual(values = PatientColors)+
        theme_bw(base_size = 14)+
        labs(title = "distribution within tumor", x = "Inter-sample distance", y = "Difference in average ATAC signal")+
        theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))
    })
    
    output$ATACcelPlot <- renderPlot({
      CEL_plot <-  ATAC_SpecificModule() %>%filter(CEL != "Border")
      ggplot(CEL_plot, aes(x = CEL, y = ATAC, color = RNAsubtype, group = CEL))+
        geom_boxplot(outlier.shape = NA)+
        geom_jitter()+
        theme_bw(base_size = 14)+
        scale_color_manual(values = RNAcolors)+
        labs(title = "contrast-enhancing lesion", x = "Contrast-enhancing",y = "Average ATAC signal" )+
        theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
        stat_compare_means(method = "t.test", label = "p.format")
    })

    ######## explore scATAC data
    
    output$scATAC_reactive <- renderUI({
      req(input$scATAC_select)
      if (input$scATAC_select == "Module"){selectizeInput(inputId = "scATAC_module", label = "Select module: ", choices =scATAC_modulelist, selected = 'R_blue')}
    })
    
    output$sc_ModuleInfo <- renderText({
      req(input$scATAC_module)
      if(input$scATAC_select == "Module"){
      if (input$scATAC_module %in% RNAModules$ModuleID){
      ModuleInfo <- RNAModuleInfo %>%filter(ModuleID == input$scATAC_module)  %>% pull(`Label/Top GSEA result (FDR)`)
      ModuleGeneN <- RNAModuleInfo %>%filter(ModuleID == input$scATAC_module)  %>% pull(`Genes`)
      paste0(ModuleGeneN, " constituent genes; ",ModuleInfo)}
      else ({
        SpecificModuleLabel <-  ATACModuleInfo %>%filter(ModuleID == input$scATAC_module)%>%pull(`Top GSEA result (FDR)`)
        SpecificModulePeaks <-  ATACModuleInfo %>%filter(ModuleID == input$scATAC_module)%>%pull(Peaks)
        SpecificModuleMotif <-  ATACModuleInfo %>%filter(ModuleID == input$scATAC_module)%>%pull(`Top Motif (p-value)`)
        paste0(SpecificModulePeaks, " constituent peaks; enriched for ", SpecificModuleLabel, " and ", SpecificModuleMotif, " motif")})
    }})

    output$exploreATAC <- renderPlot({
      if (input$scATAC_select == "Sample"){
        ggplot(scATAC_overview, aes(x=UMAP_1, y = UMAP_2, color = Sample))+
        geom_point(size = 1, shape = 16, alpha = 0.7)+
        facet_wrap(.~Patient, scales = "free")+
        theme_bw(base_size = 14)+
        scale_color_brewer(palette = "Set1", direction = -1)+
        theme(strip.background = element_rect(colour="white", fill="white"))+
        theme(plot.title = element_text(face="bold"))+
        theme(strip.text.y = element_text( size = 14))+
        labs(title = "Samples")}
      else if (input$scATAC_select == "Neoplastic status"){
        ggplot(scATAC_overview, aes(x=UMAP_1, y = UMAP_2, color = Neoplastic))+
          geom_point(size = 1, shape = 16, alpha = 0.7)+
          facet_wrap(~Patient, scales = "free")+
          theme_bw(base_size = 14)+
          theme(strip.background = element_rect(colour="white", fill="white"))+
          theme(strip.text.y = element_text( size = 14))+
          theme(plot.title = element_text(face="bold"))+
          scale_color_manual(values = c("azure4", "deeppink3"))+
          labs(title = "Neoplastic status")}
      else{
        req(input$scATAC_module)
        scATAC_module <- scATAC_all %>%filter(ModuleID == input$scATAC_module)%>%left_join(scATAC_overview)%>%arrange(Score)
        Type_label <- scATAC_module %>% distinct(Type) %>% pull(Type)
        ggplot(scATAC_module, aes(x=UMAP_1, y = UMAP_2, color = Score))+
        geom_point(size = 1, shape = 16, alpha = 0.7)+
        facet_wrap(.~Patient, scales = "free")+
        labs(title = Type_label)+
        theme_bw(base_size = 14)+
        theme(plot.title = element_text(face="bold"))+
        scale_color_viridis_c(name ="Score", option = "inferno", direction = -1)+
        theme(strip.background = element_rect(colour="white", fill="white"))+
        theme(strip.text.y = element_text( size = 14, face = "italic"))}
    })

    output$sc_barcharts <- renderPlot({
      req(input$scATAC_module)
      if(input$scATAC_select == "Module"){
      scATAC_module <- scATAC_all %>%filter(ModuleID == input$scATAC_module)%>%left_join(scATAC_overview)%>%arrange(Score)
      ggplot(scATAC_module, aes(x= Neoplastic, y = Score, color = Neoplastic))+
        geom_jitter()+
        geom_violin(color = "black", fill = NA)+
        theme_bw(base_size = 11)+
        scale_color_manual(values = c("azure4", "deeppink3"))+
        facet_grid(.~Patient, scales = "free")+
        theme(strip.background = element_rect(colour="white", fill="white"))+
        stat_compare_means(method = "t.test", label = "p.format", vjust = 1)}
    })
    
    ##About page
    output$welcomeimage <- renderImage({list(src = 'www/Welcome.png', height = 600, align="left")}, deleteFile = F)
  
}

shinyApp(ui = ui, server = server)

