##Figure1###

#initialize
library(pacman)
p_load(ggplot2, tidyverse, RColorBrewer, ggpubr, openxlsx, ggrepel, ComplexHeatmap, circlize, janitor)
dataPath <- '~/Dropbox/Postdoc/Papers/ATAC Paper/'
PatientColors = c('P455' = "#ff7500", 'P475' = "#ae7000", "P498" = "#44cef6","P500" = "#1bd1a5", "P503" = "#FFD92F", "P519" = "#8d4bbb","P521" = "#ff0097", "P524" = "#BEBEBE", "P529" = "#0B0B45", "P530" = "#FF0000")
RNAcolors <- c('Classical' = "#3e5063", 'Mesenchymal' = '#0433ff', 'Proneural' = '#ff40ff', 'Neural' = '#ff9300')
RNAcolors_heatmap <- list(Subtype = c('Classical' = "#3e5063", 'Mesenchymal' = '#0433ff', 'Proneural' = '#ff40ff', 'Neural' = '#ff9300'))

#read in files
InputFile <-  read.xlsx(paste0(dataPath, 'Supplement/SupplementaryTableS1.xlsx'), startRow = 2)
SampleInfo <- InputFile%>%dplyr::select(Patient, Sample, Purity, Dist, RNAsubtype, ID)%>%mutate(Sample = as.factor(Sample))
PyCloneClusters <- read.xlsx(paste0(dataPath, 'Supplement/SupplementaryTableS2.xlsx'), sheet = 2)
FACETSResults <- read.xlsx(paste0(dataPath, 'Supplement/SupplementaryTableS2.xlsx'), sheet = 4)
RNA <- read.xlsx(paste0(dataPath, 'Supplement/SupplementaryTableS4.xlsx'), sheet = 1)
Signatures <- read.xlsx(paste0(dataPath, 'Supplement/SupplementaryTableS4.xlsx'), sheet = 2)
RNAModules <- read.xlsx(paste0(dataPath, 'Supplement/SupplementaryTableS4.xlsx'), fillMergedCells = TRUE, sheet = 3) %>%row_to_names(1)%>%mutate(ModuleID = paste0("R_", `RNA module (R_)`), Genes = as.numeric(Genes), PurityR = as.numeric(PurityR), DistR = as.numeric(DistR))%>%relocate(ModuleID)
RNAModuleGenes <- read.xlsx(paste0(dataPath, 'Supplement/SupplementaryTableS4.xlsx'),sheet = 4)
scATAC_geneactivity <- read.xlsx(paste0(dataPath, 'Supplement/SupplementaryTableS5.xlsx'))
ATACModules <- read.xlsx(paste0(dataPath, 'Supplement/SupplementaryTableS6.xlsx'), fillMergedCells = TRUE, sheet = 1) %>%row_to_names(1)%>%mutate(Peaks = as.numeric(Peaks), PurityR = as.numeric(PurityR), DistR = as.numeric(DistR))
scATAC_peaks <- read.xlsx(paste0(dataPath, 'Supplement/SupplementaryTableS7.xlsx'), sheet = 1)
scATAC_amps <- read.xlsx(paste0(dataPath, 'Supplement/SupplementaryTableS7.xlsx'), sheet = 2)
AdultFetalBrain <- read.xlsx(paste0(dataPath, 'Supplement/SupplementaryTableS7.xlsx'), sheet = 3)
PeakCorrelations <- read.xlsx(paste0(dataPath, 'Supplement/SupplementaryTableS6.xlsx'), sheet = 2)
GeneCorrelations <- read.xlsx(paste0(dataPath, 'Supplement/SupplementaryTableS4.xlsx'), sheet = 5)

##Fig1C: Pairwise distances between PyClone-defined clusters for P529 and P530##

pairwiseDistances <- data.frame(patient=character(), sample1 = numeric(), sample2 = numeric(), distance.L=numeric(),matchingpyclone = logical(), stringsAsFactors=F)
pairwiseDist2pointsLPS <- function(p1, p2){#p1 and p2 are both vectors, in which [1] is x, [2] is y, and [3] is z
  return(sqrt( ((p2[1]-p1[1])^2) + ((p2[2]-p1[2])^2) + (((p2[3]-p1[3]))^2) ))
}
for (p in unique(InputFile$Patient)){
  localP <- InputFile[which(InputFile$Patient==p),]
  allSampleCombinations <- combn(localP$Sample,2)
  for (c in 1:ncol(allSampleCombinations)){
    a <- allSampleCombinations[1,c]
    b <- allSampleCombinations[2,c]
    a.coords <- localP[which(localP$Sample==a),c('L','P','S')]
    a.pyclone <- localP[which(localP$Sample==a),'PyClone']
    b.coords <- localP[which(localP$Sample==b),c('L','P','S')]
    b.pyclone <- localP[which(localP$Sample==b),'PyClone']
    distance <- pairwiseDist2pointsLPS(a.coords,b.coords)
    matchingpyclone <- ifelse(a.pyclone == b.pyclone, TRUE, FALSE)
    toBind <- data.frame(c(patient=p, sample1 = a, sample2 = b, distance=distance, matchingpyclone = matchingpyclone))
    pairwiseDistances <- bind_rows(pairwiseDistances, toBind)
  }
}

colnames(pairwiseDistances) <- c('Patient','FirstSample','SecondSample', 'Distance', 'MatchingPyClone')
pyclone_only <- pairwiseDistances %>% filter(!is.na(MatchingPyClone))%>%mutate(MatchingPyClone = factor(MatchingPyClone, levels = c(TRUE, FALSE)))
ggplot(pyclone_only, aes(x=MatchingPyClone, y=Distance, fill = MatchingPyClone))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_brewer(palette="Set3")+
  theme_bw(base_size = 12)+
  facet_grid(Patient~., scales = "free")+
  labs(y="Distance between sample pairs (mm)", x = "Shared subclone")+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.y = element_text( size = 12))+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90))+
  stat_compare_means(method = "t.test", vjust = 1)

## Fig 1H: Plot purity versus CEL
CEL <-InputFile %>% filter(!is.na(CEL))%>% filter(CEL != "Border")
CEL$MRI <- "CEL"
CEL$MRI[CEL$CEL == "No"] <- "NE"
ggplot(CEL, aes(x = CEL, y = Purity, color=Patient, group = CEL))+
  geom_boxplot()+
  geom_jitter()+
  labs(x = "Contrast-enhancing (CE)")+
  theme_bw(base_size = 11)+
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
  scale_color_manual(values = PatientColors)+
  stat_compare_means(method = "t.test")+
  theme(legend.position = "none")+
  ylim(-0.05,1.05)

##Fig 1J: plot purity versus distance
ggplot(InputFile, aes(Dist, Purity, color=Patient, group = "none"))+
  geom_point(position=position_jitter(h=0.05, w=0.05))+
  labs(x = "Distance", y = "Purity")+
  theme_bw(base_size = 11)+
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+ 
  geom_smooth(method='lm', se = FALSE)+
  stat_cor(method = "pearson", size =3)+
  scale_color_manual(values = PatientColors)+
  xlim(0,1.06)+ylim(0,1)

##Fig S1A: PyClone clonal clusters
P530_clusters <- PyCloneClusters %>% filter(Patient == "P530")
ggplot(P530_clusters, aes(x=sample, y=cellular_prevalence, group = Pyclone.cluster))+
  geom_line(aes(color = Pyclone.cluster), lwd=0.7)+
  theme_bw(base_size = 11)+
  labs(title = paste0("Clonal clusters for Patient 530"), y = "cellular prevalence")+
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
  #scale_color_brewer(palette="Set1", direction = -1)+
  facet_wrap(Pyclone.cluster ~., nrow = 3)+
  theme(legend.position = "none", strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.y = element_text( size = 12, angle = 0))+ylim(0,1)

P529_clusters <- PyCloneClusters %>% filter(Patient == "P529")
ggplot(P530_clusters, aes(x=sample, y=cellular_prevalence, group = Pyclone.cluster))+
  geom_line(aes(color = Pyclone.cluster), lwd=0.7)+
  theme_bw(base_size = 11)+
  labs(title = paste0("Clonal clusters for Patient 529"), y = "cellular prevalence")+
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
  #scale_color_brewer(palette="Set1", direction = -1)+
  facet_wrap(Pyclone.cluster ~., nrow = 3)+
  theme(legend.position = "none", strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.y = element_text( size = 12, angle = 0))+ylim(0,1)

##Fig 1F (FACETS copy number for amplification events)
FACETSResults$Gene[FACETSResults$chrom == 7 & FACETSResults$start<55324313 & FACETSResults$end>55086714] <- "EGFR"
FACETSResults$Gene[FACETSResults$chrom == 4 & FACETSResults$start<55164414 & FACETSResults$end>55095264] <- "PDGFRA"
FACETSResults$Gene[FACETSResults$chrom == 1 & FACETSResults$start<204542871 & FACETSResults$end>204485511] <- "MDM4"
FACETSResults$Gene[FACETSResults$chrom == 2 & FACETSResults$start<16087129 & FACETSResults$end>16080683] <- "MYCN"
FACETSResults$Gene[FACETSResults$chrom == 9 & FACETSResults$start<21975132 & FACETSResults$end>21967751] <- "CDKN2A"
FACETSResults$Gene[FACETSResults$chrom == 10 & FACETSResults$start<89728532 & FACETSResults$end>89623195] <- "PTEN"

AmplificationsOnly <- FACETSResults %>% filter(tcn.em >9) %>% select(Patient, Gene) %>% filter(!is.na(Gene)) %>% distinct() %>%left_join(FACETSResults)%>% group_by(ID, Gene)%>%summarize(CN = mean(tcn.em))%>%left_join(InputFile)
ggplot(AmplificationsOnly, aes(x=Patient, y = CN,  fill = Patient))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  scale_fill_manual(values = PatientColors)+
  facet_grid(.~Gene, scales = "free", space = "free")+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.y = element_text( size = 12, angle = 0))+
  theme(legend.position = "none")+
  labs(y = "Copy Number (CN)")+
  geom_hline(yintercept = 2, color = "red",linetype='dotted')

##Fig 1G (EGFR and PDGFRA CN versus gene expression)
GenesOnly <- FACETSResults %>% filter(!is.na(Gene))%>% group_by(ID, Gene)%>%summarize(CN = mean(tcn.em))%>%left_join(InputFile)
RNA_SpecificGenes <- RNA %>% filter(Gene %in% GenesOnly$Gene)%>% pivot_longer(cols = !Gene, names_to = "ID", values_to = "RNA")
GenesOnly_RNA <- GenesOnly %>% left_join(RNA_SpecificGenes)

EGFR_P521 <-  GenesOnly_RNA %>% filter(Gene == "EGFR" & Patient == "P521") %>% mutate(Sample = as.factor(Sample))
ggplot(EGFR_P521, aes(x=CN, y = RNA,  label = Sample, color = Sample, group = "none"))+
  theme_bw()+
  stat_cor(method = "pearson", size =4)+
  geom_smooth(method='lm', se = F, aes(color = "grey", alpha = 0.5))+ 
  geom_jitter()+
  geom_text_repel()+
  scale_color_brewer(palette = "Set1")+
  theme(legend.position = "none")+
  #xlim(0, 25)+ylim(3,10)+
  ylab ("Expression (log2tpm)")+
  ggtitle("EGFR (P521)")

PDGFRA_P521 <-  GenesOnly_RNA %>% filter(Gene == "PDGFRA" & Patient == "P521") %>% mutate(Sample = as.factor(Sample))
ggplot(PDGFRA_P521, aes(x=CN, y = RNA,  label = Sample, color = Sample, group = "none"))+
  theme_bw()+
  stat_cor(method = "pearson", size =4)+
  geom_smooth(method='lm', se = F, aes(color = "grey", alpha = 0.5))+ 
  geom_jitter()+
  geom_text_repel()+
  scale_color_brewer(palette = "Set1")+
  theme(legend.position = "none")+
  xlim(0, 25)+ylim(3,12)+
  ylab ("Expression (log2tpm)")+
  ggtitle("PDGFRA (P521)")

##Fig S2C: CDKN2A and PTEN CN versus expression in P530

PTEN_P530 <-  GenesOnly_RNA %>% filter(Gene == "PTEN" & Patient == "P530") %>% mutate(Sample = as.numeric(Sample)) 
PTEN_P530$Region <- "Temporal"
PTEN_P530$Region[PTEN_P530$Sample>9] <-  "Frontal"
ggplot(PTEN_P530, aes(x=CN, y = RNA,  label = Sample, color = Region, group = "none"))+
  theme_bw()+
  stat_cor(method = "pearson", size =4)+
  geom_point()+
  geom_text_repel()+
  #theme(legend.position = "none")+
  ylab ("Expression (log2tpm)")+
  xlab("Copy number (CN)")+
  ggtitle("PTEN (P530)")+
  scale_color_manual(values = c( "deepskyblue4","deeppink4"))

CDKN2A_P530 <-  GenesOnly_RNA %>% filter(Gene == "CDKN2A" & Patient == "P530") %>% mutate(Sample = as.factor(Sample))
CDKN2A_P530$Region <- "Temporal"
CDKN2A_P530$Region[PTEN_P530$Sample>9] <-  "Frontal"
ggplot(CDKN2A_P530, aes(x=CN, y = RNA,  label = Sample, color = Region, group = "none"))+
  theme_bw()+
  stat_cor(method = "pearson", size =4)+
  geom_point()+
  geom_text_repel()+
  #theme(legend.position = "none")+
  ylab ("Expression (log2tpm)")+
  xlab("Copy number (CN)")+
  ggtitle("CDKN2A (P530)")+
  scale_color_manual(values = c( "deepskyblue4","deeppink4"))

##Fig 2G, 2J, S2C: Gene expression comparisons between P530 frontal and temporal regions
RNA_SpecificGenes <- RNA %>% filter(Gene == "SHLD2"|Gene == "PTEN"|Gene == "RNLS") %>% pivot_longer(cols = !Gene, names_to = "Sample", values_to = "RNA") %>% separate(Sample, c("Patient", "Sample")) %>% mutate(Sample = as.factor(as.numeric(Sample))) %>%mutate(Gene = gsub("SHLD2", "FAM35A", Gene))%>%mutate(Gene = factor(Gene, levels = c("PTEN", "FAM35A", "RNLS")))
RNA_SpecificGenes <- RNA %>% filter(Gene == "KLHL9"|Gene == "MTAP") %>% pivot_longer(cols = !Gene, names_to = "Sample", values_to = "RNA") %>% separate(Sample, c("Patient", "Sample")) %>% mutate(Sample = as.factor(as.numeric(Sample))) %>%mutate(Gene = factor(Gene, levels = c("MTAP", "KLHL9")))
RNA_SpecificGenes <- RNA %>% filter(Gene == "STAG3"|Gene == "CNPY4") %>% pivot_longer(cols = !Gene, names_to = "Sample", values_to = "RNA") %>% separate(Sample, c("Patient", "Sample")) %>% mutate(Sample = as.factor(as.numeric(Sample)))

P530_plot <- RNA_SpecificGenes %>%filter(Patient == "P530")%>%mutate(Sample = as.numeric(Sample), Region = "Frontal")
P530_plot$Region[P530_plot$Sample<10] <- "Temporal"
P530_plot <- P530_plot %>% mutate(Region = factor(Region, levels = c("Temporal", "Frontal")))
ggplot(P530_plot, aes(x=Region, y = RNA, fill = Region))+
  geom_boxplot()+
  geom_jitter()+
  scale_fill_manual(values =c("#C371A5", "#4B91AD"))+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(lineheight=.8, hjust = 0.5))+
  facet_wrap(Gene ~., scales = "free")+
  stat_compare_means(method = "t.test", vjust = -0.5)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  labs(y = "Gene expression (log2tpm)")


##Fig 3A: Verhaak heatmap
SignaturesList <- Signatures%>%pivot_longer(cols = everything(), names_to = "Signature", values_to = "Gene")%>% na.omit()%>%filter(str_detect(Signature, "Verhaak"))%>%left_join(RNA, by = "Gene")%>%na.omit()%>% pivot_longer(cols = !c(Gene, Signature), names_to = "Sample", values_to = "RNA") %>% separate(Sample, c("Patient", "Sample")) %>% 
  group_by(Patient, Sample,Signature) %>% summarise(Expression = mean(RNA))%>% mutate(Sample = as.numeric(Sample)) %>% arrange(Patient, Sample)
SignaturesList$ID <- paste0(SignaturesList$Patient, "_", SignaturesList$Sample)
SignaturesAvgExpression <- SignaturesList %>% pivot_wider(id_cols = Signature, names_from = ID, values_from = Expression) %>% as.data.frame()%>% 
  separate(Signature, into = c("Author", "Journal", "Year", "Label"), extra = "merge")
SampleCols <- SignaturesAvgExpression%>% dplyr::select(starts_with("P4")|starts_with("P5"))%>% as.matrix()
rownames(SampleCols) <- SignaturesAvgExpression$Label

SampleAnnot <- InputFile %>% filter(ID %in% colnames(SampleCols))%>% mutate(Sample = as.factor(Sample))
SampleAnnot$Dist[SampleAnnot$ID == "P521_5"] <- 1.2

sample_annotation = HeatmapAnnotation(Sample = anno_text(SampleAnnot$Sample, rot = 0, just = "top", gp = gpar(fontsize = 8)), Purity = anno_barplot(SampleAnnot$Purity, gp = gpar(fill = SampleAnnot$Color), axis_param = list(side = "right")),annotation_name_side = "left", annotation_name_rot = 0,gap = unit(2, "mm"),
                                       Dist = anno_points(SampleAnnot$Dist, ylim = c(0, 1.2), gp = gpar(col = "black"), axis_param = list(side = "right",at = c(0,1), labels = c("Centroid", "Periphery"))))
VerhaakHeatmap <- Heatmap(SampleCols, name = "Avg exp.", show_column_names = F,show_row_names = T, 
                          clustering_method_rows = 'ward.D2', cluster_rows =F, cluster_columns =F, show_row_dend = FALSE,
                          row_title_rot = 0, row_names_side = "left", column_split = SampleAnnot$Patient, bottom_annotation = sample_annotation,
                          column_names_gp = gpar(fontsize = 8),row_names_gp = gpar(fontsize = 8), cluster_column_slices = F, 
                          col = rev(brewer.pal(n = 11, name = "RdYlBu")))
draw(VerhaakHeatmap, merge_legend = T)


##Fig 3B: Verhaak RNA subtypes by purity and distance
SubtypePlot <- InputFile %>% filter(!is.na(RNAsubtype))
ggplot(SubtypePlot, aes(x = RNAsubtype, y = Purity, fill = RNAsubtype))+
  geom_boxplot()+
  theme_bw(base_size = 11)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab(NULL)+
  theme(legend.position = "none")+
  scale_fill_manual(values = RNAcolors)

ggplot(SubtypePlot, aes(x = RNAsubtype, y = Dist, fill = RNAsubtype))+
  geom_boxplot()+
  theme_bw(base_size = 11)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab(NULL)+
  theme(legend.position = "none")+
  scale_fill_manual(values = RNAcolors)

##Fig3C: RNA module heatmap
ModuleAnnot <- RNAModules 
colors_module = as.character(ModuleAnnot$`RNA module (R_)`)
ModuleRowAnnot <- rowAnnotation(annotation_name_side = "top", 
                                Module = anno_text(ModuleAnnot$ModuleID, gp = gpar(fontsize = 7.5, fontface = "italic")),
                                Genes = anno_barplot(ModuleAnnot$Genes, gp = gpar(fill = colors_module, fontsize = 8), axis_param = list(side = "top")), 
                                Label = anno_text(ModuleAnnot$`Label/Top GSEA result (FDR)`, gp = gpar(fontsize = 7.5)),
                                PurityR = anno_barplot(ModuleAnnot$PurityR, width = unit(1.5, "cm"), gp = gpar(fill = colors_module), axis_param = list(side = "top")),
                                DistR = anno_barplot(ModuleAnnot$DistR, width = unit(1.5, "cm"), gp = gpar(fill = colors_module), axis_param = list(side = "top")))
SampleCols <- ModuleAnnot %>% dplyr::select(starts_with("P4")|starts_with("P5"))%>% mutate_if(is.character,as.numeric)%>% as.matrix()
rownames(SampleCols) <- ModuleAnnot$ModuleID

sample_annotation = HeatmapAnnotation(Sample = anno_text(SampleAnnot$Sample, rot = 0, just = "top", gp = gpar(fontsize = 7)))
sample_annotation2 = HeatmapAnnotation(Purity = anno_barplot(SampleAnnot$Purity, gp = gpar(fill = SampleAnnot$Color), axis_param = list(side = "right")),annotation_name_side = "left", annotation_name_rot = 0,gap = unit(2, "mm"),
                                       Dist = anno_points(SampleAnnot$Dist, ylim = c(0, 1.2), gp = gpar(col = "black"), axis_param = list(side = "right",at = c(0,1), labels = c("Centroid", "Periphery"))),
                                       Subtype = SampleAnnot$RNAsubtype, col = RNAcolors_heatmap)

Heatmap <- Heatmap(SampleCols, name = "Avg Expression",  show_row_names = F, show_column_names = F, right_annotation = ModuleRowAnnot, top_annotation = sample_annotation, bottom_annotation = sample_annotation2, column_split = SampleAnnot$Patient, 
                   clustering_method_columns = 'ward.D2', clustering_method_rows = 'ward.D2',
                   column_names_gp = gpar(fontsize = 8),row_names_gp = gpar(fontsize = 8), cluster_rows = T, cluster_column_slices = F, cluster_columns = T, col = rev(brewer.pal(n = 11, name = "RdYlBu")))
draw(Heatmap, merge_legend = T, heatmap_legend_side = "bottom")

##Fig 3D/E - R_brown (oligodendrocyte) and R_blue (microglia) modules by purity 
RNASeq_All_Pivot <- RNA %>%pivot_longer(cols = -Gene, names_to = "ID", values_to = "RNA")
RNASeq_All_Ordered <- InputFile %>% inner_join(RNASeq_All_Pivot)%>%select(ID, Patient, Sample, Purity, Dist, RNAsubtype, Gene, RNA)
RNAModuleGenesBySample <- RNAModuleGenes%>% select(Module, Gene) %>% inner_join(RNASeq_All_Ordered)
TurquoiseBrown <- RNAModuleGenesBySample  %>%filter(Module == "turquoise" |Module == "brown") %>%group_by(Sample, Module, Purity, RNAsubtype) %>%summarize(RNA = mean(RNA))%>%mutate(ModuleID = paste0("R_", Module))
ggplot(TurquoiseBrown, aes(x = Purity, y = RNA, color = RNAsubtype, group = "none"))+
  geom_point()+
  facet_grid(. ~ModuleID, scales = "free", space = "free")+
  scale_color_manual(values = RNAcolors)+
  theme_bw(base_size = 11)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.x = element_text(face = "italic"))+
  theme(legend.position = "bottom")+
  geom_smooth(method='lm', se = F)+
  stat_cor(method = "pearson", size =3)+
  labs(x = "Sample purity", y = "Average expression")

##Fig4G - Modules with interpatient heterogeneity between P475 and remaining patients
InterpatientModules <- RNAModuleGenesBySample %>%filter(Module == "ivory" |Module == "maroon"|Module == "brown2") %>%group_by(Patient, Sample, Module, Dist, Purity, RNAsubtype)%>%summarize(RNA = mean(RNA)) %>%mutate(ModuleID = paste0("R_", Module))%>%mutate(ModuleID = factor(ModuleID, levels = c("R_brown2", "R_maroon","R_ivory")))%>%mutate(Group = "Other")
ggplot(InterpatientModules, aes(x = Patient, y = RNA, fill = Patient))+
  geom_boxplot()+
  facet_wrap(.~ModuleID, scales = "free")+
  scale_fill_manual(values = PatientColors)+
  theme_bw(base_size = 11)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.x = element_text(face = "italic"))+
  theme(legend.position = "none")+
  labs(x = "Patient", y = "Average expression")+
  theme(axis.text.x = element_text(angle = 90))

##Fig S4C - Expression versus purity for modules in P475 versus remaining patients
PurityFigs <- RNAModuleGenesBySample %>%filter(Module == "ivory" |Module == "maroon"|Module == "brown2"|Module == "turquoise") %>%group_by(Patient, Sample, Module, Dist, Purity, RNAsubtype)%>%summarize(RNA = mean(RNA)) %>%mutate(ModuleID = paste0("R_", Module))%>%mutate(ModuleID = factor(ModuleID, levels = c("R_turquoise", "R_brown2", "R_maroon","R_ivory")))%>%mutate(Group = "Other")
PurityFigs$Group[PurityFigs$Patient == "P475"] <- "P475"
ggplot(PurityFigs, aes(x = Purity, y = RNA, color = Patient, group = Group))+
  geom_point()+
  facet_grid(.~ModuleID, scales = "free")+
  scale_color_manual(values = PatientColors)+
  theme_bw(base_size = 11)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.x = element_text(face = "italic"))+
  theme(legend.position = "bottom")+
  geom_smooth(method='lm', se = F)+
  stat_cor(method = "pearson", size =3)+
  labs(x = "Sample purity", y = "Average expression")

##Fig 4H/I Expression versus purity for individual genes showing differences between P475 and remaining patients
SpecificGene <- "PTPRZ1"
SpecificGene <- "PTN"
SpecificGene <- "ETV1"
SpecificGene <- "NKX2-1"
RNA_SpecificGene <- RNA %>% filter(Gene == SpecificGene) %>% pivot_longer(cols = !Gene, names_to = "Sample", values_to = "RNA") %>% separate(Sample, c("Patient", "Sample")) %>% mutate(Sample = as.factor(as.numeric(Sample))) %>%left_join(SampleInfo)

P475vRest <-RNA_SpecificGene %>%mutate(Group = "Other")
P475vRest$Group[P475vRest$Patient == "P475"] <- "P475"
ggplot(P475vRest, aes(x=Purity, y = RNA, color = Patient, group = Group))+
  geom_point()+
  labs(title = SpecificGene, x= "Sample purity", y = "Gene Expression")+
  theme_bw(base_size = 11)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.y = element_text( size = 12, angle = 0))+
  theme(plot.title = element_text(lineheight=.8, face="italic",hjust = 0.5))+
  geom_smooth(method='lm', se = F)+
  stat_cor(method = "pearson", size =3)+
  scale_color_manual(values = PatientColors)


##Fig 5A & Fig 5G- RNA modules by RNA subtype
RNAModuleExpressionBySample <- RNAModuleGenesBySample %>% group_by(Module, ID, Patient, Sample, RNAsubtype) %>% summarize(RNA = mean(RNA))
Verhaak <- RNAModuleExpressionBySample %>%filter(Module == "orangered3" |Module == "brown"|Module == "plum")%>%mutate(ModuleID = paste0("R_", Module)) ##Fig5A
Verhaak <- RNAModuleExpressionBySample %>%filter(Module == "plum2" |Module == "midnightblue"|Module == "plum3"|Module == "darkred")%>%mutate(ModuleID = paste0("R_", Module)) %>%mutate(ModuleID = factor(ModuleID, levels = c("R_midnightblue", "R_plum2", "R_plum3", "R_darkred"))) ##Fig 5G
ggplot(Verhaak, aes(x = RNAsubtype, y = RNA, fill = RNAsubtype))+
  geom_boxplot()+
  facet_wrap(.~ModuleID, scales = "free")+
  scale_fill_manual(values = RNAcolors)+
  theme_bw(base_size = 11)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.x = element_text(size = 11, face = "italic"))+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(legend.position = "bottom")+
  labs(y = "Gene expression",  x= NULL)+ theme(legend.title = element_blank())


##Fig 5B - Neural samples by RNA module
NeuralList <- c("P455_1", "P455_2", "P455_10", "P503_1","P503_2", "P524_1", "P524_7")
Neural <- RNAModuleGenesBySample  %>%filter(Module == "orangered3" |Module == "brown"|Module == "plum")%>%filter(ID %in% NeuralList)%>%mutate(Sample = as.factor(Sample), ModuleID = paste0("R_", Module))
Neural %>% select(ID) %>% distinct()
ModuleFill <- as.list(Neural$Module)
names(ModuleFill) <- as.list(Neural$ModuleID)
ggplot(Neural, aes(x = Sample, y = RNA, fill = ModuleID))+
  geom_boxplot()+
  facet_grid(.~Patient, scales = "free", space = "free")+
  scale_fill_manual(values = ModuleFill)+
  theme_bw(base_size = 11)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  ylab("Gene expression")+ theme(legend.title = element_blank(), legend.text = element_text(face = "italic"))


#Fig 5C - RNA modules enriched towards the periphery 
Periph <- RNAModuleGenesBySample  %>%filter(Module == "darkseagreen4" |Module == "lightcoral"|Module == "pink") %>%group_by(Sample, Module, Dist, RNAsubtype) %>%summarize(RNA = mean(RNA))%>%mutate(ModuleID = paste0("R_", Module))
ggplot(Periph , aes(x = Dist, y = RNA, color = RNAsubtype, group = "none"))+
  geom_point()+
  facet_wrap(. ~ModuleID, scales = "free")+
  scale_color_manual(values = RNAcolors)+
  theme_bw(base_size = 10)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.x = element_text(face = "italic"))+
  theme(legend.position = "bottom")+ theme(legend.title = element_blank())+
  geom_smooth(method='lm', se = F)+
  stat_cor(method = "pearson", size =3)+
  labs(x = "Distance from centroid", y = "Average expression")

##Fig 5H - Correlation matrix for microenvironmental RNA modules 
Microenvironment <- ModuleAnnot %>% filter(`RNA module (R_)` == "blue" | `RNA module (R_)` == "midnightblue" | `RNA module (R_)` == "greenyellow"| `RNA module (R_)` == "midnightblue"| `RNA module (R_)` == "plum2"| `RNA module (R_)` == "plum3"| `RNA module (R_)` == "darkred"| `RNA module (R_)` == "brown"| `RNA module (R_)` == "plum"| `RNA module (R_)` == "orangered3")
MicroenvironmentCols <- Microenvironment %>% dplyr::select(starts_with("P4")|starts_with("P5"))%>% mutate_if(is.character,as.numeric)%>% as.matrix()
rownames(MicroenvironmentCols) <- Microenvironment$ModuleID
cor_matrix <- round(cor(t(MicroenvironmentCols)), 3) 
ModuleMatrixHeatmap <- Heatmap(cor_matrix, name = "Correlation", column_names_gp = gpar(fontsize = 10, fontface = "italic"),row_names_gp = gpar(fontsize = 10,fontface = "italic"), col = colorRamp2(c(-1,0,1), c('blue', 'white', 'red')))
draw(ModuleMatrixHeatmap)

##Fig 5I - scATAC cell distributions by sample and neoplastic status
cell_annotations_count <- scATAC_geneactivity %>% group_by(Patient, Sample) %>% summarize(Count = n())
ggplot(cell_annotations_count, aes(x=Patient, y = Count, fill = Sample))+
  geom_col(colour = "black")+
  geom_text(aes(label = Sample),position = "stack", vjust = 1.5, size = 3.5)+
  theme_bw(base_size = 10)+
  scale_fill_brewer(palette = "Set1", direction = -1)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(legend.position = "none")+
  ylab("Number of cells")

ggplot(scATAC_geneactivity, aes(x=UMAP_1, y = UMAP_2, color = Sample))+
  geom_point(size = 1, shape = 16, alpha = 0.7)+
  facet_grid(.~Patient, scales = "free")+
  theme_bw(base_size = 11)+
  scale_color_brewer(palette = "Set1", direction = -1)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.y = element_text( size = 12))+
  theme(legend.position = "none")+
  labs(title = "Samples")

ggplot(scATAC_geneactivity, aes(x=Patient, fill = Neoplastic))+
  geom_bar(colour = "black")+
  theme_bw(base_size = 11)+
  #scale_fill_brewer(palette = "Set1", direction = -1)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(legend.position = "bottom")+
  scale_fill_manual(values = c("azure4", "deeppink3"))+
  ylab("Number of cells")

ggplot(scATAC_geneactivity, aes(x=UMAP_1, y = UMAP_2, color = Neoplastic))+
  geom_point(size = 1, shape = 16, alpha = 0.7)+
  facet_grid(~Patient, scales = "free")+
  theme_bw(base_size = 11)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.y = element_text( size = 12))+
  scale_color_manual(values = c("azure4", "deeppink3"))+
  theme(legend.position = "none")+
  labs(title = "Neoplastic")

#Fig S5D: Neoplastic and non-neoplastic cell count by sample
ggplot(scATAC_geneactivity, aes(x=Sample, fill = Neoplastic))+
  geom_bar(colour = "black")+
  theme_bw(base_size = 11)+
  #scale_fill_brewer(palette = "Set1", direction = -1)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(legend.position = "bottom")+
  facet_grid(.~Patient, scales = "free", space = "free")+
  scale_fill_manual(values = c("azure4", "deeppink3"))+
  ylab("Number of cells")

#Fig S5E: Neoplastic % in scATAC versus purity from bulk WES
neoplastic_percent <- scATAC_geneactivity %>%group_by(ID) %>% summarize(Neoplastic = mean(Neoplastic))%>%left_join(SampleInfo)
ggplot(neoplastic_percent, aes(x=Neoplastic, y = Purity, group = "none", label = Sample))+
  geom_smooth(method='lm', se = F)+
  geom_point()+
  theme_bw(base_size = 11)+
  xlim(0,1)+ylim(0,1)+
  xlab("% neoplastic (scATAC)")+ ylab("Purity (bulk exome)")+
  stat_cor(method = "pearson", size =3)+
  theme(legend.position = "bottom")

# Fig 5J/K - scATAC gene activity scores for RNA modules
scATAC_geneactivity_longer <- scATAC_geneactivity%>% pivot_longer(cols = 8:ncol(scATAC_geneactivity), names_to = "ModuleID", values_to = "Score") %>%arrange(Score)
scATAC_gene_select <- scATAC_geneactivity_longer %>%filter(ModuleID == "R_blue" | ModuleID == "R_midnightblue"| ModuleID == "R_brown"|ModuleID == "R_maroon"| ModuleID == "R_darkred"| ModuleID == "R_plum2"| ModuleID == "R_plum3")%>%arrange(Score)
scATAC_gene_select$ModuleID <- factor(scATAC_gene_select$ModuleID, levels = c('R_maroon', 'R_brown', 'R_blue', 'R_midnightblue', 'R_plum2', 'R_plum3', 'R_darkred'))

ggplot(scATAC_gene_select, aes(x=UMAP_1, y = UMAP_2, color = Score))+
  geom_point(size = 1, shape = 16, alpha = 0.7)+
  facet_grid(ModuleID~Patient, scales = "free")+
  theme_bw(base_size = 11)+
  scale_color_viridis_c(option = "inferno", direction = -1, limits = c(0, 0.003))+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.y = element_text( size = 12, face = "italic"))+
  #theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
  #theme(legend.position = "none")+
  labs(title = "Gene activity score")

ggplot(scATAC_gene_select, aes(x= Neoplastic, y = Score, color = Neoplastic))+
  geom_jitter()+
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA)+
  theme_bw(base_size = 11)+
  scale_color_manual(values = c("azure4", "deeppink3"))+
  facet_grid(.~ModuleID, scales = "free", space = "free")+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.x = element_text( size = 12, angle = 0, face = "italic"))+
  theme(axis.text.x = element_text(angle = 90))+
  theme(legend.position = "none")+
  stat_compare_means(method = "t.test", label.y = 0.004)

#Fig S5F - Correlation matrix for scATAC immune/mesenchymal modules
matrix_geneactivity <- scATAC_geneactivity %>% select(R_blue, R_midnightblue, R_greenyellow, R_darkred, R_plum2, R_plum3)%>%as.matrix()
cor_matrix <- round(cor(matrix_geneactivity), 3) 
MatrixHeatmap <- Heatmap(cor_matrix, name = "Correlation", clustering_method_columns = 'ward.D2', clustering_method_rows = 'ward.D2', column_names_gp = gpar(fontsize = 8, fontface = "italic"),row_names_gp = gpar(fontsize = 8, fontface = "italic"), col = rev(brewer.pal(n = 11, name = "RdYlBu")))
draw(MatrixHeatmap)

##Fig6C - ATAC & linkage module heatmap
SampleCols <- ATACModules%>% dplyr::select(starts_with("P4")|starts_with("P5"))%>% mutate_if(is.character,as.numeric)%>% as.matrix()
SampleAnnot <- InputFile %>% filter(ID %in% colnames(SampleCols))%>% mutate(Sample = as.factor(Sample))
SampleAnnot$Dist[SampleAnnot$ID == "P521_5"] <- 1.2
sample_annotation_top = HeatmapAnnotation(Sample = anno_text(SampleAnnot$Sample, rot = 0, just = "top", gp = gpar(fontsize = 8)))
sample_annotation_bottom = HeatmapAnnotation(Purity = anno_barplot(SampleAnnot$Purity, gp = gpar(fill = SampleAnnot$Color), axis_param = list(side = "right")),annotation_name_side = "left", annotation_name_rot = 0,gap = unit(2, "mm"),
                                             Dist = anno_points(SampleAnnot$Dist, ylim = c(0, 1.2), gp = gpar(col = "black"), axis_param = list(side = "right",at = c(0,1), labels = c("Centroid", "Periphery"))),
                                             Subtype = SampleAnnot$RNAsubtype, col = RNAcolors_heatmap)
ModuleAnnot <- ATACModules
colors_module = as.character(ModuleAnnot$Module)
ModuleRowAnnot <- rowAnnotation(annotation_name_side = "top",
                                Module = anno_text(ModuleAnnot$ModuleID,gp = gpar(fontsize = 7.5, fontface = "italic")),
                                Peaks = anno_barplot(ModuleAnnot$Peaks, width = unit(1.6, "cm"),gp = gpar(fill = colors_module, fontsize = 8), axis_param = list(side = "top")), 
                                Motif = anno_text(ModuleAnnot$`Top Motif (p-value)`, gp = gpar(fontsize = 7.5)),
                                TopGene = anno_text(ModuleAnnot$`Top Linked Gene (# peaks)`, gp = gpar(fontsize = 7.5)),
                                PurityR = anno_barplot(ModuleAnnot$PurityR, width = unit(1.5, "cm"), gp = gpar(fill = colors_module), axis_param = list(side = "top")),
                                  DistR = anno_barplot(ModuleAnnot$DistR, width = unit(1.5, "cm"), gp = gpar(fill = colors_module), axis_param = list(side = "top")))


PeaksBySampleHeatmap <- Heatmap(SampleCols, name = "Average ATAC signal",  bottom_annotation = sample_annotation_bottom, cluster_columns = T, cluster_column_slices = F,
                                cluster_row_slices = F, right_annotation = ModuleRowAnnot, show_row_names = F,top_annotation = sample_annotation_top,
                                clustering_method_columns = 'ward.D2', clustering_method_rows = 'ward.D2',
                                column_split = SampleAnnot$Patient, 
                                row_split = ModuleAnnot$ModuleSet,row_title = c("ATAC Modules", "Linkage Modules"),
                                show_column_names = F,
                                column_names_gp = gpar(fontsize = 8),row_names_gp = gpar(fontsize = 8), 
                                col = rev(brewer.pal(n = 11, name = "RdYlBu")))
draw(PeaksBySampleHeatmap, heatmap_legend_side = "bottom", merge_legend = T)

## Fig6F - scATAC read-in-peak UMAPs for modules associated with microenvironment
scATAC_P521 <-  scATAC_peaks %>%filter(Patient == "P521") %>% pivot_longer(8:ncol(scATAC_peaks), names_to = "Module", values_to = "Score")%>%filter(Module == "A_thistle"|Module == "A_lavenderblush2"|Module == "A_magenta4")%>%arrange(Score)
ggplot(scATAC_P521, aes(x=UMAP_1, y = UMAP_2, color = Score))+
  geom_point(size = 1, shape = 16, alpha = 0.7)+
  facet_grid(.~Module, scales = "free")+
  theme_bw(base_size = 11)+
  scale_color_viridis_c(option = "inferno", direction = -1)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.x = element_text( size = 12, face = "italic"))

## Fig6F - Comparing ATAC modules to adult and fetal brain published scATAC datasets
AdultFetalBrain_longer <- AdultFetalBrain %>% pivot_longer(3:ncol(AdultFetalBrain), names_to = "Module", values_to = "value")
AdultOnly <- AdultFetalBrain_longer %>% filter(Module == "A_thistle1"|Module == "A_honeydew"|Module == "L_mediumpurple1"|Module == "A_lavenderblush2"|Module == "A_coral"|Module == "L_coral"|Module == "A_magenta4"|Module == "A_skyblue4"|Module == "L_plum2"|Module == "A_plum3"|Module == "A_yellowgreen")%>%arrange(desc(value))
topcells <- AdultOnly%>%filter(Set != "Fetal") %>%filter(Cell == "Oligodendrocyte"|Cell == "Microglia"|Cell == "Glutaminergic Neuron 2"|Cell == "Oligodendrocyte Precursor"|Cell == "Astrocyte 1")%>%pivot_wider(id_cols = c(Set, Cell), names_from = Module, values_from = value)%>%arrange(Cell, Set)%>%mutate(Cell = gsub("1", "", Cell), Cell = gsub("2", "", Cell), Cell = gsub("Glutaminergic", "", Cell), Cell = gsub("Oligodendrocyte Precursor", "OPC", Cell))
ColAnnot <- HeatmapAnnotation(Set =topcells$Set, Cell = anno_text(topcells$Cell))
heatmap <- Heatmap(t(as.matrix(topcells[,3:ncol(topcells)])), bottom_annotation = ColAnnot, cluster_columns = F,row_names_gp = gpar(fontface = "italic"), name = "Value",
                   col = rev(brewer.pal(n = 11, name = "RdYlBu")))
draw(heatmap)

FetalOnly <- AdultFetalBrain_longer %>% filter(Module == "A_coral"|Module == "L_coral"|Module == "A_lavenderblush2"|Module == "L_navajowhite1"|Module == "A_darkseagreen3"|Module == "A_thistle"|Module == "A_plum"|Module == "L_darkseagreen4"|Module == "A_indianred4"|Module == "A_maroon"|Module == "L_coral2")%>%arrange(desc(value))
topcells <- FetalOnly%>%filter(Set == "Fetal") %>%filter(Cell != "Microglia")%>%pivot_wider(id_cols = c(Set, Cell), names_from = Module, values_from = value)%>%arrange(Cell, Set)
ColAnnot <- HeatmapAnnotation(Set =topcells$Set, Cell = anno_text(topcells$Cell))
heatmap <- Heatmap(t(as.matrix(topcells[,3:ncol(topcells)])), bottom_annotation = ColAnnot,row_names_gp = gpar(fontface = "italic"), name = "Value",                  col = rev(brewer.pal(n = 11, name = "RdYlBu")))
draw(heatmap)


##Fig6I - NEUROD1 gene expression by patient
SpecificGene = "NEUROD1"
RNA_SpecificGene <- RNA %>% filter(Gene == SpecificGene) %>% pivot_longer(cols = !Gene, names_to = "Sample", values_to = "RNA") %>% separate(Sample, c("Patient", "Sample")) %>% mutate(Sample = as.factor(as.numeric(Sample))) 
ggplot(RNA_SpecificGene, aes(x=Patient, y = RNA, fill = Patient))+
  geom_boxplot()+
  geom_point()+
  scale_fill_manual(values =PatientColors)+
  labs(title = paste0("Gene expression: ", SpecificGene), x = NULL, y = "NEUROD1 gene expression")+
  theme_bw(base_size = 11)+
  theme(legend.position = "none")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

##Fig6K - SOX9 versus SOX10 gene expression in P521
SOX9 <- RNA %>% filter(Gene == "SOX9") %>% pivot_longer(cols = !Gene, names_to = "Sample", values_to = "SOX9") %>% separate(Sample, c("Patient", "Sample")) %>% mutate(Sample = as.factor(as.numeric(Sample))) %>%select(-Gene)
SOX10 <- RNA %>% filter(Gene == "SOX10") %>% pivot_longer(cols = !Gene, names_to = "Sample", values_to = "SOX10") %>% separate(Sample, c("Patient", "Sample")) %>% mutate(Sample = as.factor(as.numeric(Sample))) %>%select(-Gene)
SOX9_10 <- left_join(SOX9, SOX10)%>%left_join(SampleInfo)%>%filter(Patient != "P475")
SOX_P521 <- SOX9_10 %>%filter(Patient == "P521")%>%pivot_longer(cols = c("SOX9", "SOX10"), names_to = "Gene", values_to = "RNA")%>%mutate(Gene = factor(Gene, levels = c("SOX9", "SOX10")))
ggplot(SOX_P521, aes(x=Gene, y = RNA, label = Sample, color = RNAsubtype))+
  theme_bw(base_size = 11)+
  geom_point()+
  geom_label_repel(size = 3)+
  labs(x=NULL,y="Gene expression")+
  scale_color_manual(values = RNAcolors)

## Fig6L - scATAC amplification calls for EGFR versus PDGFRA
P521_amps <- scATAC_amps %>%filter(Patient == "P521")%>%left_join(scATAC_peaks)
P521_amps$Amp[P521_amps$egfr == 1] <- "EGFR"
P521_amps$Amp[P521_amps$pdgfr == 1] <- "PDGFRA"
ggplot(P521_amps, aes(x=UMAP_1, y = UMAP_2, color = Amp))+
  geom_point(size = 1, shape = 16, alpha = 0.7)+
  facet_grid(~Patient, scales = "free")+
  theme_bw(base_size = 11)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.y = element_text( size = 12))+
  scale_color_manual(values = c("EGFR" = "#BBAED1", "PDGFRA" = "#90C786"))

## Fig6M - scATAC read-in-peak UMAPs for modules with intratumoral heterogeneity in P521
scATAC_P521 <-  scATAC_peaks %>%filter(Patient == "P521") %>% pivot_longer(8:ncol(scATAC_peaks), names_to = "Module", values_to = "Score")%>%filter(Module == "A_plum3"|Module == "A_yellowgreen")%>%arrange(Score)
ggplot(scATAC_P521, aes(x=UMAP_1, y = UMAP_2, color = Score))+
  geom_point(size = 1, shape = 16, alpha = 0.7)+
  facet_grid(.~Module, scales = "free")+
  theme_bw(base_size = 11)+
  scale_color_viridis_c(option = "inferno", direction = -1)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.x = element_text( size = 12, face = "italic"))

##Fig 6N: scATAC read-in-peak score by EGFR versus PDGFRA amplification in P521
P521_plot <- P521_amps%>%filter(!is.na(Amp))%>%select(Amp, A_plum3, A_yellowgreen, A_coral) %>%pivot_longer(cols = -Amp, names_to = "Module", values_to = "Score")%>%mutate(Module = factor(Module, levels = c("A_plum3", "A_yellowgreen", "A_coral")))
ggplot(P521_plot, aes(x = Amp, y = Score, fill = Amp))+
  geom_boxplot(outlier.shape = NA)+
  facet_grid(.~Module)+
  theme_bw(base_size = 11)+
  stat_compare_means(method = "t.test", size = 3, label = "p.format")+
  scale_fill_manual(values = c("EGFR" = "#BBAED1", "PDGFRA" = "#90C786"))+
  theme(axis.text.x = element_text(angle = 90))+
  theme(strip.text.x = element_text(face = "italic", angle = 0))+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  labs(x = "Amplification")+ylim(0,4)

##Fig 7A - Association of L_salmon module with EGFR amplification
ATACModuleFigs <- ATACModules %>% dplyr::select(ModuleID|starts_with("P4")|starts_with("P5"))%>% pivot_longer(cols = -ModuleID, names_to = "ID", values_to ="ATACsignal")%>%mutate(ATACsignal = as.numeric(ATACsignal))%>%left_join(SampleInfo)
salmon4 <- ATACModuleFigs %>%filter(ModuleID == "L_salmon4")
EGFR_CN <-  FACETSResults %>%filter(Gene == "EGFR")%>%select(Patient, ID, tcn.em)%>%rename(EGFR_CN = tcn.em)
EGFR_salmon4 <-   salmon4 %>% left_join(EGFR_CN)
ggplot(EGFR_salmon4, aes(x = EGFR_CN, y = ATACsignal, color = Patient, group = "none"))+
  geom_jitter()+
  theme_bw(base_size = 11)+
  scale_color_manual(values = PatientColors)+
  labs(x = "EGFR total copy number", y = "L_salmon4 average ATAC signal")+ 
  geom_smooth(method='lm', se = F)+
  stat_cor(method = "pearson", size =3)

##Fig 7B
EGFR <- scATAC_peaks%>%filter(Patient == "P521"|Patient == "P529") %>% select(Cell, UMAP_1, UMAP_2, Neoplastic, L_salmon4)%>%left_join(scATAC_amps) %>%mutate(egfr = as.logical(egfr))%>%filter(Neoplastic == T)

ggplot(EGFR, aes(x = egfr, y = L_salmon4, fill = egfr))+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(Patient ~ .)+  theme_bw(base_size = 11)+
  stat_compare_means(method = "t.test", size = 3, label = "p.format")+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  stat_n_text(size = 3)+
  labs(x= "EGFR amplification", y = "L_salmon4 scATAC score")+
  scale_fill_brewer(palette = "Set1")

##Fig 7E = ATAC modules by distance from centroid
DistFig<- ATACModuleFigs  %>%filter(ModuleID == "L_mediumorchid" |ModuleID == "A_lavenderblush3" |ModuleID == "L_lightcyan1"|ModuleID == "L_thistle"|ModuleID == "L_orangered3") %>%mutate(ModuleID = factor(ModuleID, levels = c("L_mediumorchid", "A_lavenderblush3", "L_lightcyan1", "L_orangered3", "L_thistle")))
ggplot(DistFig , aes(x = Dist, y = ATACsignal, color = RNAsubtype, group = "none"))+
  geom_point()+
  facet_wrap(. ~ModuleID, scales = "free", nrow = 1)+
  scale_color_manual(values = RNAcolors)+
  theme_bw(base_size = 10)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.x = element_text(face = "italic"))+
  theme(legend.position = "bottom")+ theme(legend.title = element_blank())+
  geom_smooth(method='lm', se = F)+
  stat_cor(method = "pearson", size =3)+
  labs(x = "Distance from centroid", y = "Average ATAC signal")

##Fig7H -- Individual peaks and genes correlation with purity and distance

PeakCorrelations_fig_AP1 <- PeakCorrelations %>% filter(ATACModule == "lavenderblush3"|LinkageModule == "lightcyan1"|LinkageModule == "orangered3"|LinkageModule == "thistle")%>%distinct(Peak, PurityR, DistR)%>%mutate(FigLabel = "AP1/Mesenchymal")
PeakCorrelations_fig_mediumorchid <- PeakCorrelations %>% filter(LinkageModule == "mediumorchid")%>%distinct(Peak, PurityR, DistR)%>%mutate(FigLabel = "Neuronal hijacking")
PeakCorrelations_fig_salmon4 <- PeakCorrelations %>% filter(LinkageModule == "salmon4")%>%distinct(Peak, PurityR, DistR)%>%mutate(FigLabel = "EGFR amplification")
PeakCorrelations_fig_navajowhite1 <- PeakCorrelations %>% filter(LinkageModule == "navajowhite1"|ATACModule == "plum3"|ATACModule == "yellowgreen"|ATACModule == "darkseagreen3"|ATACModule == "thistle"|ATACModule == "plum")%>%distinct(Peak, PurityR, DistR)%>%mutate(FigLabel = "Neurodevelopmental")
PeakCorrelations_fig_oligodendrocyte <- PeakCorrelations %>% filter(ATACModule == "lavenderblush2"|ATACModule == "coral"|LinkageModule == "coral")%>%distinct(Peak, PurityR, DistR)%>%mutate(FigLabel = "Oligodendrocyte")
PeakCorrelations_fig_neuron <- PeakCorrelations %>% filter(ATACModule == "thistle1"|ATACModule == "honeydew"|LinkageModule == "mediumpurple1")%>%distinct(Peak, PurityR, DistR)%>%mutate(FigLabel = "Neuron")
PeakCorrelations_fig_microglia <- PeakCorrelations %>% filter(ATACModule == "magenta4"|ATACModule == "skyblue4"|LinkageModule == "plum2")%>%distinct(Peak, PurityR, DistR)%>%mutate(FigLabel = "Microglia")

Correlations_summary <- ATACModules
Correlations_summary$FigLabel[Correlations_summary$ModuleID == "A_lavenderblush3" | Correlations_summary$ModuleID == "L_lightcyan1"| Correlations_summary$ModuleID == "L_orangered3" | Correlations_summary$ModuleID == "L_thistle"] <- "AP1/Mesenchymal"
Correlations_summary$FigLabel[Correlations_summary$ModuleID == "L_mediumorchid"] <- "Neuronal hijacking"
Correlations_summary$FigLabel[Correlations_summary$ModuleID == "L_navajowhite1"|Correlations_summary$ModuleID == "A_plum3"|Correlations_summary$ModuleID == "A_yellowgreen"|Correlations_summary$ModuleID == "A_darkseagreen3" | Correlations_summary$ModuleID == "A_thistle"| Correlations_summary$ModuleID == "A_plum"] <- "Neurodevelopmental"
Correlations_summary$FigLabel[Correlations_summary$ModuleID == "L_salmon4"] <- "EGFR amplification"
Correlations_summary$FigLabel[Correlations_summary$ModuleID == "A_lavenderblush2" | Correlations_summary$ModuleID == "A_coral"| Correlations_summary$ModuleID == "L_coral"] <- "Oligodendrocyte"
Correlations_summary$FigLabel[Correlations_summary$ModuleID == "A_thistle1" | Correlations_summary$ModuleID == "A_honeydew"| Correlations_summary$ModuleID == "L_mediumpurple1"] <- "Neuron"
Correlations_summary$FigLabel[Correlations_summary$ModuleID == "A_magenta4" | Correlations_summary$ModuleID == "A_skyblue4"| Correlations_summary$ModuleID == "L_plum2"] <- "Microglia"
Correlations_summary <- Correlations_summary %>% filter(!is.na(FigLabel))%>%relocate(FigLabel)

Number_peaks <- Correlations_summary %>% group_by(FigLabel)%>% summarize(Peaks = sum(Peaks))
PeakCorrelations_fig <- bind_rows(PeakCorrelations_fig_AP1, PeakCorrelations_fig_mediumorchid, PeakCorrelations_fig_salmon4, PeakCorrelations_fig_navajowhite1, PeakCorrelations_fig_oligodendrocyte, PeakCorrelations_fig_neuron, PeakCorrelations_fig_microglia)%>%left_join(Number_peaks)%>%arrange(desc(Peaks))

GeneCorrelations$FigLabel[GeneCorrelations$RNAModule == "brown"] <- "Oligodendrocyte"
GeneCorrelations$FigLabel[GeneCorrelations$RNAModule == "plum"] <- "Astrocyte"
GeneCorrelations$FigLabel[GeneCorrelations$RNAModule == "orangered3"] <- "Neuron"
GeneCorrelations$FigLabel[GeneCorrelations$RNAModule == "blue"] <- "Microglia"
GeneCorrelations$FigLabel[GeneCorrelations$RNAModule == "greenyellow"] <- "Interferon response"
GeneCorrelations$FigLabel[GeneCorrelations$RNAModule == "midnightblue"|GeneCorrelations$RNAModule == "plum2"|GeneCorrelations$RNAModule == "plum3"|GeneCorrelations$RNAModule == "darkred"] <- "AP1/Mesenchymal"
GeneCorrelations$FigLabel[GeneCorrelations$RNAModule == "pink"|GeneCorrelations$RNAModule == "darkseagreen4"|GeneCorrelations$RNAModule == "lightcoral"] <- "Neuronal hijacking"
GeneCorrelations$FigLabel[GeneCorrelations$RNAModule == "maroon"|GeneCorrelations$RNAModule == "turquoise"|GeneCorrelations$RNAModule == "brown2"|GeneCorrelations$RNAModule == "ivory"] <- "Neurodevelopmental"

GeneCorrelations <- GeneCorrelations %>%filter(!is.na(FigLabel))%>%arrange(desc(RNAModule))%>%rename(`RNA module (R_)` = RNAModule)
RNAModuleSummary <- GeneCorrelations %>%distinct(`RNA module (R_)`, FigLabel)%>%left_join(RNAModules)%>%mutate(ModuleID = paste0("R_", `RNA module (R_)`))%>%mutate(PurityR = as.numeric(PurityR), DistR = as.numeric(DistR))

FigLabelColors <- c("AP1/Mesenchymal","Oligodendrocyte","Neuronal hijacking","Neuron", "Neurodevelopmental","Microglia", "EGFR amplification","Astrocyte", "Interferon response")
Colors = c("#e3342f", "#f66d9b", "#9561e2", "#ff781f", "#4dc0b5", "#3490dc", "#e1ad01", "#811453", "#8B8000")
names(Colors) <- FigLabelColors

ggplot(PeakCorrelations_fig, aes(x=PurityR, y = DistR, color = FigLabel, group = "none"))+
  geom_point(aes(alpha = 0.8))+
  labs(x= "Purity R", y= "Distance R")+
  theme_bw(base_size = 11)+
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+ 
  geom_smooth(method='lm', se = FALSE)+
  stat_cor(method = "pearson", size =3)+
  geom_label_repel(data = Correlations_summary, aes(PurityR, DistR, label = ModuleID), size = 3, fontface = "italic")+
  scale_color_manual(values = Colors, breaks=c("Oligodendrocyte", "Astrocyte", "Neuron", "Microglia", "Neurodevelopmental", "EGFR amplification", "Neuronal hijacking", "AP1/Mesenchymal", "Interferon response"))+
  xlim(-1,1)+ylim(-0.6,0.6)

ggplot(GeneCorrelations, aes(x=Gene_PurityR, y = Gene_DistR, color = FigLabel, group = "none"))+
  geom_point(aes(alpha = 0.8))+
  labs(x= "Purity R", y= "Distance R")+
  theme_bw(base_size = 11)+
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+ 
  geom_smooth(method='lm', se = FALSE)+
  stat_cor(method = "pearson", size =3)+
  geom_label_repel(data = RNAModuleSummary, aes(PurityR, DistR, label = ModuleID), size = 3, fontface = "italic")+
  scale_color_manual(values = Colors, breaks=c("Oligodendrocyte", "Astrocyte", "Neuron", "Microglia", "Neurodevelopmental", "EGFR amplification", "Neuronal hijacking", "AP1/Mesenchymal", "Interferon response"))+
  xlim(-1,1)+ylim(-0.6,0.6)



