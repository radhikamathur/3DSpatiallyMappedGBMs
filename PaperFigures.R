## R code to analyze data and generate figures for the paper: "Glioblastoma evolution and heterogeneity from a 3D whole-tumor perspective" by Mathur et al.

#install and load R packages
library(pacman)
p_load(ggplot2, tidyverse, RColorBrewer, ggpubr, openxlsx, ggrepel, ComplexHeatmap, circlize, janitor, viridis, EnvStats, VennDiagram)

#read in data from Supplementary Tables (change 'dataPath' to the correct folder location)

dataPath <- '~/Dropbox/Postdoc/Papers/ATAC Paper/Revision/Supplementary Tables/'
InputFile <-  read.xlsx(paste0(dataPath, 'TableS1 - Patient and sample information.xlsx'),sheet =3, startRow = 2, )
SampleInfo <- InputFile%>%dplyr::select(Patient, Sample, Purity, Dist, RNAsubtype, ID)%>%mutate(Sample = as.factor(Sample))
PyCloneClusters <- read.xlsx(paste0(dataPath, 'TableS3 - Genomic characterization.xlsx'), sheet = 3)
FACETSResults <- read.xlsx(paste0(dataPath, 'TableS3 - Genomic characterization.xlsx'), sheet = 5)
RNA <- read.xlsx(paste0(dataPath, 'TableS5 - Transcriptome analysis.xlsx'), sheet = 2)
Signatures <- read.xlsx(paste0(dataPath, 'TableS5 - Transcriptome analysis.xlsx'), sheet = 5)
RNAModules <- read.xlsx(paste0(dataPath, 'TableS5 - Transcriptome analysis.xlsx'), fillMergedCells = TRUE, sheet = 3) %>%row_to_names(1)%>%mutate(ModuleID = paste0("R_", `RNA module (R_)`), Genes = as.numeric(Genes), PurityR = as.numeric(PurityR), DistR = as.numeric(DistR))%>%relocate(ModuleID)
RNAModuleGenes <- read.xlsx(paste0(dataPath, 'TableS5 - Transcriptome analysis.xlsx'),sheet = 4)
scATAC_geneactivity <- read.xlsx(paste0(dataPath, 'TableS7 - Single-cell analysis.xlsx'), sheet = 2)
ATACLinkages <- read.xlsx(paste0(dataPath, 'TableS6 - Chromatin landscape analysis.xlsx'), sheet = 2)
ATACModules <- read.xlsx(paste0(dataPath, 'TableS6 - Chromatin landscape analysis.xlsx'), fillMergedCells = TRUE, sheet = 3) %>%row_to_names(1)%>%mutate(Peaks = as.numeric(Peaks), PurityR = as.numeric(PurityR), DistR = as.numeric(DistR))
scATAC_peaks <- read.xlsx(paste0(dataPath, 'TableS7 - Single-cell analysis.xlsx'), sheet = 3)
scATAC_amps <- read.xlsx(paste0(dataPath, 'TableS7 - Single-cell analysis.xlsx'), sheet = 4)
AdultFetalBrain <- read.xlsx(paste0(dataPath, 'TableS7 - Single-cell analysis.xlsx'), sheet = 5)
PeakCorrelations <- read.xlsx(paste0(dataPath, 'TableS6 - Chromatin landscape analysis.xlsx'), sheet = 6)
GeneCorrelations <- read.xlsx(paste0(dataPath, 'TableS5 - Transcriptome analysis.xlsx'), sheet = 6)

# save all data as RData file for easier loading
save(InputFile, SampleInfo, PyCloneClusters, FACETSResults, RNA,Signatures, RNAModules, RNAModuleGenes,scATAC_geneactivity,ATACLinkages, ATACModules,scATAC_peaks, scATAC_amps, AdultFetalBrain,PeakCorrelations, GeneCorrelations, file = "~/Dropbox/Postdoc/Papers/ATAC Paper/Other/FigureData.RData")
load("~/Dropbox/Postdoc/Papers/ATAC Paper/Other/FigureData.RData")

# set color scales for figures
PatientColors = c('P455' = "#ff7500", 'P475' = "#ae7000", "P498" = "#44cef6","P500" = "#1bd1a5", "P503" = "#FFD92F", "P519" = "#8d4bbb","P521" = "#ff0097", "P524" = "#BEBEBE", "P529" = "#7b7b96", "P530" = "#FF0000")
RNAcolors <- c('Classical' = "#3e5063", 'Mesenchymal' = '#0433ff', 'Proneural' = '#ff40ff', 'Neural' = '#ff9300')
RNAcolors_heatmap <- list(Subtype = c('Classical' = "#3e5063", 'Mesenchymal' = '#0433ff', 'Proneural' = '#ff40ff', 'Neural' = '#ff9300'))

##Fig 1A - Count the number of samples with each data type and generate summary plot
DataTypes_RNA <- InputFile %>% filter(`RNA-Seq` == "Yes")%>%group_by(Patient) %>% count(name = "RNA-Seq")
DataTypes_Exome <- InputFile %>% filter(`Exome-Seq` == "Yes")%>%group_by(Patient) %>% count(name = "Exome-Seq")
DataTypes_ATAC <- InputFile %>% filter(`ATAC-Seq` == "Yes")%>%group_by(Patient) %>% count(name = "ATAC-Seq")
DataTypes_HiC <- InputFile %>% filter(`HiC` == "Yes")%>%group_by(Patient) %>% count(name = "Hi-C")
DataTypes_snATAC <- InputFile %>% filter(`snATAC` == "Yes")%>%group_by(Patient) %>% count(name = "snATAC")

DataTypes <- reduce(list(DataTypes_RNA, DataTypes_Exome, DataTypes_ATAC, DataTypes_HiC, DataTypes_snATAC), dplyr::left_join, by = 'Patient')%>%mutate(Patient = factor(Patient))
DataTypes$`Hi-C`[DataTypes$Patient == "P455"|DataTypes$Patient == "P475"|DataTypes$Patient == "P498"|DataTypes$Patient == "P503"|DataTypes$Patient == "P519"|DataTypes$Patient == "P521"] <- 1

DataTypes_long <- DataTypes %>%pivot_longer(cols = -Patient, names_to = "Data", values_to = "Samples")%>%mutate(Patient = factor(Patient, levels = rev(levels(DataTypes$Patient))), `Data Type` = factor(Data, levels = c("Exome-Seq", 'RNA-Seq', 'ATAC-Seq','snATAC','Hi-C'))) %>%filter(Patient != "P565")
DataTypes_n <- DataTypes_long %>%filter(!is.na(Samples)) %>%group_by(`Data Type`)%>%summarize(Total = sum(Samples))
ggplot()+
  theme_bw(base_size = 12)+
  geom_point(data=DataTypes_long, aes(x = `Data Type`, y = Patient, color = `Data Type`, size = Samples))+
  geom_text(data=DataTypes_long, aes(x = `Data Type`, y = Patient, label = Samples), color = "black", size = 4)+
  labs(x="Samples by data type")+scale_x_discrete(position = "top")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  scale_color_manual(values = c("#977DC3", "#70AD47", "#ED7D31", "#009FEA", "#FFC001"))+
  geom_text(data = DataTypes_n, aes(x = `Data Type`, y = 0, vjust = 1.5, label = Total))+
  coord_cartesian(clip = "off")+
  theme(legend.position = "bottom") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## Fig 1B: Plot purity versus CEL (after filtering out samples on the border)
CEL <-InputFile %>% filter(!is.na(CEL))%>% filter(CEL != "Border")
ggplot(CEL, aes(x = CEL, y = Purity, color=Patient, group = CEL))+
  geom_boxplot(show.legend = F)+
  geom_jitter()+
  labs(x = "Contrast-enhancing (CE)", y = "Purity (ψ)")+
  theme_bw(base_size = 11)+
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
  scale_color_manual(values = PatientColors)+
  stat_compare_means(method = "t.test", show.legend = F)+
  ylim(-0.05,1.05)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

CEL_count <- CEL %>%group_by(CEL) %>% count()

##Fig 1C: plot purity versus relative distance from centroid 
ggplot(InputFile, aes(Dist, Purity, color=Patient, group = "none"))+
  geom_point()+
  labs(x = "Relative distance from centroid (d)", y = "Purity (ψ)")+
  theme_bw(base_size = 11)+
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+ 
  geom_smooth(method='lm', se = FALSE)+
  stat_cor(method = "pearson", size =3, show.legend = F)+
  scale_color_manual(values = PatientColors)+
  xlim(0,1.06)+ylim(0,1)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##Fig1H: Generate pairwise inter-sample distance matrix and plot for genetic subclones from P529 and P530

patients_distances <- InputFile %>%filter(Patient != "P565") %>% pull(Patient)

pairwiseDistances <- data.frame(patient=character(), sample1 = numeric(), sample2 = numeric(), distance.L=numeric(),matchingpyclone = logical(), stringsAsFactors=F)

pairwiseDist2pointsLPS <- function(p1, p2){#p1 and p2 are both vectors, in which [1] is x, [2] is y, and [3] is z
  return(sqrt( ((p2[1]-p1[1])^2) + ((p2[2]-p1[2])^2) + (((p2[3]-p1[3]))^2) ))
}
for (p in unique(patients_distances)){
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

distance_table <- pairwiseDistances %>% mutate(FirstSample = paste0(Patient, "_", FirstSample), SecondSample = paste0(Patient, "_",SecondSample)) %>% select(-Patient, -MatchingPyClone) %>%filter(!is.na(Distance))
distance_matrix <- distance_table %>% pivot_wider(id_cols = "FirstSample", names_from = "SecondSample", values_from = "Distance") 

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
  stat_compare_means(method = "t.test", vjust = 1)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

avg_distance <- pyclone_only %>% group_by(Patient, MatchingPyClone) %>% summarize(Distance = mean(Distance))


##Fig S1F/H: Plot cellular prevalence of PyClone-defined clonal clusters across samples from P529 and P530
P530_clusters <- PyCloneClusters %>% filter(Patient == "P530")
ggplot(P530_clusters, aes(x=sample, y=cellular_prevalence, group = Pyclone.cluster))+
  geom_line(aes(color = Pyclone.cluster), lwd=0.7)+
  theme_bw(base_size = 11)+
  labs(title = paste0("Clonal clusters for Patient 530"), y = "cellular prevalence")+
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
  facet_wrap(Pyclone.cluster ~., nrow = 3)+
  theme(legend.position = "none", strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.y = element_text( size = 12, angle = 0))+ylim(0,1)

P529_clusters <- PyCloneClusters %>% filter(Patient == "P529")
ggplot(P530_clusters, aes(x=sample, y=cellular_prevalence, group = Pyclone.cluster))+
  geom_line(aes(color = Pyclone.cluster), lwd=0.7)+
  theme_bw(base_size = 11)+
  labs(title = paste0("Clonal clusters for Patient 529"), y = "cellular prevalence")+
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+
  facet_wrap(Pyclone.cluster ~., nrow = 3)+
  theme(legend.position = "none", strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.y = element_text( size = 12, angle = 0))+ylim(0,1)

## Fig S1A - Calculate inter-sample differences in purity and generate plot comparing with inter-sample distances

pairwiseMatrix1 <- data.frame(patient=character(), sample1 = numeric(), sample2 = numeric(), distance.L=numeric(), purity_difference = numeric(), stringsAsFactors=F)

for (p in unique(patients_distances)){
  localP <- InputFile[which(InputFile$Patient==p),]
  allSampleCombinations <- combn(localP$Sample,2)
  for (c in 1:ncol(allSampleCombinations)){
    a <- allSampleCombinations[1,c]
    b <- allSampleCombinations[2,c]
    a.coords <- localP[which(localP$Sample==a),c('L','P','S')]
    a.purity <- localP[which(localP$Sample==a),'Purity']
    b.coords <- localP[which(localP$Sample==b),c('L','P','S')]
    b.purity <- localP[which(localP$Sample==b),'Purity']
    distance <- pairwiseDist2pointsLPS(a.coords,b.coords)
    purity_difference <- a.purity - b.purity
    toBind <- data.frame(c(patient=p, sample1 = a, sample2 = b, distance=distance, purity_difference = purity_difference))
    pairwiseMatrix1 <- bind_rows(pairwiseMatrix1, toBind)
  }
}
colnames(pairwiseMatrix1) <- c('Patient','FirstSample','SecondSample', 'Distance', 'PurityDifference')
pairwiseMatrix1$Discard <- F
pairwiseMatrix1$Discard[pairwiseMatrix1$Patient == "P530" & pairwiseMatrix1$FirstSample<10 & pairwiseMatrix1$SecondSample>9] <- T
pairwiseMatrix1 <- pairwiseMatrix1 %>%filter(Discard == F)
purity_distance_regression <- lm(abs(PurityDifference) ~ Distance, data = pairwiseMatrix1)
summary(purity_distance_regression)
pairwiseMatrix1 %>% ggplot(aes(x=Distance, y = abs(PurityDifference), color = Patient, group = "none"))+
  geom_point()+
  geom_smooth(method='lm', se = F)+
  stat_cor(method = "pearson", size =3, label.x.npc = "center")+
  scale_color_manual(values = PatientColors)+
  theme_bw(base_size = 11)+
  xlab(bquote(Inter-sample~distance~ (d[ab])))+ylab(bquote(Purity~difference~(abs(ψ[b] - ψ[a]))))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



##Fig 2A: Plot FACETS-derived copy number estimates for oncogene ampification events
FACETSResults$Gene[FACETSResults$chrom == 7 & FACETSResults$start<55324313 & FACETSResults$end>55086714] <- "EGFR"
FACETSResults$Gene[FACETSResults$chrom == 4 & FACETSResults$start<55164414 & FACETSResults$end>55095264] <- "PDGFRA"
FACETSResults$Gene[FACETSResults$chrom == 1 & FACETSResults$start<204542871 & FACETSResults$end>204485511] <- "MDM4"
FACETSResults$Gene[FACETSResults$chrom == 2 & FACETSResults$start<16087129 & FACETSResults$end>16080683] <- "MYCN"
FACETSResults$Gene[FACETSResults$chrom == 9 & FACETSResults$start<21975132 & FACETSResults$end>21967751] <- "CDKN2A"
FACETSResults$Gene[FACETSResults$chrom == 10 & FACETSResults$start<89728532 & FACETSResults$end>89623195] <- "PTEN"

AmplificationsOnly <- FACETSResults %>% filter(tcn.em >9) %>% select(Patient, Gene) %>% filter(!is.na(Gene)) %>% distinct() %>%left_join(FACETSResults)%>% group_by(ID, Gene)%>%summarize(CN = mean(tcn.em))%>%left_join(InputFile)%>%mutate(LowPurity = "No")
AmplificationsOnly$LowPurity[AmplificationsOnly$Purity<0.5] <- "Yes"
ggplot(AmplificationsOnly, aes(x=Patient, y = CN,  fill = Patient))+
  theme_bw()+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(shape = LowPurity))+
  scale_fill_manual(values = PatientColors)+
  facet_grid(.~Gene, scales = "free", space = "free")+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.y = element_text( size = 16, angle = 0))+
  theme(legend.position = "none")+
  labs(y = "Copy Number (CN)")+
  geom_hline(yintercept = 2, color = "red",linetype='dotted')+
  scale_shape_manual(values = c("Yes" = 8, "No" = 16))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##Fig2H/K: ATAC-Seq signal and gene expression comparisons between P530 frontal and temporal regions
ATAC_SpecificPeak <- ATACLinkages %>%filter(Gene == "MTAP" & Type == "Promoter")%>% pivot_longer(cols = 4:ncol(ATACLinkages), names_to = "ID", values_to = "ATAC")%>%left_join(SampleInfo)
P530_plot <- ATAC_SpecificPeak %>%filter(Patient == "P530")%>%mutate(Sample = as.numeric(Sample), Region = "Frontal")
P530_plot$Region[P530_plot$Sample<10] <- "Temporal"
P530_plot <- P530_plot %>% mutate(Region = factor(Region, levels = c("Temporal", "Frontal")))
ggplot(P530_plot, aes(x=Region, y = ATAC, fill = Region))+
  geom_boxplot()+
  geom_jitter()+
  scale_fill_manual(values =c("#C371A5", "#4B91AD"))+
  theme_bw(base_size = 12)+
  theme(plot.title = element_text(lineheight=.8, hjust = 0.5))+
  stat_compare_means(method = "t.test", vjust = -0.5, label = "p.format")+
  labs(title = "MTAP promoter", y = "ATAC signal")+
  theme(title = element_text(size = 12),axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ATAC_SpecificPeak <- ATACLinkages %>%filter(Gene == "PTEN" & Type == "Promoter")%>% pivot_longer(cols = 4:ncol(ATACLinkages), names_to = "ID", values_to = "ATAC")%>%left_join(SampleInfo)
P530_plot <- ATAC_SpecificPeak %>%filter(Patient == "P530")%>%mutate(Sample = as.numeric(Sample), Region = "Frontal")
P530_plot$Region[P530_plot$Sample<10] <- "Temporal"
P530_plot <- P530_plot %>% mutate(Region = factor(Region, levels = c("Temporal", "Frontal")))
ggplot(P530_plot, aes(x=Region, y = ATAC, fill = Region))+
  geom_boxplot()+
  geom_jitter()+
  scale_fill_manual(values =c("#C371A5", "#4B91AD"))+
  theme_bw(base_size = 12)+
  theme(plot.title = element_text(lineheight=.8, hjust = 0.5))+
  stat_compare_means(method = "t.test", vjust = -0.5, label = "p.format")+
  labs(title = "PTEN promoter", y = "ATAC signal")+
  theme(title = element_text(size = 12),axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

RNA_SpecificGenes <- RNA %>% filter(Gene == "SHLD2"|Gene == "PTEN"|Gene == "RNLS") %>% pivot_longer(cols = !Gene, names_to = "Sample", values_to = "RNA") %>% separate(Sample, c("Patient", "Sample")) %>% mutate(Sample = as.factor(as.numeric(Sample))) %>%mutate(Gene = gsub("SHLD2", "FAM35A", Gene))%>%mutate(Gene = factor(Gene, levels = c("PTEN", "FAM35A", "RNLS")))
RNA_SpecificGenes <- RNA %>% filter(Gene == "KLHL9"|Gene == "MTAP") %>% pivot_longer(cols = !Gene, names_to = "Sample", values_to = "RNA") %>% separate(Sample, c("Patient", "Sample")) %>% mutate(Sample = as.factor(as.numeric(Sample))) %>%mutate(Gene = factor(Gene, levels = c("MTAP", "KLHL9")))

P530_plot <- RNA_SpecificGenes %>%filter(Patient == "P530")%>%mutate(Sample = as.numeric(Sample), Region = "Frontal")
P530_plot$Region[P530_plot$Sample<10] <- "Temporal"
P530_plot <- P530_plot %>% mutate(Region = factor(Region, levels = c("Temporal", "Frontal")))
ggplot(P530_plot, aes(x=Region, y = RNA, fill = Region))+
  geom_boxplot()+
  geom_jitter()+
  scale_fill_manual(values =c("#C371A5", "#4B91AD"))+
  theme_bw(base_size = 12)+
  facet_wrap(Gene ~., scales = "free")+
  stat_compare_means(method = "t.test", vjust = -0.5,  label = "p.format")+
  theme(strip.text = element_text(size = 12),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  labs(y = "Gene expression", x = NULL)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##Fig S2C: EGFR and PDGFRA copy-number versus gene expression
GenesOnly <- FACETSResults %>% filter(!is.na(Gene))%>% group_by(ID, Gene)%>%summarize(CN = mean(tcn.em))%>%left_join(InputFile)
RNA_SpecificGenes <- RNA %>% filter(Gene %in% GenesOnly$Gene)%>% pivot_longer(cols = !Gene, names_to = "ID", values_to = "RNA")
GenesOnly_RNA <- GenesOnly %>% left_join(RNA_SpecificGenes)%>%mutate(LowPurity = "No")
GenesOnly_RNA$LowPurity[GenesOnly_RNA$Purity<0.5] <- "Yes"

EGFR_P521 <-  GenesOnly_RNA %>% filter(Gene == "EGFR" & Patient == "P521") %>% mutate(Sample = as.factor(Sample))
ggplot(EGFR_P521, aes(x=CN, y = RNA,  label = Sample, color = Sample, group = "none"))+
  theme_bw()+
  stat_cor(method = "pearson", size =4, label.y.npc = "top")+
  geom_smooth(method='lm', se = F, aes(color = "grey", alpha = 0.5))+ 
  geom_point(aes(shape = LowPurity))+
  scale_shape_manual(values = c("Yes" = 8, "No" = 16))+
  geom_text_repel()+
  scale_color_brewer(palette = "Set1")+
  theme(legend.position = "none")+
  #xlim(0, 25)+ylim(3,10)+
  ylab ("Expression (log2tpm+1)")+
  ggtitle("EGFR (P521)")+
  geom_vline(xintercept = 2, color = "red",linetype='dotted')

PDGFRA_P521 <-  GenesOnly_RNA %>% filter(Gene == "PDGFRA" & Patient == "P521") %>% mutate(Sample = as.factor(Sample))
ggplot(PDGFRA_P521, aes(x=CN, y = RNA,  label = Sample, color = Sample, group = "none"))+
  theme_bw()+
  stat_cor(method = "pearson", size =4)+
  geom_smooth(method='lm', se = F, aes(color = "grey", alpha = 0.5))+ 
  geom_point(aes(shape = LowPurity))+
  scale_shape_manual(values = c("Yes" = 8, "No" = 16))+
  geom_text_repel()+
  scale_color_brewer(palette = "Set1")+
  theme(legend.position = "none")+
  xlim(0, 25)+ylim(3,12)+
  ylab ("Expression (log2tpm+1)")+
  ggtitle("PDGFRA (P521)")+
  geom_vline(xintercept = 2, color = "red",linetype='dotted')

##Fig S2E: CDKN2A and PTEN copy-number versus expression in P530

PTEN_P530 <-  GenesOnly_RNA %>% filter(Gene == "PTEN" & Patient == "P530") %>% mutate(Sample = as.numeric(Sample)) 
PTEN_P530$Region <- "Temporal"
PTEN_P530$Region[PTEN_P530$Sample>9] <-  "Frontal"
ggplot(PTEN_P530, aes(x=CN, y = RNA,  label = Sample, color = Region, group = "none"))+
  theme_bw()+
  #stat_cor(method = "pearson", size =4)+
  geom_text_repel(show.legend = FALSE)+
  geom_point(aes(shape = LowPurity))+
  scale_shape_manual(values = c("Yes" = 8, "No" = 16))+
  ylab ("Expression (log2tpm+1)")+
  xlab("Copy number (CN)")+
  ggtitle("PTEN (P530)")+
  scale_color_manual(values = c( "deepskyblue4","deeppink4"))

CDKN2A_P530 <-  GenesOnly_RNA %>% filter(Gene == "CDKN2A" & Patient == "P530") %>% mutate(Sample = as.factor(Sample))
CDKN2A_P530$Region <- "Temporal"
CDKN2A_P530$Region[PTEN_P530$Sample>9] <-  "Frontal"
ggplot(CDKN2A_P530, aes(x=CN, y = RNA,  label = Sample, color = Region, group = "none"))+
  theme_bw()+
  #stat_cor(method = "pearson", size =4)+
  geom_point(aes(shape = LowPurity))+
  scale_shape_manual(values = c("Yes" = 8, "No" = 16))+
  geom_text_repel(show.legend = FALSE)+
  #theme(legend.position = "none")+
  ylab ("Expression (log2tpm+1)")+
  xlab("Copy number (CN)")+
  ggtitle("CDKN2A (P530)")+
  scale_color_manual(values = c( "deepskyblue4","deeppink4"))

##Fig 4A: Heatmap of average gene expression across samples for subtypes defined by Verhaak et al. Cancer Cell 2010
SignaturesList <- Signatures%>%pivot_longer(cols = everything(), names_to = "Signature", values_to = "Gene")%>% na.omit()%>%filter(str_detect(Signature, "Verhaak"))%>%left_join(RNA, by = "Gene")%>%na.omit()%>% pivot_longer(cols = !c(Gene, Signature), names_to = "Sample", values_to = "RNA") %>% separate(Sample, c("Patient", "Sample")) %>% 
  group_by(Patient, Sample,Signature) %>% summarise(Expression = mean(RNA))%>% mutate(Sample = as.numeric(Sample)) %>% arrange(Patient, Sample)
SignaturesList$ID <- paste0(SignaturesList$Patient, "_", SignaturesList$Sample)
SignaturesAvgExpression <- SignaturesList %>% pivot_wider(id_cols = Signature, names_from = ID, values_from = Expression) %>% as.data.frame()%>% 
  separate(Signature, into = c("Author", "Journal", "Year", "Label"), extra = "merge")
SampleCols <- SignaturesAvgExpression%>% dplyr::select(starts_with("P4")|starts_with("P5"))%>% select(-starts_with("P56"))%>%as.matrix()
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
draw(VerhaakHeatmap)


##Fig 4B: Subtypes by purity and relative distance from centroid
SubtypePlot <- InputFile %>% filter(!is.na(RNAsubtype))
ggplot(SubtypePlot, aes(x = RNAsubtype, y = Purity, fill = RNAsubtype))+
  geom_boxplot()+
  theme_bw(base_size = 11)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab(NULL)+ylab("Purity (ψ)" )+
  theme(legend.position = "none")+
  scale_fill_manual(values = RNAcolors)+
  stat_compare_means(method = "anova", size = 3, label.x = 1.5)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(SubtypePlot, aes(x = RNAsubtype, y = Dist, fill = RNAsubtype))+
  geom_boxplot()+
  theme_bw(base_size = 11)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab(NULL)+ylab("Relative distance (d)")+
  theme(legend.position = "none")+
  scale_fill_manual(values = RNAcolors)+
  stat_compare_means(method = "anova", size = 3, label.y = 0.25, label.x = 1.5)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##Fig4C: RNA module heatmap
ModuleAnnot <- RNAModules 
colors_module = as.character(ModuleAnnot$`RNA module (R_)`)
ModuleRowAnnot <- rowAnnotation(annotation_name_side = "top", 
                                Module = anno_text(ModuleAnnot$ModuleID, gp = gpar(fontsize = 7.5, fontface = "italic")),
                                Genes = anno_barplot(ModuleAnnot$Genes, gp = gpar(fill = colors_module, fontsize = 8), axis_param = list(side = "top")), 
                                Label = anno_text(ModuleAnnot$`Label/Top GSEA result (FDR)`, gp = gpar(fontsize = 7.5)),
                                PurityR = anno_barplot(ModuleAnnot$PurityR, width = unit(1.5, "cm"), gp = gpar(fill = colors_module), axis_param = list(side = "top")),
                                DistR = anno_barplot(ModuleAnnot$DistR, width = unit(1.5, "cm"), gp = gpar(fill = colors_module), axis_param = list(side = "top")))

PatientCols <- ModuleAnnot %>% dplyr::select(ModuleID, starts_with("P4")|starts_with("P5"))%>% pivot_longer(cols = -ModuleID,names_to = "ID", values_to = "RNA") %>%separate(ID, into = "Patient") %>%filter(Patient != "P565")%>%mutate(RNA = as.numeric(RNA))
CohortVariance <- PatientCols %>% select(-Patient) %>%group_by(ModuleID) %>%summarize(var = var(RNA))%>%column_to_rownames("ModuleID") %>% as.matrix()

SampleCols <- ModuleAnnot %>% dplyr::select(starts_with("P4")|starts_with("P5"))%>% select(-starts_with("P56"))%>% mutate_if(is.character,as.numeric)%>% as.matrix()
rownames(SampleCols) <- ModuleAnnot$ModuleID

sample_annotation = HeatmapAnnotation(Sample = anno_text(SampleAnnot$Sample, rot = 0, just = "top", gp = gpar(fontsize = 7)))
sample_annotation2 = HeatmapAnnotation(Purity = anno_barplot(SampleAnnot$Purity, gp = gpar(fill = SampleAnnot$Color), axis_param = list(side = "right")),annotation_name_side = "left", annotation_name_rot = 0,gap = unit(2, "mm"),
                                       Dist = anno_points(SampleAnnot$Dist, ylim = c(0, 1.2), gp = gpar(col = "black"), axis_param = list(side = "right",at = c(0,1), labels = c("Centroid", "Periphery"))),
                                       Subtype = SampleAnnot$RNAsubtype, col = RNAcolors_heatmap)

Heatmap <- Heatmap(SampleCols, name = "Average Expression",  show_row_names = F, show_column_names = F,  top_annotation = sample_annotation, bottom_annotation = sample_annotation2, column_split = SampleAnnot$Patient, 
                   clustering_method_columns = 'ward.D2', clustering_method_rows = 'ward.D2',
                   column_names_gp = gpar(fontsize = 8),row_names_gp = gpar(fontsize = 8), cluster_rows = T, cluster_column_slices = F, cluster_columns = T, col = rev(brewer.pal(n = 11, name = "RdYlBu")))
HeatmapCohortVariance <- Heatmap(CohortVariance, name = "Variance", show_row_names = F, show_column_names = F, right_annotation = ModuleRowAnnot, col=rev(brewer.pal(11,"RdBu")))
draw(Heatmap + HeatmapCohortVariance, merge_legend = T, heatmap_legend_side = "bottom")

##Fig 4D,4G - RNA modules by subtype
RNASeq_All_Pivot <- RNA %>%pivot_longer(cols = -Gene, names_to = "ID", values_to = "RNA")
RNASeq_All_Ordered <- RNASeq_All_Pivot %>%left_join(SampleInfo)
RNAModuleGenesBySample <- RNAModuleGenes%>% select(Module, Gene) %>% inner_join(RNASeq_All_Ordered)
RNAModuleExpressionBySample <- RNAModuleGenesBySample %>% group_by(Module, ID, Patient, Sample, RNAsubtype) %>% summarize(RNA = mean(RNA))

Verhaak <- RNAModuleExpressionBySample %>%filter(Module == "orangered3" |Module == "brown"|Module == "plum")%>%mutate(ModuleID = paste0("R_", Module)) %>%mutate(ModuleID = factor(ModuleID, levels = c("R_orangered3", "R_plum","R_brown")))
Verhaak <- RNAModuleExpressionBySample %>%filter(Module == "plum2" |Module == "midnightblue"|Module == "plum3"|Module == "darkred")%>%mutate(ModuleID = paste0("R_", Module)) %>%mutate(ModuleID = factor(ModuleID, levels = c("R_darkred","R_plum3", "R_plum2","R_midnightblue", "R_blue"))) ##Fig 5G
Verhaak <- RNAModuleExpressionBySample %>%filter(Module == "darkseagreen4" |Module == "lightcoral"|Module == "pink") %>%mutate(ModuleID = paste0("R_", Module)) 

ggplot(Verhaak, aes(x = RNAsubtype, y = RNA, fill = RNAsubtype))+
  geom_boxplot()+
  facet_wrap(.~ModuleID, scales = "free", nrow= 1)+
  scale_fill_manual(values = RNAcolors)+
  theme_bw(base_size = 11)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.x = element_text(size = 11, face = "italic"))+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(legend.position = "bottom")+
  labs(y = "Average expression",  x= NULL)+ theme(legend.title = element_blank())+
  stat_compare_means(method = "anova", size = 3.5, hjust = -0.1, vjust = 0.5, label = "p.format")

ggplot(Verhaak, aes(x = ModuleID, y = RNA, color = RNAsubtype))+
  geom_jitter(size = 1, width = 0.2)+
  theme_bw(base_size = 11)+
  scale_color_manual(values = RNAcolors)+
  theme(axis.text.y = element_text(face = "italic"))+
  xlab(NULL)+ylab("Average expression")+coord_flip()+
  theme(legend.position = "bottom")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#Fig 4E/S4C - RNA modules enriched towards the periphery 
Periph <- RNAModuleGenesBySample  %>%filter(Module == "lightcoral") %>%group_by(Patient, Sample, Module, Dist, RNAsubtype) %>%summarize(RNA = mean(RNA))%>%mutate(ModuleID = paste0("R_", Module))
Periph <- RNAModuleGenesBySample  %>%filter(Module == "pink"|Module == "darkseagreen4") %>%group_by(Patient, Sample, Module, Dist, RNAsubtype) %>%summarize(RNA = mean(RNA))%>%mutate(ModuleID = paste0("R_", Module))

ggplot(Periph , aes(x = Dist, y = RNA, color = RNAsubtype, group = "none"))+
  geom_point()+
  facet_wrap(.~ModuleID, scales = "free")+
  scale_color_manual(values = RNAcolors)+
  theme_bw(base_size = 10)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.x = element_text(face = "italic"))+
  theme(legend.position = "bottom")+ theme(legend.title = element_blank())+
  geom_smooth(method='lm', se = F)+
  stat_cor(method = "pearson", size =3)+
  labs(x = "Relative distance (d)", y = "Average expression")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##Fig 4I - Association of R_darkred with CE regions

CE_list <- InputFile %>%select(ID, CEL) %>% filter(!is.na(CEL))%>% filter(CEL != "Border")
Darkred <- RNAModuleExpressionBySample %>%filter(Module == "darkred") %>%inner_join(CE_list) 
ggplot(Darkred, aes(x = CEL, y = RNA, color = RNAsubtype, group = CEL))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter()+
  theme_bw(base_size = 10)+
  scale_color_manual(values = RNAcolors)+
  labs(title = "R_darkred", x = "Contrast-enhancing",y = "Average expression" )+
  theme(plot.title = element_text(face = "italic", hjust = 0.4, size = 11))+
  theme(legend.position = "none")+
  stat_compare_means(method = "t.test", label = "p.format", label.y = 5)+ ylim(2,5.2)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##Fig S4B - Neural samples by RNA module

NeuralList <- c("P455_1", "P455_2", "P455_10", "P503_1","P503_2", "P524_1", "P524_7")
Neural <- RNAModuleGenesBySample  %>%filter(Module == "orangered3" |Module == "brown"|Module == "plum")%>%filter(ID %in% NeuralList)%>%mutate(Sample = as.factor(Sample), ModuleID = paste0("R_", Module))
Neural %>% select(ID) %>% distinct()
ModuleFill <- as.list(Neural$Module)
names(ModuleFill) <- as.list(Neural$ModuleID)
ggplot(Neural, aes(x = Sample, y = RNA, fill = ModuleID))+
  geom_boxplot(outlier.shape = NA)+
  facet_grid(.~Patient, scales = "free", space = "free")+
  scale_fill_manual(values = ModuleFill)+
  theme_bw(base_size = 11)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  ylab("Average expression")+ theme(legend.title = element_blank(), legend.text = element_text(face = "italic"))+ylim(0,8)+
  theme(legend.position = "bottom")+xlab(NULL)

##Fig S4I - Correlation matrix for microenvironmental RNA modules 
Microenvironment <- ModuleAnnot %>% filter(`RNA module (R_)` == "blue" | `RNA module (R_)` == "midnightblue" | `RNA module (R_)` == "greenyellow"| `RNA module (R_)` == "midnightblue"| `RNA module (R_)` == "plum2"| `RNA module (R_)` == "plum3"| `RNA module (R_)` == "darkred"| `RNA module (R_)` == "brown"| `RNA module (R_)` == "plum"| `RNA module (R_)` == "orangered3")
MicroenvironmentCols <- Microenvironment %>% dplyr::select(starts_with("P4")|starts_with("P5"))%>% mutate_if(is.character,as.numeric)%>% as.matrix()
rownames(MicroenvironmentCols) <- Microenvironment$ModuleID
cor_matrix <- round(cor(t(MicroenvironmentCols)), 3) 
ModuleMatrixHeatmap <- Heatmap(cor_matrix, name = "Correlation", column_names_gp = gpar(fontsize = 10, fontface = "italic"),row_names_gp = gpar(fontsize = 10,fontface = "italic"), col = colorRamp2(c(-1,0,1), c('blue', 'white', 'red')))
draw(ModuleMatrixHeatmap)


# Fig 5A - R_turquoise module versus purity

PurityFigs <- RNAModuleGenesBySample %>%filter(Module == "turquoise") %>%group_by(Patient, Sample, Module, Dist, Purity, RNAsubtype)%>%summarize(RNA = mean(RNA)) %>%mutate(ModuleID = paste0("R_", Module))%>%mutate(ModuleID = factor(ModuleID, levels = c("R_turquoise", "R_brown2", "R_maroon","R_ivory")))%>%mutate(Group = "Other")
ggplot(PurityFigs, aes(x = Purity, y = RNA, color = Patient, group = "none"))+
  geom_point()+
  scale_color_manual(values = PatientColors)+
  theme_bw(base_size = 11)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.y = element_text(face = "italic"))+
  geom_smooth(method='lm', se = F)+
  stat_cor(method = "pearson", size =4, show.legend = F)+
  labs(x = "Sample purity", y = "Average expression", title = "R_turquoise")+
  theme(plot.title = element_text(lineheight=.8, face="italic",hjust = 0.5))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


##Fig 5C - Differential module enrichment in GBMs with primitive neuronal component (plotting against purity)
PurityFigs <- RNAModuleGenesBySample %>%filter(Module == "ivory" |Module == "maroon"|Module == "brown2") %>%group_by(Patient, Sample, Module, Dist, Purity, RNAsubtype)%>%summarize(RNA = mean(RNA)) %>%mutate(ModuleID = paste0("R_", Module))%>%mutate(ModuleID = factor(ModuleID, levels = c("R_turquoise", "R_brown2", "R_maroon","R_ivory")))%>%mutate(Group = "Other")
PurityFigs$Group[PurityFigs$Patient == "P475"] <- "P475"
ggplot(PurityFigs, aes(x = Purity, y = RNA, color = Patient, group = Group))+
  geom_point()+
  facet_wrap(.~ModuleID, scales = "free")+
  scale_color_manual(values = PatientColors)+
  theme_bw(base_size = 11)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.x = element_text(face = "italic"))+
  geom_smooth(method='lm', se = F)+
  stat_cor(method = "pearson", size =3, show.legend = F)+
  labs(x = "Sample purity", y = "Average expression")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##FigS5C - Differential module enrichment in GBMs with primitive neuronal component (with P565)
InterpatientModules <- RNAModuleGenesBySample %>%filter(Module == "ivory" |Module == "maroon"|Module == "brown2") %>%group_by(Patient, Sample, Module, Dist, Purity, RNAsubtype)%>%summarize(RNA = mean(RNA)) %>%mutate(ModuleID = paste0("R_", Module))%>%mutate(ModuleID = factor(ModuleID, levels = c("R_ivory", "R_turquoise", "R_brown2", "R_maroon")))%>%mutate(Group = "Other")%>%mutate(Patient = replace_na(Patient, "P565"))
ggplot(InterpatientModules, aes(x = Patient, y = RNA, fill = Patient))+
  geom_boxplot()+
  facet_wrap(.~ModuleID, scales = "free", nrow = 1)+
  scale_fill_manual(values = PatientColors)+
  theme_bw(base_size = 11)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.x = element_text(face = "italic"))+
  theme(legend.position = "none")+
  labs(x = "Patient", y = "Average expression")+
  theme(axis.text.x = element_text(angle = 90))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


##Fig 6A - scATAC cell distributions by sample and neoplastic status
cell_annotations_count <- scATAC_geneactivity %>% group_by(Patient, Sample) %>% summarize(Count = n())
ggplot(cell_annotations_count, aes(x=Patient, y = Count, fill = Sample))+
  geom_col(colour = "black")+
  geom_text(aes(label = Sample),position = "stack", vjust = 1.3, size = 3.5)+
  theme_bw(base_size = 12)+
  scale_fill_brewer(palette = "Set1", direction = -1)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(legend.position = "none")+
  ylab("Number of cells")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(scATAC_geneactivity, aes(x=UMAP_1, y = UMAP_2, color = Sample))+
  geom_point(size = 1, shape = 16, alpha = 0.7)+
  facet_wrap(.~Patient, scales = "free", nrow = 1)+
  theme_bw(base_size = 11)+
  scale_color_brewer(palette = "Set1", direction = -1)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.y = element_text( size = 12))+
  theme(legend.position = "none")+
  labs(title = "Samples")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(scATAC_geneactivity, aes(x=Patient, fill = Neoplastic))+
  geom_bar(colour = "black")+
  theme_bw(base_size = 12)+
  #scale_fill_brewer(palette = "Set1", direction = -1)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(legend.position = "bottom")+
  scale_fill_manual(values = c("azure4", "deeppink3"))+
  ylab("Number of cells")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(scATAC_geneactivity, aes(x=UMAP_1, y = UMAP_2, color = Neoplastic))+
  geom_point(size = 1, shape = 16, alpha = 0.7)+
  facet_wrap(~Patient, scales = "free", nrow =1)+
  theme_bw(base_size = 11)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.y = element_text( size = 12))+
  scale_color_manual(values = c("azure4", "deeppink3"))+
  theme(legend.position = "none")+
  labs(title = "Neoplastic")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#Fig S6A: Neoplastic % in scATAC versus purity from bulk WES
neoplastic_percent <- scATAC_geneactivity %>%group_by(ID) %>% summarize(Neoplastic = mean(Neoplastic))%>%left_join(SampleInfo)
ggplot(neoplastic_percent, aes(x=Neoplastic, y = Purity, group = "none", label = Sample))+
  geom_smooth(method='lm', se = F)+
  geom_point()+
  theme_bw(base_size = 11)+
  xlim(0,1)+ylim(0,1)+
  xlab("% neoplastic (scATAC)")+ ylab("Purity (bulk exome)")+
  stat_cor(method = "pearson", size =3)+
  theme(legend.position = "bottom")

# Fig S6B - scATAC gene activity scores for RNA modules
scATAC_geneactivity_longer <- scATAC_geneactivity%>% pivot_longer(cols = 8:ncol(scATAC_geneactivity), names_to = "ModuleID", values_to = "Score") %>%arrange(Score)
scATAC_gene_select <- scATAC_geneactivity_longer %>%filter(ModuleID == "R_blue" | ModuleID == "R_midnightblue"| ModuleID == "R_brown"|ModuleID == "R_maroon")%>%arrange(Score)
scATAC_gene_select$ModuleID <- factor(scATAC_gene_select$ModuleID, levels = c('R_maroon', 'R_brown', 'R_blue', 'R_midnightblue'))

ggplot(scATAC_gene_select, aes(x=UMAP_1, y = UMAP_2, color = Score))+
  geom_point(size = 1, shape = 16, alpha = 0.7)+
  facet_wrap(ModuleID~Patient, scales = "free")+
  theme_bw(base_size = 11)+
  scale_color_viridis_c(option = "inferno", direction = -1, limits = c(0, 0.003))+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.x = element_blank())+
  labs(x = NULL, y = NULL)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

ggplot(scATAC_gene_select, aes(x= Neoplastic, y = Score, color = Neoplastic))+
  geom_jitter()+
  geom_violin(color = "black", fill = NA)+
  theme_bw(base_size = 11)+
  scale_color_manual(values = c("azure4", "deeppink3"))+
  facet_grid(ModuleID~., scales = "free")+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.y = element_text( size = 12, face = "italic"))+
  theme(legend.position = "none")+
  stat_compare_means(method = "t.test", label = "p.format", vjust = 1)

scATAC_gene_select <- scATAC_geneactivity_longer %>%filter(ModuleID == "R_darkred")%>%arrange(Score)
ggplot(scATAC_gene_select, aes(x= Neoplastic, y = Score, color = Neoplastic))+
  geom_jitter()+
  geom_violin(color = "black", fill = NA)+
  theme_bw(base_size = 11)+
  scale_color_manual(values = c("azure4", "deeppink3"))+
  facet_grid(ModuleID~., scales = "free")+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.y = element_text( size = 12, face = "italic"))+
  theme(legend.position = "none")+
  stat_compare_means(method = "t.test", label = "p.format", vjust = 1)



##Fig S6C - single-cell correlation matrices for microenvironmental RNA modules 
matrix_geneactivity <- scATAC_geneactivity %>% select(R_blue, R_midnightblue, R_greenyellow, R_darkred, R_plum, R_plum2, R_plum3, R_brown, R_orangered3)%>%as.matrix()
cor_matrix <- round(cor(matrix_geneactivity), 3) 
MatrixHeatmap <- Heatmap(cor_matrix, name = "Correlation",  column_names_gp = gpar(fontsize = 10, fontface = "italic"),row_names_gp = gpar(fontsize = 10,fontface = "italic"), col = colorRamp2(c(-1,0,1), c('blue', 'white', 'red')))
draw(MatrixHeatmap)

##Fig5G - ATAC & linkage module heatmap
rownames(ATACModules) <- NULL
SampleCols <- ATACModules%>% column_to_rownames(var = "ModuleID")%>% dplyr::select(starts_with("P4")|starts_with("P5"))%>% mutate_if(is.character,as.numeric)%>% as.matrix()
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
                                #TopGene = anno_text(ModuleAnnot$`Top Linked Gene (# peaks)`, gp = gpar(fontsize = 7.5)),
                                PurityR = anno_barplot(ModuleAnnot$PurityR, width = unit(1.5, "cm"), gp = gpar(fill = colors_module), axis_param = list(side = "top")),
                                  DistR = anno_barplot(ModuleAnnot$DistR, width = unit(1.5, "cm"), gp = gpar(fill = colors_module), axis_param = list(side = "top")))


PeaksBySampleHeatmap <- Heatmap(SampleCols, name = "Average ATAC signal",  bottom_annotation = sample_annotation_bottom, cluster_columns = T, cluster_column_slices = F,
                                cluster_row_slices = F, show_row_names = F,top_annotation = sample_annotation_top,
                                clustering_method_columns = 'ward.D2', clustering_method_rows = 'ward.D2',
                                column_split = SampleAnnot$Patient, 
                                row_split = ModuleAnnot$ModuleSet,row_title = c("ATAC Modules", "Linkage Modules"),
                                show_column_names = F,
                                column_names_gp = gpar(fontsize = 8),row_names_gp = gpar(fontsize = 8), 
                                col = rev(brewer.pal(n = 11, name = "RdYlBu")))

ATACVariance <-  ATACModules%>% dplyr::select(ModuleID, starts_with("P4")|starts_with("P5")) %>%pivot_longer(cols = -ModuleID, names_to = "ID", values_to = "ATAC") %>%group_by(ModuleID) %>%summarize(Variance = var(ATAC))%>% column_to_rownames(var = "ModuleID")%>%as.matrix()
HeatmapATACVariance <- Heatmap(ATACVariance, name = "Variance", show_row_names = F, show_column_names = F, right_annotation = ModuleRowAnnot, col=rev(brewer.pal(11,"RdBu")))


draw(PeaksBySampleHeatmap + HeatmapATACVariance, heatmap_legend_side = "bottom", merge_legend = T)


## Fig6D - Comparing ATAC modules to adult and fetal brain published scATAC datasets
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

## Fig 5H - L_navajowhite1 expression in P475 v. rest
PurityFig <-ATACModules %>%filter(ModuleID == "L_navajowhite1")%>%pivot_longer(cols = 12:ncol(ATACModules), names_to = "ID", values_to = "ATAC")%>%select(ID, ATAC) %>%left_join(SampleInfo)%>%mutate(Group = "Other", ATAC = as.numeric(ATAC))
  PurityFig$Group[PurityFig$Patient == "P475"] <- "P475"
ggplot(PurityFig, aes(x = Purity, y = ATAC, color = Patient, group = Group))+
  geom_point()+
  scale_color_manual(values = PatientColors)+
  theme_bw(base_size = 11)+
  theme(legend.position = "right", plot.title = element_text(face = "italic"))+
  geom_smooth(method='lm', se = F)+
  stat_cor(method = "pearson", size =3, show.legend = F)+
  labs(title = "L_navajowhite", x = "Sample purity", y = "Average ATAC signal")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


##Fig5J - NEUROD1 gene expression by patient
SpecificGene = "NEUROD1"
RNA_SpecificGene <- RNA %>% filter(Gene == SpecificGene) %>% pivot_longer(cols = !Gene, names_to = "Sample", values_to = "RNA") %>% separate(Sample, c("Patient", "Sample")) %>% mutate(Sample = as.factor(as.numeric(Sample))) 
ggplot(RNA_SpecificGene, aes(x=Patient, y = RNA, fill = Patient))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter()+
  scale_fill_manual(values =PatientColors)+
  labs(x = NULL, y = paste0(SpecificGene, " gene expression"))+
  theme_bw(base_size = 11)+
  theme(legend.position = "none")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##FigS6G - SOX9 versus SOX10 gene expression in P521
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

## Fig6E - scATAC amplification calls for EGFR versus PDGFRA
P521_amps <- scATAC_amps %>%filter(Patient == "P521")%>%left_join(scATAC_peaks)
P521_amps$Amp[P521_amps$egfr == 1] <- "EGFR"
P521_amps$Amp[P521_amps$pdgfr == 1] <- "PDGFRA"
ggplot(P521_amps, aes(x=UMAP_1, y = UMAP_2, color = Amp))+
  geom_point(size = 1, shape = 16, alpha = 0.7)+
  facet_grid(~Patient, scales = "free")+
  theme_bw(base_size = 11)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.y = element_text( size = 12))+
  scale_color_manual(values = c("EGFR" = "#BBAED1", "PDGFRA" = "#90C786"))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

## scATAC read-in-peak UMAPs for modules with intratumoral heterogeneity in P521
scATAC_P521 <-  scATAC_peaks %>%filter(Patient == "P521") %>% pivot_longer(8:ncol(scATAC_peaks), names_to = "Module", values_to = "Score")%>%filter(Module == "A_plum3"|Module == "A_yellowgreen")%>%arrange(Score)
ggplot(scATAC_P521, aes(x=UMAP_1, y = UMAP_2, color = Score))+
  geom_point(size = 1, shape = 16, alpha = 0.7)+
  facet_grid(.~Module, scales = "free")+
  theme_bw(base_size = 11)+
  scale_color_viridis_c(option = "inferno", direction = -1)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.x = element_text( size = 12, face = "italic"))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


##Fig 6F: scATAC read-in-peak score by EGFR versus PDGFRA amplification in P521
P521_plot <- P521_amps%>%filter(!is.na(Amp))%>%select(Amp, A_plum3, A_yellowgreen, A_coral) %>%pivot_longer(cols = -Amp, names_to = "Module", values_to = "Score")%>%mutate(Module = factor(Module, levels = c("A_plum3", "A_yellowgreen", "A_coral")))
ggplot(P521_plot, aes(x = Amp, y = Score, fill = Amp))+
  geom_violin(width = 1.2)+
  facet_grid(.~Module)+
  theme_bw(base_size = 11)+
  stat_compare_means(method = "t.test", size = 3, label = "p.format")+
  scale_fill_manual(values = c("EGFR" = "#BBAED1", "PDGFRA" = "#90C786"))+
  theme(axis.text.x = element_blank())+
  theme(strip.text.x = element_text(face = "italic", angle = 0))+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  labs(x = "Amplification")+ylim(0,4)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

##Fig 6H - Association of L_salmon module with EGFR amplification
ATACModuleFigs <- ATACModules %>% dplyr::select(ModuleID|starts_with("P4")|starts_with("P5"))%>% pivot_longer(cols = -ModuleID, names_to = "ID", values_to ="ATACsignal")%>%mutate(ATACsignal = as.numeric(ATACsignal))%>%left_join(SampleInfo)
salmon4 <- ATACModuleFigs %>%filter(ModuleID == "L_salmon4")
EGFR_CN <-  FACETSResults %>%filter(Gene == "EGFR")%>%select(Patient, ID, tcn.em)%>%rename(EGFR_CN = tcn.em)
EGFR_salmon4 <-   salmon4 %>% left_join(EGFR_CN)
#PatientColors_bulk <- PatientColors[names(PatientColors) != "P524" &names(PatientColors) != "P529"]
ggplot(EGFR_salmon4, aes(x = EGFR_CN, y = ATACsignal, color = Patient, group = "none"))+
  geom_jitter()+
  theme_bw(base_size = 11)+
  scale_color_manual(values = PatientColors)+
  labs(title = "L_salmon", x = "EGFR total copy number", y = "Average ATAC signal")+ 
  geom_smooth(method='lm', se = F)+
  stat_cor(method = "pearson", size =3, show.legend = F, label.y = 3)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(face = "italic"))+
  theme(plot.title = element_text(hjust = 0.5))

##Fig 6I - scATAC-level association of L_salmon module with EGFR amplification (in P521 and P529 with intratumoral heterogeneity)
EGFR <- scATAC_peaks%>%filter(Patient == "P521"|Patient == "P529") %>% select(Cell, UMAP_1, UMAP_2, Neoplastic, L_salmon4)%>%left_join(scATAC_amps) %>%mutate(egfr = as.logical(egfr))%>%filter(Neoplastic == T)
ggplot(EGFR, aes(x = egfr, y = L_salmon4, fill = egfr))+
  geom_violin(outlier.shape = NA)+
  facet_wrap(Patient ~ .)+  theme_bw(base_size = 12)+
  stat_compare_means(method = "t.test", size = 3.5, label = "p.format")+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  labs(title = "L_salmon", x= "EGFR amplification", y = "Score")+
  scale_fill_brewer(palette = "Set1", direction = -1)+
  theme(legend.position = "none")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(face = "italic"))+
  theme(plot.title = element_text(hjust = 0.5))


##ELOVL2 and NOVA1 enhancer ATAC signal verus expression colored by EGFR_CN

ELOVL2_enhancer <- ATACLinkages %>%filter(Peak == "chr6:11035803-11036303")%>%select(-Gene, -Type, -Peak)%>% pivot_longer(cols = everything(), names_to = "ID", values_to = "ATAC")
ELOVL2_gene <- RNA %>%filter(Gene == "ELOVL2")%>% dplyr::select(-Gene)%>% pivot_longer(cols = everything(), names_to = "ID", values_to = "RNA")
MergedSignal <- left_join(ELOVL2_enhancer,ELOVL2_gene, by = "ID")%>%left_join(EGFR_CN)
ggplot(MergedSignal, aes(x=ATAC, y = RNA, color = EGFR_CN, group = "none"))+
  geom_point()+
  theme_bw(base_size=11)+ theme(plot.title = element_text(face = "bold"))+ 
  geom_smooth(method='lm')+  
  stat_cor(method = "pearson", size =3)+
  labs(x = paste0("ATAC signal at ELOVL2 enhancer"), y = paste0("Expression of ELOVL2"))+
  scale_color_viridis(direction = -1)

NOVA1_enhancer <- ATACLinkages %>%filter(Peak == "chr14:27063972-27064472")%>%select(-Gene, -Type, -Peak)%>% pivot_longer(cols = everything(), names_to = "ID", values_to = "ATAC")
NOVA1_gene <- RNA %>%filter(Gene == "NOVA1")%>% dplyr::select(-Gene)%>% pivot_longer(cols = everything(), names_to = "ID", values_to = "RNA")
MergedSignal <- left_join(NOVA1_enhancer,NOVA1_gene, by = "ID")%>%left_join(EGFR_CN)
ggplot(MergedSignal, aes(x=ATAC, y = RNA, color = EGFR_CN, group = "none"))+
  geom_point()+
  theme_bw(base_size=11)+ theme(plot.title = element_text(face = "bold"))+ 
  geom_smooth(method='lm')+  
  stat_cor(method = "pearson", size =3)+
  labs(x = paste0("ATAC signal at NOVA1 enhancer"), y = paste0("Expression of NOVA1"))+
  scale_color_viridis(direction = -1)

##Fig S6I = ATAC modules by distance from centroid
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
  stat_cor(method = "pearson", size =3, show.legend =F)+
  labs(x = "Distance from centroid", y = "Average ATAC signal")

#Fig 6K - scATAC enrichment of AP1 modules
scATAC_AP1 <-  scATAC_peaks  %>% pivot_longer(8:ncol(scATAC_peaks), names_to = "Module", values_to = "Score")%>%filter(Module == "L_orangered3")%>%filter(Patient == "P519" | Patient == "P529")%>%arrange(Score)
ggplot(scATAC_AP1, aes(x=UMAP_1, y = UMAP_2, color = Score))+
  geom_point(size = 1, shape = 16, alpha = 0.8)+
  facet_wrap(Patient~., scales = "free", nrow = 1)+
  theme_bw(base_size = 11)+
  scale_color_viridis_c(option = "inferno", direction = -1)+
  labs(title = "L_orangered3")+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.x = element_text( size = 12))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  theme(plot.title = element_text(face = "italic"))+
  theme(plot.title = element_text(hjust = 0.5))


##Fig7A -- Individual peaks and genes correlation with purity and distance

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
GeneCorrelations$FigLabel[GeneCorrelations$RNAModule == "red"] <- "Choroid/cilia"
GeneCorrelations$FigLabel[GeneCorrelations$RNAModule == "midnightblue"|GeneCorrelations$RNAModule == "plum2"|GeneCorrelations$RNAModule == "plum3"|GeneCorrelations$RNAModule == "darkred"] <- "AP1/Mesenchymal"
GeneCorrelations$FigLabel[GeneCorrelations$RNAModule == "pink"|GeneCorrelations$RNAModule == "darkseagreen4"|GeneCorrelations$RNAModule == "lightcoral"] <- "Neuronal hijacking"
GeneCorrelations$FigLabel[GeneCorrelations$RNAModule == "maroon"|GeneCorrelations$RNAModule == "turquoise"|GeneCorrelations$RNAModule == "brown2"|GeneCorrelations$RNAModule == "ivory"] <- "Neurodevelopmental"

GeneCorrelations <- GeneCorrelations %>%filter(!is.na(FigLabel))%>%arrange(desc(RNAModule))%>%rename(`RNA module (R_)` = RNAModule)
RNAModuleSummary <- GeneCorrelations %>%distinct(`RNA module (R_)`, FigLabel)%>%left_join(RNAModules)%>%mutate(ModuleID = paste0("R_", `RNA module (R_)`))%>%mutate(PurityR = as.numeric(PurityR), DistR = as.numeric(DistR))

FigLabelColors <- c("AP1/Mesenchymal","Oligodendrocyte","Neuronal hijacking","Neuron", "Neurodevelopmental","Microglia", "EGFR amplification","Astrocyte", "Interferon response", "Choroid/cilia")
Colors = c("#e3342f", "#f66d9b", "#9561e2", "#ff781f", "#4dc0b5", "#3490dc", "#e1ad01", "#811453", "#8B8000", "#d88888")
names(Colors) <- FigLabelColors

ggplot(PeakCorrelations_fig, aes(x=PurityR, y = DistR, color = FigLabel, group = "none"))+
  geom_point(aes(alpha = 0.8))+
  labs(x= "Purity R", y= "Distance R")+
  theme_bw(base_size = 11)+
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+ 
  #geom_smooth(method='lm', se = FALSE)+
  #stat_cor(method = "pearson", size =3, show.legend = F)+
  geom_label_repel(data = Correlations_summary, aes(PurityR, DistR, label = ModuleID), size = 3, fontface = "italic", show.legend = F)+
  scale_color_manual(values = Colors, breaks=c("Oligodendrocyte", "Astrocyte", "Neuron", "Microglia", "Neurodevelopmental", "EGFR amplification", "Neuronal hijacking", "AP1/Mesenchymal", "Interferon response"))+
  xlim(-1,1)+ylim(-0.6,0.6)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(GeneCorrelations, aes(x=Gene_PurityR, y = Gene_DistR, color = FigLabel, group = "none"))+
  geom_point(aes(alpha = 0.8))+
  labs(x= "Purity R", y= "Distance R")+
  theme_bw(base_size = 11)+
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+ 
  #geom_smooth(method='lm', se = FALSE)+
  #stat_cor(method = "pearson", size =3, show.legend = F)+
  geom_label_repel(data = RNAModuleSummary, aes(PurityR, DistR, label = ModuleID), size = 3, fontface = "italic", show.legend = F)+
  scale_color_manual(values = Colors, breaks=c("Oligodendrocyte", "Neuron", "Microglia", "Neurodevelopmental", "Neuronal hijacking", "AP1/Mesenchymal", "EGFR amplification", "Astrocyte", "Interferon response", "Choroid/cilia"))+
  xlim(-1,1)+ylim(-0.6,0.6)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##Fig 7B: scATAC figure annotated
Correlations_summary_select <- Correlations_summary %>% select(FigLabel, ModuleID)%>%filter(ModuleID != "A_plum3"& ModuleID != "A_plum"& ModuleID != "A_darkseagreen3"& ModuleID != "A_thistle"& ModuleID != "A_yellowgreen")
scATAC_summary <- scATAC_peaks %>% pivot_longer(8:ncol(scATAC_peaks), names_to = "ModuleID", values_to = "Score")%>%left_join(Correlations_summary_select)%>%filter(!is.na(FigLabel))
scATAC_summary_ID <- scATAC_summary%>%arrange(desc(Score))%>% distinct(Cell, .keep_all = T)%>%arrange(Score)
scATAC_summary_ID_patient <- scATAC_summary_ID%>%group_by(ID, ModuleID, FigLabel) %>%count()
ggplot(scATAC_summary_ID, aes(x=UMAP_1, y = UMAP_2, color = FigLabel))+
  geom_point(size = 1, shape = 16, alpha = 0.8)+
  facet_wrap(.~Patient, scales = "free", nrow = 1)+
  theme_bw(base_size = 11)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  scale_color_manual(values = Colors)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## Fig 6C: scATAC non_neoplastic
scATAC_nonneoplastic <- scATAC_summary_ID
scATAC_nonneoplastic$FigLabel[scATAC_nonneoplastic$FigLabel != "Oligodendrocyte" &scATAC_nonneoplastic$FigLabel != "Neuron" &scATAC_nonneoplastic$FigLabel != "Microglia"] <-"Other"
ggplot(scATAC_nonneoplastic, aes(x=UMAP_1, y = UMAP_2, color = FigLabel))+
  geom_point(size = 1, shape = 16)+
  facet_wrap(.~Patient, scales = "free", nrow = 2)+
  theme_bw(base_size = 11)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.x = element_text( size = 12))+
  scale_color_manual(values = c("#ED7B30", "#2E6CCF", "#48A0A3", "grey"))+
  theme(legend.position = "bottom")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# spatial figure for graphical abstract
ggplot(PeakCorrelations_fig, aes(x=PurityR, y = DistR, color = FigLabel,  group = "none"))+
  geom_point(aes(alpha = 0.8))+
  labs(x= "Purity", y= "Distance")+
  theme_bw(base_size = 11)+
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust = 0.5))+ 
  #geom_smooth(method='lm', se = FALSE)+
  #stat_cor(method = "pearson", size =3, show.legend = F)+
  #geom_label_repel(data = Correlations_summary, aes(PurityR, DistR, label = FigLabel), size = 3, fontface = "italic", show.legend = F)+
  scale_color_manual(values = Colors, breaks=c("Oligodendrocyte", "Astrocyte", "Neuron", "Microglia", "Neurodevelopmental", "EGFR amplification", "Neuronal hijacking", "AP1/Mesenchymal", "Interferon response"))+
  xlim(-1,1)+ylim(-0.6,0.6)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


## scATAC figures for graphical abstract
scATAC_P529 <- scATAC_geneactivity %>%filter(Patient == "P529")%>%mutate(Clone = "Red")
scATAC_P529$Clone[scATAC_P529$Sample == "10"] <- "Blue"
ggplot(scATAC_P529, aes(x=UMAP_1, y = UMAP_2, color = Clone))+
  geom_point(size = 1, shape = 16, alpha = 0.7)+
  theme_bw(base_size = 11)+
  scale_color_manual(values = c("#91B2BE", "#ECA29E"))+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(strip.text.y = element_text( size = 12))+
  theme(legend.position = "none")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

scATAC_summary_ID_P529 <- scATAC_summary_ID %>%filter(Patient == "P529")
ggplot(scATAC_summary_ID_P529, aes(x=UMAP_1, y = UMAP_2, color = FigLabel))+
  geom_point(size = 1, shape = 16, alpha = 0.8)+
  facet_wrap(.~Patient, scales = "free")+
  theme_bw(base_size = 11)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  scale_color_manual(values = Colors)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())



#RNA Module / Signature overlap
AllSignatures <- Signatures%>%pivot_longer(cols = everything(), names_to = "Signature", values_to = "Gene")%>% na.omit()
AllSignatures_RNAModule <- AllSignatures %>%left_join(RNAModuleGenes)
Red_choroid <- AllSignatures_RNAModule %>%filter(Module == "red")%>%filter(str_detect(Signature, "Choroid"))
turquoise_ipc <- AllSignatures_RNAModule %>%filter(Module == "turquoise")%>%filter(str_detect(Signature, "IPC.div2"))


##venn diagrams
Signatures <- Signatures %>%pivot_longer(cols = everything(), names_to = "Signature", values_to = "Gene") %>% separate(Signature, into = c("Data", "Journal", "Year", "Label"), extra = "merge") 

Turquoise <- RNAModuleGenes %>% filter(Module == "turquoise")%>%select(Gene) %>% distinct()
IPCdiv2 <- Signatures %>% filter(Label == "IPC-div2_upreg")%>%select(Gene)%>%filter(!is.na(Gene)) %>% distinct()
venn.diagram(x=list(Turquoise$Gene, IPCdiv2$Gene), category.names = c("R_turquoise", "IPC-div2"), filename = 'turquoise_IPC_venn.png', 
             output = T, fontfamily ="Helvetica", fill = c("turquoise", "grey"), cat.pos = 180, cat.fontface = "italic", cat.fontfamily = "Helvetica", height = 1200, width = 1200, disable.logging = TRUE)

Brown <- RNAModuleGenes %>% filter(Module == "brown")%>%select(Gene) %>% distinct()
Oligodendrocyte <- Signatures %>% filter(Label == "OL_myel")%>%select(Gene)%>%filter(!is.na(Gene)) %>% distinct()
venn.diagram(x=list(Brown$Gene, Oligodendrocyte$Gene), category.names = c("R_brown", "Oligodendrocyte"), filename = 'brown_oligo_venn.png', 
             output = T,fontfamily ="Helvetica", fill = c("brown", "orange2"), cat.pos = 180, cat.fontface = "italic", cat.fontfamily = "Helvetica", height = 1100, width = 1100, disable.logging = TRUE)

Midnightblue <- RNAModuleGenes %>% filter(Module == "midnightblue")%>%select(Gene) %>% distinct()
Mesenchymal <- Signatures %>% filter(Label == "Mesenchymal")%>%select(Gene)%>%filter(!is.na(Gene)) %>% distinct()
venn.diagram(x=list(Midnightblue$Gene, Mesenchymal$Gene), category.names = c("", ""), filename = 'Midnightblue_mesenchymal_venn_labeled.png', 
             output = T, fontfamily ="Helvetica", fill = c("darkblue", "blue"),  cat.cex = 0.5, cat.fontface = "italic", cat.fontfamily = "Helvetica", height = 1100, width = 1100, disable.logging = TRUE)

