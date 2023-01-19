##Figure1###

#initialize
library(pacman)
p_load(ggplot2, tidyverse, RColorBrewer, ggpubr, openxlsx, usethis)
dataPath <- '~/Dropbox/Postdoc/Papers/ATAC Paper/'
InputFile <-  read.xlsx(paste0(dataPath, 'Supplement/SupplementaryTableS1.xlsx'), startRow = 2)
PatientColors = c('P455' = "#ff7500", 'P475' = "#ae7000", "P498" = "#44cef6","P500" = "#1bd1a5", "P503" = "#FFD92F", "P519" = "#8d4bbb","P521" = "#ff0097", "P524" = "#BEBEBE", "P529" = "#0B0B45", "P530" = "#FF0000")
PyCloneClusters <- read.xlsx(paste0(dataPath, 'Supplement/SupplementaryTableS2.xlsx'), sheet = 2)
FACETSResults <- read.xlsx(paste0(dataPath, 'Supplement/SupplementaryTableS2.xlsx'), sheet = 4)

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



FACETSResults$Gene[FACETSResults$chrom == 9 & FACETSResults$start<21975132 & FACETSResults$end>21967751] <- "CDKN2A"
FACETSResults$Gene[FACETSResults$chrom == 10 & FACETSResults$start<89728532 & FACETSResults$end>89623195] <- "PTEN"
FACETSResults$Gene[FACETSResults$chrom == 13 & FACETSResults$start<49056026 & FACETSResults$end>48877883] <- "RB1"
FACETSResults$Gene[FACETSResults$chrom == 17 & FACETSResults$start<7590868 & FACETSResults$end>7571720] <- "TP53"

