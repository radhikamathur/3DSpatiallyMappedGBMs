##Figure1###

#initialize
library(pacman)
p_load(ggplot2, tidyverse, RColorBrewer, ggpubr, openxlsx)
dataPath <- '~/Dropbox/Postdoc/Papers/ATAC Paper/'

##Fig1C: Pairwise distances between PyClone-defined clusters for P529 and P530##

InputFile <-  read.xlsx(paste0(dataPath, 'Supplement/SupplementaryTableS1.xlsx'), startRow = 2)
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


