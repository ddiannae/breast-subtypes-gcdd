library(dplyr)

setwd("/mnt/ddisk/transpipeline-data/breast-data/subtypes")
conds <- c("healthy", "luma", "lumb", "her2", "basal")

for(cond in conds) {
  net <- read.delim(file = paste0("networks/aracne-cluster/1e8/13319/", cond, "-13319-annotated.tsv"), 
                    sep = "\t", stringsAsFactors = F)
  net <- net %>% mutate(Intra = Target.Chr == Source.Chr, 
                        Distance = ifelse(Intra, ifelse(net$Target.Start > net$Source.Start, net$Target.Start - net$Source.Start,
                                                                                                                     net$Source.Start - net$Target.Start), 0) )
  net <- net %>% select(Source, MI, Target, Distance, Intra)
  net2 <- net
  colnames(net2) <- c("Target", "MI", "Source", "Distance", "Intra")
  net2 <- net2 %>% select(Source, MI, Target, Distance, Intra)
  net3 <- rbind(net, net2)
  maxMI <- max(net3$MI)
  net3$MI <- net3$MI/maxMI
  mean.vals <- net3 %>% group_by(Source) %>% summarise(meanMI = mean(MI), sumDistance = sum(Distance), sumIntra = sum(Intra), meanDistance = sumDistance/sumIntra) 
  write.table(mean.vals, file = paste0("networks/aracne-cluster/1e8/13319/", cond, "-meanNormMI-distance-byGene.tsv"),
              sep = " \t", row.names = F, col.names = T, quote = F)
}  
  