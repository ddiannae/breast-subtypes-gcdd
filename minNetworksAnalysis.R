library(ggplot2)
library(dplyr)
library(data.table)

all.cond <- lapply(conds, function(cond){
  network <- read.delim(paste0("networks/aracne-cluster/1e8/13319/", cond,"-13319-annotated.tsv"), 
                        stringsAsFactors = F)
  network <- network %>% mutate(Distance = ifelse(Target.Start > Source.Start, Target.Start - Source.Start,
                                                  Source.Start - Target.Start))
  
  colnames(annot) <-  c("Source", "Source.Chr", "Source.Band")
  network <- merge(annot, network, all.y = T)
  colnames(annot) <-  c("Target", "Target.Chr", "Target.Band")
  network <- merge(annot, network, all.y = T)

  network$InteractionType <- ifelse(network$Inter, "Trans", "Inter-Cytoband")
  network$InteractionType <- ifelse(network$Source.Band == network$Target.Band, "Intra-Cytoband", network$InteractionType)
  network$InteractionType <- factor(network$InteractionType, levels = c("Intra-Cytoband", "Inter-Cytoband", "Trans"))
  return(network)
})
all.cond <- rbindlist(all.cond)
all.cond <- all.cond[!is.na(all.cond$InteractionType), ]
all.cond <- all.cond %>% group_by(Cond) %>% mutate(max.byCond.MI = max(MI))
all.cond <- all.cond %>% mutate(norm.MI = MI / max.byCond.MI)
all.cond.summ.norm <- all.cond %>% group_by(Cond, InteractionType) %>% summarise(mean.MI = mean(norm.MI),
                                                                      sd.MI = sd(norm.MI), median.mi = median(norm.MI), max.MI = max(norm.MI), min.MI = min(norm.MI),
                                                                            n = n())
all.cond.summ.norm$Cond <- factor(all.cond.summ.norm$Cond, levels = conds) 
all.cond.summ.norm <- all.cond.summ.norm %>% arrange(Cond)

ggplot(all.cond.summ, aes(x = Cond, y = mean.MI, group = InteractionType)) +   # group = id is important!
  geom_path(aes(size = n, color =InteractionType, y = mean.MI, group = InteractionType),
            alpha = 0.5,
            lineend = 'round', linejoin = 'round')  +
  scale_size(breaks = NULL, range = c(1, 5))

#ggplot(network, aes(x = Distance, y = MI, color = InteractionType)) +
#  geom_point() + 
#  scale_color_manual(name = "Interaction", values =  c("#1696AC", "#73299e", "#C85200"))
#dev.off()

#ggplot(network, aes(x = MI, color = InteractionType)) +
#  geom_density() + 
#  scale_color_manual(name = "Interaction", values =  c("#1696AC", "#73299e", "#C85200"))

ggplot(all.cond[all.cond$Cond == "healthy", ], aes(x = InteractionType, y = MI, fill = InteractionType)) +
  geom_boxplot() + 
  scale_fill_manual(name = "Interaction", values =  c("#1696AC", "#73299e", "#C85200"))




# inter.by.chr <- network %>% group_by(chr) %>% tally()
# inter.by.chr <- merge(inter.by.chr, genes.by.chr, by = "chr")
# inter.by.chr$nnorm <- inter.by.chr$n / inter.by.chr$total
# inter.by.chr <- inter.by.chr[match(chrs, inter.by.chr$chr), ]
# inter.by.chr$chr <- factor(inter.by.chr$chr, levels = inter.by.chr$chr)
# 
# png(paste0(cond, "_norm_inter_by_chr.png"), width = 1600, height = 1200)
# ggplot(inter.by.chr, aes(x = chr, fill = chr, y = nnorm)) + 
#   geom_bar(stat="identity") + ylab("Norm interactions") +  xlab("Chr") +
#   theme(legend.position = "none", axis.text=element_text(size=25), 
#         axis.title=element_text(size=25)) 
# dev.off() 
# 
# png(paste0(cond, "_inter_by_chr.png"), width = 1600, height = 1200)
# ggplot(inter.by.chr, aes(x = chr, fill = chr, y = n)) + 
#   geom_bar(stat="identity") + ylab("Interactions") + 
#   theme(legend.position = "none", axis.text=element_text(size=25), 
#         axis.title=element_text(size=25)) 
# dev.off() 
# write.table(genes.by.chr, file = paste0(cond, "_all_interactions_by_chr.tsv"), quote = F, row.names = F, sep = "\t")
