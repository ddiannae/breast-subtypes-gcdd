library(ggplot2)
library(ggthemes)
library(dplyr)
library(tidyr)
library(readr)

setwd("/media/ddisk/transpipeline-data/breast-data/subtypes")

subtypes.pal <- c("#C7C7C7", "#FF9E4A", "#67BF5C", "#729ECE", "#ED665D")

annot <- read.delim("/media/ddisk/transpipeline-data/biomarts/Biomart_EnsemblG94_GRCh38_p12_karyotype.txt",
                    col.names = c("ensemblID", "chr", "start", "end", "band"), stringsAsFactors = F)
annot <- annot[!duplicated(annot$ensemblID), c("ensemblID", "chr", "band")]
conds <- c("basal", "luma", "lumb",  "her2", "healthy")
cond <- conds[1]
all.nets <- lapply(conds, function(cond) {
  network <-read_tsv(file = paste0("networks/aracne-cluster/1e8/13319/", cond, "-13319-norm-mi.tsv"),
                     col_names = c("source", "MI", "target", "row_number"), skip = 1)
  colnames(annot) <-  c("source", "source_chr", "source_band")
  network <- network %>%  inner_join(annot, by = "source")
  colnames(annot) <-  c("target", "target_chr", "target_band")
  network <- network %>%  inner_join(annot, by = "target")
  network$inter <- ifelse(network$source_chr == network$target_chr, F, T)
  network$inter_type <- ifelse(network$inter, "Trans", "Inter-Cytoband")
  network$inter_type <- ifelse(network$source_band == network$target_band & !network$inter,
                                "Intra-Cytoband", network$inter_type)
  
  # ggplot(network, aes(x=MI, color=inter_type)) +
  #   geom_density() + 
  #   scale_color_manual(name = "Interaction", values =  c("#1696AC", "#73299e", "#C85200"))
  # 
  # all.summ <- network %>% group_by(inter_type) %>% summarise(mean_mi = mean(MI), 
  #                                                             median_mi = median(MI), 
  #                                                             sd_mi = sd(MI), 
  #                                                             min_mi = min(MI),
  #                                                             max_mi = max(MI),
  #                                                             sum_mi = sum(MI),
  #                                                             n_interactions = n())
  # intraMI <- network %>% filter(inter == F)
  # intra.summ <- intraMI %>% group_by(source_chr, inter_type) %>%
  #   summarise(mean_mi = mean(MI), 
  #             median_mi = median(MI), 
  #             sd_mi = sd(MI), 
  #             min_mi = min(MI),
  #             max_mi = max(MI),
  #             sum_mi = sum(MI),
  #             n_interactions = n())
  # write_tsv(all.summ, path = paste0("networks/aracne-cluster/1e8/13319/MI_summary_", cond, ".tsv"))
  # write_tsv(intra.summ, path = paste0("networks/aracne-cluster/1e8/13319/MI_intra_summary_", cond, ".tsv"))
  network$cond <- cond
  return(network)
})

### Los boxplots no dicen nada
# net8<- net[net$Target.Chr == "8" | net$Source == "8", 
#             c("Source", "Target", "MI", "InteractionType", "Cond")]
# net8$Chr <- "Chr 8"
# net17 <-  net[net$Target.Chr == "17" | net$Source == "17", 
#               c("Source", "Target", "MI", "InteractionType", "Cond")]
# net17$Chr <- "Chr 17"
# net.8.17 <- rbind(net8, net17)
# net.8.17$Chr <- factor(net.8.17$Chr, levels = c("Chr 8", "Chr 17"))
# net.8.17$InteractionType <- factor(net.8.17$InteractionType, levels = c("Intra-Cytoband", "Inter-Cytoband", "Trans"))
# 
# 
# p <- ggplot(net.8.17, aes(x = InteractionType, y = MI,  fill = InteractionType)) + 
#   geom_boxplot() +
#   theme_few() +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         strip.text.x = element_text(size = 16),
#         axis.title.x = element_blank()) + 
#   geom_text(stat="count", aes(label=paste0("(",..count.., ")")), 
#             y=min(net.8.17$MI) - 0.005, size = 3) +
#   facet_wrap(vars(Cond, Chr), nrow = 1) + 
#   ylab("Mutual Information") +
#   scale_fill_manual(name = "Interaction", values =  c("#1696AC", "#73299e", "#C85200"))
# 
# chrs <- as.character(c(1:22, "X"))
# chromosomes.pal <- c("#D909D1", "#0492EE", "#5DA0CB", "#106F35", "#5BD2AE", "#199F41", 
#                      "#FE0F43", "#00FFCC", "#F495C5", "#E1BF5D", "#5F166F", "#088ACA",
#                      "#41CFE0", "#0F0A71", "#FFFF99", "#B06645", "#800092", "#B925AE",
#                      "#B1B719", "#CB97E8", "#130B9E", "#E12B29", "#79A5B9")
# 
# names(chromosomes.pal) <- c("22","11","12","13","14","15","16","17","18","19","1" ,"2" ,"3" ,"4" ,"5" ,
#                             "6" ,"7" ,"X" ,"8" ,"9" ,"20","10","21")
# chromosomes.pal <- chromosomes.pal[chrs]
# chromosomes.pal <- chromosomes.pal[c("8", "17", "8", "17")]
# g <- ggplot_gtable(ggplot_build(p))
# strips <- which(grepl('strip-', g$layout$name))
# for (i in 1:4) {
#   k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[2]]$childrenOrder))
#   g$grobs[[strips[i]]]$grobs[[2]]$children[[k]]$gp$fill <- chromosomes.pal[i]
# }
# plot(g)

# cis.counts <- net %>% group_by(Cond, Source.Chr, InteractionType) %>% 
#   summarise(n = n()) %>% filter(InteractionType %in% c("Inter-Cytoband", "Intra-Cytoband"))
# 
# net %>% group_by(Cond, InteractionType) %>% tally()
# intravsinter <- cis.counts %>% spread(InteractionType, n) %>%
#   filter(`Intra-Cytoband` > `Inter-Cytoband`) 
# intravsinter %>% group_by(Cond) %>% tally()

all.nets <- bind_rows(all.nets)
all.nets$cond <- factor(all.nets$cond, levels = c("healthy", "luma", "lumb", "her2",  "basal"))
levels(all.nets$cond) <- c("Healthy", "LumA", "LumB", "Her2", "Basal")
all.nets$inter_type <- factor(all.nets$inter_type, levels = c("Intra-Cytoband", "Inter-Cytoband", "Trans"))

p <- ggplot(all.nets, aes(x = inter_type, y = MI,  fill = inter_type)) + 
  geom_boxplot() +
  theme_few(base_size = 30) +
  facet_wrap(~cond, nrow = 1) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(),
        axis.title.x = element_blank()) +
  ylab("Mutual Information") +
  scale_fill_manual(name = "Interaction", values =  c("#1696AC", "#73299e", "#C85200"))

g <- ggplot_gtable(ggplot_build(p))
strips <- which(grepl('strip-', g$layout$name))
for (i in 1:5) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- subtypes.pal[i]
}
plot(g)

net_counts <- all.nets %>% group_by(cond, inter_type) %>% tally()
p <- ggplot(net_counts, aes(x = inter_type, y = n,  fill = inter_type)) + 
  geom_bar(stat = "identity") +
  theme_few(base_size = 30) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(),
        axis.title.x = element_blank()) + 
  facet_wrap(~cond, nrow = 1) + 
  ylab("Number of interactions") +
  scale_fill_manual(name = "Interaction", 
                    values =  c("#1696AC", "#73299e", "#C85200")) +
  geom_text(stat="count", aes(label = n), 
            y=-250, size = 6) 

g <- ggplot_gtable(ggplot_build(p))
strips <- which(grepl('strip-', g$layout$name))
for (i in 1:5) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- subtypes.pal[i]
}
plot(g)

all.summ <- lapply(conds, function(cond) {
  summ <- read_tsv(file = paste0("networks/aracne-cluster/1e8/13319/MI_summary_", cond, ".tsv"))
  #summ <- pivot_longer(summ, cols = colnames(summ)[2:ncol(summ)], names_to = "statistic")
  summ$cond <- cond
  return(summ)
})

all.summ <- bind_rows(all.summ)
all.summ$cond <- factor(all.summ$cond, levels = c("healthy", "luma", "lumb", "her2",  "basal"))
levels(all.summ$cond) <- c("Healthy", "LumA", "LumB", "Her2", "Basal")
all.summ$inter_type <- factor(all.summ$inter_type, 
                              levels = c("Intra-Cytoband", "Inter-Cytoband", "Trans"))

p <- ggplot(all.summ %>% select(sum_mi, sd_mi, cond, inter_type), 
            aes(x = inter_type, y = sum_mi,  fill = inter_type)) + 
  geom_bar(stat = "identity") +
 # geom_errorbar(aes(ymin = mean_mi - sd_mi, ymax = mean_mi + sd_mi), width=.2,
#                position = position_dodge(.9), color = "black") +
  theme_few(base_size = 30) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 16),
        axis.title.x = element_blank()) + 
  facet_wrap(~cond, nrow = 1) + 
  ylab("MI") +
  scale_fill_manual(name = "Interaction", values =  c("#1696AC", "#73299e", "#C85200"))

subtypes.pal <- c("#C7C7C7", "#FF9E4A", "#67BF5C", "#729ECE", "#ED665D")
g <- ggplot_gtable(ggplot_build(p))
strips <- which(grepl('strip-', g$layout$name))
for (i in 1:5) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- subtypes.pal[i]
}
plot(g)

net.counts.by.chr <- net %>% filter(Inter == F)
net.counts.by.chr <- net.counts.by.chr %>% group_by(InteractionType, Source.Chr, Cond) %>% tally()
net.counts.by.chr %>% filter(Cond %in% c("basal", "healthy")) %>% arrange(Cond, Source.Chr) %>% print(n = Inf)
