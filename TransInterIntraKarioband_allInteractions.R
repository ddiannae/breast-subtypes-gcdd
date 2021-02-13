setwd("/media/ddisk/transpipeline-data/breast-data/subtypes")
library(plyr)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggthemes)
annot <- read.delim("/media/ddisk/transpipeline-data/biomarts/Biomart_EnsemblG94_GRCh38_p12_karyotype.txt",
                    col.names = c("EnsemblID", "Chr", "Start", "End", "Band"), stringsAsFactors = F)
annot <- annot[!duplicated(annot$EnsemblID), c("EnsemblID", "Chr", "Band")]
conds <- c("basal", "luma", "lumb",  "her2", "healthy")
cond <- conds[1]
lapply(conds, function(cond) {
  MImatrix <-read_tsv(file = paste0("networks/aracne-cluster/1/", cond, "_normMI.tsv"))
  genes <- colnames(MImatrix)
  MImatrix$source <- genes
  MImatrix <- pivot_longer(MImatrix, genes, values_drop_na = TRUE,
                           names_to = "target", values_to = "MI")
  
  MImatrix <- MImatrix %>% filter(source != target)
  colnames(annot) <-  c("source", "source_chr", "source_band")
  MImatrix <- MImatrix %>%  inner_join(annot, by = "source")
  colnames(annot) <-  c("target", "target_chr", "target_band")
  MImatrix <- MImatrix %>%  inner_join(annot, by = "target")
  MImatrix$inter <- ifelse(MImatrix$source_chr == MImatrix$target_chr, F, T)
  MImatrix$inter_type <- ifelse(MImatrix$inter, "Trans", "Inter-Cytoband")
  MImatrix$inter_type <- ifelse(MImatrix$source_band == MImatrix$target_band & !MImatrix$inter,
                                "Intra-Cytoband", MImatrix$inter_type)
  
  ggplot(MImatrix, aes(x=MI, color=inter_type)) +
    geom_density() + 
    scale_color_manual(name = "Interaction", values =  c("#1696AC", "#73299e", "#C85200"))
  
  all.summ <- MImatrix %>% group_by(inter_type) %>% summarise(mean_mi = mean(MI), 
                                                              median_mi = median(MI), 
                                                              sd_mi = sd(MI), 
                                                              min_mi = min(MI),
                                                              max_mi = max(MI),
                                                              sum_mi = sum(MI),
                                                              n_interactions = n())
  intraMI <- MImatrix %>% filter(inter == F)
  intra.summ <- intraMI %>% group_by(source_chr, inter_type) %>%
    summarise(mean_mi = mean(MI), 
              median_mi = median(MI), 
              sd_mi = sd(MI), 
              min_mi = min(MI),
              max_mi = max(MI),
              sum_mi = sum(MI),
              n_interactions = n())
  write_tsv(all.summ, path = paste0("networks/aracne-cluster/1/MI_summary_", cond, ".tsv"))
  write_tsv(intra.summ, path = paste0("networks/aracne-cluster/1/MI_intra_summary_", cond, ".tsv"))
})

all.summ <- lapply(conds, function(cond) {
  summ <- read_tsv(file = paste0("networks/aracne-cluster/1/MI_summary_", cond, ".tsv"))
  #summ <- pivot_longer(summ, cols = colnames(summ)[2:ncol(summ)], names_to = "statistic")
  summ$cond <- cond
  return(summ)
})

all.summ <- bind_rows(all.summ)
all.summ$cond <- factor(all.summ$cond, levels = c("healthy", "luma", "lumb", "her2",  "basal"))
levels(all.summ$cond) <- c("Healthy", "LumA", "LumB", "Her2", "Basal")
all.summ$inter_type <- factor(all.summ$inter_type, 
                              levels = c("Intra-Cytoband", "Inter-Cytoband", "Trans"))

p <- ggplot(all.summ %>% select(mean_mi, sd_mi, cond, inter_type), 
            aes(x = inter_type, y = mean_mi,  fill = inter_type)) + 
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean_mi - sd_mi, ymax = mean_mi + sd_mi), width=.2,
                position = position_dodge(.9), color = "black") +
  theme_few(base_size = 20) +
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

