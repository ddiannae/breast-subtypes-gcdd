library(dplyr)
library(tidyr)

setwd("/media/ddisk/transpipeline-data/breast-data/subtypes/")

conds <- c("healthy", "luma", "lumb", "basal", "her2")
chrs <- c(as.character(1:22), "X")
n.tests <- 1000

for(cond in conds) {
  cat("Working with condition", cond, "\n")
  interactions <- read.csv(file = paste0("networks/network-tables/", cond, "-norm-mi-interactions.tsv"),
                           sep = "\t", stringsAsFactors = F)
  vertices <- read.csv(file = paste0("networks/network-tables/", cond, "-vertices.tsv"),
                       sep = "\t", stringsAsFactors = F)
  interactions <- interactions %>% filter(InteractionType != "Trans")
  genes.in.cis <- union(interactions$Source, interactions$Target)
  
  all.chrs <- lapply(chrs, function(chr) {
    cat("\tWorking with chromosome", chr, "\n")
    genes.in.chr <- vertices %>% filter(Chr == chr, EnsemblID %in% genes.in.cis)
    interactions.in.chr <- bind_rows(interactions %>% semi_join(genes.in.chr, 
                                                                by = c("Source" = "EnsemblID")), 
                                     interactions %>% semi_join(genes.in.chr, 
                                                                by = c("Target" = "EnsemblID")))
    n.inter <- nrow(interactions.in.chr)
    all.intrak.vals <- parallel::mclapply(X = 1:n.tests, mc.cores = 6, function(k) {
      sinter <- lapply(1:n.inter, function(i){
        rows <- sample_n(genes.in.chr, 2, replace = F)
        return(tibble(Source = rows[1, "EnsemblID"],
                      Target = rows[2, "EnsemblID"], 
                      InteractionType = ifelse(rows[1, "Band"] == rows[2, "Band"], 
                                               "Intra-Cytoband", "Inter-Cytoband")))
      })
      sinter <- bind_rows(sinter)
      n.intrak <- sinter %>% filter(InteractionType == "Intra-Cytoband") %>% nrow()
      n.interk <- sinter %>% filter(InteractionType == "Inter-Cytoband") %>% nrow()
      return(n.intrak/n.interk)
    })
    all.intrak.vals <- unlist(all.intrak.vals)
    all.intrak.vals <- all.intrak.vals
    n.intrak <- interactions.in.chr %>% filter(InteractionType == "Intra-Cytoband" ) %>% nrow()
    n.interk <- interactions.in.chr %>% filter(InteractionType == "Inter-Cytoband" ) %>% nrow()
    return(list(chr = chr, random.intrak = all.intrak.vals, 
                real.intrak = n.intrak/n.interk, n = n.inter))
  })
  
  random.intraks <- bind_cols(lapply(all.chrs, "[[", "random.intrak"))
  colnames(random.intraks) <- paste0("chr", lapply(all.chrs, "[[", "chr"))
  write.table(random.intraks, file = paste0("networks/nullintrak/", cond, "-random-intrak-proportion.tsv"),
              sep = "\t", col.names = T, row.names = F, quote = F)
  real.intraks <-  bind_rows(lapply(all.chrs, "[", c(1, 3, 4)))
  write.table(real.intraks, file = paste0("networks/nullintrak/", cond, "-real-intrak-proportion.tsv"),
              sep = "\t", col.names = T, row.names = F, quote = F)
}


### now get the p-values and 
for(cond in conds) {
  random.intra <- read.csv(file = paste0("networks/nullintrak/", cond, "-random-intrak-proportion.tsv"), sep = "\t")
  real.intra <- read.csv(file = paste0("networks/nullintrak/", cond, "-real-intrak-proportion.tsv"), sep = "\t")
  pvals <- lapply(1:ncol(random.intra), function(j) {
    random.vals <- random.intra[[j]]/real.intra[j, 3]
    real.mu <- real.intra[j,2]/real.intra[j,3]
    random.sd <- sd(random.vals)
    random.mean <- mean(random.vals)
    zscore <- (real.mu - random.mean)/random.sd
    ttest <- t.test(x = random.vals, mu =  real.mu)
    return(list(pval = ttest$p.value, zscore = zscore))
  })
  real.intra$pval <- unlist(lapply(pvals, "[[", "pval"))
  real.intra$zscore <- unlist(lapply(pvals, "[[", "zscore"))
  write.table(real.intra, file = paste0("networks/nullintrak/", cond, "-real-intrak-proportion.tsv"),
              sep = "\t", col.names = T, row.names = F, quote = F)
}

## add number of chromosomes

for(cond in conds) {
  vertices <- read.csv(file = paste0("networks/network-tables/", cond, "-vertices.tsv"),
                       sep = "\t", stringsAsFactors = F)
  interactions <- read.csv(file = paste0("networks/network-tables/", cond, "-norm-mi-interactions.tsv"),
                           sep = "\t", stringsAsFactors = F)
  genes.in.cis <- union(interactions$Source, interactions$Target)
  real.intra <- read.csv(file = paste0("networks/nullintrak/", cond, "-real-intrak-proportion.tsv"), 
                         sep = "\t", stringsAsFactors = F)
  all.chrs <- lapply(chrs, function(chr) {
    genes.in.chr <- vertices %>% filter(Chr == chr, EnsemblID %in% genes.in.cis) %>% nrow()
    return(list(genes = genes.in.chr, chr = chr))
  })
  all.chrs <- bind_rows(all.chrs)
  real.intra <- real.intra %>% inner_join(all.chrs, by = "chr")
  write.table(real.intra, file = paste0("networks/nullintrak/", cond, "-real-intrak-proportion.tsv"),
              sep = "\t", col.names = T, row.names = F, quote = F)
}

            

