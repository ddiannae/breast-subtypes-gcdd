setwd("/mnt/antares/mnt/ddisk/transpipeline-data/breast-data/subtypes/networks/")

annot <- read.delim(file = paste0("../../../biomarts/Biomart_EnsemblG94_GRCh38_p12_NCBI.txt"), 
                  sep = "\t", stringsAsFactors = F)
annot <- annot[, c("Gene.stable.ID", "HGNC.symbol")]
colnames(annot) <- c("ensemblID", "symbol")
conds <- c("healthy")

all.clusters <- lapply(conds, function(cond) {
  enrichments <- read.delim(file = paste0("enrichment/", cond, "_enrichment.tsv"), 
                    sep = "\t", stringsAsFactors = F)
  clusters <- read.delim(file = paste0("enrichment/", cond, "_clusters.tsv"), 
                            sep = "\t", stringsAsFactors = F)
  clusters$cond <- cond
  return(merge(enrichments, clusters, by = "com_id"))
})

all.clusters <- plyr::ldply(all.clusters)
all.clusters <- all.clusters[, c("com_id", "name", "observed_hits",  "chr", "pg_gene", "type", "cond")]
all.clusters <- merge(all.clusters, annot, by.x = "pg_gene", by.y = "ensemblID")

#all.clusters$cond <- factor(all.clusters$cond, levels = c("luma", "lumb", "her2", "basal"))
#levels(all.clusters$cond) <- c("LumA", "LumB", "Her2", "Basal")
all.clusters$cond <- "Healthy"
all.clusters$type <- factor(all.clusters$type, levels = c("inter", "intra"))
levels(all.clusters$type) <- c("Trans", "Cis")
#write.table(all.clusters[all.clusters$type == "inter", ], file = paste0("enrichment/inter_cluster.tsv"), sep = "\t", quote = F, row.names = F)
#write.table(all.clusters[all.clusters$type == "intra", ], file = paste0("enrichment/intra_clusters.tsv"), sep = "\t", quote = F, row.names = F)
write.table(all.clusters, file = paste0("enrichment/healthy.tsv"), sep = "\t", quote = F, row.names = F)

      