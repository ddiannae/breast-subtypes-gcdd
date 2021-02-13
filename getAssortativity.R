library(igraph)
library(dplyr)
setwd("/mnt/ddisk/transpipeline-data/breast-data/subtypes")

conds <- c("luma", "lumb", "basal", "her2")
lapply(conds, function(cond){
    
  net <- read.delim(file = paste0("networks/chucho-tables/", cond, "-interactions.tsv"), 
                                         sep = "\t", stringsAsFactors = F)
  vertices <- read.delim(file = paste0("networks/chucho-tables/", cond, "-vertices.tsv"), 
                         sep = "\t", stringsAsFactors = F)
  membership <- read.delim(file = paste0("networks/enrichment/", cond, "-membership.tsv"), 
                            sep = "\t", stringsAsFactors = F, col.names = c("com_id", "EnsemblID"))
  vertices <- vertices %>% inner_join(membership, by = "EnsemblID" )
  vertices$up <- ifelse(vertices$LFC >0, T, F)
  network <- graph_from_data_frame(net, directed=F, vertices = vertices)
  assortativity.nominal(network, types = as.numeric(as.factor(V(network)$up)), directed = F)
  
  all.aval <- lapply(unique(vertices$com_id), function(cid){
    sgraph <- induced.subgraph(network, V(network)[V(network)$com_id == cid])
    return(list(com_id = cid, assortativity_deg = assortativity.nominal(sgraph, types = as.numeric(as.factor(V(sgraph)$up)), directed = F)))
  })
  all.aval <- bind_rows(all.aval)
  all.aval[is.na(all.aval$assortativity_deg), "assortativity_deg"]  <- as.numeric(1)
  
  clusters.info <- read.delim(file = paste0("networks/enrichment/", cond, "-clusters.tsv"), 
                                            sep = "\t", stringsAsFactors = F)
  clusters.info <- clusters.info %>% inner_join(all.aval, by = "com_id")
  write.table(clusters.info, file = paste0("networks/enrichment/", cond, "-clusters.tsv"), sep ="\t",  row.names = F, quote = F)

})
