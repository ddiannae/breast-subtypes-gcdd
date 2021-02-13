library(Hmisc)
setwd("/home/dianae/Workspace/breast-cancer-subtypes/")
source("enrichment.R")

getComInfo <- function(cmembership, network){
  comp.info <- lapply(unique(cmembership), function(idc){
    mem <- names(cmembership[cmembership == idc])
    typec <- "intra"
    com <- induced.subgraph(network, mem)
    if(mean(E(com)$Inter) > 0.5) {
      typec <- "inter"
    }
    prs <- page.rank(com)
    chrs <- table(V(com)$Chr)
    return(data.frame(type = typec, com_id = idc,
                      pg_gene = names(which.max(prs$vector))[1], 
                      chr = names(which.max(chrs))[1], length = length(V(com))))
  })
  comp.info <- plyr::compact(comp.info)
  comp.info <- plyr::ldply(comp.info)
  return(comp.info)
}

getAllVertices <- function(network) {
  all.vertices <- network[, c("Target", "Target.Start", "Target.End", "Target.Chr")]
  colnames(all.vertices) <- c("Source", "Source.Start", "Source.End", "Source.Chr")
  all.vertices  <- rbind(all.vertices, network[, c("Source", "Source.Start", "Source.End", "Source.Chr")])
  colnames(all.vertices) <- c("EnsemblID", "Start", "End", "Chr")
  return(all.vertices[!duplicated(all.vertices$EnsemblID), ])
}

setwd("/mnt/ddisk/transpipeline-data/breast-data/subtypes")
#GO.BP <- getGODataset("BP")
#save(GO.BP, file = "rdata/GS_BP_12062019.RData", compress = "xz")
load("rdata/GS_BP_12062019.RData")

cond <- "luma"
net <- read.delim(file = paste0("networks/aracne-cluster/1e8/13319/", cond, "-13319-annotated.tsv"), 
                                       sep = "\t", stringsAsFactors = F)

net <- graph_from_data_frame(net[, c("Source", "Target", "MI", "Inter")], 
                                 directed=F, vertices = getAllVertices(net))

net.components <- components(net)
net.com.mem <- net.components$membership[net.components$membership !=  which.max(net.components$csize)]
luma.components.enrich <- getGOEnrichments(net.com.mem, GO.BP, V(net)$name)
luma.components.info <- getComInfo(net.com.mem, net)
luma.components.enrich$com_id <- paste0("component_", luma.components.enrich$com_id)
luma.components.info$com_id <- paste0("component_", luma.components.info$com_id)

big.component <-  names(which(net.components$membership == which.max(net.components$csize)))
big.component <- induced.subgraph(net, big.component)
big.comp.comm <- cluster_infomap(graph = big.component, nb.trials = 10)
names(big.comp.comm$membership) <- big.comp.comm$names
luma.commmunities.enrich <- getGOEnrichments(big.comp.comm$membership, GO.BP, V(net)$name)
luma.communities.info <- getComInfo(big.comp.comm$membership, net)
luma.commmunities.enrich$com_id <- paste0("community_", luma.commmunities.enrich$com_id)
luma.communities.info$com_id <- paste0("community_", luma.communities.info$com_id)

luma.erichments <- rbind(luma.components.enrich, luma.commmunities.enrich)
luma.erichments$name <-  capitalize(luma.erichments$name)
luma.erichments$gene_hits <- unlist(luma.erichments$gene_hits)

write.table(luma.erichments, file = paste0("networks/enrichment/luma_enrichment.tsv"), sep ="\t",  row.names = F, quote = F)
luma.clusters <- rbind(luma.components.info, luma.communities.info)
write.table(luma.clusters, file = paste0("networks/enrichment/luma_clusters.tsv"), sep ="\t",  row.names = F, quote = F)

luma.mem <- data.frame(community = paste0("component_", net.com.mem), gene = names(net.com.mem))
luma.mem <- rbind(luma.mem, data.frame(community = paste0("community_", big.comp.comm$membership), gene = names(big.comp.comm$membership) ))
write.table(luma.mem, file = paste0("networks/enrichment/luma_membership.tsv"), sep ="\t",  row.names = F, quote = F)

conds <- c("her2", "basal")
for(cond in conds){
  net <- read.delim(file = paste0("networks/aracne-cluster/1e8/13319/", cond, "-13319-annotated.tsv"), 
                    sep = "\t", stringsAsFactors = F)
  
  net <- graph_from_data_frame(net[, c("Source", "Target", "MI", "Inter")], 
                               directed=F, vertices = getAllVertices(net))
  net.components <- components(net)
  enrichments <- getGOEnrichments(net.components$membership, GO.BP, V(net)$name)
  enrichments$name <-  capitalize(enrichments$name)
  enrichments$gene_hits <- unlist(enrichments$gene_hits)
  write.table(enrichments, file = paste0("networks/enrichment/", cond, "-enrichment.tsv"), sep ="\t",  row.names = F, quote = F)
  
  components <- getComInfo(net.components$membership, net)
  write.table(components, file = paste0("networks/enrichment/", cond, "-clusters.tsv"), sep ="\t",  row.names = F, quote = F)
  
  membership <- data.frame(community = net.components$membership, gene = names(net.components$membership))
  write.table(membership, file = paste0("networks/enrichment/", cond, "-membership.tsv"), sep ="\t",  row.names = F, quote = F)
}

cond <- "healthy"
net <- read.delim(file = paste0("networks/aracne-cluster/1e8/13319/", cond, "-13319-annotated.tsv"), 
                  sep = "\t", stringsAsFactors = F)

net <- graph_from_data_frame(net[, c("Source", "Target", "MI", "Inter")], 
                             directed=F, vertices = getAllVertices(net))

communities <- cluster_infomap(graph = net, nb.trials = 10)
names(communities$membership) <- communities$names
healthy.commmunities.enrich <- getGOEnrichments(communities$membership, GO.BP, V(net)$name)
healthy.commmunities.enrich$name <-  capitalize(healthy.commmunities.enrich$name)
healthy.commmunities.enrich$gene_hits <- unlist(healthy.commmunities.enrich$gene_hits)
write.table(healthy.commmunities.enrich, file = paste0("networks/enrichment/", cond, "_enrichment.tsv"), sep ="\t",  row.names = F, quote = F)

healthy.communities.info <- getComInfo(communities$membership, net)
write.table(healthy.communities.info, file = paste0("networks/enrichment/", cond, "_clusters.tsv"), sep ="\t",  row.names = F, quote = F)

membership <- data.frame(community = communities$membership, gene = names(communities$membership))
write.table(membership, file = paste0("networks/enrichment/", cond, "_membership.tsv"), sep ="\t",  row.names = F, quote = F)
