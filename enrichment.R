library(HTSanalyzeR)
library(org.Hs.eg.db)
library(GO.db)

getGODataset <- function(dataSet) {
  toEnsemblID <- function(GO_list){
    return(mapIds(org.Hs.eg.db, GO_list , 'ENSEMBL', 'ENTREZID'))
  }
  GS.GO <- GOGeneSets(species="Hs", ontologies=c(dataSet))                    
  GS.GO <- lapply(GS.GO, toEnsemblID)
}

## It should work for components and communities, so we're working with the membership vector of
## both objects.
## https://igraph.org/r/doc/communities.html
## https://igraph.org/r/doc/components.html

getGOEnrichments <- function(cmembership, GS.GO, universe) {
  enrichments <- lapply(unique(cmembership), function(idc){
    mem <- names(cmembership[cmembership == idc])
    if(length(mem) > 5){
      results.GS.GO <- multiHyperGeoTest(collectionOfGeneSets = GS.GO,
                                         universe = universe,
                                         hits = mem , minGeneSetSize = 10,
                                         pAdjustMethod = "BH")
      colnames(results.GS.GO) <- c("universe_size", "gene_set_size", "total_hits", "expected_hits",
                                      "observed_hits", "p_val", "adj_pval")
      
      results.GS.GO <- as.data.frame(results.GS.GO)
      results.GS.GO <- results.GS.GO[results.GS.GO$adj_pval < 0.001,]
      results.GS.GO <- results.GS.GO[, -1]
      
      if(nrow(results.GS.GO) > 0) {
        results.GS.GO$go_term <- row.names(results.GS.GO)
        rownames(results.GS.GO) <- Term(rownames(results.GS.GO))
        results.GS.GO$name <- rownames(results.GS.GO)
        results.GS.GO$com_id <- idc
        results.GS.GO$gene_hits <- lapply(lapply(GS.GO[results.GS.GO$go_term], intersect,  mem), paste, collapse = ", ")
        return(results.GS.GO)
      }
    }
    return(NULL)
  })
  enrichments <- plyr::compact(enrichments)
  enrichments <- plyr::ldply(enrichments)
  return(enrichments)
}
