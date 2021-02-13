setwd("/mnt/ddisk/transpipeline-data/breast-data/subtypes")
library(dplyr)

annot <- read.delim("/mnt/ddisk/transpipeline-data/biomarts/Biomart_EnsemblG94_GRCh38_p12_karyotype.txt",
                    col.names = c("EnsemblID", "Chr", "Start", "End", "Band"), stringsAsFactors = F)
annot <- annot[!duplicated(annot$EnsemblID), ]
conds <- c("basal", "luma", "lumb",  "her2", "healthy")
net <- lapply(conds, function(cond){
  network <- read.delim(paste0("networks/aracne-cluster/1e8/13319/", cond, "-13319-norm-mi.tsv"), 
                        stringsAsFactors = F)
  network$Cond <- cond
  return(network)
})
net <- bind_rows(net)
colnames(annot) <-  c("Source", "Source.Chr", "Source.Start", "Source.End", "Source.Band")
net <- merge(annot, net)
colnames(annot) <-  c("Target", "Target.Chr", "Target.Start", "Target.End",  "Target.Band")
net <- merge(annot, net)
net$Inter <- ifelse(net$Source.Chr == net$Target.Chr, F, T)
net$InteractionType <- ifelse(net$Inter, "Trans", "Inter-Cytoband")
net$InteractionType <- ifelse(net$Source.Band == net$Target.Band & !net$Inter, "Intra-Cytoband", net$InteractionType)

annot.symbol <- read.delim("/mnt/ddisk/transpipeline-data/biomarts/Biomart_EnsemblG94_GRCh38_p12.txt",
                           colClasses = c("character",  "NULL",  "NULL",  "NULL", "NULL", "NULL",  "character" ), 
                           stringsAsFactors = F, 
                           col.names = c("EnsemblID",  "NULL",  "NULL",  "NULL", "NULL", "NULL",  "Symbol"))

annot.symbol <- annot.symbol[!duplicated(annot.symbol$EnsemblID), ]
targets <- net %>% select(Target, Target.Chr, Target.Start, Target.End, Target.Band)
sources <- net %>% select(Source, Source.Chr, Source.Start, Source.End, Source.Band)
colnames(targets) <- c("EnsemblID", "Chr", "Start", "End", "Band")
colnames(sources)  <- c("EnsemblID", "Chr", "Start", "End", "Band")
vertices <- rbind(targets, sources)
vertices <- vertices[!duplicated(vertices$EnsemblID), ]
vertices <- vertices %>% inner_join(annot.symbol, by = "EnsemblID")

housekeeping <- read.delim("/mnt/ddisk/transpipeline-data/annotations/housekeeping.txt", header = F, 
                           colClasses = c("character", "NULL"))
housekeeping <- housekeeping[, 1, drop= T]
housekeeping <- stringr::str_trim(housekeeping)
housekeeping <- tolower(housekeeping)
vertices$isHK <- ifelse(tolower(vertices$Symbol) %in% housekeeping, T, F)
vertices %>% filter(isHK == T)

tfs <- read.delim("/mnt/ddisk/transpipeline-data/annotations/tfs_2018.txt", header = F, col.names = c("EnsemblID"))
vertices$isTF <- ifelse(vertices$EnsemblID %in% tfs$EnsemblID, T, F)

kariotype <- read.delim("/mnt/ddisk/transpipeline-data/annotations/chromosome.band.hg38.txt", header = T, stringsAsFactors = F)
names(kariotype)[1] <- "chr"
tbands <- kariotype %>% mutate(arm = unlist(lapply(as.vector(strsplit(name, split = "")), "[[", 1)))
tbands <- tbands %>% group_by(chr, arm) %>% summarise(t = max(name))
tbands <- tbands %>% group_by(chr) %>% summarise(t1 = min(t), t2 = max(t))
tbands$chr <- stringr::str_replace(tbands$chr, "chr", "")
vertices <- vertices %>% inner_join(tbands, by = c("Chr" = "chr"))
vertices$isExtreme <- vertices$Band == vertices$t1 | vertices$Band == vertices$t2
vertices$isExtreme <- as.integer(vertices$isExtreme)
vertices$isHK <- as.integer(vertices$isHK)
vertices$isTF <- as.integer(vertices$isTF)

conds <- c("luma", "lumb", "her2", "basal")
names(conds) <- c("LumA", "LumB", "Her2", "Basal")

for(i in 1:4){
  net.cond <- net%>% filter(Cond == conds[i]) %>% select(Source, Target, MI, InteractionType) %>%
    arrange(desc(MI)) 
  net.cond$RowNumber <- 1:nrow(net.cond)
  vertices.cond <- vertices %>% 
    filter(vertices$EnsemblID %in% net.cond$Source | vertices$EnsemblID %in% net.cond$Target) 
  deg.cond <- read.delim(paste0("deg/ebayes-", names(conds)[i], ".tsv"),  header = T, stringsAsFactors = F, 
                         colClasses = c("character", "NULL", "numeric", "numeric", "NULL", "NULL"), 
                         col.names = c( "EnsemblID", "NULL",  "LFC", "P-val","NULL", "NULL"))
  vertices.cond <- vertices.cond %>% inner_join(deg.cond, by = "EnsemblID")
  vertices.cond <- vertices.cond %>% select(-t1, -t2)
  write.table(vertices.cond, file = paste0("networks/network-tables/", conds[i], "-vertices.tsv") , 
            quote = F, row.names = F, col.names = T, sep = "\t")
 write.table(net.cond, file = paste0("networks/network-tables/", conds[i], "-norm-mi-interactions.tsv"),
            quote = F, row.names = F, col.names = T, sep = "\t")
}

net.cond <- net%>% filter(Cond == "healthy") %>% select(Source, Target, MI, InteractionType) %>%
  arrange(desc(MI)) 
net.cond$RowNumber <- 1:nrow(net.cond)
vertices.cond <- vertices %>% 
  filter(vertices$EnsemblID %in% net.cond$Source | vertices$EnsemblID %in% net.cond$Target) 
vertices.cond <- vertices.cond %>% select(-t1, -t2)
write.table(vertices.cond, file = paste0("networks/network-tables/healthy-vertices.tsv") , 
            quote = F, row.names = F, col.names = T, sep = "\t")
write.table(net.cond, file = paste0("networks/network-tables/healthy-norm-mi-interactions.tsv"),
            quote = F, row.names = F, col.names = T, sep = "\t")

