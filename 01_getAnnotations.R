setwd("/home/diana/Workspace/INMEGEN/subtipos_mama")

annot <- read.delim("Ensembl93_GRCh38.12.txt", stringsAsFactors = T)
head(annot)
colnames(annot) <- c("Chr", "start", "end", "symbol")
## Remove null symbols
which(is.na(annot$symbol) | annot$symbol == "")
annot <- annot[!is.na(annot$symbol), ]
annot <- annot[annot$symbol != "", ]
## Remove non conventional chromosomes
levels(annot$Chr)
annot<-annot[annot$Chr%in%c(as.character(1:22), "X", "Y"),]
annot$Chr<-droplevels(annot$Chr)
dim(annot)
# [1] 37398     4
length(which(duplicated(annot$symbol)))
##[1] 36. No son tantos
head.dup <- head(which(duplicated(annot$symbol)))
symbols.dup <- annot[head.dup, "symbol"]
dup.symbols <- annot[annot$symbol %in% symbols.dup, ]
dup.symbols
annot <- data.frame(annot, stringsAsFactors = F)
dim(annot)
## Save annotation File
save(annot, file="annot.RData", compress="xz")


### Reading gene files
basal <-read.delim("Basal_genes.txt", col.names = c("symbol"), stringsAsFactors = F)
head(basal)
which(is.na(basal$symbol))
which(basal$symbol == "")

# All gene lists are the same. They got removed
# healthy <-read.delim("Healthy_genes.txt", col.names = c("symbol"), stringsAsFactors = F)
# head(healthy)
# her2 <-read.delim("Her2_genes.txt", col.names = c("symbol"), stringsAsFactors = F)
# lumA <-read.delim("LumA_genes.txt", col.names = c("symbol"), stringsAsFactors = F)
# lumB <-read.delim("LumB_genes.txt", col.names = c("symbol"), stringsAsFactors = F)
# normal <-read.delim("Normal_genes.txt", col.names = c("symbol"), stringsAsFactors = F)

mergeWithAnnotation <- function(annot, df, name) {
  merged.annot <- merge(df, annot, by="symbol", all = F, all.x = T)
  cat("Initial annot for ", name, " ")
  print(dim(df))
  # [1] 15281     1
  cat("Merged data frame dimensions ")
  print(dim(merged.annot))
  # [1] 15288     4
  not.in.biomart <- merged.annot[which(is.na(merged.annot$Chr)), "symbol"]
  cat("Length of not found symbols ")
  print(length(not.in.biomart))
  #[1] 459. Tiramos
  write.table(not.in.biomart, file = paste(name, "not_found.txt", sep = "_"), 
              row.names = F, quote = F, col.names = F)
  merged.annot <- merged.annot[!(merged.annot$symbol %in% not.in.biomart), ]
  cat("Merged data frame dimensions without not found ")
  print(dim(merged.annot))
  #[1] 14829     4
  ## No tenemos forma de decidir qué duplicados quitar.
  merged.annot <- merged.annot[-which(duplicated(merged.annot$symbol)), ]
  cat("Merged data frame dimensions without duplicates ")
  print(dim(merged.annot))
  #[1] 14810     4
  rownames(merged.annot) <- merged.annot$symbol
  return(merged.annot) 
}

genes.annot  <- mergeWithAnnotation(annot, basal, "genes")
save(genes.annot, file="genes.annot.RData", compress="xz")

#Initial annot for  genes  [1] 15281     1
#Merged data frame dimensions [1] 15288     4
#Length of not found symbols [1] 459
#Merged data frame dimensions without not found [1] 14829     4
#Merged data frame dimensions without duplicates [1] 14810     4

