require(pbcmc)
library("BiocParallel")
library("ggplot2")
library("reshape2")
library ("png") 
library("grid")
library("gridExtra")
library("parallel")

library("NOISeq")

setwd("/mnt/ddisk/transpipeline-data/breast-data/subtypes")
load("/mnt/ddisk/transpipeline-data/breast-data/rdata/Length.full_GC.full_Between.tmm_Norm_cpm10.RData")

#### Get Gene Symbols
annot <- read.delim("/mnt/ddisk/transpipeline-data/biomarts/Biomart_EnsemblG94_GRCh38_p12_NCBI.txt", header = T, stringsAsFactors = F)
genes <- merge(annot, norm.data.cpm10$Annot, by.x = "Gene.stable.ID", by.y = "EnsemblID")
genes <- genes[!duplicated(genes$Gene.stable.ID), c("Gene.stable.ID", "HGNC.symbol", "NCBI.gene.ID")]
names(genes)<-c("probe", "NCBI.gene.symbol", "EntrezGene.ID")
rownames(genes) <- genes$probe

#### Get File Info
norm.data.cpm10$Targets$FileName = paste(substring(norm.data.cpm10$Targets$File,20, length(norm.data.cpm10$Targets$File)), "gz", sep=".")
filenames.cancer <- norm.data.cpm10$Targets[norm.data.cpm10$Targets$Group == "T", "FileName"]
#write.table(filenames.cancer, file = "/mnt/ddisk/transpipeline-data/annotations/filenames-breast-cancer.2019-01-09.tsv", sep="\t", quote=T, row.names=FALSE)
filenames.healthy <- norm.data.cpm10$Targets[norm.data.cpm10$Targets$Group == "N", "FileName"]
norm.data.cpm10$Targets[norm.data.cpm10$Targets$Group == "N", "FileName"]substring(filenames.healthy, 2, length(filenames.healthy))
#write.table(filenames.healthy, file = "/mnt/ddisk/transpipeline-data/annotations/filenames-breast-healthy.2019-01-09.tsv", sep="\t", quote=T, row.names=FALSE)

exp.caseids = read.csv("/mnt/ddisk/transpipeline-data/breast-data-cnvs/exp-caseids-cancer.csv", header = T, col.names = c("FileName", "SubmitterID", "CaseID", "FileID"), stringsAsFactors = F)
fileNames.cancer <- data.subtypes.arsyn$Targets["Group" == "T", "fileName"]

M <- norm.data.cpm10$M[genes$probe, norm.data.cpm10$Targets$Group == "T"]
pam50obj <- PAM50(exprs=log2(1+ M), annotation=genes)

pam50obj <- filtrate(pam50obj, verbose = T)
pam50obj <- classify(pam50obj, std="median", verbose = T)
head(classification(pam50obj))
parameters(pam50obj)

### Troubles with the original permutate function
pam50obj <- permutate(pam50obj, BPPARAM = MulticoreParam(workers = 7, progressbar = TRUE),  pCutoff = 0.03,  corCutoff = 0.08)
# pCutoff= 0.01 -> 34 samples Assigned
# pCutoff= 0.05 -> 134 samples Assigned
# pCutoff= 0.1, corCutoff = 0.05 -> 207 samples Assigned
# pCutoff = 0.1,  corCutoff = 0.01 -> 225 samples Assigned

# pCutoff = 0.25,  corCutoff = 0.01 -> [1] 415    
subtypes <- permutation(pam50obj)$subtype
dim(subtypes[subtypes$Permuted == "Assigned", ])
dim(subtypes[subtypes$Permuted == "Assigned" & subtypes$Subtype == "LumA", ])
# [1] 197   5
dim(subtypes[subtypes$Permuted == "Assigned" & subtypes$Subtype == "LumB", ])
# [1] 171  5
dim(subtypes[subtypes$Permuted == "Assigned" & subtypes$Subtype == "Basal", ])
# [1] 216  5
dim(subtypes[subtypes$Permuted == "Assigned" & subtypes$Subtype == "Her2", ])
# [1] 88
dim(subtypes[subtypes$Permuted == "Ambiguous", ])
# [1] 116  5

save(pam50obj, file = "pam50obj_0.30p_0.01c.RData", compress="xz")
assigned.subtypes <- subtypes[subtypes$Permuted == "Assigned", ]
assigned.subtypes$ID <- rownames(assigned.subtypes)
assigned.subtypes <- assigned.subtypes[, c("ID", "Subtype")]
healthy.subtypes <- norm.data.cpm10$Targets[norm.data.cpm10$Targets$Group == "N", "ID", drop = FALSE]
healthy.subtypes$Subtype <- "Healthy"
assigned.subtypes <- rbind(assigned.subtypes, healthy.subtypes)
#load( "pam50obj_0.30p_0.01c.RData")
targets.merged <- merge(exp.caseids, norm.data.cpm10$Targets, by = "FileName" )
targets.merged <- merge(assigned.subtypes, targets.merged, by = "ID" )
rownames(targets.merged) <- targets.merged$ID
targets.merged <- targets.merged[, c("ID","Subtype","FileName","SubmitterID","CaseID")]
targets.merged <- targets.merged[targets.merged$Subtype != "Normal", ]
colnames(targets.merged)

write.table(targets.merged, file = "caseid-filename-subtype.tsv", sep="\t", row.names = F, col.names = T, quote = F)

PLOTSDIR <- "/mnt/ddisk/transpipeline-data/breast-data/subtypes/plots"
#targets.merged <- read.delim("caseid-filename-subtype.tsv", sep="\t", header = T )

## Filtrar solo tumores con subtipos
targets.merged <- merge(norm.data.cpm10$Targets, targets.merged, by= "ID", all.x = T)
rownames(targets.merged) <- targets.merged$ID
targets.merged[targets.merged$Group == "N", "Subtype"] <- "Healthy"
targets.merged <- targets.merged[!is.na(targets.merged$Subtype), ]
data.subtypes <- list(M=norm.data.cpm10$M[, rownames(targets.merged)], 
                      Annot=norm.data.cpm10$Annot[, ], 
                      Targets=targets.merged)

### PCA with cpm (filter but no ARSyn)
pca.results <- pca.results <- PCA.GENES(t(log2(1 + data.subtypes$M)))
  
## Variance explained by each component
pdf(file=paste(PLOTSDIR, "subtypes_corrected_cpm10_PCAVariance.pdf", sep="/"), 
      width = 4*2, height = 4*2)
barplot(pca.results$var.exp[,1], xlab = "PC", ylab = "Explained variance")
dev.off()
cat("PCA variance norm plot generated.\n")
  
## Loading plot
pdf(file=paste(PLOTSDIR, "subtypes_corrected_cpm10_PCALoading.pdf", sep="/"), 
      width = 4*2, height = 4*2)
  plot(pca.results$loadings[,1:2], col = 1, pch = 20, cex = 0.5,
       xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
       ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
       main = "PCA loadings",
       xlim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1),
       ylim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1))  
  dev.off()
  cat("PCA loading norm plot generated.\n")
  
  ## Score plot
  mycol <- as.character(data.subtypes$Targets$Subtype)
  mycol[mycol == 'Basal'] <- "black"
  mycol[mycol == 'Her2'] <- "red2"
  mycol[mycol == 'LumA'] <- "green"
  mycol[mycol == 'LumB'] <- "blue"
  mycol[mycol == 'Healthy'] <- "pink"
  
  pdf(file=paste(PLOTSDIR,  "subtypes_corrected_cpm10_PCAScore.pdf", sep="/"), 
      width = 5*2, height = 5)
  par(mfrow = c(1,2))
  
  # PC1 & PC2
  rango <- diff(range(pca.results$scores[,1:2]))
  plot(pca.results$scores[,1:2], col = "white",
       xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
       ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
       main = "PCA scores",
       xlim = range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1),
       ylim = range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1))
  points(pca.results$scores[,1], pca.results$scores[,2], col = mycol, cex = 1.5)  
  legend("topright", c("Basal", "Her2", "LumA", "LumB", "Healthy"), col = c("black", "red2", "green", "blue", "pink"), ncol = 5, pch = 1)
  
  # PC1 & PC3
  rango2 = diff(range(pca.results$scores[,c(1,3)]))
  plot(pca.results$scores[,c(1,3)], col = "white",
       xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
       ylab = paste("PC 3 ", round(pca.results$var.exp[3,1]*100,0), "%", sep = ""),
       main = "PCA scores",
       xlim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1),
       ylim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1))
  points(pca.results$scores[,1], pca.results$scores[,3], col = mycol, cex = 1.5)
  legend("topright", c("Basal", "Her2", "LumA", "LumB", "Healthy"), col = c("black", "red2", "green", "blue", "pink"), ncol = 5, pch = 1)
  dev.off()
  cat("PCA scores norm plot generated.\n")

## ARSyN to reduce batch effect
##########################################
  data.subtypes$Targets$Subtype <- as.factor(data.subtypes$Targets$Subtype)
  mydata <- NOISeq::readData(
    data =data.subtypes$M, 
    factors = data.subtypes$Targets[, "Subtype",drop = FALSE])
  
  cat("Performing ARSyN for batch correction")
  myARSyN <- ARSyNseq(mydata, norm = "n", logtransf = F)
  pca.dat <- dat(myARSyN, type = "PCA", logtransf = F)
  pca.results <- pca.dat@dat$result
  
  ## Variance explained by each component
  pdf(file=paste(PLOTSDIR, "subtypes_corrected_cpm10_arsyn_PCAVariance.pdf", sep="/"),
      width = 4*2, height = 4*2)
  barplot(pca.results$var.exp[,1], xlab = "PC", ylab = "Explained variance", ylim = c(0,0.4))
  dev.off()
  cat("PCA variance arsyn plot generated.\n")
  
  
  ## Loading plot
  pdf(file=paste(PLOTSDIR, "subtypes_corrected_cpm10_arsyn_PCALoading.pdf", sep="/"), 
      width = 4*2, height = 4*2)
  plot(pca.results$loadings[,1:2], col = 1, pch = 20, cex = 0.5,
       xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
       ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
       main = "PCA loadings",
       xlim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1),
       ylim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1))  
  dev.off()
  cat("PCA loading arsyn plot generated.\n")
  
  ## Score plot
  mycol <- as.character(data.subtypes$Targets$Subtype)
  mycol[mycol == 'Basal'] <- "black"
  mycol[mycol == 'Her2'] <- "red2"
  mycol[mycol == 'LumA'] <- "green"
  mycol[mycol == 'LumB'] <- "blue"
  mycol[mycol == 'Healthy'] <- "pink"
  
  pdf(file=paste(PLOTSDIR, "subtypes_corrected_cpm10_arsyn_PCAScoreARSyN.pdf", sep="/"), 
      width = 5*2, height = 5)
  par(mfrow = c(1,2))
  
  # PC1 & PC2
  rango <- diff(range(pca.results$scores[,1:2]))
  plot(pca.results$scores[,1:2], col = "white",
       xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
       ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
       main = "PCA scores",
       xlim = range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1),
       ylim = range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1))
  points(pca.results$scores[,1], pca.results$scores[,2], col = mycol, cex = 1.5)  
  legend("topright", c("Basal", "Her2", "LumA", "LumB", "Healthy"), col = c("black", "red2", "green", "blue", "pink"), ncol = 5, pch = 1)
  
  # PC1 & PC3
  rango2 = diff(range(pca.results$scores[,c(1,3)]))
  plot(pca.results$scores[,c(1,3)], col = "white",
       xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
       ylab = paste("PC 3 ", round(pca.results$var.exp[3,1]*100,0), "%", sep = ""),
       main = "PCA scores",
       xlim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1),
       ylim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1))
  points(pca.results$scores[,1], pca.results$scores[,3], col = mycol, cex = 1.5)
  legend("topright", c("Basal", "Her2", "LumA", "LumB", "Healthy"), col = c("black", "red2", "green", "blue", "pink"), ncol = 5, pch = 1)
  dev.off()
  cat("PCA scores arsyn plot generated.\n")
  
  pl<-ggplot(data=melt(log(assayData(myARSyN)$exprs+1)), aes(x=value, group=Var2, colour=Var2))+geom_density(show.legend = F)
  png(file=paste(PLOTSDIR, "corrected_cpm10_arsyn_densitylog.png",  sep="/"),
      width = 2048, height = 1024, pointsize = 20)
  print(pl)
  dev.off()
  
  cat('ARSyN data. Final dimensions: ', paste(dim(assayData(myARSyN)$exprs), collapse=", "), '.\n')
  
  ##Saving everything
  data.subtypes.arsyn <- list(M = assayData(myARSyN)$exprs, Annot = data.subtypes$Annot, 
                                Targets = data.subtypes$Targets)
  
  stopifnot(nrow(data.subtypes.arsyn$M) == nrow(data.subtypes.arsyn$Annot))
  stopifnot(all(row.names(data.subtypes.arsyn$M) == row.names(data.subtypes.arsyn$Annot)))
  
  
  save(data.subtypes.arsyn, file="data_subtypes_arsyn.RData", compress="xz")
  
  cat("Generating data matrices with arsyn for Aracne\n")
  ## Data matrices for Aracne
  ## ALL = healthy | cancer
  M <- as.data.frame(data.subtypes.arsyn$M)
  M <- cbind(gene=as.character(data.subtypes.arsyn$Annot$EnsemblID), M)

  #normal samples
  basal <- as.data.frame(data.subtypes.arsyn$M[,data.subtypes.arsyn$Targets$Subtype == "Basal"])
  basal <- cbind(gene=as.character(data.subtypes.arsyn$Annot$EnsemblID), basal)
  
  luma <- as.data.frame(data.subtypes.arsyn$M[,data.subtypes.arsyn$Targets$Subtype == "LumA"])
  luma <- cbind(gene=as.character(data.subtypes.arsyn$Annot$EnsemblID), luma)
  
  lumb <- as.data.frame(data.subtypes.arsyn$M[,data.subtypes.arsyn$Targets$Subtype == "LumB"])
  lumb <- cbind(gene=as.character(data.subtypes.arsyn$Annot$EnsemblID), lumb)
  
  her2 <- as.data.frame(data.subtypes.arsyn$M[,data.subtypes.arsyn$Targets$Subtype == "Her2"])
  her2 <- cbind(gene=as.character(data.subtypes.arsyn$Annot$EnsemblID), her2)
  
  healthy <- as.data.frame(data.subtypes.arsyn$M[,data.subtypes.arsyn$Targets$Subtype == "Healthy"])
  healthy <- cbind(gene=as.character(data.subtypes.arsyn$Annot$EnsemblID), healthy)
  
  symbols <-as.character(norm.data.cpm10$Annot$EnsemblID)
  
  write.table(basal, file = "Basal_cpm10_arsyn.tsv", sep="\t", quote=FALSE, row.names=FALSE)
  write.table(luma, file = "LumA_cpm10_arsyn.tsv", sep="\t", quote=FALSE, row.names=FALSE)
  write.table(lumb, file = "LumB_cpm10_arsyn.tsv", sep="\t", quote=FALSE, row.names=FALSE)
  write.table(her2, file = "Her2_cpm10_arsyn.tsv", sep="\t", quote=FALSE, row.names=FALSE)
  write.table(healthy, file = "Healthy_cpm10_arsyn.tsv", sep="\t", quote=FALSE, row.names=FALSE)
  
  save(data.subtypes.arsyn, file="data_subtypes_arsyn.RData", 
       compress="xz")
