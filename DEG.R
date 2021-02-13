#################################################################################
## Proyecto: Mecanismos de desregulación del metabolismo energético asociados a 
##           cáncer de mama
################################################################################
## 3) DIFFERENTIAL EXPRESSION
## Autor: Cristóbal Fresno
## Fecha: 2019/02/05
## Modificado por : Erandi Serrano
## https://github.com/CSB-IG/tcgarnaseqbc/blob/master/DifGenes.R
#################################################################################
#options(width=120)

# Instalar las bibliotecas necesarias
# NOIseq
#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#    BiocManager::install("NOISeq", version = "3.8")
# EDAseq
#if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#      BiocManager::install("EDASeq", version = "3.8")
#      BiocManager::install("Glimma", version = "3.8")
library("NOISeq")
library("EDASeq")
library("Glimma")

#Insatalar limma
#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#   BiocManager::install("limma", version = "3.8")
library("limma")
library("edgeR")
library("ggplot2")

setwd("/Volumes/nawal/PHD-projects.bkp/CANCER-SUBTIPOS/subtipos.prj")
DATADIR <- 'pipeline/data/'
DEGDATA <- paste(DATADIR, "deg", sep = "")
dir.create(DEGDATA)
RDATA <- paste(DATADIR, "rdata", sep = "")
f <- paste(RDATA, "Length.full_GC.full_Between.tmm_Norm_cpm10_arsyn_all.tsv", sep="/")

# tname <- 'tBa' # idetificacion de los tumores
# tname <- 'tLA' # idetificacion de los tumores
# tname <- 'tLB' # idetificacion de los tumores
tname <- 'tH2' # idetificacion de los tumores


Matrix <-  read.delim(f, stringsAsFactors = FALSE)
Normal <- Matrix[,grep('^Ctr',names(Matrix))]
Tumor <- Matrix[,grep(paste("^", tname, sep = ""),names(Matrix))]
M <- cbind(Normal,Tumor)
rownames(M) <-  Matrix[,1]
#M1 <- log2(M) # NO ES LO MISMO QUE CPM

#### Correcciones y acomodado de los Data Frames
Group <- factor(substr(colnames(M), start = 1, stop = 3))
targets <- data.frame(ID = colnames(M), Group = Group)
Tumor <- list(M = cpm(M, log = TRUE), Targets = targets)

######### DISENIO
# ¿Para el diseño correcto ya no hay interceptor "Group-1"?
# El resultado es que muchos genes son diferencialmente expresados
design <- model.matrix(~1+Group, data = Tumor$Targets)
rownames(design) <- Tumor$Targets$ID
M <- Tumor$M
row.names(M) <- Matrix[,1]

gname <- paste("Group",tname,sep = "") # NOmbre del grupo en el experimento, por ejemplo: GrouptLA

######### AJUSTE
fit <- lmFit(M, design)
head(fit$coefficients)
fit2 <- eBayes(fit)
fit2$fdr<-apply(fit2$"p.value", 2, p.adjust, method="fdr")
ggplot(as.data.frame(fit$coefficients), aes(x=gname))+geom_density()
nrow(M)

#####  Buscando los genes diferenciales cortados por los alphas
alphas <- c(0.05, 10^(-2:-10), 10^(seq(-20, -100, by=-10)))
degCount<-sapply(alphas, function(alpha){
  table(fit2$fdr[,gname]<alpha)
})
colnames(degCount)<-alphas
# plot(degCount["TRUE",])
degCount<-as.data.frame(t(degCount))
degCount$alpha<-alphas
degCount$Noise<-alphas*nrow(M)
degCount




##### [1] Select genes diferenciales por el metodo de EBAYES
##### Ajustar para cada subtipo y caso del grupo
# Pval <- 1e-50 # BASAL 10% de los genes totales(~13000) son TRUE
# Pval <- 1e-30 # LumA 
# Pval <- 1e-40 # LumB
Pval <- 1e-30 # Her2

result <- decideTests(fit2, adjust.method="fdr", p.value=Pval)
pdf(file=paste(DEGDATA, "01-venndiagram.pdf", sep="/"))
vennDiagram(result)
dev.off()
genesFULL<-cbind(
  Coef=fit2$fdr[fit2$fdr[,gname]<Pval,gname],
  Pvalue=fit2$p.value[fit2$fdr[,gname]<Pval,gname]
)
fname <- paste("Ebayes-",tname,".tsv",sep = "")
write.table(genesFULL, file = paste(DEGDATA, fname, sep="/"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE )
#write.table(fit2, file = paste(DEGDATA, "DifExpMatrix.tsv", sep="/"), sep = "\t", row.names = F)




##### [2] Select genes diferenciales por el metodo de TREAT
lfc<-seq(0, 3, by=0.5)
B <- 5
fitTreat<-lapply(lfc, function(x){
  ##Ajustando el modelo
  aux<-treat(fit, lfc=x)
  aux$fdr<-apply(aux$p.value, 2, p.adjust, method="fdr")
  ##Estadistico B
  p<-1-aux$fdr[, gname] ##Probabilidad de salir diferencialºº
  aux$B<-log(p/(1-p))
  ##Diferenciales
  aux$dif<-data.frame(
    alpha=alphas,
    lfc=rep(x, length(alphas)),
    Dif=t(sapply(alphas, function(alpha){
      table(factor(aux$fdr[,gname]<alpha & (aux$B>5), levels=c(TRUE, FALSE)))
    })),
    Noise=alphas*nrow(M)
  )
  return(aux)
})
names(fitTreat)<-paste("lfc", lfc, sep="")
difTreat<-do.call(rbind, lapply(fitTreat, function(x){x$dif}))
difTreat

# Pval <- 1e-50             # Basal pvalue = 1e-50
# LFCsel <- fitTreat$lfc0   # Basal lfc = 0.0
# Pval <- 1e-02             # LumA pvalue = 1e-02
# LFCsel <- fitTreat$lfc0.5 # LumA lfc = 0.5
# Pval <- 1e-40             # LumB pvalue = 1e-40
# LFCsel <- fitTreat$lfc0   # LumB lfc = 0.0
Pval <- 1e-08             # LumB pvalue = 1e-40
LFCsel <- fitTreat$lfc0.5   # LumB lfc = 0.0

LFCsel$deg<-LFCsel$fdr[,gname, drop=FALSE]<Pval & (LFCsel$B>5)
table(LFCsel$deg)

#Glimma plot
dt <- decideTests(LFCsel) #### NO ENTENDI
#La ocupamos despues para pintar las redes y por eso la guardamos...
write.table(dt, file = "dt.txt", quote = FALSE, sep = " ", row.names = TRUE, col.names = FALSE )
summary(dt)
# glMDPlot(path = DEGDATA, fitTreat$lfc0, counts = M, groups = Group, status = dt)

#Salida
p<-1-LFCsel$fdr[, gname]
genesFULL<-cbind(
  Coef=LFCsel$coefficients,
  Diff=LFCsel$deg[, gname],
  p.value=LFCsel$p.value[, gname],
  FDR=LFCsel$fdr[, gname],
  B=log(p/(1-p)),
  Exp=LFCsel$coefficients%*%t(unique(design))
)

head(genesFULL)
fname <- paste("Treat-",tname,".tsv",sep = "")
write.table(genesFULL, file =  paste(DEGDATA, fname, sep="/"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE )



