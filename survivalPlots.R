setwd("/mnt/ddisk/transpipeline-data/breast-data/subtypes/")
library(survival)
library(survminer)

subtypes.pancancer <- read.delim("brca_tcga_pan_can_atlas_2018_clinical_data.tsv", stringsAsFactors = F)
files.pbcmc <- read.delim("id-filenames-submitterid-caseid-breast-cancer.2019-01-09.tsv", stringsAsFactors = F)
subtypes.pbcmc <- read.delim("id-subtype.tsv", stringsAsFactors = F)
subtypes.pbcmc <- merge(files.pbcmc, subtypes.pbcmc, by = "ID")
subtypes.pbcmc <- subtypes.pbcmc[, c("ID", "SubmitterID", "Subtype")]

subtypes.pancancer$Subtype <- stringr::str_remove(subtypes.pancancer$Subtype, "BRCA_") 
subtypes.all <- merge(subtypes.pancancer, subtypes.pbcmc, by.x = "Patient.ID", by.y = "SubmitterID")
subtypes.all$Diff <- subtypes.all$Subtype.x != subtypes.all$Subtype.y
dim(subtypes.all[subtypes.all$Diff == T & subtypes.all$Subtype.x != "", ])
subtypes.all <- subset(subtypes.all, select = -c( Patient.ID, Diff))
colnames(subtypes.all) <- c("SampleID", "Status", "SubtypeCBio", "LastCommunication", "SurvivalMonths", "ID",  "Subtype")
subtypes.all$SurvivalDays <- subtypes.all$SurvivalMonths * 30

subtypes.diff <- subtypes.all[subtypes.all$Diff == T, c("Subtype.x", "Subtype.y", "Diff")]

subtypes.all <- subtypes.all[subtypes.all$SurvivalDays != 0 & !is.na(subtypes.all$SurvivalDays), ]

survival.time <- ifelse(is.na(subtypes.all$LastCommunication) | subtypes.all$SurvivalDays > subtypes.all$LastCommunication,
                                         as.integer(subtypes.all$SurvivalDays), subtypes.all$LastCommunication)
surv_object <- Surv(time = survival.time, event = subtypes.all$Status == "DECEASED")
subtypes.all$Subtype <- factor(subtypes.all$Subtype, levels = c("LumA", "LumB", "Her2", "Basal"))
fit1 <- survfit(surv_object ~ subtypes.all$Subtype, data = subtypes.all)
ggsurvplot(fit1, pval = F, xlim = c(0, 1000), ylim = c(0.8, 1.0),
           palette = c( "#ED665D", "#729ECE",  "#FF9E4A", "#67BF5C"),
           break.x.by = 200, censor = T, legend.title = "", 
           legend.labs = c("Basal", "Her2", "LumA", "LumB"), 
           ggtheme = theme_survminer(base_size = 20),
           font.x = c(25),
           font.y = c(25),
           font.tickslab = c(20),
           font.legend = c(20)
           ) 

+ theme(
             #legend.spacing.y = unit(0.0001, "cm"),
             legend.title = element_text(size = 12),
             legend.text = element_text(size = 20),
             axis.title=element_text(size=16)
           )
