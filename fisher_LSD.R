library("reshape2")
library("ggplot2")
library(ggthemes)

setwd("/mnt/ddisk/transpipeline-data/breast-data/subtypes")
subtypes.pal <- c("#C7C7C7", "#FF9E4A", "#67BF5C", "#729ECE", "#ED665D")

####### Esta
datos <- data.frame(
    Subtype = c("Basal", "Her2", "LumB", "LumA", "Healthy"),
    LogMI = c(-2.79, -2.88, -3.00, -3.22, -3.36),
    SE = rep(3.1E-03, 5),
    Group = c("E", "D", "C", "B", "A"),
    D = rep("1e+6", 5)
  )

datos$Subtype <- factor(
  as.character(datos$Subtype),
  levels = c("Healthy", "LumA", "LumB", "Her2", "Basal")
)
datos

p <- ggplot(
  data = datos,
  aes(
    x = Subtype,
    y = exp(LogMI),
    ymin = exp(LogMI-SE),
    ymax = exp(LogMI),
    label = Group,
    fill = Subtype
  )
) +
  geom_bar(stat = "identity") +
  geom_errorbar()+
  geom_text(aes(y=0.075), size = 10) +
  ylim(c(0, 0.080)) +
  ylab("Mutual Information") +
  xlab("") +
  theme_few(base_size = 30) +
  theme(legend.position = "none", axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        plot.title = element_text(size = 25)) +
  labs(title = "Fisher's Least Significant Difference") +
  scale_fill_manual(values = subtypes.pal)
p

base <- 1e6
log(base)
log(base)^2
log(base)^3
log(base)^4






