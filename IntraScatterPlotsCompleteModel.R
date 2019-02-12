suppressMessages(library("car"))
suppressMessages(library("ggplot2"))
suppressMessages(library("reshape"))
suppressMessages(library("cowplot"))
suppressMessages(library("grid"))
suppressMessages(library("scales"))
suppressMessages(library("plyr"))

setwd("/mnt/ddisk/transpipeline-data/breast-data/subtypes")

datos <- read.delim(file = "intra/intra-fixed-bin-size/100/all.txt", header = T)
names(datos) <- c("Bin", "Distance", "MI", "Subtype", "Chr")
datos$logD <- log(datos$Distance)
datos$logD2 <- log(datos$Distance)^2
datos$logD3 <- log(datos$Distance)^3
datos$logD4 <- log(datos$Distance)^4
datos$logD5 <- log(datos$Distance)^5
datos$logMI <- log(datos$MI)
datos$Subtype <- relevel(datos$Subtype, ref = "healthy")
levels(datos$Subtype) <- c("Healthy", "Basal", "Her2", "LumA", "LumB")
datos$Chr <- relevel(datos$Chr, ref = "21")
head(datos)

model <- lm(
  data = datos,
  logMI ~ Subtype*Chr*logD+Subtype*Chr*logD2+Subtype*Chr*logD3+Subtype*Chr*logD4
)
#par(mfrow = c(4,1))
plot(model)
#summary(model)
datos$fitted <- exp(fitted(model))
datos$Chr <- factor(datos$Chr, levels = c(as.character(1:22), "X"))
datos$Subtype <-  factor(datos$Subtype, levels = c("Healthy", "LumA", "LumB", "Her2", "Basal"))
  
## Supplementario todas las figuras
p <- ggplot(
  data = datos,
  aes(
    x = Distance/1e6,
    y = MI
   # group = Chr
  )
) +
  geom_point() +
  geom_line(
    aes(
      x = Distance/1e6,
      y = fitted,
      color = "pink"
    ), size = 1.2
  ) +
  facet_grid(Chr ~ Subtype, scales = "free_y") +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 30),
    strip.text.y  = element_text(size = 24),
    axis.title=element_text(size=30,face="bold")
  ) + xlab("Distance (Million bp)") + ylab("Mutual Information")

p


### All y Cromosoma 1
datos_chr1 <- datos[datos$Chr == "1", ]
datos_chr1$Type <- "Chr 1 interactions"

datos_all <- datos
datos_all$Type <- "All cis- interactions"
datos_all <- rbind(datos_chr1, datos_all)
datos_chr1_d <- datos_chr1
datos_chr1_d$Type <- "All cis- interactions"
datos_chr1_d <-rbind(datos_chr1, datos_chr1_d)

#change this to get the correct lines
datos_chr1_d <- datos_chr1_d[order(datos_chr1_d$Type, datos_chr1_d$Subtype), ]

p <- ggplot(
  data = datos_all,
  aes(
    x = Distance/1e6,
    y = MI
  )
) +
  geom_point(
    color = "gray38"
  ) +
  facet_grid(
    Type ~ Subtype
    #Subtype ~ Type
  ) +   geom_line(
    data = datos_chr1_d,
    aes(
      x = datos_chr1_d$Distance/1e6,
      y = datos_chr1_d$fitted,
      color = "pink"
    ), size = 1.2
  ) +
  theme_bw(base_size = 18) +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 24),
    strip.text.y = element_text(size = 14),
    axis.title=element_text(size=24)
  ) + xlab("Distance (Million bp)") +
  ylab("Mutual Information") +
  scale_x_continuous(breaks = seq(0, 300, by = 100))
p
#### solo Chromosoma 1
datos_chr1 <- datos_chr1[order(datos_chr1$Type, datos_chr1$Subtype), ]
p <- ggplot(
  data = datos_chr1,
  aes(
    x = Distance/1e6,
    y = MI
  )
) +
  geom_point(
    color = "gray38"
  ) +
  facet_grid(
    Type ~ Subtype
    #Subtype ~ Type
  ) +   geom_line(
    data = datos_chr1,
    aes(
      x = datos_chr1$Distance/1e6,
      y = datos_chr1$fitted,
      color = "pink"
    ), size = 1.2
  ) +
  theme_bw(base_size = 34) +
  theme(
    legend.position = "none",
    strip.text.x = element_blank(),
    strip.text.y = element_text(size = 32),
    axis.title=element_text(size=34)
  ) + xlab("Distance (Million bp)") +
  ylab("Mutual Information") +
  scale_x_continuous(breaks = seq(0, 300, by = 100))
p

### Figura 4

p2 <- ggplot(
  data = subset(
    datos,
    Subtype %in% c("Healthy", "Basal") &
      Chr == "1"  
  ),
  aes(
    x = Distance/1e6,
    y = MI,
    shape = Subtype, 
    color=Subtype
  ) 
) + geom_point() + 
  scale_colour_grey(start = 0.4, end = 0.7) + 
  geom_line(
    data = subset(
      datos,
      Subtype %in% c("Basal") &
        Chr == "1"  
    ),
    aes(
      x = Distance/1e6,
      y = fitted
    ),
    size = 1.2,
    color = hue_pal()(1)
  ) + 
  geom_line(
    data = subset(
      datos,
      Subtype %in% c("Healthy") &
        Chr == "1"  
    ),
    aes(
      x = Distance/1e6,
      y = fitted
    ),
    size = 1.2,
    color = hue_pal()(2)[2]
  ) +
  theme_classic(base_size = 22) +
  theme(
    legend.position = "bottom"
  ) +
  xlab("Distance (Million bp)") +
  ylab("Mutual Information") + 
  guides(shape = guide_legend(override.aes = list(size = 5)))
  
pzoom <- p2
pzoom <- pzoom +
  xlim(0, 20) +
  ylim(0.020, 0.08) +
  theme_classic(base_size = 16) +
  theme(legend.position = "none", axis.title.x=element_blank(),
        axis.title.y=element_blank(),) 
vp <- viewport(
  width = 0.48, 
  height = 0.38, 
  x = 0.48,
  y = 0.58, 
  just = c("left", "bottom")
)
p2; print(pzoom, vp = vp)

