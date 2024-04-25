library(tidyverse)
library(svglite)
theme_set(theme_bw())

h<-read.csv("/Users/ascarpa/Downloads/dmel_twocenturies_local/RM_TEs/RM-readable.tsv", header=TRUE, sep = "\t")
# rm SW pid contig qstart qend strand te rstart rend score strain  len year continent lat lon
# names(h)<-c("te","chr","strand","start","end","qstart","qend","div","fraglen","telen")
names(h)<-c("rm","SW","div", "chr", "qstart","qend", "strand", "te","rstart","rend","score","strain","fraglen", "year", "continent", "lat", "lon")
h$te <- sub("Transib_Riccardo", "Transib1", h$te)
h$te <- sub("DMIFACA", "I-ele", h$te)
h$te <- sub("PPI251", "P-ele", h$te)
h$te <- sub("Micropia", "MLE", h$te)


h$telen <- ifelse(h$te == "Blood", 7410,
                  ifelse(h$te == "Tirant", 8526,
                         ifelse(h$te == "I-ele", 5371,
                                ifelse(h$te == "412", 7567,
                                       ifelse(h$te == "Hobo", 2959,
                                              ifelse(h$te == "Opus", 7521,
                                                     ifelse(h$te == "Transib1", 3030,
                                                            ifelse(h$te == "MLE", 5360,
                                                                   ifelse(h$te == "Spoink", 5216,
                                                                          ifelse(h$te == "Souslik", 5275,
                                                                                 ifelse(h$te == "P-ele", 2907, NA)))))))))))






h$lenfraction <- h$fraglen / h$telen
tp<-subset(h,fraglen>300 & div<21)
tp$fragtype="canonical"
tp[tp$div>1.5,]$fragtype="degraded"

# Subset for strain and TEs
tp <- subset(tp, strain == "TOM008")
tp <- subset(tp, te == "Spoink" | te == "MLE" | te == "Souslik" | te == "Transib1")

tp$fragtype <- factor(tp$fragtype, levels=c("canonical","degraded"))


unique(tp$chr)
tp <- tp %>% 
  mutate(chr = ifelse(chr == "CM034739.1", "X", chr)) %>% 
  mutate(chr = ifelse(chr == "CM034740.1", "2L", chr)) %>% 
  mutate(chr = ifelse(chr == "CM034741.1", "2R", chr)) %>% 
  mutate(chr = ifelse(chr == "CM034742.1", "3L", chr)) %>% 
  mutate(chr = ifelse(chr == "CM034743.1", "3R", chr)) %>% 
  mutate(chr = ifelse(chr == "CM034744.1", "4", chr))

custom_te_order <- c("Spoink", "MLE", "Souslik", "Transib1")
tp$te <- factor(tp$te, levels = custom_te_order, ordered = TRUE)

p<-ggplot(tp,aes(x = qstart, y = 0, xend = qstart, yend = lenfraction,color=div))+
  geom_segment(size=1)+facet_grid(te+fragtype~chr,scales="free_x",space="free_x")+
  ggtitle("TOM008 (2016)") +
  xlab("position")+ylab("length fraction of TE insertion [0-1]")+
  scale_colour_gradientn(colours = c("#0C6B37", "#FBB324", "#BC2023"))+
  scale_x_continuous(breaks=c(0,10000000,20000000,30000000),
                     labels=c("0","10m","20m","30m"))+
  scale_y_continuous(breaks=c(0,0.5,1.0))+
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = "right")

plot(p)

ggsave("/Users/ascarpa/Downloads/dmel_twocenturies_local/RM_TEs/images/chr.png", plot = p, width = 18, height = 10, dpi = 300)
