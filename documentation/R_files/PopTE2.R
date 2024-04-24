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
tp$fragtype="consensus"
tp[tp$div>1.5,]$fragtype="diverged"
tp <- subset(tp, strain == "TOM008")

tp$fragtype <- factor(tp$fragtype, levels=c("consensus","diverged"))
tp<-tp[order(tp$te, tp$fragtype, tp$chr, tp$qstart, decreasing = FALSE),]


if (tp$strain[1] == "CanS") {
  tp <- rbind(tp, mock_df)
}

tp$xp <- 0
tp$qstart <- as.integer(tp$qstart)
tp$telen <- as.integer(tp$telen)

c <- 0

#v2
for (i in 1:nrow(tp)) {
  if (i == 1){
    c <- 1
  } else {
    if (tp$te[i] != tp$te[i - 1]){
      c <- 1
    }
    else{
      if(tp$fragtype[i] != tp$fragtype[i - 1]) {
        c <- 1
      }
      else {
        if (tp$chr[i] == tp$chr[i - 1]){
          if ((tp$qstart[i] - tp$qstart[i-1]) < tp$telen[i]){
            c <- c - 1
          }
        }
      }
    }
  }
  tp$xp[i] <- c
  c <- c +1
}


tp$xp=-1*tp$xp

#Sorting from longest maintaining the pairs
#1 retain th information of the pairs
tp <- mutate(tp, infopairs = paste(te, fragtype, xp, sep = ""))
#2 sort
tp<-tp[order(tp$te, tp$fragtype, tp$fraglen,decreasing=TRUE),]
#3 move the sorter of the pair in new df called tp_pair
v_infoTE <- c()
#4 make 2 df in s keep the second of the pair
tp_pair_l <- tp[0, ]
tp_pair_s <- tp[0, ]

for (i in 1:nrow(tp)) {
  if (tp$infopairs[i] %in% v_infoTE) {
    tp_pair_s <- rbind(tp_pair_s, tp[i, ])
  } else {
    tp_pair_l <- rbind(tp_pair_l, tp[i, ])
    v_infoTE <- c(v_infoTE, tp$infopairs[i])
  }
}

#5 id for visualization
tp_pair_l<-tp_pair_l %>% group_by(te,fragtype) %>% mutate(id = row_number())
tp_pair_l$id=-1*tp_pair_l$id
#6 Give the id to tp_infopair_s PROBLEM IN THE MERGE
#info_id <- merge(tp_pair_s, tp_pair_l, by = "infopairs", all.x = TRUE)[, c("infopairs", "id")]

for (i in 1:nrow(tp_pair_s)) {
  for (k in 1:nrow(tp_pair_l)) {
    if (tp_pair_s$infopairs[i] == tp_pair_l$infopairs[k]) {
      tp_pair_s$id[i] <- tp_pair_l$id[k]}
  }
}

#7 Combine the two df back
new_tp <- rbind(tp_pair_l, tp_pair_s)


new_tp %>% filter(fragtype == "consensus") %>% summarise(sum(as.numeric(fraglen)))
new_tp %>% group_by(fragtype) %>% summarise(sum(as.numeric(fraglen)))

#Too many small highly diverged fragments in 412 removed for readability
new_tp <- subset(new_tp, id > -51)
# Order TEs
new_tp$te <- factor(new_tp$te, levels = c("Blood", "412", "Opus", "Tirant", "I-ele", "Hobo", "P-ele", "Spoink", "MLE", "Souslik", "Transib1"))

p<-ggplot(new_tp,aes(x = rstart, y = id, xend = rend, yend = id,color=div))+
  geom_segment(size=1)+facet_grid(fragtype~te,scales="free",space="free")+
  xlab("position")+ylab("TE fragment")+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 16),
        legend.position = "nonw")+
  scale_colour_gradientn(colours = c("#0C6B37", "#FBB324", "#BC2023"))

plot(p)

# bed file for diverged
diverged_tp <- subset(new_tp, fragtype == "diverged")
diverged_tp <- diverged_tp %>%
  mutate(name = paste0(te, "_", row_number()))

columns_bed <- c("chr", "qstart", "qend", "strand", "name")
diverged_tp <- diverged_tp[, columns_bed, drop = FALSE]

# names
colnames(diverged_tp) <- c("chrom", "start", "end", "strand", "name")

# order
diverged_tp <- diverged_tp %>%
  select("chrom", "start", "end", "name", "strand")


# Write dataframe to a BED file without headers
write.table(diverged_tp, "/Users/ascarpa/Downloads/dmel_twocenturies_local/Supp_fig/PopTE2/diverged_tp.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
