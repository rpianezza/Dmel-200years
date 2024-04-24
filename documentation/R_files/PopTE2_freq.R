library(tidyverse)
theme_set(theme_bw())

df<-read.csv("/Volumes/Element/Backup/Two_centuries_dmel/PopTE2/map_trimmed/filter/Pool_trimmed_filter44.teinsertions", header=FALSE, sep = "\t")
colnames(df)<-c("ID", "chr", "postion", "strand", "TE_family", "order", "support", "comment", "TE_frequency")


df %>% 
  filter(TE_family == "Transib1") %>% 
  arrange(desc(TE_frequency))

TE_order <- c("Blood_diverged","Opus_diverged" ,"412_diverged","Tirant_diverged","I-ele_diverged","Hobo_diverged",
              "Spoink_diverged", "MLE_diverged","Souslik_diverged","Transib1_diverged",
              "Blood","Opus" ,"412","Tirant","I-ele","Hobo" ,"P-ele","Spoink", "MLE","Souslik","Transib1",)

df$TE_family <- factor(df$TE_family, levels = TE_order)

df <- df %>% mutate(TE_type = ifelse(grepl("diverged", TE_family, ignore.case = TRUE), "diverged", "canonical"))


ggplot(df, aes(x = TE_family, y = TE_frequency, fill = TE_type)) +
  geom_point() +
  geom_boxplot() +
  labs(x = "TE Family", y = "TE Frequency") +  
  ggtitle("TE frequency by TE Family") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  scale_fill_manual(values = c("diverged" = "#F8766D", "canonical" = "#619CFF"))+
  labs(fill = "Insertion type")



