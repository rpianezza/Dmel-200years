library(tidyverse)
library(ggpubr)
theme_set(theme_bw())

df<-read.csv("/Volumes/Element/Backup/dmel_twocenturies_local/Revision/final_USA_UK/output_sh_3_Ukraine/Ukraine-min-count3.teinsertions_fixed", header=FALSE, sep = "\t")
df2<-read.csv("/Volumes/Element/Backup/dmel_twocenturies_local/Revision/final_USA_UK/output_sh_3_UK/UK-min-count3.teinsertions_fixed", header=FALSE, sep = "\t")
df3<-read.csv("/Volumes/Element/Backup/dmel_twocenturies_local/Revision/final_USA_UK/output_sh_3/USA-min-count3.teinsertions_fixed", header=FALSE, sep = "\t")

colnames(df)<-c("ID", "chr", "postion", "strand", "TE_family", "order", "support", "comment", "TE_frequency")
colnames(df2)<-c("ID", "chr", "postion", "strand", "TE_family", "order", "support", "comment", "TE_frequency")
colnames(df3)<-c("ID", "chr", "postion", "strand", "TE_family", "order", "support", "comment", "TE_frequency")

df <- df %>%
  mutate(
    TE_family = str_replace(TE_family, "MLE", "McLE"),
    TE_family = str_replace(TE_family, "MLE_diverged", "McLE_diverged")
  )

df2 <- df2 %>%
  mutate(
    TE_family = str_replace(TE_family, "MLE", "McLE"),
    TE_family = str_replace(TE_family, "MLE_diverged", "McLE_diverged")
  )

df3 <- df3 %>%
  mutate(
    TE_family = str_replace(TE_family, "MLE", "McLE"),
    TE_family = str_replace(TE_family, "MLE_diverged", "McLE_diverged")
  )

TE_order <- c("Blood","Opus" ,"412","Tirant","I-ele","Hobo" ,"P-ele","Spoink", "McLE","Souslik","Transib1",
              "Blood_diverged","Opus_diverged" ,"412_diverged","Tirant_diverged","I-ele_diverged","Hobo_diverged",
              "Spoink_diverged", "McLE_diverged","Souslik_diverged","Transib1_diverged")

df$TE_family <- factor(df$TE_family, levels = TE_order)
df2$TE_family <- factor(df2$TE_family, levels = TE_order)
df3$TE_family <- factor(df3$TE_family, levels = TE_order)

df <- df %>% mutate(TE_type = ifelse(grepl("diverged", TE_family, ignore.case = TRUE), "diverged", "canonical"))
df2 <- df2 %>% mutate(TE_type = ifelse(grepl("diverged", TE_family, ignore.case = TRUE), "diverged", "canonical"))
df3 <- df3 %>% mutate(TE_type = ifelse(grepl("diverged", TE_family, ignore.case = TRUE), "diverged", "canonical"))


df <- df %>%
  group_by(TE_family) %>%
  filter(n() >= 3) %>%
  ungroup()

df2 <- df2 %>%
  group_by(TE_family) %>%
  filter(n() >= 3) %>%
  ungroup()

df3 <- df3 %>%
  group_by(TE_family) %>%
  filter(n() >= 3) %>%
  ungroup()


g1 <- ggplot(df, aes(x = TE_family, y = TE_frequency, fill = TE_type)) +
  geom_point() +
  geom_boxplot() +
  labs(x = NULL, y = "TE Frequency") +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values = c("diverged" = "#F8766D", "canonical" = "#619CFF"))+
  labs(fill = "Insertion type")


g2 <- ggplot(df2, aes(x = TE_family, y = TE_frequency, fill = TE_type)) +
  geom_point() +
  geom_boxplot() +
  labs(x = NULL, y = "TE Frequency") +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values = c("diverged" = "#F8766D", "canonical" = "#619CFF"))+
  labs(fill = "Insertion type")

g3 <- ggplot(df3, aes(x = TE_family, y = TE_frequency, fill = TE_type)) +
  geom_point() +
  geom_boxplot() +
  labs(x = "TE Family", y = "TE Frequency") +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  scale_fill_manual(values = c("diverged" = "#F8766D", "canonical" = "#619CFF"))+
  labs(fill = "Insertion type")

g <- ggarrange(
  g1, g2, g3,
  ncol = 1, nrow = 3,
  labels = c("A", "B", "C"),
  align = "hv"
)

plot(g)
