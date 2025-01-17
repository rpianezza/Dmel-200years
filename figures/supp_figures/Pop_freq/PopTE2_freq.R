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

TE_order <- c("Blood","Opus" ,"412","Tirant","I-ele","Hobo" ,"P-ele","Spoink", "MLE","Souslik","Transib1",
              "Blood_diverged","Opus_diverged" ,"412_diverged","Tirant_diverged","I-ele_diverged","Hobo_diverged",
              "Spoink_diverged", "MLE_diverged","Souslik_diverged","Transib1_diverged")

df$TE_family <- factor(df$TE_family, levels = TE_order)
df2$TE_family <- factor(df2$TE_family, levels = TE_order)
df3$TE_family <- factor(df3$TE_family, levels = TE_order)

df <- df %>% mutate(TE_type = ifelse(grepl("diverged", TE_family, ignore.case = TRUE), "diverged", "canonical"))
df2 <- df2 %>% mutate(TE_type = ifelse(grepl("diverged", TE_family, ignore.case = TRUE), "diverged", "canonical"))
df3 <- df3 %>% mutate(TE_type = ifelse(grepl("diverged", TE_family, ignore.case = TRUE), "diverged", "canonical"))


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

g3 <- ggplot(df2, aes(x = TE_family, y = TE_frequency, fill = TE_type)) +
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

#java -jar popte2-v1.10.03.jar se2pe --fastq1 trimmed/SRR3939083_1_clean.fastq.gz --fastq2 trimmed/SRR3939083_2_clean.fastq.gz --bam1 map/USA_2014_1.sam --bam2 map/USA_2014_2.sam --sort --output output/USA.sort.bam 
#java -jar popte2-v1.10.03.jar ppileup --bam output/USA.sort.bam --map-qual 15 --hier ../../Supp_fig/PopTE2/names_TEs.hier --output output/USA.ppileup.gz
#java -jar popte2-v1.10.03.jar identifySignatures --ppileup output/USA.ppileup.gz --mode separate --output output/USA3.signatures --min-count 3
#java -jar popte2-v1.10.03.jar frequency --ppileup output/UK.ppileup.gz --signature output/UK3.signatures --output output/UK3.freqsig 
#java -jar popte2-v1.10.03.jar filterSignatures --input output/USA3.freqsig --output output/USA3.filter.freqsig --max-otherte-count 2 --max-structvar-count 2
#java -jar popte2-v1.10.03.jar pairupSignatures --signature output/USA3.filter.freqsig --ref-genome ../../Supp_fig/PopTE2/Dmel_TE_merged.fa --hier ../../Supp_fig/PopTE2/names_TEs.hier --min-distance -200 --max-distance 300 --output output/USA3.filter.teinsertions

#Remove the commas the german java puts to separate the thousands: sed 's/,//g' output_final/USA-min-count3.freqsig > output_final/USA-min-count3_fixed.freqsig

#java -jar popte2-v1.10.03.jar se2pe --fastq1 trimmed/SRR3939083_1_clean.fastq.gz --fastq2 trimmed/SRR3939083_2_clean.fastq.gz --bam1 map/USA_2014_1.sam --bam2 map/USA_2014_2.sam --sort --output output/USA.sort.bam 
#java -jar popte2-v1.10.03.jar ppileup --bam output/USA.sort.bam --map-qual 15 --hier ../../Supp_fig/PopTE2/names_TEs.hier --output output/USA.ppileup.gz
#java -jar popte2-v1.10.03.jar identifySignatures --ppileup output_final/USA.ppileup.gz --mode separate --output output_final/USA-min-count3.signatures --min-count 3
#java -jar popte2-v1.10.03.jar frequency --ppileup output_final/USA.ppileup.gz --signature output_final/USA-min-count3.signatures --output output_final/USA-min-count3.freqsig
#sed 's/,//g' output_final/USA-min-count3.freqsig > output_final/USA-min-count3_fixed.freqsig
#java -jar popte2-v1.10.03.jar filterSignatures --input output_final/USA-min-count3_fixed.freqsig --output output_final/USA-min-count3.filter.freqsig  --max-otherte-count 2 --max-structvar-count 2
#sed 's/,//g' output_final/USA-min-count3.filter.freqsig > output_final/USA-min-count3_fixed.filter.freqsig
#java -jar popte2-v1.10.03.jar pairupSignatures --signature output_final/USA-min-count3_fixed.filter.freqsig --ref-genome ../../Supp_fig/PopTE2/Dmel_TE_merged.fa --hier names_TEs.hier --min-distance -200 --max-distance 300 --output output_final/USA-min-count-3.teinsertions
#sed 's/\([0-9]\),\([0-9]\)/\1.\2/g' USA-min-count-3.teinsertions > USA-min-count-3.teinsertions_fixed
