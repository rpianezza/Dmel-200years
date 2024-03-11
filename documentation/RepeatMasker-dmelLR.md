Drosophila melanogaster - TE invasions in 49 LR assemblies
================

``` r
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(svglite))
theme_set(theme_bw())
```

``` r
dmel_lr_meta <- read_tsv("/Volumes/Storage/dmel-full-story/longread-metadata.txt") %>% select(strain, year, continent, lat, lon)
```

    ## Rows: 49 Columns: 7
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (4): GCA, strain, continent, location
    ## dbl (3): year, lat, lon
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
RM <- read_tsv("/Volumes/Storage/dmel-full-story/RepeatMasker/RM-longreads/merged.clean.sum", col_names = c("rm","SW","pid","contig","qstart","qend","strand","te","rstart","rend","score","strain")) %>% mutate(len = ifelse(qstart>qend, qstart-qend, qend-qstart), strain = gsub("D.mel.","",strain)) %>% inner_join(dmel_lr_meta, by="strain") %>% filter(te!="Shellder") %>% distinct()
```

    ## Rows: 296691 Columns: 12
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (5): rm, contig, strand, te, strain
    ## dbl (7): SW, pid, qstart, qend, rstart, rend, score
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
write_tsv(RM, "/Volumes/Storage/GitHub/Dmel-200years/data/RM-readable.tsv")
```

``` r
RM_full_len <- RM %>% filter(!(te %in% c("Hobo","PPI251","Transib_Riccardo"))) %>% select(strain, te, SW, pid, len, score, year) %>% filter(score > 0.8, pid < 1)
RM_DNA_TE <- RM %>% filter(te %in% c("Hobo","Transib_Riccardo","PPI251")) %>% select(strain, te, SW, pid, len, score, year) %>% filter(score > 0.5, pid < 1)
RM_P <- RM %>% filter(te == "PPI251") %>% select(strain, te, SW, pid, len, score, year) %>% filter(score > 0.25, pid < 1)
RM_count <- bind_rows(RM_full_len, RM_DNA_TE, RM_P)

strains <- RM %>% select(strain, year) %>% distinct()
te <- RM %>% select(te) %>% distinct()
combinations <- expand.grid(strain = strains$strain, te = te$te)
repeated_combinations <- combinations[rep(1:nrow(combinations), each = 11), ]
dummy <- as_tibble(repeated_combinations) %>% distinct() %>% inner_join(strains, by="strain")

RM_plottable <- RM_count %>% group_by(strain, te, year) %>% summarise(insertions = n()) %>% right_join(dummy, by=c("strain","te","year")) %>% mutate(insertions = ifelse(is.na(insertions), -1, insertions)) %>% mutate(presence = ifelse(insertions>1, "present", "absent"), year = ifelse(is.na(year), 1975, year), te=ifelse(te=="Transib_Riccardo","Transib1",te))
```

    ## `summarise()` has grouped output by 'strain', 'te'. You can override using the
    ## `.groups` argument.

``` r
RM_plottable$te <- factor(RM_plottable$te, levels = c("Blood", "Opus", "412", "Tirant", "DMIFACA", "Hobo", "PPI251", "Spoink", "Micropia", "Souslik", "Transib1"))

insertions_plot <- ggplot(RM_plottable, aes(x=reorder(strain,year), y=insertions, fill=presence))+
  geom_bar(stat = "identity")+
    labs(x="", y="copynumber", fill="TE")+
    facet_wrap(~te, ncol=1)+
scale_fill_manual(values = c("violet", "darkblue"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4), legend.position = "top", legend.key.size = unit(0.2, "cm"))

ggsave("/Volumes/Storage/dmel-full-story/figures/LR-08.png", insertions_plot, height = 20)
```

    ## Saving 7 x 20 in image

``` r
to_map <- RM_plottable %>% inner_join(dmel_lr_meta, by=c("strain","year")) %>% ungroup()
to_map$presence <- factor(to_map$presence, levels = c("absent", "present"))

world_map <- map_data("world")
world_map <- subset(world_map, region != "Antarctica")

tomap_sous <- to_map %>% filter(te == "Souslik") %>% mutate(year = ifelse(year > 2003, ">2003", "<2003"))
tomap_tra <- to_map %>% filter(te == "Transib1") %>% mutate(year = ifelse(year > 2011, ">2011", "<2011"))
tomap_mic <- to_map %>% filter(te == "Micropia") %>% mutate(year = ifelse(year > 1995, ">1995", "<1995"))

(lr_map_sous <- ggplot() +
    geom_map(data = world_map, map = world_map,
             aes(long, lat, map_id = region),
             color = "white", fill = "cornsilk3", linewidth = 0) +
   geom_point(data = tomap_sous, aes(x = lon, y = lat, color = presence), size = 4, position = position_jitter(width = 1, height = 1), alpha = 0.5)) +
   scale_colour_manual(values = c("violet", "darkblue")) +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), legend.position = "bottom") +
  facet_wrap(~year)+
  labs(color = "Souslik")
```

    ## Warning in geom_map(data = world_map, map = world_map, aes(long, lat, map_id =
    ## region), : Ignoring unknown aesthetics: x and y

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

![](RepeatMasker-dmelLR_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
(lr_map_tra <- ggplot() +
    geom_map(data = world_map, map = world_map,
             aes(long, lat, map_id = region),
             color = "white", fill = "cornsilk3", linewidth = 0) +
   geom_point(data = tomap_tra, aes(x = lon, y = lat, color = presence), size = 4, position = position_jitter(width = 1, height = 1), alpha = 0.5)) +
   scale_colour_manual(values = c("violet", "darkblue")) +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), legend.position = "bottom") +
  facet_wrap(~year)+
  labs(color = "Transib")
```

    ## Warning in geom_map(data = world_map, map = world_map, aes(long, lat, map_id = region), : Ignoring unknown aesthetics: x and y
    ## Removed 1 rows containing missing values (`geom_point()`).

![](RepeatMasker-dmelLR_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
(lr_map_mic <- ggplot() +
    geom_map(data = world_map, map = world_map,
             aes(long, lat, map_id = region),
             color = "white", fill = "cornsilk3", linewidth = 0) +
   geom_point(data = tomap_mic, aes(x = lon, y = lat, color = presence), size = 4, position = position_jitter(width = 1, height = 1), alpha = 0.5)) +
   scale_colour_manual(values = c("violet", "darkblue")) +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), legend.position = "bottom") +
  facet_wrap(~year)+
  labs(color = "Micropia")
```

    ## Warning in geom_map(data = world_map, map = world_map, aes(long, lat, map_id = region), : Ignoring unknown aesthetics: x and y
    ## Removed 1 rows containing missing values (`geom_point()`).

![](RepeatMasker-dmelLR_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->
