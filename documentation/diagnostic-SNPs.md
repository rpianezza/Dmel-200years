Drosophila melanogaster - Define diagnostic SNPs for all TEs
================

``` r
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(svglite))
theme_set(theme_bw())
```

## File import

Read metadata and copynumber estimates

``` r
(metadata <- read_tsv("/Volumes/Storage/dmel-full-story/metadata.tsv"))
```

    ## Rows: 585 Columns: 9
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (6): Sample, strain, publication, study, study_id, location
    ## dbl (3): year, lat, lon
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

    ## # A tibble: 585 × 9
    ##    Sample      strain   publication    study study_id  year location   lat   lon
    ##    <chr>       <chr>    <chr>          <chr> <chr>    <dbl> <chr>    <dbl> <dbl>
    ##  1 SRR23876563 H10      https://doi.o… PRJN… Shpak20…  1815 Lund, S…    55    13
    ##  2 SRR23876564 H4       https://doi.o… PRJN… Shpak20…  1815 Lund, S…    55    13
    ##  3 SRR23876566 H6       https://doi.o… PRJN… Shpak20…  1815 Smaland…    60    15
    ##  4 SRR23876582 H13      https://doi.o… PRJN… Shpak20…  1815 Lund, S…    55    13
    ##  5 SRR23876583 H12      https://doi.o… PRJN… Shpak20…  1815 Lund, S…    55    13
    ##  6 SRR23876584 H11      https://doi.o… PRJN… Shpak20…  1815 Lund, S…    55    13
    ##  7 SRR23876562 H9       https://doi.o… PRJN… Shpak20…  1850 Passau,…    49    13
    ##  8 SRR23876569 H25      https://doi.o… PRJN… Shpak20…  1850 Passau,…    49    13
    ##  9 SRR23876565 H5       https://doi.o… PRJN… Shpak20…  1875 Zealand…    55    12
    ## 10 ERR6474638  Oregon-R https://doi.o… PRJN… Burny20…  1925 <NA>        NA    NA
    ## # ℹ 575 more rows

Which TEs do we need to check for diagnostic SNPs? Only re-invading TEs,
thus:

- Blood
- Tirant
- I-element
- Souslik
- Transib1

Read SNPs info from DeviaTE.

``` r
read_deviate <- function(path, meta){
  transib <- read_table(path, col_names = c("TE","Sample","pos","ref","A","C","G","T","cov","phys_cov","hq_cov","snp","refsnp","int_del","int_del_freq","trunc_left","trunc_rigth","ins","del","annotation")) %>% select(TE, Sample, pos, A, C, G, T, cov, hq_cov, snp) %>% mutate(Sample = gsub(".fastq.sort.bam","",Sample)) %>% inner_join(meta, by="Sample") %>% select(-study, -study_id, -publication)
}

blood <- read_deviate("/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Blood.deviaTE.txt", metadata)
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   TE = col_character(),
    ##   Sample = col_character(),
    ##   ref = col_character(),
    ##   snp = col_logical(),
    ##   refsnp = col_logical(),
    ##   int_del = col_logical(),
    ##   int_del_freq = col_logical(),
    ##   trunc_left = col_logical(),
    ##   ins = col_logical(),
    ##   del = col_logical(),
    ##   annotation = col_logical()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

    ## Warning: 9827 parsing failures.
    ##  row        col           expected actual                                                                              file
    ## 2799 trunc_left 1/0/T/F/TRUE/FALSE  2.975 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Blood.deviaTE.txt'
    ## 3857 trunc_left 1/0/T/F/TRUE/FALSE  4.27  '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Blood.deviaTE.txt'
    ## 4145 trunc_left 1/0/T/F/TRUE/FALSE  1.708 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Blood.deviaTE.txt'
    ## 4163 trunc_left 1/0/T/F/TRUE/FALSE  1.956 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Blood.deviaTE.txt'
    ## 4290 trunc_left 1/0/T/F/TRUE/FALSE  4.215 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Blood.deviaTE.txt'
    ## .... .......... .................. ...... .................................................................................
    ## See problems(...) for more details.

``` r
tirant <- read_deviate("/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Tirant.deviaTE.txt", metadata)
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   TE = col_character(),
    ##   Sample = col_character(),
    ##   ref = col_character(),
    ##   snp = col_logical(),
    ##   refsnp = col_logical(),
    ##   int_del = col_logical(),
    ##   int_del_freq = col_logical(),
    ##   ins = col_logical(),
    ##   del = col_logical(),
    ##   annotation = col_logical()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

    ## Warning: 339 parsing failures.
    ##   row col           expected          actual                                                                               file
    ##  8114 ins 1/0/T/F/TRUE/FALSE 8114:8120:0.083 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Tirant.deviaTE.txt'
    ## 26174 del 1/0/T/F/TRUE/FALSE 596:598:1.29    '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Tirant.deviaTE.txt'
    ## 26730 ins 1/0/T/F/TRUE/FALSE 1152:1154:1.075 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Tirant.deviaTE.txt'
    ## 52317 ins 1/0/T/F/TRUE/FALSE 1161:1164:1.602 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Tirant.deviaTE.txt'
    ## 77329 del 1/0/T/F/TRUE/FALSE 595:598:0.676   '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Tirant.deviaTE.txt'
    ## ..... ... .................. ............... ..................................................................................
    ## See problems(...) for more details.

``` r
Iele <- read_deviate("/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/I-ele.deviaTE.txt", metadata)
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   TE = col_character(),
    ##   Sample = col_character(),
    ##   ref = col_character(),
    ##   snp = col_logical(),
    ##   refsnp = col_logical(),
    ##   int_del = col_logical(),
    ##   int_del_freq = col_logical(),
    ##   ins = col_character(),
    ##   del = col_character(),
    ##   annotation = col_character()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

``` r
souslik <- read_deviate("/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Souslik.deviaTE.txt", metadata)
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   TE = col_character(),
    ##   Sample = col_character(),
    ##   ref = col_character(),
    ##   snp = col_logical(),
    ##   refsnp = col_logical(),
    ##   int_del = col_logical(),
    ##   int_del_freq = col_logical(),
    ##   ins = col_logical(),
    ##   del = col_logical(),
    ##   annotation = col_logical()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

    ## Warning: 923 parsing failures.
    ##   row col           expected          actual                                                                                file
    ##  4216 del 1/0/T/F/TRUE/FALSE 4216:4218:0.386 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Souslik.deviaTE.txt'
    ##  4555 del 1/0/T/F/TRUE/FALSE 4555:4557:1.074 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Souslik.deviaTE.txt'
    ##  9491 del 1/0/T/F/TRUE/FALSE 4216:4218:0.537 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Souslik.deviaTE.txt'
    ##  9830 del 1/0/T/F/TRUE/FALSE 4555:4557:1.285 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Souslik.deviaTE.txt'
    ## 14766 del 1/0/T/F/TRUE/FALSE 4216:4218:0.471 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Souslik.deviaTE.txt'
    ## ..... ... .................. ............... ...................................................................................
    ## See problems(...) for more details.

``` r
transib <- read_deviate("/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Transib1.deviaTE.txt", metadata)
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   TE = col_character(),
    ##   Sample = col_character(),
    ##   ref = col_character(),
    ##   snp = col_logical(),
    ##   refsnp = col_logical(),
    ##   int_del = col_logical(),
    ##   int_del_freq = col_logical(),
    ##   ins = col_character(),
    ##   del = col_logical(),
    ##   annotation = col_logical()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

    ## Warning: 1380 parsing failures.
    ##  row col           expected                          actual                                                                                 file
    ## 2602 del 1/0/T/F/TRUE/FALSE 2602:2609:0.386                 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Transib1.deviaTE.txt'
    ## 2697 del 1/0/T/F/TRUE/FALSE 2697:2699:0.496                 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Transib1.deviaTE.txt'
    ## 2917 del 1/0/T/F/TRUE/FALSE 2917:2919:0.441,2917:2920:0.441 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Transib1.deviaTE.txt'
    ## 4175 del 1/0/T/F/TRUE/FALSE 1145:1151:0.345                 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Transib1.deviaTE.txt'
    ## 5632 del 1/0/T/F/TRUE/FALSE 2602:2609:0.192                 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Transib1.deviaTE.txt'
    ## .... ... .................. ............................... ....................................................................................
    ## See problems(...) for more details.

``` r
te412 <- read_deviate("/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/412.deviaTE.txt", metadata)
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   Sample = col_character(),
    ##   ref = col_character(),
    ##   snp = col_logical(),
    ##   refsnp = col_logical(),
    ##   int_del = col_logical(),
    ##   int_del_freq = col_logical(),
    ##   ins = col_logical(),
    ##   del = col_logical(),
    ##   annotation = col_character()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

    ## Warning: 278 parsing failures.
    ##    row col           expected          actual                                                                            file
    ##  30620 ins 1/0/T/F/TRUE/FALSE 352:354:5.051   '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/412.deviaTE.txt'
    ##  37673 ins 1/0/T/F/TRUE/FALSE 7405:7407:4.473 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/412.deviaTE.txt'
    ##  45770 ins 1/0/T/F/TRUE/FALSE 368:370:3.524   '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/412.deviaTE.txt'
    ## 113359 ins 1/0/T/F/TRUE/FALSE 7421:7423:1.854 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/412.deviaTE.txt'
    ## 113873 ins 1/0/T/F/TRUE/FALSE 368:370:2.26    '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/412.deviaTE.txt'
    ## ...... ... .................. ............... ...............................................................................
    ## See problems(...) for more details.

``` r
micropia <- read_deviate("/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Micropia.deviaTE.txt", metadata)
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_logical(),
    ##   TE = col_character(),
    ##   Sample = col_character(),
    ##   pos = col_double(),
    ##   ref = col_character(),
    ##   A = col_double(),
    ##   C = col_double(),
    ##   G = col_double(),
    ##   T = col_double(),
    ##   cov = col_double(),
    ##   phys_cov = col_double(),
    ##   hq_cov = col_double()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

    ## Warning: 8900 parsing failures.
    ##  row         col           expected actual                                                                                 file
    ## 2050 trunc_left  1/0/T/F/TRUE/FALSE  0.523 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Micropia.deviaTE.txt'
    ## 2078 trunc_rigth 1/0/T/F/TRUE/FALSE  0.028 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Micropia.deviaTE.txt'
    ## 2092 trunc_left  1/0/T/F/TRUE/FALSE  0.331 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Micropia.deviaTE.txt'
    ## 2098 trunc_left  1/0/T/F/TRUE/FALSE  0.028 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Micropia.deviaTE.txt'
    ## 2117 trunc_rigth 1/0/T/F/TRUE/FALSE  0.496 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Micropia.deviaTE.txt'
    ## .... ........... .................. ...... ....................................................................................
    ## See problems(...) for more details.

``` r
spoink <- read_deviate("/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Spoink.deviaTE.txt", metadata)
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_logical(),
    ##   TE = col_character(),
    ##   Sample = col_character(),
    ##   pos = col_double(),
    ##   ref = col_character(),
    ##   A = col_double(),
    ##   C = col_double(),
    ##   G = col_double(),
    ##   T = col_double(),
    ##   cov = col_double(),
    ##   phys_cov = col_double(),
    ##   hq_cov = col_double()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

    ## Warning: 23267 parsing failures.
    ##  row         col           expected actual                                                                               file
    ## 1530 trunc_left  1/0/T/F/TRUE/FALSE  0.055 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Spoink.deviaTE.txt'
    ## 1998 trunc_left  1/0/T/F/TRUE/FALSE  0.083 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Spoink.deviaTE.txt'
    ## 2061 trunc_rigth 1/0/T/F/TRUE/FALSE  0.165 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Spoink.deviaTE.txt'
    ## 2114 trunc_left  1/0/T/F/TRUE/FALSE  1.157 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Spoink.deviaTE.txt'
    ## 2190 trunc_rigth 1/0/T/F/TRUE/FALSE  0.441 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Spoink.deviaTE.txt'
    ## .... ........... .................. ...... ..................................................................................
    ## See problems(...) for more details.

``` r
opus <- read_deviate("/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Opus.deviaTE.txt", metadata)
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   TE = col_character(),
    ##   Sample = col_character(),
    ##   ref = col_character(),
    ##   snp = col_logical(),
    ##   refsnp = col_logical(),
    ##   int_del = col_logical(),
    ##   int_del_freq = col_logical(),
    ##   ins = col_logical(),
    ##   del = col_character(),
    ##   annotation = col_logical()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

    ## Warning: 14 parsing failures.
    ##     row col           expected          actual                                                                             file
    ## 2580813 ins 1/0/T/F/TRUE/FALSE 1110:1112:1.571 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Opus.deviaTE.txt'
    ## 3542657 ins 1/0/T/F/TRUE/FALSE 266:272:3.453   '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Opus.deviaTE.txt'
    ## 3549660 ins 1/0/T/F/TRUE/FALSE 7269:7275:3.275 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Opus.deviaTE.txt'
    ## 3550178 ins 1/0/T/F/TRUE/FALSE 266:272:6.576   '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Opus.deviaTE.txt'
    ## 3557181 ins 1/0/T/F/TRUE/FALSE 7269:7275:6.576 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Opus.deviaTE.txt'
    ## ....... ... .................. ............... ................................................................................
    ## See problems(...) for more details.

``` r
hobo <- read_deviate("/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Hobo.deviaTE.txt", metadata)
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   TE = col_character(),
    ##   Sample = col_character(),
    ##   ref = col_character(),
    ##   snp = col_logical(),
    ##   refsnp = col_logical(),
    ##   int_del = col_logical(),
    ##   int_del_freq = col_logical(),
    ##   ins = col_character(),
    ##   del = col_logical(),
    ##   annotation = col_logical()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

    ## Warning: 14 parsing failures.
    ##   row col           expected          actual                                                                             file
    ##  1184 del 1/0/T/F/TRUE/FALSE 1184:1189:0.468 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Hobo.deviaTE.txt'
    ##  8625 del 1/0/T/F/TRUE/FALSE 2707:2712:3.925 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Hobo.deviaTE.txt'
    ## 11584 del 1/0/T/F/TRUE/FALSE 2707:2712:3.01  '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Hobo.deviaTE.txt'
    ## 14543 del 1/0/T/F/TRUE/FALSE 2707:2712:5.772 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Hobo.deviaTE.txt'
    ## 17502 del 1/0/T/F/TRUE/FALSE 2707:2712:5.591 '/Volumes/Storage/dmel-full-story/diagnostic_snps/deviaTE_SNPs/Hobo.deviaTE.txt'
    ## ..... ... .................. ............... ................................................................................
    ## See problems(...) for more details.

## Define diagnostic SNPs

``` r
find_diagnostic <- function(old_sample, new_sample, deviaTE, output){
  
snps_new <- deviaTE %>% filter(snp==TRUE, Sample==new_sample) %>%
   mutate(max = pmax(A, C, G, T, na.rm = TRUE), new_allele = case_when(max==A ~ "A", max==C ~ "C", max==T ~ "T", max==G ~ "G"), new_maf = max/cov) %>% select(TE, pos, new_allele, new_maf)

snps_old <- deviaTE %>% filter(Sample==old_sample, pos %in% snps_new$pos) %>%
    mutate(max = pmax(A, C, G, T, na.rm = TRUE), old_allele = case_when(max==A ~ "A", max==C ~ "C", max==T ~ "T", max==G ~ "G"), old_maf = max/cov) %>% select(TE, pos, old_allele, old_maf)

diagnostic_snps <- inner_join(snps_new, snps_old, by=c("TE","pos")) %>% filter(new_allele!=old_allele, old_maf > 0.95)

toplot <- deviaTE %>% 
  inner_join(diagnostic_snps, by="pos") %>%
  filter(pos %in% diagnostic_snps$pos) %>%
  pivot_longer(c(A, C, G, T)) %>%
  rename(allele_cov = value, base = name) %>%
  mutate(prop = allele_cov / cov)

unique_pos <- unique(toplot$pos)

# If there are more than 20 unique "pos" values, sample 20 of them randomly for visualization
if (length(unique_pos) > 20) {
  set.seed(123)
  sampled_pos <- sample(unique_pos, 20)
  toplot <- toplot %>% filter(pos %in% sampled_pos)
}

plot_snps <- ggplot(toplot, aes(x = reorder(Sample, year), y = prop, fill = base)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ pos, ncol = 1) +
  labs(x = "year", y = "allele frequency", title = "diagnostic SNPs") +
  theme(legend.position = "top", axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(output, plot_snps, height = 20)

diagnostic_snps
}
```

``` r
(blood_diagnostic <- find_diagnostic("SRR23876563", "SRR8494428", blood, "/Volumes/Storage/dmel-full-story/diagnostic_snps/blood.png"))
```

    ## Saving 7 x 20 in image

    ## Warning: Removed 4 rows containing missing values (`position_stack()`).

    ## # A tibble: 66 × 6
    ##    TE      pos new_allele new_maf old_allele old_maf
    ##    <chr> <dbl> <chr>        <dbl> <chr>        <dbl>
    ##  1 Blood    77 G            0.794 A                1
    ##  2 Blood   117 T            0.788 C                1
    ##  3 Blood   133 C            0.771 T                1
    ##  4 Blood   229 G            0.852 A                1
    ##  5 Blood   425 T            0.848 C                1
    ##  6 Blood   456 T            0.844 C                1
    ##  7 Blood   471 G            0.861 A                1
    ##  8 Blood   472 G            0.857 A                1
    ##  9 Blood   479 T            0.854 G                1
    ## 10 Blood   494 G            0.837 A                1
    ## # ℹ 56 more rows

``` r
(tirant_diagnostic <- find_diagnostic("ERR6474638", "SRR8494428", tirant, "/Volumes/Storage/dmel-full-story/diagnostic_snps/tirant.png"))
```

    ## Saving 7 x 20 in image

    ## Warning: Removed 384 rows containing missing values (`position_stack()`).

    ## # A tibble: 90 × 6
    ##    TE       pos new_allele new_maf old_allele old_maf
    ##    <chr>  <dbl> <chr>        <dbl> <chr>        <dbl>
    ##  1 Tirant    26 C            0.755 T                1
    ##  2 Tirant  1225 C            0.792 T                1
    ##  3 Tirant  2362 C            0.800 T                1
    ##  4 Tirant  2390 T            0.755 C                1
    ##  5 Tirant  2395 C            0.737 T                1
    ##  6 Tirant  2403 A            0.748 G                1
    ##  7 Tirant  2413 G            0.741 T                1
    ##  8 Tirant  2421 T            0.750 A                1
    ##  9 Tirant  2425 C            0.740 T                1
    ## 10 Tirant  2446 T            0.748 A                1
    ## # ℹ 80 more rows

``` r
(Iele_diagnostic <- find_diagnostic("ERR6474638", "SRR8494428", Iele, "/Volumes/Storage/dmel-full-story/diagnostic_snps/Iele.png"))
```

    ## Saving 7 x 20 in image

    ## Warning: Removed 24 rows containing missing values (`position_stack()`).

    ## # A tibble: 98 × 6
    ##    TE        pos new_allele new_maf old_allele old_maf
    ##    <chr>   <dbl> <chr>        <dbl> <chr>        <dbl>
    ##  1 DMIFACA    37 T            0.527 A            0.997
    ##  2 DMIFACA   131 T            0.531 C            0.966
    ##  3 DMIFACA   401 T            0.588 C            0.956
    ##  4 DMIFACA   448 A            0.504 T            0.998
    ##  5 DMIFACA  1115 C            0.589 T            1    
    ##  6 DMIFACA  1133 A            0.552 C            0.954
    ##  7 DMIFACA  1156 A            0.513 C            0.974
    ##  8 DMIFACA  1198 A            0.515 G            0.998
    ##  9 DMIFACA  1361 T            0.509 G            0.970
    ## 10 DMIFACA  1403 A            0.619 G            0.964
    ## # ℹ 88 more rows

``` r
(souslik_diagnostic <- find_diagnostic("ERR6474638", "SRR8494428", souslik, "/Volumes/Storage/dmel-full-story/diagnostic_snps/souslik.png"))
```

    ## Saving 7 x 20 in image

    ## Warning: Removed 136 rows containing missing values (`position_stack()`).

    ## # A tibble: 181 × 6
    ##    TE        pos new_allele new_maf old_allele old_maf
    ##    <chr>   <dbl> <chr>        <dbl> <chr>        <dbl>
    ##  1 Souslik   244 A            0.800 G                1
    ##  2 Souslik   310 C            0.797 G                1
    ##  3 Souslik   328 G            0.803 T                1
    ##  4 Souslik   329 A            0.800 G                1
    ##  5 Souslik   333 A            0.803 G                1
    ##  6 Souslik   402 A            0.815 G                1
    ##  7 Souslik   415 G            0.786 T                1
    ##  8 Souslik   445 C            0.758 A                1
    ##  9 Souslik   450 T            0.758 G                1
    ## 10 Souslik   464 A            0.766 G                1
    ## # ℹ 171 more rows

``` r
(transib_diagnostic <- find_diagnostic("ERR6474638", "SRR8494428", transib, "/Volumes/Storage/dmel-full-story/diagnostic_snps/transib.png"))
```

    ## Saving 7 x 20 in image

    ## Warning: Removed 96 rows containing missing values (`position_stack()`).

    ## # A tibble: 21 × 6
    ##    TE         pos new_allele new_maf old_allele old_maf
    ##    <chr>    <dbl> <chr>        <dbl> <chr>        <dbl>
    ##  1 Transib1   182 A            0.853 C            1    
    ##  2 Transib1  1119 C            0.878 A            1    
    ##  3 Transib1  1124 T            0.882 A            1    
    ##  4 Transib1  1430 A            0.824 G            1    
    ##  5 Transib1  1438 A            0.832 G            1    
    ##  6 Transib1  1554 A            0.829 T            1    
    ##  7 Transib1  1806 T            0.869 C            1    
    ##  8 Transib1  1848 G            0.815 A            0.988
    ##  9 Transib1  2168 T            0.852 G            1    
    ## 10 Transib1  2288 G            0.853 A            1    
    ## # ℹ 11 more rows

``` r
(te412_diagnostic <- find_diagnostic("SRR23876563", "SRR8494428", te412, "/Volumes/Storage/dmel-full-story/diagnostic_snps/412.png"))
```

    ## Saving 7 x 20 in image

    ## Warning: Removed 84 rows containing missing values (`position_stack()`).

    ## # A tibble: 30 × 6
    ##       TE   pos new_allele new_maf old_allele old_maf
    ##    <dbl> <dbl> <chr>        <dbl> <chr>        <dbl>
    ##  1   412  2608 A            0.871 G                1
    ##  2   412  2629 G            0.868 T                1
    ##  3   412  2660 C            0.885 T                1
    ##  4   412  2685 C            0.898 A                1
    ##  5   412  2711 A            0.884 T                1
    ##  6   412  3575 C            0.891 T                1
    ##  7   412  3788 C            0.890 T                1
    ##  8   412  3861 C            0.892 T                1
    ##  9   412  3912 C            0.884 A                1
    ## 10   412  4034 A            0.886 G                1
    ## # ℹ 20 more rows

``` r
#(opus_diagnostic <- find_diagnostic("SRR23876563", "SRR8494428", opus, "/Volumes/Storage/dmel-full-story/diagnostic_snps/opus.png"))
(micropia_diagnostic <- find_diagnostic("ERR6474638", "SRR8494428", micropia, "/Volumes/Storage/dmel-full-story/diagnostic_snps/micropia.png"))
```

    ## Saving 7 x 20 in image

    ## Warning: Removed 144 rows containing missing values (`position_stack()`).

    ## # A tibble: 14 × 6
    ##    TE         pos new_allele new_maf old_allele old_maf
    ##    <chr>    <dbl> <chr>        <dbl> <chr>        <dbl>
    ##  1 Micropia  2066 T            0.881 C            0.986
    ##  2 Micropia  2075 T            0.878 C            0.986
    ##  3 Micropia  2085 T            0.832 C            1    
    ##  4 Micropia  2087 C            0.834 T            1    
    ##  5 Micropia  2090 T            0.829 A            1    
    ##  6 Micropia  2096 T            0.782 G            1    
    ##  7 Micropia  2117 T            0.803 C            1    
    ##  8 Micropia  2120 G            0.806 A            1    
    ##  9 Micropia  2123 T            0.806 C            1    
    ## 10 Micropia  2131 A            0.784 T            1    
    ## 11 Micropia  2132 C            0.786 T            1    
    ## 12 Micropia  2153 C            0.836 T            1    
    ## 13 Micropia  2156 A            0.836 C            0.992
    ## 14 Micropia  2171 T            0.899 A            1

``` r
(spoink_diagnostic <- find_diagnostic("ERR6474638", "SRR8494428", spoink, "/Volumes/Storage/dmel-full-story/diagnostic_snps/spoink.png"))
```

    ## Saving 7 x 20 in image

    ## Warning: Removed 320 rows containing missing values (`position_stack()`).

    ## # A tibble: 97 × 6
    ##    TE       pos new_allele new_maf old_allele old_maf
    ##    <chr>  <dbl> <chr>        <dbl> <chr>        <dbl>
    ##  1 Spoink  2130 G            0.786 T            1    
    ##  2 Spoink  2133 C            0.783 T            1    
    ##  3 Spoink  2139 T            0.758 C            1    
    ##  4 Spoink  2145 C            0.751 T            1    
    ##  5 Spoink  2163 A            0.754 G            1    
    ##  6 Spoink  2181 T            0.850 C            1    
    ##  7 Spoink  2183 A            0.851 G            1    
    ##  8 Spoink  2232 T            0.889 C            1    
    ##  9 Spoink  2238 T            0.857 A            1    
    ## 10 Spoink  2243 C            0.838 A            0.986
    ## # ℹ 87 more rows

``` r
(hobo_diagnostic <- find_diagnostic("ERR6474638", "SRR8494428", hobo, "/Volumes/Storage/dmel-full-story/diagnostic_snps/hobo.png"))
```

    ## Saving 7 x 20 in image

    ## Warning: Removed 28 rows containing missing values (`position_stack()`).

    ## # A tibble: 15 × 6
    ##    TE      pos new_allele new_maf old_allele old_maf
    ##    <chr> <dbl> <chr>        <dbl> <chr>        <dbl>
    ##  1 Hobo    359 C            0.831 G            1    
    ##  2 Hobo    471 G            0.861 A            1    
    ##  3 Hobo    509 A            0.840 T            1    
    ##  4 Hobo    736 T            0.874 A            1    
    ##  5 Hobo    845 T            0.898 A            0.996
    ##  6 Hobo    978 A            0.801 G            1    
    ##  7 Hobo   1011 C            0.695 G            1    
    ##  8 Hobo   1266 T            0.672 A            1    
    ##  9 Hobo   1458 C            0.662 G            1    
    ## 10 Hobo   2129 A            0.820 T            1    
    ## 11 Hobo   2138 T            0.804 A            1    
    ## 12 Hobo   2171 A            0.781 G            1    
    ## 13 Hobo   2210 G            0.788 A            1    
    ## 14 Hobo   2261 C            0.815 T            1    
    ## 15 Hobo   2382 A            0.764 T            0.993

``` r
calculate_oldness <- function(deviaTE, diagnostic){
  
  prop <- deviaTE %>% filter(pos %in% diagnostic$pos) %>% pivot_longer(c(A, C, G, T)) %>% rename(allele_cov=value, new_allele=name) %>% mutate(prop = allele_cov/cov) %>% inner_join(diagnostic, by=c("TE", "pos", "new_allele"))
  
  (by_sample <- prop %>% group_by(TE, Sample) %>% summarise(new_alleles = mean(prop)) %>% mutate(new_alleles = if_else(is.nan(new_alleles), 0, new_alleles)) %>% ungroup())
}

blood_oldness <- calculate_oldness(blood, blood_diagnostic)
```

    ## `summarise()` has grouped output by 'TE'. You can override using the `.groups`
    ## argument.

``` r
tirant_oldness <- calculate_oldness(tirant, tirant_diagnostic)
```

    ## `summarise()` has grouped output by 'TE'. You can override using the `.groups`
    ## argument.

``` r
Iele_oldness <- calculate_oldness(Iele, Iele_diagnostic)
```

    ## `summarise()` has grouped output by 'TE'. You can override using the `.groups`
    ## argument.

``` r
souslik_oldness <- calculate_oldness(souslik, souslik_diagnostic)
```

    ## `summarise()` has grouped output by 'TE'. You can override using the `.groups`
    ## argument.

``` r
transib_oldness <- calculate_oldness(transib, transib_diagnostic)
```

    ## `summarise()` has grouped output by 'TE'. You can override using the `.groups`
    ## argument.

``` r
te412_oldness <- calculate_oldness(te412, te412_diagnostic) %>% mutate(TE = as.character(TE))
```

    ## `summarise()` has grouped output by 'TE'. You can override using the `.groups`
    ## argument.

``` r
spoink_oldness <- calculate_oldness(spoink, spoink_diagnostic)
```

    ## `summarise()` has grouped output by 'TE'. You can override using the `.groups`
    ## argument.

``` r
micropia_oldness <- calculate_oldness(micropia, micropia_diagnostic)
```

    ## `summarise()` has grouped output by 'TE'. You can override using the `.groups`
    ## argument.

``` r
hobo_oldness <- calculate_oldness(hobo, hobo_diagnostic)
```

    ## `summarise()` has grouped output by 'TE'. You can override using the `.groups`
    ## argument.

``` r
oldness <- bind_rows(blood_oldness, tirant_oldness, Iele_oldness, souslik_oldness, transib_oldness, te412_oldness, spoink_oldness, micropia_oldness, hobo_oldness)
```

``` r
(invaders <- read_csv("/Volumes/Storage/dmel-full-story/analysis/csv/new-TEs-v4.csv", show_col_types = FALSE) %>% filter(Sample!="Sample") %>% type_convert())
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   Sample = col_character(),
    ##   TE = col_character(),
    ##   All_reads = col_double(),
    ##   HQ_reads = col_double()
    ## )

    ## # A tibble: 6,996 × 4
    ##    Sample     TE       All_reads HQ_reads
    ##    <chr>      <chr>        <dbl>    <dbl>
    ##  1 ERR6474638 412          52.3     42.6 
    ##  2 ERR6474638 Blood        28.6     25.1 
    ##  3 ERR6474638 DMIFACA      17.9     13.3 
    ##  4 ERR6474638 Hobo          4.06     2.13
    ##  5 ERR6474638 Micropia      0.1      0   
    ##  6 ERR6474638 Opus         22.7     18.8 
    ##  7 ERR6474638 PPI251        0        0   
    ##  8 ERR6474638 Shellder      0.13     0.01
    ##  9 ERR6474638 Souslik       3.15     2.1 
    ## 10 ERR6474638 Spoink        0.83     0.1 
    ## # ℹ 6,986 more rows

``` r
(final_dataset <- left_join(invaders, oldness, by=c("Sample","TE")) %>% mutate(reinvasion = ifelse(TE %in% c("PPI251","Opus"), FALSE, TRUE), TE=case_when(TE=="DMIFACA"~"I-element", TE=="PPI251" ~ "P-element", TRUE ~ TE)) %>% inner_join(metadata, by="Sample") %>% filter(TE!="Shellder"))
```

    ## # A tibble: 6,413 × 14
    ##    Sample     TE    All_reads HQ_reads new_alleles reinvasion strain publication
    ##    <chr>      <chr>     <dbl>    <dbl>       <dbl> <lgl>      <chr>  <chr>      
    ##  1 ERR6474638 412       52.3     42.6     0.956    TRUE       Orego… https://do…
    ##  2 ERR6474638 Blood     28.6     25.1     0.863    TRUE       Orego… https://do…
    ##  3 ERR6474638 I-el…     17.9     13.3     0.00531  TRUE       Orego… https://do…
    ##  4 ERR6474638 Hobo       4.06     2.13    0.000458 TRUE       Orego… https://do…
    ##  5 ERR6474638 Micr…      0.1      0       0.00266  TRUE       Orego… https://do…
    ##  6 ERR6474638 Opus      22.7     18.8    NA        FALSE      Orego… https://do…
    ##  7 ERR6474638 P-el…      0        0      NA        FALSE      Orego… https://do…
    ##  8 ERR6474638 Sous…      3.15     2.1     0.00231  TRUE       Orego… https://do…
    ##  9 ERR6474638 Spoi…      0.83     0.1     0.000187 TRUE       Orego… https://do…
    ## 10 ERR6474638 Tira…      0.48     0.08    0.00107  TRUE       Orego… https://do…
    ## # ℹ 6,403 more rows
    ## # ℹ 6 more variables: study <chr>, study_id <chr>, year <dbl>, location <chr>,
    ## #   lat <dbl>, lon <dbl>

``` r
write_tsv(final_dataset, "/Volumes/Storage/dmel-full-story/dataset.tsv")
```
