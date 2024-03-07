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
metadata <- read_tsv("/Volumes/Storage/dmel-full-story/metadata.tsv")
```

    ## Rows: 585 Columns: 9
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (6): Sample, strain, publication, study, study_id, location
    ## dbl (3): year, lat, lon
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

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

# If there are more than 20 unique "pos" values, sample 20 of them randomly
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
