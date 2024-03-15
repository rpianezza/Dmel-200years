Drosophila melanogaster - TE trees annotation
================

``` r
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ape))
theme_set(theme_bw())
```

From the RepeatMasker output folder (containing all the ori.out files),
I defragment the output and made it readable in R, concatenating all the
files and adding the column “species” based on the name of each file.

    for i in *ori.out; do python ../rm-defragmenter.py --rm $i --fai ../../ref/new-TEs-v4.fasta.fai --dist 250 > $i.def; done

    for file in *def; do
        filename=$(basename "$file" .fa.ori.out.def)
        awk -v filename="$filename" '{print filename "\t" $0}' "$file"
    done > /Volumes/Storage/dmel-full-story/RepeatMasker/output/all.fa.ori.out.def.tsv

Read the resulting file, add info for each TE (len in bp) and filter the
RM output by keeping only full length insertions with less than 10%
divergence from the consensus.

``` r
defragmented <- read_tsv("/Volumes/Storage/dmel-full-story/RepeatMasker/output/all.fa.ori.out.def.tsv", col_names = c("species","te","contig","strand","qstart","qend","pdiv","rstart","rend"))
```

    ## Rows: 1193036 Columns: 9
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (4): species, te, contig, strand
    ## dbl (5): qstart, qend, pdiv, rstart, rend
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
te_info <- read_tsv("/Volumes/Storage/dmel-full-story/ref/new-TEs-v4.fasta.fai", col_names = c("te","full_len","x","y","z")) %>% select(te, full_len)
```

    ## Rows: 12 Columns: 5
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): te
    ## dbl (4): full_len, x, y, z
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
full_len <- defragmented %>% inner_join(te_info, by="te") %>% mutate(len = qend-qstart, fraction = len/full_len, strand = ifelse(strand=="C", "-", strand)) %>% filter(pdiv < 10, fraction > 0.8, fraction < 1.25)

to_bed <- full_len %>% select(contig, qstart, qend, te, species, strand)

species_list <- to_bed %>% select(species) %>% distinct() %>% pull()
te_list <- to_bed %>% select(te) %>% distinct() %>% pull()
```

For each TE, create a subfolder and extract a bed file for each species
(each RM output file) in the TE subfolder.

``` r
for (t in te_list) {
  dir.create(paste0("/Volumes/Storage/dmel-full-story/TE-trees/trees-RM/bed/", t), showWarnings = FALSE)
}

for (s in species_list) {
  spec <- to_bed %>% filter(species == s)
  for (t in te_list) {
    transposon <- spec %>% filter(te == t)
    if (length(transposon$qstart)>1){
    path <- paste0("/Volumes/Storage/dmel-full-story/TE-trees/trees-RM/bed/", t, "/", s, ".bed")
    #write_tsv(transposon, path, col_names=FALSE)
    }
  }
}
```

Extract fasta from bed, rename fasta sequences according to the species
name, concatenate the fasta, perform MSA and convert to NEXUS.

    for bed_file in *.bed; do
        fasta_file="/Volumes/Storage/data/assemblies-droso/${bed_file%.bed}.fa"
        output_file="/Volumes/Storage/dmel-full-story/TE-trees/trees-RM/bed/PPI251/$(basename "$bed_file" .bed).fasta"
        bedtools getfasta -s -fi "$fasta_file" -bed "$bed_file" -fo "$output_file"
    done

    for file in *.fasta; do python /Volumes/Storage/dmel-full-story/TE-trees/rename-insertions.py "$file" "${file%.fasta}.renamed.fasta"; done

    cat *renamed.fasta > PPI251.fasta

    muscle -in PPI251.fasta -out PPI251.MSA

    seqret -auto -sequence PPI251.MSA -outseq PPI251.nexus -osformat nexus

The NEXUS file can be uploaded on BEAUti to get an XML file, which can
be then used as input for BEAST which will generate a .trees file. The
.trees file can be given as input to TreeAnnotator, which will extract a
consensus tree.

Here, I read the resulting trees, one for each of the 11 TEs. I also
read the colorcode file, containing info about how each species should
be colored in the tree.

For the Spoink tree, I had to modify the labels to format them correctly
(I used a previous version of the renaming script and thus the labels
were different).

``` r
blood <- read.nexus("/Volumes/Storage/dmel-full-story/TE-trees/trees-LTR/Blood.tree")
te412 <- read.nexus("/Volumes/Storage/dmel-full-story/TE-trees/trees-LTR/412-412.tree")
opus <- read.nexus("/Volumes/Storage/dmel-full-story/TE-trees/trees-LTR/Opus-Opus.tree")
tirant <- read.nexus("/Volumes/Storage/dmel-full-story/TE-trees/trees-LTR/Tirant.tree")
iele <- read.nexus("/Volumes/Storage/dmel-full-story/TE-trees/trees-RM/Iele-DMIFACA.tree")
hobo <- read.nexus("/Volumes/Storage/dmel-full-story/TE-trees/trees-RM/Hobo-Hobo.tree")
pele <- read.nexus("/Volumes/Storage/dmel-full-story/TE-trees/trees-RM/Pele-PPI251.tree")
spoink <- read.nexus("/Volumes/EXT-RICCARDO/DoubleTrouble/TE-trees/101+Petrov-refined/Spoink-refined.tree")
micropia <- read.nexus("/Volumes/Storage/dmel-full-story/TE-trees/trees-LTR/Micropia.tree")
souslik <- read.nexus("/Volumes/Storage/dmel-full-story/TE-trees/trees-LTR/Souslik.tree")
transib <- read.nexus("/Volumes/Storage/dmel-full-story/TE-trees/trees-RM/Transib1-Transib1.tree")

colorcode <- read_tsv("/Volumes/Storage/dmel-full-story/TE-trees/color-code-R-groups.txt") #%>% mutate(RGB = ifelse(group=="D.melanogaster group", "yellow", RGB))
```

    ## Rows: 51 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (4): species, label, group, RGB
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
colorcode_micropia <- read_tsv("/Volumes/Storage/dmel-full-story/TE-trees/color-code-R-group-micropia.txt")
```

    ## Rows: 51 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (4): species, label, group, RGB
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
labels_spoink <- spoink$tip.label %>% as_tibble() %>% mutate(value = substr(value, 2, nchar(value) - 1), value = ifelse(substr(value, 1, 1) != "D", paste0("D.", value), value), value = substr(value, 1, 5))
spoink$tip.label <- labels_spoink$value
```

Function to adjust the label of each tree, assign them to the respective
color as specified in the colorcode file and finally plot the tree, by
saving the image to the given path in pdf.

``` r
label_tree <- function(tree, path, cc) {
  labels <- tree$tip.label %>% as_tibble() %>% rename(label = value) %>% mutate(label = ifelse(substr(label, 2, 2) != ".", paste0(substr(label, 1, 1), ".", substr(label, 2, nchar(label))), label), label = substr(label, 1, 5)) %>% separate(label, into = c("genus", "spe"), sep="\\.") %>% mutate(spe = tolower(spe), label = paste0(genus, ".", spe)) %>% left_join(cc, by="label")

tree$tip.label <- labels$species

tip_label_colors <- labels$RGB

pdf(path, width = 21, height = 21)
(plot <- plot(tree, cex = 2, type = "fan", tip.color = tip_label_colors, node.pos = 1))
dev.off()

return(plot)
}
```

Run the function on all of our 11 trees.

``` r
blood_tree <- label_tree(blood, "/Volumes/Storage/dmel-full-story/TE-trees/blood.pdf", colorcode)
opus_tree <- label_tree(opus, "/Volumes/Storage/dmel-full-story/TE-trees/opus.pdf", colorcode)
te412_tree <- label_tree(te412, "/Volumes/Storage/dmel-full-story/TE-trees/412.pdf", colorcode)
tirant_tree <- label_tree(tirant, "/Volumes/Storage/dmel-full-story/TE-trees/tirant.pdf", colorcode)
iele_tree <- label_tree(iele, "/Volumes/Storage/dmel-full-story/TE-trees/iele.pdf", colorcode)
hobo_tree <- label_tree(hobo, "/Volumes/Storage/dmel-full-story/TE-trees/hobo.pdf", colorcode)
pele_tree <- label_tree(pele, "/Volumes/Storage/dmel-full-story/TE-trees/pele.pdf", colorcode)
spoink_tree <- label_tree(spoink, "/Volumes/Storage/dmel-full-story/TE-trees/spoink.pdf", colorcode)
micropia_tree <- label_tree(micropia, "/Volumes/Storage/dmel-full-story/TE-trees/micropia.pdf", colorcode_micropia)
souslik_tree <- label_tree(souslik, "/Volumes/Storage/dmel-full-story/TE-trees/souslik.pdf", colorcode)
transib_tree <- label_tree(transib, "/Volumes/Storage/dmel-full-story/TE-trees/transib.pdf", colorcode)
```
