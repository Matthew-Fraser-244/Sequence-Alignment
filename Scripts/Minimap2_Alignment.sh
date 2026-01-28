# Minimap2 Alignment
```{bash}
# Initialise Conda
source ~/.bashrc
conda init

# Specify WD for Files
cd /Users/matthewfraser/Desktop/Matthew_data/Sequence_Alignment/Example_Data

# If needed, use SAMtools to extract different sequences into own fasta
# samtools faidx .fasta "" > .fasta

# If needed, use seqtk to convert the existing fasta into the reverse complement
# seqtk seq -r .fasta > .fasta

# Use cat command to combine multiple fasta files into a single new fasta for alignment
# cat .fasta .fasta > .fasta

# Create a path to the FASTA file
fasta_file="ANNOTATED_SANCTUARY.fasta"

# Create PAF file with Alignments
# Change fasta file name and output paf file
echo "Starting minimap2 alignment"
minimap2 -X -N 5 -p 0.5 -c ANNOTATED_SANCTUARY.fasta ANNOTATED_SANCTUARY.fasta > ANNOTATED_SANCTUARY.paf

```


```{r}
# Load required R Libraries
library(gggenomes)
library(ggplot2)
library(tidyverse)
library(ggnewscale)
library(dplyr)

# Specify WD for Files, Ensure Same As Before
setwd("/Users/matthewfraser/Desktop/Matthew_data/Sequence_Alignment/Example_Data")

# Read in minimap2 links from specified .paf file
# Filter out alignments less than 2200bp for clarity
minimap2_links <- read_paf("/Users/matthewfraser/Desktop/Matthew_data/Sequence_Alignment/Example_Data/ANNOTATED_SANCTUARY.paf")
minimap2_links <-minimap2_links %>% filter(map_length>=2200)

# If Alignment Involves GFF files (Eg. for Gene Annotations), processes and filters data into dataframe for use in sequence alignment
filelist<-list.files(pattern = "\\.gff")
gfflist <-lapply(filelist, function(file) {
  read.delim(file, comment.char = "#", header = FALSE)})

# If needed, can specify gff files manually and furhter specify sequences within file
# gff_CS10 <- read.delim("~/Desktop/Matthew_data/Sequence_Alignment/Example_Data/ANNOTATED_SANCTUARY.gff", header = FALSE, comment.char = "#", sep = "\t", stringsAsFactors = FALSE) %>%
#  filter(V1 %in% c("CS10_Chr08_Sanc"))

# Combine the gffs into one and assign column and gene names
# Make sure to specify whether using gfflist or individual gff files in list()
comb_gff<-do.call(rbind, gfflist)
# comb_gff<-do.call(rbind,  list())
colnames(comb_gff) <- c("chr","source","type","start","end","dot","strand","dot2","Feature")

# Extract annotation names that in the "Features" column of comb_gff but not given distinct annotation types
comb_gff <- comb_gff %>%
 mutate(Name = str_extract(Feature, "Name=[^;]+"),
         Name = str_replace(Name, "Name=", "")) %>%
  mutate(type = ifelse(type == "gene" & !is.na(Name), Name, type)) %>%
  select(-Name)

# Specify gene annotations present within "type" column of comb_gff dataframe for plotting
# If some of the annotations have slight name variations, can allow these to be recognised as single type
# Can specify priority, to dictate which of any overlapping annotations are shown over others
genes<- comb_gff %>%
  transmute(
    seq_id=chr,
    start= start,
    end = end,
    strand= strand,
    feat_id = type)%>%
filter(feat_id %in% c(
    "YR", "HET", "ClassII", "NLR", "Patatin-PLP", "TPR-repeat",
    "ClassI", "DUF3723", "Secondary Metabolite", "DUF3435", "ToxhAT", "primer_bind", "gene", "HP", "IR", "TIR", "HaT_Transposase", "ToxA", "FabD/lysophospholipase-like protein", "Tc5 Transposase", "TPR_12", "Patatin", "TOXA; gene", "TOXA%3B gene"
  )) %>%
mutate(
    feat_id = ifelse(feat_id == "Toxin_ToxA", "ToxA", feat_id),
    feat_id = ifelse(feat_id == "ToxhAT", "ToxTA", feat_id),
    feat_id = ifelse(feat_id == "TOXA; gene", "ToxA", feat_id),
    feat_id = ifelse(feat_id == "TOXA%3B gene", "ToxA", feat_id),
    feat_id = ifelse(feat_id == "Tc5", "Tc5 Transposase", feat_id),
    feat_id = ifelse(feat_id == "hp", "HP", feat_id),
    feat_id = ifelse(feat_id == "ir", "IR", feat_id),
    feat_id = ifelse(feat_id == "Patatin", "Patatin-PLP", feat_id),
    feat_id = ifelse(feat_id == "TPR_12", "TPR-repeat", feat_id),
    priority = case_when(
      feat_id == "IR" ~ 99,
      feat_id == "TIR" ~ 90,
      feat_id == "ToxA" ~ 85,
      feat_id == "gene" ~ 1,  
      TRUE ~ 50)
  ) %>%
  arrange(priority)

# Read the FASTA file into a seq object, ensure name is changed
# Inlcude only sequences filtered that link with both Sanctuary and Horizon
Sanc_seq <- read_seqs("/Users/matthewfraser/Desktop/Matthew_data/Sequence_Alignment/Example_Data/ANNOTATED_SANCTUARY.fasta")

```

# Minimap2 Plotting: Aligning Multiple WGS
```{r}
# For Minimap2 is the same format as the plotting for LASTZ except there is no perc_id parameter so cannot create links for % identity.

# Create bins for each sequence and put into single dataframe
# This allows to edit order, names, and group multiple sequences in the same bin for adjacent alignment against a single sequence.
Sanc_bins <- tibble(
  seq_id = c(
  "CS10_Chr08_Sanc",
  "WAI2411_Chr02_Sanc",
  "WAI2431_Chr07_Sanc",
  "WAI2432_tig045_chr08_chr07_Sanc",
  "WAI3285_Chr02_Sanc",
  "WAI3295_tig122_Sanc",
  "WAI3382_Chr06_Sanc",
  "WAI3384_Chr16_Sanc",
  "WAI3398_Chr04_Sanc")
) %>%
  mutate(
    bin_id = case_when(
      seq_id == "CS10_Chr08_Sanc"                  ~ "Sanctuary (CS10)",
      seq_id == "WAI2411_Chr02_Sanc"              ~ "Sanctuary (WAI2411)",
      seq_id == "WAI2431_Chr07_Sanc"              ~ "Sanctuary (WAI2431)",
      seq_id == "WAI2432_tig045_chr08_chr07_Sanc"  ~ "Sanctuary (WAI2432)",
      seq_id == "WAI3285_Chr02_Sanc"              ~ "Sanctuary (WAI3285)",
      seq_id == "WAI3295_tig122_Sanc"              ~ "Sanctuary (WAI3295)",
      seq_id == "WAI3382_Chr06_Sanc"              ~ "Sanctuary (WAI3382)",
      seq_id == "WAI3384_Chr16_Sanc"              ~ "Sanctuary (WAI3384)",
      seq_id == "WAI3398_Chr04_Sanc"              ~ "Sanctuary (WAI3398)",
      TRUE ~ seq_id 
    )
  ) %>%
  left_join(Sanc_seq %>% select(seq_id, length), by = "seq_id")

# Arrange bins for specific plotting order
Sanc_bins_ordered1 <- Sanc_bins %>%
  mutate(bin_id = factor(bin_id, levels = 
c("Sanctuary (CS10)", "Sanctuary (WAI2411)", "Sanctuary (WAI2431)", "Sanctuary (WAI2432)", "Sanctuary (WAI3285)", "Sanctuary (WAI3295)", "Sanctuary (WAI3382)", "Sanctuary (WAI3384)", "Sanctuary (WAI3398)"
)
)) %>% filter(!is.na(bin_id)) %>%
  arrange(bin_id)

# Create gggenomes plot using Minimap2 links based on strand direction
minimap2_plot <- gggenomes(seqs = Sanc_bins_ordered1, links = minimap2_links, genes = genes) + geom_seq() + geom_bin_label(expand_left = .4, size = 4) +
  geom_link(aes(fill = strand), color = NA, alpha = 0.2, linewidth = 1) +
  new_scale_fill() +
  geom_gene(aes(fill = feat_id), color = NA, size = 3) +
  scale_fill_manual(values = c(
  "YR" = "maroon","HET" = "skyblue","NLR" = "purple","Patatin-PLP" = "yellow","TPR-repeat" = "darkgreen","ClassII" = "blue","ClassI" = "orange","ToxTA" = "red","Secondary Metabolite" = "#f0e68c","DUF3723" = "#cd653f","DUF3435" = "#29AB87","SDR" = "limegreen","gene" = "grey", "TIR" = "#D1E231", "IR"= "salmon", "HP"= "darkblue", "HaT_Transposase" = "#FE6633", "ToxA" = "turquoise", "FabD/lysophospholipase-like protein" = "beige", "Tc5 Transposase" = "#B8B0FF")) 
print(minimap2_plot)

```
