# Including AT-Richness
```{bash}
# Set WD and activate Conda
cd /Users/matthewfraser/Desktop/Matthew_data/Sequence_Alignment/Example_Data

source ~/.bashrc
conda init
conda activate bedtools

# Specify input FASTA file
FASTA_INPUT="ANNOTATED_SANCTUARY.fasta"

# Extract base name for output files
NAME=$(basename "$FASTA_INPUT")
NAME="${NAME%.*}"

# Index the FASTA
samtools faidx "$FASTA_INPUT"

# Create 1 kb windows across all sequences
bedtools makewindows -g "${FASTA_INPUT}.fai" -w 1000 > "${NAME}_1kb_windows.bed"

# Calculate AT-richness for each window
bedtools nuc -fi "$FASTA_INPUT" -bed "${NAME}_1kb_windows.bed" | \
awk 'BEGIN{OFS="\t"} NR>1 {print $1, $2, $3, $4}' > "${NAME}_at_richness.bed"

```

```{r}
# Read in the AT-richness file into a feature object
AT_richness <- read_feats("/Users/matthewfraser/Desktop/Matthew_data/Sequence_Alignment/Example_Data/ANNOTATED_SANCTUARY_at_richness.bed",
                        col_names = c("seq_id", "start", "end", "score"))%>% 
                        mutate(score =as.numeric(score), feat_id='AT_cont')


# If using zoomed sequences add the following adjustment for AT richness coordinates
AT_richness_zoom <- AT_richness %>%
  left_join(Zoomed_Sanc, by = "seq_id") %>%
  mutate(
    start = if_else(!is.na(zoom_start), start - zoom_start, start),
    end   = if_else(!is.na(zoom_start), end   - zoom_start, end)
  ) %>%
  filter(
    (is.na(zoom_length)) |
    (start >= 0 & end <= zoom_length)
  ) %>%
  select(-zoom_start, -zoom_end, -zoom_length)

# Add the following to the gggenomes string in ggplot
gggenomes(feats =list(at=AT_richness_zoom))

```
