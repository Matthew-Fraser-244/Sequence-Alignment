# Sequence-Alignment
A collection of scripts for aligning and plotting DNA sequences using Lastz, Minimap2 and gggenomes. 

These can produce high-quality synteny maps between sequences, showcasing links in percent identity and strand direction alignment. Annotations present within gff files can also be shown. 

These scripts can be utilised in some of the following analysis scenarios:
* Plotting one or multiple aligned whole-genome sequences to another.
* Creating multiple plots of a single subject sequence to mutliple different query sequences.
* Plotting a specific zoomed-in region of one or multiple aligned whole-genome sequences to another. 

Additional script has been included should you want to determine the AT-richness distribution for the sequences in the alignments and then include this in the plot. An R notebook containing all scripts is also included. 

Example fasta and gff files have been provided from the following publication: https://doi.org/10.1128/mbio.01371-25
The scripts have been written using this data to align the giant mobile <i>Starship</i> elements from several fungal isolates together.
