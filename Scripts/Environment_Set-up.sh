```{bash}
# Initialise Conda
source ~/.bashrc
conda init

# Create conda environment called GenomeAlign
conda create -n GenomeAlign
conda activate GenomeAlign

# Install LASTZ, Minimap2, SAMtools and seqtk from Bioconda
conda install -c Bioconda lastz minimap2 samtools seqtk

conda deactivate

# Create conda environment for bedtools
conda create -n bedtools
conda activate bedtools

conda install -c Bioconda bedtools

```
