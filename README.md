# Extract rDNA sequences from fungal genomes

This utility combines the power of Barrnap — a fast and accurate tool to identify the location of ribosomal RNA genes — with ITSx — software that extracts rDNA sequences from genomic FASTA files — to quickly extract rDNA sequences from fungal genomes. 

For more information on these software:
* [Barrnap (BAsic Rapid Ribosomal RNA Predictor)](https://github.com/tseemann/barrnap)
* [ITSx](https://microbiology.se/software/itsx/)

This script was originally built by [Fantin Mesny](https://github.com/fantin-mesny/Extract-ITS-sequences-from-a-fungal-genome), patched by (Pepijn Kooij)[https://github.com/pwkooij/Extract-ITS-sequences-from-a-fungal-genome], and refactored and further modified by [Alden Dirks](https://github.com/aldendirks/Extract-ITS-sequences-from-a-fungal-genome). 

## Installation

Install dependencies with conda and clone the GitHub repository. Add the script to your `PATH` for easier execution. 

```
conda create -n extractITS_env barrnap biopython itsx pandas
conda activate extractITS_env
git clone https://github.com/aldendirks/Extract-ITS-sequences-from-a-fungal-genome.git
cd Extract-ITS-sequences-from-a-fungal-genome
chmod +x extractITS.py
```

## Usage

```
extractITS.py -i genome.fasta -o ./output/ --which all --prefix species1
```
