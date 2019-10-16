# Pipeline for gene annotation, clustering and extract data for ML

## Prerequisites
- Python 3.7
- Prokka
- RGI
- Roary
- KMC
## Pipeline
![Diagram](https://github.com/quangnuce/bioinformatics/blob/master/pipeline.png)
## Usage
Make a folder contains all DNA sequences in fasta format (.fasta, .fna)
### Run gene annotation and clustering 
Parameter:
- Fasta folder: contains all fasta files
- Output folder: stores all file exported from running pipeline
```bash
python pipeline.py Example Example/Out
```
### Export data for ML
Parameter:
- Roary output folder: make sure it contains gene_presence_absence.csv
- RGI output folder
- K-mer dump folder: output after run kmc and kmc_dump
- K-mer index file: list of all k-mers with their indexes. Kmer dump files and k-mer index files must have same k.
- ML output folder: hold all exported files can be use in ML
```bash
python combine2MLinput.py Example/Out/ROARY_OUT Example/Out/RGI_OUT Example/Out/KMC_OUT Kleb.ArrInds Example/Out/ML_OUT
```
