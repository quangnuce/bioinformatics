# Pipeline for gene annotation, clustering and extract data for ML

## Prerequisites
- Python 3.6
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
- CARD data: path to card.json (RGI). Get from here: https://card.mcmaster.ca/latest/data
```bash
python pipeline.py Example/fna Example/OUT card.json
```
Output:
- Prokka output folder: protein annotation (.faa, .ffn, .gff)
- RGI output folder: AMR annotation (.json,.txt)
- ROARY output: gene cluster (pangenome)
### Export pangenome data (gene cluster) for ML 
Parameter:
- Prokka output folder
- Roary output folder: make sure it contains gene_presence_absence.csv
- RGI output folder
- ML output folder: hold all exported files can be use in ML
```bash
python combine2MLInput_GeneCluster.py Example/OUT/PROKKA_OUT Example/OUT/ROARY_OUT Example/OUT/RGI_OUTPUT Example/OUT/ML_INPUT
```
### Export kmer data  for ML 
Parameter:
- Fasta folder: contains all fasta files
- ML output folder
- k (10-20)
- Prokka output folder
- RGI output folder
```bash
python combile2MLInput_kmer.py Example/fna Example/OUT/ML_INPUT 10 Example/OUT/PROKKA_OUT Example/OUT/RGI_OUTPUT
```
Output:
- KMC output folder: .kmers files
- whole_genome_kmer folder: kmers for whole genome
- amr_genome_kmer folder: kmers for DNA sequences extracted from  AMR genes
- non_amr_gennome_kmer folder: kmers for DNA sequences extracted from non-AMR genes
