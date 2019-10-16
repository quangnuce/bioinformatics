# Pipeline for gene annotation, clustering and extract data for ML

## Prerequisites
- Python 3.7
- Prokka
- RGI
- Roary
- KMC

## Usage
Make a folder contains all DNA sequences in fasta format (.fasta, .fna)
### Run gene annotation and clustering 
```bash
python pipeline.py <fasta folder> <output folder>
```
### Export data for ML

```bash
python combine2MLinput.py <roary output folder> <rgi output folder> <kmer-index file> <ML data output folder>
```
