# Benchmarking gene embeddings from sequence, expression, network, and text models for functional prediction tasks
This repository contains code and datasets for benchmarking 38 gene embedding methods across individual gene attributes, paired gene interactions, and gene set relationships.

## Organization
This repo is organized into several sections, part of which is stored on zenodo.
- `data`: contains datasets and metadata used for benchmarking
  - `embeddings`: preprocessed embeddings for genes from various methods (on zenodo)
  - `gmt`: gene set files used for benchmarking 
  - `matched_pairs`: files used to map one annotation to another
  - `obo`: ontology files for hierarchical biological relationships
  - `slim_sets`: subsets of annotation terms
  - `embed_meta.csv`: metadata file detailing the embedding methods, their training input type, algorithm, and dimension.
- `src`: contains the code used for preprocessing, summarizing, and benchmarking embeddings across our functional prediction tasks
  - `gene_level_benchmark`: benchmarking disease gene prediction (OMIM) and gene function prediction (GO).
  - `gene_pair_benchmark`: benchmarking genetic interaction (e.g., SL/NG) and transcription factor target (TF) prediction.
  - `gene_set_benchmark`: benchmarking matching pathways (GO/KEGG) and disease/tissue (OMIM/Brenda).
  - `preprocess`: preprocessing embeddings
  - `summary.ipynb`: summarizes tested embeddings



## Running the scripts
We recommend using conda for installing all necessary packages. Once conda is installed 
get started by creating and activating the virtual environment.

 ```bash
 conda env create -f env.yml
 conda activate gene_embed_benchmark 
 ```
 To run the gene set benchmarks, first dowload the ANDES tool:
```bash 
git clone https://github.com/ylaboratory/ANDES.git
```
Make sure to navigate to the appropriate directory and follow any additional instructions provided in the [ANDES repository](https://github.com/ylaboratory/ANDES) for setting up and running the tool.

## Citation
