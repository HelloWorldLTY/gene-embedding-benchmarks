# Benchmarking gene embeddings from sequence, expression, network, and text models for functional prediction tasks
This repository contains code and datasets for benchmarking gene embedding methods across individual gene attributes, paired gene interactions, and gene set relationships.

## About
Gene embeddings have emerged as transformative tools in computational biology, enabling the efficient translation of complex biological datasets into compact vector representations. This study presents a comprehensive benchmark by evaluating 38 classic and state-of-the-art gene embedding methods across a spectrum of functional prediction tasks. These embeddings, derived from data sources such as amino acid sequences, gene expression profiles, protein-protein interaction networks, and biomedical literature, are assessed for their performance in predicting individual gene attributes, paired gene interactions, and gene set relationships. Our analysis reveals that biomedical literature-based embeddings consistently excel in general predictive tasks, amino acid sequence embeddings outperform in functional and genetic interaction predictions, and gene expression embeddings are particularly well-suited for disease-related tasks. Importantly, we find that the type of training data has a greater influence on performance than the specific embedding construction method, with embedding dimensionality having only minimal impact. By elucidating the strengths and limitations of various gene embeddings, this work provides guidance for selecting and successfully leveraging gene embeddings for downstream biological prediction tasks.

## Organization
This repo is organized into several sections, part of which is stored on zenodo.
- `bin`: contains binaries and intermediate files from the benchmarking experiments, which includes the fold and holdout splits that we used in our tests
- `data`: contains datasets and metadata used for benchmarking
  - `embeddings`: preprocessed embeddings for genes from various methods (on zenodo)
  - `gmt`: gene set files used for benchmarking 
  - `matched_pairs`: files used to map one annotation to another
  - `obo`: ontology files for hierarchical biological relationships
  - `paired_gene_interaction_data`: files used for benchmarking paired genetic interactions (on zenodo)
  - `slim_sets`: subsets of annotation terms
  - `embed_meta.csv`: metadata file detailing the embedding methods, their training input type, algorithm, and dimension.
- `results`: contains the results of the gene level and gene pair benchmarking experiments
  - `andes_results`: contains the scores from the gene set benchmarks (on zenodo)
- `src`: contains the code used for preprocessing, summarizing, and benchmarking embeddings across our functional prediction tasks
  - `gene_level_benchmark`: code used for benchmarking disease gene prediction (OMIM) and gene function prediction (GO).
  - `gene_pair_benchmark`: code used for benchmarking genetic interaction (e.g., SL/NG) and transcription factor target (TF) prediction.
  - `gene_set_benchmark`: code used for benchmarking matching pathways (GO/KEGG) and disease/tissue (OMIM/Brenda).
  - `preprocess`: code used for preprocessing embeddings
  - `summary.py`: code used for summarizing the tested embeddings


## Running the scripts
We recommend using conda for installing all necessary packages. Once conda is installed, get started by creating and activating the virtual environment.

 ```bash
 conda env create -f env.yml
 conda activate gene_embed_benchmark 
 ```
 To run the gene set benchmarks, first download the ANDES tool:
```bash 
git clone https://github.com/ylaboratory/ANDES.git
```
Make sure to navigate to the appropriate directory and follow any additional instructions provided in the [ANDES repository](https://github.com/ylaboratory/ANDES) for setting up and running the tool.

## Citation
