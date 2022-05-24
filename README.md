# scHarmonization
Harmonize annotations and other information for an integrated single cell dataset 

# Overview

## Basic harmonization:

1. Prepare a seurat object by finding highly variable genes and exporting and clean object to anndata.
2. Run scvi with manually selected parameters
3. Add scvi results to seurat object and caculate SNN and UMAP using Seurat.
4. Run initial leiden clustering for preliminary analysis (increase resolution until a user specified number of target clusters and use that clustering level)
5. Find Marker genes for preliminary clusters using batch-aware (or basic) wilcoxon test.
6. (Not in slurm) Manually curate object by removing unwanted cells (e.g. additional Doublets clusters or low-quality cell clusters that only formed after integration ) and re-run SNN and UMAP.

## Advanced Harmonization
1. Run leiden clustering over use-specified resolutions (wide range)
2. (Not in slurm) Manual select a subset of cluster levels and export to a table.
3. Construct mrtree on pre-selected clustering levels.
4. Calculate marker genes
5. TODO

## Parameter overview

TODO: explain major params

TODO: add json example

### Old Step outline:

- Prepare input dataset (export to anndata etc.)
- Run scVI to integrate the merged dataset using use-specfied parameters (use scIntegration to find a set of optimal parameters)
- Create harmonized seurat
  - Read results from scVI and create a new object
  - Calculate a SNN (Seurat) on integrated result
  - Run UMAP (and save model in object)
- Run clustering on SNN of integrated result
  - general option: run one round / repeat until a certain numer of clusters is reached
  - [HypoMap]: Run multiple runs and combine to hierarchical tree using [mrtree](https://github.com/pengminshi/MRtree)
- Run marker detection using batch-aware approach (optionally exclude features)
  - [HypoMap]: Run marker detection to sibling clusters in hierarchical tree using batch-aware approach (optionally exclude features) 
  - [HypoMap]: Prune clustering result based on markers to siblings (if no differences --> merge siblings)
  - [HypoMap]: Re-Run marker detection on pruned clusters using batch-aware approach
- Run [scMRMR](https://github.sf.mpg.de/lsteuernagel/scMRMR) to find a small set of descriptive markers for each cluster
  - [HypoMap]: Run [scMRMR](https://github.sf.mpg.de/lsteuernagel/scMRMR) using a combination of global and siblings markers to obtain most informative genes
  - [HypoMap]: Annotate clusters in hierarchical tree using most informative gene(s)
- [HypoMap]: Run [scCoco](https://github.sf.mpg.de/lsteuernagel/scCoco) to infer likely locations of each cluster
  - [HypoMap]: Combine with Dataset origin to prioritize voxels
  - [HypoMap]:Break region annotation down to one most likely subregion/nucleus
- Export final object
  - [Manual]: Add a manual script for finalizing the object (e.g. pruning metadata etc.)
