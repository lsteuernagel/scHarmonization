# scHarmonization
Harmonize annotations and other information for an integrated single cell dataset 


# Overview

Steps:

- Run scVI to integrate the merged dataset using use-spcfied parameters (use scIntegration to find a set of optimal parameters)
- Create a SNN (Seurat) on integrated result
  - Run UMAP (and save model)
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
