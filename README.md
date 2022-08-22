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
1. Run leiden clustering over user-specified resolutions (wide range)
2. (Not in slurm) Manual select a subset of cluster levels and export to a table.
3. Construct mrtree on pre-selected clustering levels.
4. Calculate marker genes
5. TODO

## Pipeline

This section includes explanations of the different scripts involved and should allow to reproduce the above steps.

- The pipeline requires a seurat object stored as an rds file that contains all Datasets/Batches that will be integrated and a metadata column specifying for each , a json with features to exclude from HVG search 

- The parameters for the pipeline are stored in json format as well (see below).

- The main scripts are [R/executeIntegration.R] and [R/executeIntegration.R] to prepare the integration objects etc. and run the hyperparameter serach with scVI. See the script for detils on the slurm pipeline and which scripts are used for each of the steps outlined above.


# Details & Parameters 

## Prepare object

The [run_script] takes the input seurat object and calculates highly variable genes etc. Then it saves the data as a h5 anndata object.
Important parameters inlcude:

-

-

## Run scvi

The [run_script] python script executes scvi model setting and training to obtain the scvi integrated low dimensional representation. We are not using the corrected counts but do export them as well. The scvi hyperparameters can be used using default values or an "optimal" set can be searched for using [scIntegration]. 
Important parameters (mostly scvi hyperparameters (see also documentation:)):
-
-

## Integrated Seurat object

Using the [mmm] script, the scvi results are added to the merged seurat and uing Seurat's FindNeighbors and RunUMAP a SNN graph and a UMAP based on the integrated embedding are calculated.
Important parameters:

## Inital clustering:

Using the python script [nnn] the integrated data is clustered into a set of preliminary clusters that can be used for inital data exploration ( when only running the first prt of the pipeline this can also be used as the main clustering). 
Important parameters :
- nClustrers: How many clusters are expected by the user (will increase resolution until this is met.)
- ....

**TBC**

