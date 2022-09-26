# scHarmonization
Harmonize annotations and other information for an integrated single cell dataset 

# Overview

The basic harmonization can be used as a straightforward approach to creating integrated single-cell datasets using a mix of scvi, seurat and scanpy solutions in a slurm workflow for HPCs.
The core scvi model is relatively robust with default parameters. For specific projects like [HypoMap](https://www.repository.cam.ac.uk/handle/1810/340518) the [scIntegration](https://github.com/lsteuernagel/scIntegration) pipeline can be used to evalute hyperparameters. 

The second part of the pipeline is used to add additional detailed annotations (optimized clustering, naming, brain region prediction) and is mostly optimized towards the HypoMap project but can be adapted for other projects as well.

## Docker image

All required packages and software is available via this [docker image] (LINK), which can also be pulled as a singularity image (which is the current default in the slurm bash scripts of this pipeline).

## Basic harmonization:
1. Prepare a seurat object by finding highly variable genes and exporting and clean object to anndata.
2. Run scvi with manually selected parameters
3. Add scvi results to seurat object and caculate SNN and UMAP using Seurat.
4. Run initial leiden clustering for preliminary analysis (increase resolution until a user specified number of target clusters and use that clustering level)
5. Find Marker genes for preliminary clusters using batch-aware (or basic) wilcoxon test.
6. (Not in slurm) Manually curate object by removing unwanted cells (e.g. additional Doublets clusters or low-quality cell clusters that only formed after integration ) and re-run SNN and UMAP.

## Advanced Harmonization

This part is optimized towards the HypoMap and mght require large chnages to properly work with other data.

1. Run leiden clustering over user-specified resolutions (wide range)
2. (Not in slurm) Manually select a subset of cluster levels and export to a table.
3. Construct mrtree on pre-selected clustering levels.
4. Calculate marker genes
5. Prune mrtree based on marker genes
6. Re-calculate marker genes after pruning
7. Annotate mrtree nodes using specific marker genes (and optionally manual labels)
8. Predict spatial locations
9. Finalize and export curated object as seurat and anndata.

## Pipeline

This section includes explanations of the different scripts involved and should allow to reproduce the above steps.

- The pipeline requires a seurat object stored as an rds file that contains all Datasets/Batches that will be integrated and a metadata column specifying for each , a json with features to exclude from HVG search 

- The parameters for the pipeline are stored in json format as well. See this file for the HypoMap: [data/parameters_harmonization_v2_8.json(data/parameters_harmonization_v2_8.json).

- The main scripts are [R/execute_scHarmonization_basic.R](R/execute_scHarmonization_basic.R) and [R/execute_scHarmonization_advanced.R](R/execute_scHarmonization_advanced.R) to prepare the integration objects etc. and run the hyperparameter serach with scVI. See the script for detils on the slurm pipeline and which scripts are used for each of the steps outlined above.

# Details: Basic harmonization: 

## Prepare object

The [R/run_scripts/prepare_harmonization.R](R/run_scripts/prepare_harmonization.R) takes the input seurat object and calculates highly variable genes etc. Then it saves the data as an h5 anndata object.
Important parameters include:

- **feature_set_size** : A vector of different sizes for HVGs

- **batch_var** : the batch ID variable in the metadata

## Run scvi

The [python/integrate_scVI_v015.py](python/integrate_scVI_v015.py) script executes scvi model setting and training to obtain the scvi integrated low dimensional representation. We are not using the corrected counts but do export them as well. The scvi hyperparameters can be used using default values or an "optimal" set can be searched for using [scIntegration](LINK). 
Important parameters are mostly scvi hyperparameters (see also the [documentation](https://docs.scvi-tools.org/en/stable/api/reference/scvi.model.SCVI.html))

## Integrated Seurat object

Using the [R/run_scripts/basic_harmonization.R](R/run_scripts/basic_harmonization.R) script, the scvi results are added to the merged seurat and uing Seurat's FindNeighbors and RunUMAP a SNN graph and a UMAP based on the integrated embedding are calculated.

## Inital clustering:

Using the python script [python/basic_leiden_clustering.py](python/basic_leiden_clustering.py) the integrated data is clustered into a set of preliminary clusters that can be used for inital data exploration ( when only running the first prt of the pipeline this can also be used as the main clustering). 
Important parameters :
- **nClustrers**: How many clusters are expected by the user (will increase resolution until this is met.)

## Initial marker detection:

Using the [R/run_scripts/basic_marker_detection.R](R/run_scripts/basic_marker_detection.R) script, markers for each of the abve clusters are detected using Seurat's FindMarkers function or an adapted version using a [stratified wilcoxon test](https://github.com/KChen-lab/stratified-tests-for-seurat).

## Manual curation based on inital harmonization

This is done by the user. See [curate_hypoMap_2.R ](curate_hypoMap_2.R ) for the HypoMap version.

# Details: Advanced harmonization: 

##  Leiden clustering 

Using the [python/basic_leiden_clustering.py](python/basic_leiden_clustering.py) script leiden clusterings over a wide range of user-specified resolutions are obtained.

## Select final cluster levels

From the leiden clusters computed in the previous step, select a subset as input for mrtree. This is not part of the slurm pipeline because it is different for each dataset and project. For HypoMap, please see: [R/curation_scripts/cluster_selection_hypoMap_2.R](R/curation_scripts/cluster_selection_hypoMap_2.R).

## Construct hierachical cluster tree

Using mrtree and the clusters selected in the previous step a hierachical cluster treeis built using the script [R/run_scripts/mrtree_construction.R](R/run_scripts/mrtree_construction.R).

## Cluster marker detection on tree

Marker genes for each cluster on each level are detected using Seurat's FindMarkers function or an adapted version using a [stratified wilcoxon test](https://github.com/KChen-lab/stratified-tests-for-seurat). The script [R/run_scripts/mrtree_marker_detection.R](R/run_scripts/mrtree_marker_detection.R) traverses the tree and finds markers vs both all clusters and the sibling clusters in the tree specifically. Most parameters for this step specifiy detection thresholds etc.

## Prune hierachical cluster tree

Using the script [R/run_scripts/mrtree_pruning.R](R/run_scripts/mrtree_pruning.R) and the sibling marker genes calculated before, the tree is pruned. I.e. if there are less than X sibling markers between sibling clusters they are merged to one node in the tree to avoid over fragmentation on the lowest tree level.
After this the cluster markers are re-run similar to the previous step on the pruned tree.

## Cluster annotation

Annotate mrtree nodes using specific marker genes (and optionally manual labels provided in the parameter file) using the script [R/run_scripts/
mrtree_annotation.R](R/run_scripts/mrtree_annotation.R).

This script pulls its parameters and manual annotations from a separate json file: [data/parameters_annotation_v2_2.json](data/parameters_annotation_v2_2.json) which is specific to HypoMap.

## Spatial predictions

Using the script [R/run_scripts/region_prediction.R](R/run_scripts/region_prediction.R) and the [scCoco package](LINK) which provides additional functions around the [cocoframer R package](https://github.com/AllenInstitute/cocoframer) that allows to query the [Allen brain atlas API](atlas.brain-map.org) for each cluster potential brain regions of origin are predicted. This step (and scripts) is optimized for HypoMap and probably does not apply seamlessly to other projects.

## Finalize and export 

The final seurat object is manually curated (outside of the pipeline) by adding all information created above and then exported as .rds and .h5ad.
For HypoMap this script was used: [R/curation_scripts/finalize_hypoMap_2.R](R/curation_scripts/finalize_hypoMap_2.R) with parameters from: [data/parameters_annotation_v2_2.json](data/parameters_annotation_v2_2.json).



