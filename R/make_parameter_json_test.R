## This script creates the json with general parameters --> make other jsons or edit this script if other params should be used
# requires some manually decisions which are added here.

param_list = list()

# must be loaded from params:
param_list$harmonization_folder_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_harmonization_test/"#"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_harmonization/"
#param_list$processed_suffix = "_seurat_processed"

# for final merged object:
param_list$merged_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/hypoMap_merged_filtered.rds"#"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/hypoMap_v2_downsampled_example.rds"#
param_list$new_name_suffix = "hypoMap_test"#"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/hypoMap_merged_filtered.rds"
param_list$additional_clustering_suffix = ""

# signature for evaluation
param_list$genes_to_exclude_file = "data/features_exclude_list2.json"

# general
param_list$n_cores = 50
param_list$id_column = "Cell_ID"
param_list$global_seed = 123456
param_list$sample_column = "Sample_ID"
param_list$batch_var = "Batch_ID"
param_list$feature_set_size = 3000
param_list$feature_set_file = paste0(param_list$harmonization_folder_path,"feature_set.json")
param_list$assay_name = "RNA"
param_list$integration_name = "scvi"

# scvi integration:
param_list$categorical_covariates = c("Dataset",param_list$batch_var)
param_list$continuous_covariates =character(0)
param_list$n_layers = 3
param_list$n_latent = 65
param_list$n_hidden = 256
param_list$dropout_rate = 0.1
param_list$max_epochs = 400
param_list$early_stopping = FALSE
param_list$dispersion = "gene"
param_list$gene_likelihood = "zinb"
param_list$use_cuda =FALSE

## general harmonization
param_list$k_param = 30
param_list$dist_type="cosine"

## initial clustering
param_list$target_clusterN_initial = 200
param_list$start_res_initial = 10
param_list$end_res_initial = 20
param_list$step_size_initial = 1
param_list$include_low_res_initial = FALSE

## full clustering
param_list$target_clusterN = 300
param_list$start_res = 1
param_list$end_res = 40
param_list$step_size = 1
param_list$include_low_res = TRUE
param_list$min_cells_valid = 3

# mrtree
param_list$clusters_for_mrtree_file = "mrtree_input_labels.txt"
param_list$use_recon_labelmat = TRUE # avoids skipping the lowest level of cluster when building the matrix !
param_list$specificity_base = 0.001
param_list$n_cores_markers = 4

# pruning:
param_list$min_cells = 5 # if beelow --> merge with neighbor
param_list$min_specificity = 0.5 # min specificity for a sibling marker to count
param_list$max_pvalue_prune = 0.01 # max pvalue for a sibling marker to count
param_list$min_sibling_markers = 3 # how many sibling markers are required to not merge
param_list$min_prune_level = 4 # highest level is 2 (1 does not exist because level is based on destination node)
param_list$start_nodes_pruning_markers = c("K2-0","K2-1") # use this when there are multiple marker tables after splitting the marker detection
param_list$old_prefix = "K"
param_list$new_prefix = "C"

# basic marker detection
param_list$basic_marker_filename = "_inital_markers"
param_list$assay_markers ="RNA"
param_list$assay_slot = "data"
param_list$test.use = "wilcox-stratified" # either "wilcox-stratified" or a basic Seurat FindMarkers test
param_list$logfc.threshold = 0.3
param_list$min.pct = 0.1
param_list$min.diff.pct = 0.05
param_list$max.cells.per.ident = 20000
param_list$min.cells.feature = 10
param_list$min.cells.group =  10
param_list$base = 2
param_list$only.pos = TRUE
param_list$batch_var = param_list$batch_var

# annotation
param_list$start_nodes_annotation_markers = c("C2-1","C2-2") # can also just be = start_nodes_pruning_markers  # use this when there are multiple marker tables after splitting the marker detection

# save
scUtils::writeList_to_JSON(param_list,filename = "data/parameters_harmonization_v2_1_test.json")
