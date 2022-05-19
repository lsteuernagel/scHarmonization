## This script creates the json with general parameters --> make other jsons or edit this script if other params should be used
# requires some manually decisions which are added here.

param_list = list()

# must be loaded from params:
param_list$harmonization_folder_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_harmonization_test/"#"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_harmonization/"
#param_list$processed_suffix = "_seurat_processed"

# for final merged object:
param_list$merged_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/hypoMap_v2_downsampled_example.rds"#"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/hypoMap_merged_filtered.rds"
param_list$new_name_suffix = "hypoMap_test"#"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/hypoMap_merged_filtered.rds"

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
param_list$n_latent = 85
param_list$n_hidden = 256
param_list$dropout_rate = 0.1
param_list$max_epochs = 50
param_list$early_stopping = FALSE
param_list$dispersion = "gene"
param_list$gene_likelihood = "zinb"
param_list$use_cuda =FALSE

## general harmonization
param_list$k_param = 30
param_list$dist_type="cosine"

## clustering
param_list$target_clusterN = 200
param_list$start_res = 8
param_list$end_res = 20

# save
scUtils::writeList_to_JSON(param_list,filename = "data/parameters_harmonization_v2_1_test.json")
