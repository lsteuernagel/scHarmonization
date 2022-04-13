## This script creates the json with general parameters --> make other jsons or edit this script if other params should be used
# requires some manually decisions which are added here.

param_list = list()

# must be loaded from params:
param_list$harmonization_folder_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_harmonization/"#"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_harmonization/"
#param_list$processed_suffix = "_seurat_processed"

# for final merged object:
param_list$merged_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/hypoMap_merged_filtered.rds"#"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/hypoMap_merged_filtered.rds"
param_list$new_name_suffix = "hypoMap_harmonized"#"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/hypoMap_merged_filtered.rds"

# signature for evaluation
param_list$genes_to_exclude_file = "data/features_exclude_list.json"

# general
param_list$n_cores = 50
param_list$id_column = "Cell_ID"
param_list$global_seed = 123456
param_list$sample_column = "Sample_ID"
param_list$batch_var = "Batch_ID"
param_list$feature_set_sizes = c(750,1000,1250,1500,2000,2500,3000)
param_list$assay_name = "RNA"

# scvi integration:
param_list$scvi_models_to_test = 30 # depends on the number of the models in scvi arg list
param_list$latent_space_sizes = c(50,65,80,95,110,140)
# TODO

# evaluation:
param_list$ntrees_mixing = 20000
param_list$sampsize_pct = 0.3333
param_list$max_for_norm = 0.01
param_list$k_param = 20
param_list$dist_type = "cosine"
## asw
param_list$subset_cells = TRUE
param_list$target_clusterN = 200
param_list$start_res_asw = 5
param_list$end_res_asw = 15

# save
scUtils::writeList_to_JSON(param_list,filename = "data/parameters_harmonization_v2_3.json")
