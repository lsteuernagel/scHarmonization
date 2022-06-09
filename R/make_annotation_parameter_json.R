## This script creates the json with general parameters --> make other jsons or edit this script if other params should be used
# requires some manually decisions which are added here.

param_list = list()

# must be loaded from params:
param_list$harmonization_folder_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_harmonization/"#"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_harmonization/"

# for final curated object:
param_list$new_name_suffix = "hypoMap_harmonized_curated"#
param_list$marker_suffix = "pruned"

# signature for exclusion
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

# annotation general
param_list$start_nodes_annotation_markers = c("C2-1","C2-2") # can also just be = start_nodes_pruning_markers  # use this when there are multiple marker tables after splitting the marker detection
param_list$min_specificity = 1.5 # ?
param_list$min_specificity_sibling_children = 2.5
param_list$limit_factor = 5
param_list$max_score_siblings_children = 20
param_list$reverse_order = TRUE

# new cluster names:
# convert as.list to store names in json
param_list$manual_names_annotation = as.list(c("C2-1" = "Neurons",
                                           "C2-2"="Non-Neurons",
                                           "C7-1" = "vGLUT",
                                           "C7-2"= "GABA",
                                           "C7-3"="Astro-Ependymal",
                                           "C7-4" = "Oligo+Precursor",
                                           "C7-5" = "Immune",
                                           "C7-6" = "ParsTuber",
                                           "C7-7" = "Vascular",
                                           "C23-1"= "vGLUT-1",
                                           "C23-2"= "vGLUT-2",
                                           "C23-3"= "vGLUT-3",
                                           "C23-4"= "vGLUT-4",
                                           "C23-5"= "vGLUT-5",
                                           "C23-6"= "vGLUT-6-Pmch",
                                           "C23-7"= "vGLUT-7-Foxb1",
                                           "C23-8"= "vGLUT-8",
                                           "C23-9"= "vGLUT-9",
                                           "C23-10"="GABA-3",
                                           "C23-11"="GABA-2",
                                           "C23-12"="GABA-1",
                                           "C23-13"="GABA-4",
                                           "C23-14"="GABA-5",
                                           "C23-15"="GABA-6-Chat",
                                           "C23-16"="Ependymal-like",
                                           "C23-17"= "Astrocytes",
                                           "C23-18"= "Oligodendrocytes",
                                           "C23-19"= "OPC",
                                           "C23-20"= "Immune",
                                           "C23-21"= "ParsTuber",
                                           "C23-22" = "Mural+Fibro",
                                           "C23-23" = "Endothelial",
                                           "C62-36" = "GABA-1.Mixed",
                                           "C62-46" = "Tanycytes",
                                           "C62-47"= "Ependymal",
                                           "C62-48"= "Hypendymal",
                                           "C62-49"= "Erythroid-like",
                                           "C62-55" = "Dividing.OPC",
                                           "C62-60"= "Fibroblasts",
                                           "C62-59"= "Mural",
                                           "C180-53"= "Unassigned",
                                           "C180-3"= "Neural stem cells"))

# file ABA data
param_list$aba_gene_to_id_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/allen_brain/ish_data_voxel/aba_gene_to_id.tsv" # can also provide NULL then scCoco will load independetly
param_list$aba_ish_matrix_file  = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/allen_brain/ish_data_voxel/aba_ish_matrix_energy.tsv"
param_list$aba_ccf_grid_annotation_file  = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/allen_brain/ish_data_voxel/aba_ccf_grid_annotation.rds"
param_list$mba_ontology_flatten_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/allen_brain/ish_data_voxel/mba_ontology_flatten.tsv"

# region prediction ABA:
param_list$cluster_column = "C180"
param_list$start_node_tree = "C2-1"
# filtering of marker genes
param_list$min_specificity = 4
param_list$max_pval_adj = 0.001
param_list$max_pct_2 = 0.3
param_list$min_ids = 3 # how many valid ids are required to calculate results !?
param_list$topn_results = 4 # report topn results back
param_list$target_structure_id = "1097" # hypothalamaus
param_list$max_ids_to_include = 10000 # how many genes in ABA enrichment (10k or Inf for all)
param_list$min_score_region_summary = 0.75
param_list$target_level = "8"
param_list$max_region_weight_value = 0.75 # all region weights (percentages) above this are set to 1

# region prediction with manual suggeseted regions per data:
param_list$suggested_region_file = "data/hypoMap_suggested_region_per_dataset.json"
param_list$min_dataset_cell_value = 2000
param_list$dataset_column = "Dataset"

# save
scUtils::writeList_to_JSON(param_list,filename = "data/parameters_annotation_v2_1.json")







