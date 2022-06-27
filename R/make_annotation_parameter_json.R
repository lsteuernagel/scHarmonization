## This script creates the json with general parameters --> make other jsons or edit this script if other params should be used
# requires some manually decisions which are added here.

param_list = list()

# must be loaded from params:
param_list$harmonization_folder_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2c_harmonization/"#"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_harmonization/"

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
                                               "C7-1" = "GLU",
                                               "C7-2"= "GABA",
                                               "C7-3"="Astro-Ependymal",
                                               "C7-4" = "Oligo+Precursor",
                                               "C7-5" = "Immune",
                                               "C7-6" = "ParsTuber",
                                               "C7-7" = "Vascular",
                                               # glu
                                               "C25-1"= "GLU-1",
                                               "C25-2"= "GLU-2",
                                               "C25-3"= "GLU-3",
                                               "C25-4"= "GLU-4",
                                               "C25-5"= "GLU-5",
                                               "C25-6"= "GLU-6", #oxt
                                               "C66-22" = "GLU-6.Caprin2",
                                               "C25-7"= "GLU-7", # Pmch
                                               "C66-24" = "GLU-7.Pmch",
                                               "C25-8"= "GLU-8", # Foxb1
                                               "C66-25" = "GLU-8.Foxb1",
                                               "C66-26" = "GLU-8.Prkch",
                                               "C25-9"= "GLU-9",
                                               "C66-27" = "GLU-9.Gnrh1",
                                               # gaba
                                               "C25-10"="GABA-2",
                                               "C25-11"="GABA-1",
                                               "C25-12"="GABA-3",
                                               "C25-13"="GABA-4",
                                               "C25-14"="GABA-5",
                                               "C66-48"="GABA-5.Meis2",
                                               "C25-15"="GABA-6", # satb2 th
                                               "C66-49"="GABA-6.Satb2",
                                               "C25-16"="GABA-7", # Chat
                                               "C66-50" = "GABA-7.Chat",
                                               # non-neuron
                                               "C25-17"= "Ependymal-like",
                                               "C25-18"= "Astrocytes",
                                               "C25-19"= "Oligodendrocytes",
                                               "C25-20"= "OPC",
                                               "C25-21"= "Immune",
                                               "C25-22"= "Immune-Trbc1",
                                               "C25-23"= "ParsTuber",
                                               "C25-24" = "Mural+Endothelial",
                                               "C25-25" = "Fibroblasts",

                                               # difficult neurons:
                                               "C185-3" = "NSCs",
                                               "C66-16" = "GLU-4.Mixed",
                                               "C66-28" = "GABA-2.Mixed",
                                               "C66-23"= "GLU-6.Mixed",
                                               #"C66-7"= "GABA-2.Mixed", # fine
                                              # "C66-14"= "GLU-3.Mixed",
                                              # "C66-31"= "GABA-2.Mixed",
                                               #"C185-53"= "Unassigned",
                                              # "C185-3"= "NSCs",

                                               # specific noin-neurons
                                               "C66-51" = "Tanycytes",
                                               "C66-52"= "Ependymal",
                                               "C66-63"= "Mural",
                                               "C66-64"= "Endothelial",
                                               # on lowe rlvele !
                                               "C185-160" = "OPC.Dividing",
                                               "C185-137"= "Erythroid-like",
                                               "C185-138"= "Hypendymal"
                                               #"C66-59"= "Mural",
                                           ))

# file ABA data
param_list$aba_gene_to_id_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/allen_brain/ish_data_voxel/aba_gene_to_id.tsv" # can also provide NULL then scCoco will load independetly
param_list$aba_ish_matrix_file  = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/allen_brain/ish_data_voxel/aba_ish_matrix_energy.tsv"
param_list$aba_ccf_grid_annotation_file  = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/allen_brain/ish_data_voxel/aba_ccf_grid_annotation.rds"
param_list$mba_ontology_flatten_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/allen_brain/ish_data_voxel/mba_ontology_flatten.tsv"

# region prediction ABA:
param_list$cluster_column = "C286"
param_list$start_node_tree = "C2-1"
# filtering of marker genes
param_list$min_specificity = 4
param_list$max_pval_adj = 0.001
param_list$max_pct_2 = 0.3
# for prediction with ABA iSH
param_list$min_ids = 3 # how many valid ids are required to calculate results !?
param_list$topn_results = 4 # report topn results back
param_list$target_structure_id = "1097" # hypothalamaus
param_list$max_ids_to_include = 10000 # how many genes in ABA enrichment (10k or Inf for all)
param_list$target_level = "8"
# region prediction with manual suggested regions per data:
param_list$suggested_region_file = "data/hypoMap_suggested_region_per_dataset.json"
param_list$min_dataset_cell_value = 2000
param_list$dataset_column = "Dataset"
param_list$max_region_weight_value = 0.75 # all region weights (percentages) above this are set to 1
# for final filtering (prediction not meeting the criteria are set to NA):
param_list$min_score_region_summary = 0.75
param_list$min_pct_of_genes = 0.3333
param_list$min_number_of_ids = 3


# save
scUtils::writeList_to_JSON(param_list,filename = "data/parameters_annotation_v2_2.json")







