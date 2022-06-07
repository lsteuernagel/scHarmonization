##########
### Load parameters and packages
##########

message(Sys.time(),": Starting mrtree pruning .." )

message(" Load parameters and packages ")

require(tidyverse)
require(Seurat)
require(Matrix)

source("R/harmonization_functions.R")
source("R/stratified_wilcoxon_functions.R")

# get params-filename from commandline
command_args<-commandArgs(TRUE)
param_file = command_args[1]
# read all parameters and filepaths
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

#test:
parameter_list = jsonlite::read_json("data/parameters_harmonization_v2_1_test.json")
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})
parameter_list$marker_suffix = "pruned"
parameter_list$new_name_suffix=paste0(parameter_list$new_name_suffix,"_curated")

# read features to excludes
features_exclude_list= unlist(jsonlite::read_json(parameter_list$genes_to_exclude_file))
#features_exclude_list = lapply(features_exclude_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# load seurat
harmonized_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,".rds"))

##########
### ...
##########

neighbors = as.Neighbor(harmonized_seurat_object@graphs$SNN_scvi)
neighbors = as.Neighbor(harmonized_seurat_object@graphs$NN_scvi)

current_dataset = "CampbellDropseq"
label_vec = rep(NA,length(harmonized_seurat_object@meta.data$Author_CellType))
label_vec[harmonized_seurat_object@meta.data$Dataset == current_dataset] = harmonized_seurat_object@meta.data$Author_CellType[harmonized_seurat_object@meta.data$Dataset == current_dataset]
table(label_vec)

prop_labels = mapscvi::propagate_labels_prob(neighbors_object = neighbors,
                                             label_vec = label_vec,
                                             with_euclidean = FALSE,
                                             add_entropy = TRUE,
                                             add_to_seurat = FALSE)
table(prop_labels$predicted)


##########
### ...
##########

neighbors = as.Neighbor(curated_seurat_object@graphs$SNN_scvi)
neighbors = as.Neighbor(curated_seurat_object@graphs$NN_scvi)

current_dataset = "CampbellDropseq"
label_vec = rep(NA,length(curated_seurat_object@meta.data$Author_CellType))
label_vec[curated_seurat_object@meta.data$Dataset == current_dataset] = curated_seurat_object@meta.data$Author_CellType[curated_seurat_object@meta.data$Dataset == current_dataset]
table(label_vec)

prop_labels = mapscvi::propagate_labels_prob(neighbors_object = neighbors,
                                             label_vec = label_vec,
                                             with_euclidean = FALSE,
                                             add_entropy = TRUE,
                                             add_to_seurat = FALSE)
table(prop_labels$predicted)

prop_labels_filtered = prop_labels[prop_labels$prediction_probability > 0.1 & prop_labels$prediction_entropy < 0.6 ,] %>% na.omit()

## maybe better to do with random forest ?

## or scANVI  with initialized ScVI !?

# https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scarches_scvi_tools.html



