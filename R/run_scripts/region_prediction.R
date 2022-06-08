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
parameter_list = jsonlite::read_json("data/parameters_harmonization_v2_4.json")
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})
parameter_list$marker_suffix = "pruned"
parameter_list$new_name_suffix=paste0(parameter_list$new_name_suffix,"_curated")

# read features to excludes
features_exclude_list= unlist(jsonlite::read_json(parameter_list$genes_to_exclude_file))
#features_exclude_list = lapply(features_exclude_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# load seurat
harmonized_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,".rds"))


##########
### Load marker genes
##########




##########
### Clean up marker genes
##########

# all_markers_filtered


# make list per cluster
markers_per_cluster = split(all_markers_filtered$gene,f = all_markers_filtered$cluster_1)
markers_weight_per_cluster = split(all_markers_filtered$specificity,f = all_markers_filtered$cluster_1)
markers_weight_per_cluster = sapply(markers_weight_per_cluster,function(x,maxv=100){x[x>maxv] = maxv; x = x/sum(x)*length(x); return(x)})

##########
### Load aba ish data
##########

aba_gene_to_id = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/allen_brain/ish_data_voxel/aba_gene_to_id.tsv",data.table = FALSE)
aba_ish_matrix = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/allen_brain/ish_data_voxel/aba_ish_matrix_energy.tsv",data.table = FALSE,header = TRUE)
aba_ccf_grid_annotation = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/allen_brain/ish_data_voxel/aba_ccf_grid_annotation.rds")
mba_ontology_flatten= data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/allen_brain/ish_data_voxel/mba_ontology_flatten.tsv",data.table = FALSE)

##########
### Run region enrichment with scCoco
##########

library(scCoco)

# run scCoco regions per geneset
hypoMap_region_annotation_full = findRegions_genesets(gene_set = markers_per_cluster,
                                                      min_ids = 5,
                                                      topn=5,
                                                      target_structure_id = "1097",
                                                      max_ids_to_include = Inf,
                                                      # gene_set_weights = markers_weight_per_cluster,
                                                      aba_gene_to_id = aba_gene_to_id,
                                                      aba_ish_matrix = aba_ish_matrix,
                                                      aba_ccf_grid_annotation = aba_ccf_grid_annotation ,
                                                      mba_ontology_flatten= mba_ontology_flatten
)

# run scCoco summarised regions  per geneset
hypoMap_region_annotation_df = summariseRegions_genesets(findRegion_result = hypoMap_region_annotation_full,min_score=0.75)
#table(hypoMap_region_annotation_df$Region)

##########
### Load manual per dataset curation
##########



hypoMap_region_annotation_df$Region_curated = hypoMap_region_annotation_df$Region
hypoMap_region_annotation_df$Region_curated[hypoMap_region_annotation_df$Region_curated %in% c("Lateral mammillary nucleus","Medial mammillary nucleus")] = "Mammillary nucleus"
hypoMap_region_annotation_df$Region_curated[hypoMap_region_annotation_df$Region_curated %in% c("Median eminence")] = "Arcuate hypothalamic nucleus"
hypoMap_region_annotation_df$Region_curated[hypoMap_region_annotation_df$Region_curated %in% c("Periventricular hypothalamic nucleus, intermediate part","Anteroventral periventricular nucleus","Periventricular hypothalamic nucleus, preoptic part","Periventricular hypothalamic nucleus, posterior part")] = "Periventricular hypothalamic nucleus"
hypoMap_region_annotation_df$Region_curated[hypoMap_region_annotation_df$Region_confidence == 0] = NA













