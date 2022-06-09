##########
### Scrip Aim
##########

# This script gathers all data from the hypoMap harmoinzation and curation, to save the most relevant files and a final seurat object (and anndata) into a sepratae folder.

##########
### Load parameters and packages
##########

message(Sys.time(),": Starting mrtree pruning .." )

message(" Load parameters and packages ")

require(tidyverse)
require(Seurat)
require(Matrix)
require(scCoco)

source("R/harmonization_functions.R")

#load parameters from annotation:
parameter_list = jsonlite::read_json("data/parameters_annotation_v2_1.json")
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# define new folder for output
parameter_list$final_output_folder = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_final/" # can also be subfolder of parameter_list$harmonization_folder_path

min_specificity_markers = 0.5

##########
### Load objects
##########

# load seurat
curated_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,".rds"))

# load mrtree clustering
mrtree_result = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_mrtree_clustering_results",".rds"))

# annotation results
annotation_result = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_annotation_result.txt"),data.table = F)
annotation_labelmat =  data.table::fread(file = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_annotation_labelmat.txt"),data.table = F)
annotation_markers_filtered =  data.table::fread(file =paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_annotation_markers_filtered.txt"),data.table = F)


# region prediction result
region_prediction = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_region_prediction.txt"),data.table = F)



##########
### Adjust annotation and add to seurat
##########

# add labels with annotation
label_df = cbind(mrtree_result$labelmat,annotation_labelmat %>% dplyr::select(-Cell_ID)) %>% as.data.frame()

## add labels to seurat:
curated_seurat_object@meta.data = cbind(curated_seurat_object@meta.data , label_df)



##########
### Load and filter marker objects
##########

# read all available marker tables if marker detection was run on multiple subsets !!
markers_comparisons_all_list=list()
markers_comparisons_siblings_list =list()
for(current_start_node in parameter_list$start_nodes_annotation_markers){
  # load markers all
  filename=paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",current_start_node,"_markers_all_",parameter_list$marker_suffix,".tsv")
  if(file.exists(filename)){
    markers_comparisons_all_list[[current_start_node]] = data.table::fread(filename,data.table = F)
  }else{
    message("Cannot find markers stored in : ",filename)
  }
  # siblings
  filename_sib=paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",current_start_node,"_markers_siblings_",parameter_list$marker_suffix,".tsv")
  if(file.exists(filename_sib)){
    markers_comparisons_siblings_list[[current_start_node]] = data.table::fread(filename_sib,data.table = F)
  }else{
    message("Cannot find markers stored in : ",filename_sib)
  }
}
markers_comparisons_all = as.data.frame(do.call(rbind,markers_comparisons_all_list))
markers_comparisons_siblings  = as.data.frame(do.call(rbind,markers_comparisons_siblings_list))
message("All markers for: ",length(unique(markers_comparisons_all$cluster_id))," clusters available")
message("Sibling markers for: ",length(unique(markers_comparisons_siblings$cluster_id))," clusters available")

# load markers for exclusion ad further filter
# load additional remove
additional_remove_genes = unlist(jsonlite::read_json(unlist(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,parameter_list$marker_suffix,"_additionally_removed_markers.json"))))
# gather some other genes that are not informative during annoation:
other_genes_remove = rownames(harmonized_seurat_object@assays$RNA@counts)[grepl("RP|-ps",rownames(harmonized_seurat_object@assays$RNA@counts))]
# make list of all genes that should be removed:
all_exclusion_genes = unique(c(features_exclude_list,additional_remove_genes,other_genes_remove))

# filter marker genes:
markers_comparisons_all = markers_comparisons_all %>% dplyr::filter(! gene %in% all_exclusion_genes &
                                                                      specificity >= min_specificity_markers &
                                                                      p_val_adj < parameter_list$max_pval_adj &
                                                                      cluster_id %in% cluster_column_ids &
                                                                      cluster_id %in% all_children_of_start_node &
                                                                      pct.2 < 0.5
)

# filter marker genes:
markers_comparisons_siblings = markers_comparisons_siblings %>% dplyr::filter(! gene %in% all_exclusion_genes &
                                                                      specificity >= min_specificity_markers
                                                                      p_val_adj < parameter_list$max_pval_adj &
                                                                      cluster_id %in% cluster_column_ids &
                                                                      cluster_id %in% all_children_of_start_node &
                                                                      pct.2 < 0.5
)


# join labels with annotation to marker tables

# save marker table is @misc slot

# store edgelist in @misc


##########
### Adjust predicted regions and summarize
##########

# add region for coloring:
hypothalamus_regions_mapping = data.table::fread("data/region_prediction_mapping.tsv",data.table = FALSE)
result_region = dplyr::left_join(result_region,hypothalamus_regions_mapping %>% dplyr::select(name,Region_curated_Color=name_summarized),by=c("Region_curated"="name"))

# read color scheme --> important !!
#color_value_vector =unlist(jsonlite::read_json("data/region_prediction_mapping_colors.json"))


# overwrite some SCN regions : ...

# overwrite Ghrh C280-121

# overwrite Ghrh C280-21 Prdm8

# overwrite -14 ?


##########
### Clean up metadata
##########


# TODO

##########
### save files
##########

# check that folder exists
system(paste0("mkdir -p ",parameter_list$final_output_folder))

# save other relevant files
# TODO


# save seurat
# TODO

# save anndata




