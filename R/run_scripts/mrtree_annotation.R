##########
### Load parameters and packages
##########

message(Sys.time(),": Starting mrtree annotation .." )

message(" Load parameters and packages ")

require(tidyverse)
require(Seurat)
require(Matrix)

source("R/harmonization_functions.R")

# get params-filename from commandline
command_args<-commandArgs(TRUE)
param_file = command_args[1]
# read all parameters and filepaths
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# test on full
# parameter_list = jsonlite::read_json("data/parameters_annotation_v2_1.json")
# parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# read features to excludes
features_exclude_list= unlist(jsonlite::read_json(parameter_list$genes_to_exclude_file))
#features_exclude_list = lapply(features_exclude_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# load seurat
harmonized_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,".rds"))

# load mrtree clustering
mrtree_result = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_mrtree_clustering_results",".rds"))

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
other_genes_remove = rownames(harmonized_seurat_object@assays$RNA@counts)[grepl("RP|Gm|Rik|-ps",rownames(harmonized_seurat_object@assays$RNA@counts))]
# make list of all genes that should be removed:
all_exclusion_genes = unique(c(features_exclude_list,additional_remove_genes,other_genes_remove))

##########
### Load parameters and packages
##########

# get key objects
edgelist = mrtree_result$edgelist
labelmat = mrtree_result$labelmat

edgelist = edgelist[,c("from","to","level")]
cluster_levels = as.data.frame(labelmat) %>% tidyr::pivot_longer(everything(),names_to = "clusterlevel", values_to = "cluster")  %>% dplyr::group_by(cluster) %>%
  dplyr::add_count(name="ncells") %>% dplyr::distinct(clusterlevel,cluster,ncells)
edgelist = dplyr::left_join(edgelist,cluster_levels,by=c("to"="cluster")) %>% dplyr::arrange(level)
all_nodes = unique(edgelist$to)

##########
### Run annotation
##########

message("Starting annotation: ")
print(parameter_list$manual_names_annotation)

annotation_results = annotate_tree(edgelist = edgelist,#[1:32,],#edgelist = edgelist[1:291,],
                                   labelmat = labelmat,
                                   markers_comparisons_all = markers_comparisons_all,
                                   markers_comparisons_siblings = markers_comparisons_siblings,
                                   preferred_genes=character(0),
                                   manual_names= parameter_list$manual_names_annotation,
                                   overwrite_with_manual = TRUE,
                                   manual_exclude_genes=all_exclusion_genes,
                                   max_pval_adj= parameter_list$max_pval_adj,
                                   min_specificity = parameter_list$min_specificity,
                                   min_specificity_sibling_children= parameter_list$min_specificity_sibling_children,
                                   scale_preferred=1,
                                   limit_factor=parameter_list$limit_factor,
                                   max_score_siblings_children = parameter_list$max_score_siblings_children,
                                   reverse_order = parameter_list$reverse_order)

message("Formating annotation results")

## get annotation
annotation_df = annotation_results$annotation_df
# update cluster names with ID
annotation_df$clean_names_withID = paste0(annotation_df$cluster_id,": ",annotation_df$clean_names)
# get labelmat and add cell ids
labelmat_updated = cbind(Cell_ID=harmonized_seurat_object@meta.data$Cell_ID,mrtree_result$labelmat)
# make a per cell cluster annotation
cell_cluster_map = labelmat_updated %>% as.data.frame() %>% tidyr::gather(-Cell_ID,key="clusterlevel",value="cluster_id")
cell_cluster_map$clusterlevel = gsub("_pruned","",cell_cluster_map$clusterlevel)
# make wide version of labelmat ( that can be added to seurat metadata)
annotation_df_wide = annotation_df %>% dplyr::left_join(cell_cluster_map,by=c("clusterlevel"="clusterlevel","cluster_id"="cluster_id")) %>%
  dplyr::select(Cell_ID,clusterlevel,clean_names_withID)  %>% tidyr::spread(key = clusterlevel,value = clean_names_withID)
colnames(annotation_df_wide)[2:ncol(annotation_df_wide)] = paste0(colnames(annotation_df_wide)[2:ncol(annotation_df_wide)],"_named")
# sort columns in mrtree annotation wide
vec_with_numbers = as.numeric(stringr::str_extract(colnames(annotation_df_wide),"[0-9]+"))
names(vec_with_numbers) = colnames(annotation_df_wide)
sorted_colnames = names(sort(vec_with_numbers,na.last = FALSE))
annotation_df_wide = annotation_df_wide[,sorted_colnames]
# ensure order makes sense
annotation_df_wide = annotation_df_wide[match(harmonized_seurat_object@meta.data$Cell_ID,annotation_df_wide$Cell_ID),]

# remove existing if running again on same object:
# harmonized_seurat_object@meta.data = harmonized_seurat_object@meta.data[,!grepl("C*\\_named",colnames(harmonized_seurat_object@meta.data))]
# # add to seurat
# harmonized_seurat_object@meta.data= dplyr::left_join(harmonized_seurat_object@meta.data,annotation_df_wide,by="Cell_ID")
# rownames(harmonized_seurat_object@meta.data) =  harmonized_seurat_object@meta.data$Cell_ID

## clean markers
descriptive_markers_df = annotation_results$descriptive_markers_df
descriptive_markers_df = dplyr::left_join(descriptive_markers_df,annotation_df[,c("cluster_id","clean_names","clean_names_withID")],by=c("cluster_id"="cluster_id"))

##########
### Save annotation
##########

data.table::fwrite(annotation_df,file = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_annotation_result.txt") ,sep="\t")

data.table::fwrite(annotation_df_wide,file = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_annotation_labelmat.txt") ,sep="\t")

data.table::fwrite(descriptive_markers_df,file = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_annotation_markers_filtered.txt") ,sep="\t")



