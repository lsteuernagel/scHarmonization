##########
### Load parameters and packages
##########

message(Sys.time(),": Starting region prediction .." )

message(" Load parameters and packages ")

require(tidyverse)
require(Seurat)
require(Matrix)
require(scCoco)

source("R/harmonization_functions.R")
#source("R/stratified_wilcoxon_functions.R")

# get params-filename from commandline
command_args<-commandArgs(TRUE)
param_file = command_args[1]
# read all parameters and filepaths
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

#test:
# parameter_list = jsonlite::read_json("data/parameters_annotation_v2_1.json")
# parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})
# parameter_list$cluster_column ="C280"

# read features to excludes
features_exclude_list= unlist(jsonlite::read_json(parameter_list$genes_to_exclude_file))

# load seurat
curated_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,".rds"))

# load mrtree clustering
mrtree_result = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_mrtree_clustering_results",".rds"))
#mrtree_result = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_harmonization/backup_marker_tree/hypoMap_harmonized_curated_pruned_mrtree_clustering_results.rds")

##########
### Load marker genes
##########

# read all available marker tables if marker detection was run on multiple subsets !!
markers_comparisons_all_list=list()
for(current_start_node in parameter_list$start_nodes_annotation_markers){
  # load markers all
  filename=paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",current_start_node,"_markers_all_",parameter_list$marker_suffix,".tsv")
  if(file.exists(filename)){
    markers_comparisons_all_list[[current_start_node]] = data.table::fread(filename,data.table = F)
  }else{
    message("Cannot find markers stored in : ",filename)
  }
}
markers_comparisons_all = as.data.frame(do.call(rbind,markers_comparisons_all_list))
message("All markers for: ",length(unique(markers_comparisons_all$cluster_id))," clusters available")

##########
### Clean up marker genes
##########

# load additional remove
#additional_remove_genes = jsonlite::read_json(unlist(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,parameter_list$marker_suffix,"_additionally_removed_markers.json")))
additional_remove_genes = jsonlite::read_json(unlist(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"raw","_additionally_removed_markers.json")))

# which clusters are in cluster_column
#cluster_column_ids = as.character(unique(curated_seurat_object@meta.data[,parameter_list$cluster_column]))
cluster_column_ids = as.character(unique(mrtree_result$labelmat[,parameter_list$cluster_column]))

# subset to nodes below: param_list$start_node_tree
all_children_of_start_node = scUtils::find_children(parameter_list$start_node_tree,mrtree_result$edgelist)


# all_markers_filtered
all_markers_filtered = markers_comparisons_all %>% dplyr::filter(! gene %in% additional_remove_genes &
                                                                   specificity >= parameter_list$min_specificity &
                                                                   p_val_adj < parameter_list$max_pval_adj &
                                                                   cluster_id %in% cluster_column_ids &
                                                                   cluster_id %in% all_children_of_start_node &
                                                                   pct.2 < parameter_list$max_pct_2
)

# manually remove some genes that are broken in ABA ISH
all_markers_filtered = all_markers_filtered[!all_markers_filtered$gene %in% c("Nts","Sim1","Th"),]

# make list per cluster
markers_per_cluster = split(all_markers_filtered$gene,f = all_markers_filtered$cluster_id)
markers_weight_per_cluster = split(all_markers_filtered$specificity,f = all_markers_filtered$cluster_id)
markers_weight_per_cluster = sapply(markers_weight_per_cluster,function(x,maxv=100){x[x>maxv] = maxv; x = x/sum(x)*length(x); return(x)})

##########
### Load aba ish data
##########

# load with filenames to parameter json
aba_gene_to_id = data.table::fread(parameter_list$aba_gene_to_id_file,data.table = FALSE)
aba_ish_matrix = data.table::fread(parameter_list$aba_ish_matrix_file,data.table = FALSE,header = TRUE)
aba_ccf_grid_annotation = readRDS(parameter_list$aba_ccf_grid_annotation_file)
mba_ontology_flatten= data.table::fread(parameter_list$mba_ontology_flatten_file,data.table = FALSE)

##########
### Run region enrichment with scCoco
##########

message(Sys.time(),": Run prediction with ABA ISH .." )

# run scCoco regions per geneset
hypoMap_region_annotation_full = findRegions_genesets(gene_set = markers_per_cluster,
                                                      min_ids = parameter_list$min_ids,
                                                      topn= parameter_list$topn_results,
                                                      target_structure_id = parameter_list$target_structure_id,
                                                      max_ids_to_include = parameter_list$max_ids_to_include,
                                                      #  gene_set_weights = markers_weight_per_cluster,
                                                      aba_gene_to_id = aba_gene_to_id,
                                                      aba_ish_matrix = aba_ish_matrix,
                                                      aba_ccf_grid_annotation = aba_ccf_grid_annotation ,
                                                      mba_ontology_flatten= mba_ontology_flatten,
                                                      target_level = as.character(parameter_list$target_level)
)

# run scCoco summarised regions  per geneset
hypoMap_region_annotation_df = summariseRegions_genesets(findRegion_result = hypoMap_region_annotation_full,min_score=parameter_list$min_score_region_summary)


##########
### Load manual per dataset curation
##########

message(Sys.time(),": Combine with per dataset region curation .." )

require(dplyr)
# load suggested_region_per_dataset
suggested_region_per_dataset = jsonlite::read_json(parameter_list$suggested_region_file)
suggested_region_per_dataset = lapply(suggested_region_per_dataset,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# get params
cluster_column = parameter_list$cluster_column
min_dataset_cell_value = parameter_list$min_dataset_cell_value
dataset_column = parameter_list$dataset_column

# update column name with dataset column
colnames(curated_seurat_object@meta.data)[colnames(curated_seurat_object@meta.data)==dataset_column] = "Dataset"

# add label mat to seurat:
#curated_seurat_object@meta.data = curated_seurat_object@meta.data[,!grepl("C[0-9]+",colnames(curated_seurat_object@meta.data))]
curated_seurat_object@meta.data = cbind(curated_seurat_object@meta.data,mrtree_result$labelmat)

# make a summary of each dataset contribution to each cluster:
total_dataset_counts = curated_seurat_object@meta.data %>% # [map_seurat_object@meta.data[,cluster_column] %in% current_node,]
  dplyr::group_by(Dataset) %>%
  dplyr::summarise(total_dataset_occ = dplyr::n())
total_dataset_counts$total_dataset_occ[total_dataset_counts$total_dataset_occ < min_dataset_cell_value] = min_dataset_cell_value
total_cluster_counts = curated_seurat_object@meta.data %>% # [map_seurat_object@meta.data[,cluster_column] %in% current_node,]
  dplyr::group_by(!!sym(cluster_column)) %>%
  dplyr::summarise(total_cluster_occ = n())
dataset_info = curated_seurat_object@meta.data %>% # [map_seurat_object@meta.data[,cluster_column] %in% current_node,]
  dplyr::group_by(Dataset,!!sym(cluster_column)) %>%
  dplyr::summarise(Dataset_occ = n()) %>% ungroup()  %>%
  dplyr::left_join(total_dataset_counts,by=c("Dataset"="Dataset")) %>%
  dplyr::left_join(total_cluster_counts,by=cluster_column) %>% ungroup() %>%
  dplyr::mutate(adjusted_counts = (Dataset_occ / total_dataset_occ) * total_cluster_occ) %>%
  dplyr::group_by(!!sym(cluster_column)) %>%
  dplyr::mutate(dataset_pct = Dataset_occ / sum(Dataset_occ), adjusted_dataset_pct = adjusted_counts / sum(adjusted_counts)) %>% dplyr::arrange(desc(adjusted_dataset_pct))#

# matrix which regions can be counted from which dataset
suggested_region_per_dataset_matrix = t(scCoco::intersect_from_list(suggested_region_per_dataset))

# matrix with datasets weights epr cluster
dataset_info_wide = dataset_info %>% dplyr::select(Dataset,cluster = !!sym(cluster_column),adjusted_dataset_pct) %>%
  tidyr::spread(key = Dataset,value = adjusted_dataset_pct) %>% as.data.frame()
rownames(dataset_info_wide) = dataset_info_wide$cluster
dataset_info_wide = as.matrix(dataset_info_wide[,2:ncol(dataset_info_wide)])
# order according to other one
dataset_info_wide = dataset_info_wide[,match(rownames(suggested_region_per_dataset_matrix),colnames(dataset_info_wide))]
dataset_info_wide[is.na(dataset_info_wide)] = 0

# multiply the two
## this matrix gives weights for different regions for each cluster based on the (asjusted) contribution of each dataset to the cluster
# the regions each dataset could originate from are manually curated above
regionweight_per_cluster = (dataset_info_wide) %*% suggested_region_per_dataset_matrix

regionweight_per_cluster = as.matrix(regionweight_per_cluster)
regionweight_per_cluster = regionweight_per_cluster + (1- parameter_list$max_region_weight_value)
regionweight_per_cluster[regionweight_per_cluster > 1] = 1


##########
### Combine Allen brain atlas prediction with manual curation per dataset
##########

## matrix from enrichemt scores:
scores_per_target_level_region_all=hypoMap_region_annotation_full$scores_per_target_level_region_all
scores_per_target_level_region_matrix = scores_per_target_level_region_all[,3:ncol(scores_per_target_level_region_all)]
rownames(scores_per_target_level_region_matrix) = scores_per_target_level_region_all$topname

# match the two
shared_regions = intersect(rownames(scores_per_target_level_region_matrix),colnames(regionweight_per_cluster))
scores_per_target_level_region_matrix = scores_per_target_level_region_matrix[rownames(scores_per_target_level_region_matrix) %in% shared_regions,]
regionweight_per_cluster = regionweight_per_cluster[,colnames(regionweight_per_cluster) %in% shared_regions]
scores_per_target_level_region_matrix = scores_per_target_level_region_matrix[match(colnames(regionweight_per_cluster),rownames(scores_per_target_level_region_matrix)),]
# also remove any non-shared clusters
shared_clusters = intersect(colnames(scores_per_target_level_region_matrix),rownames(regionweight_per_cluster))
scores_per_target_level_region_matrix2 = scores_per_target_level_region_matrix[,colnames(scores_per_target_level_region_matrix) %in% shared_clusters]
regionweight_per_cluster = regionweight_per_cluster[rownames(regionweight_per_cluster) %in% shared_clusters,]
scores_per_target_level_region_matrix = scores_per_target_level_region_matrix[,match(rownames(regionweight_per_cluster),colnames(scores_per_target_level_region_matrix))]

# multiply with weights from dataset origin
regionweight_per_cluster_transposed = as.data.frame(t(regionweight_per_cluster))
# version 4: just multiply --> after adjusting regionweight_per_cluster_transposed (via regionweight_per_cluster) with parameter_list$max_region_weight_value
dataset_weight_enrichment_per_cluster =as.data.frame(scores_per_target_level_region_matrix*regionweight_per_cluster_transposed )

# copied from function:
# find top region and other regions
min_score = 0.5
all_clusters = colnames(dataset_weight_enrichment_per_cluster)
result_region_list = list()
for(cluster in all_clusters){
  result_vec = vector()
  result_vec["cluster"] = cluster
  if(cluster %in% colnames(dataset_weight_enrichment_per_cluster)){
    # if(max(dataset_weight_enrichment_per_cluster[,cluster])[1]>min_score){
    # - Region (Target level)
    result_vec["Region"] = rownames(dataset_weight_enrichment_per_cluster)[which(dataset_weight_enrichment_per_cluster[,cluster] == max(dataset_weight_enrichment_per_cluster[,cluster]))[1]]
    # - Most likely region (all levels)
    result_vec["Region_score"] = max(dataset_weight_enrichment_per_cluster[,cluster])[1]
    # - Other likely regions (all levels, cutoff ?)
    top_5_scores = sort(dataset_weight_enrichment_per_cluster[,cluster],decreasing = TRUE)[1:5]
    top_5_scores = top_5_scores[top_5_scores>min_score]
    top_5_scores = dataset_weight_enrichment_per_cluster[dataset_weight_enrichment_per_cluster[,cluster] %in% top_5_scores,c(1,2,which(colnames(dataset_weight_enrichment_per_cluster)==cluster))]
    top_5_scores = top_5_scores[order(top_5_scores[,cluster],decreasing=TRUE),]
    result_vec["Region_other"] = possible_regions = paste0(rownames(top_5_scores)[-1],collapse = " | ")
    # }else{
    #   result_vec =c(result_vec,rep(NA,4))
    #   names(result_vec)[2:5] = c("Region","Region_specific","Region_specific_enrichment","Region_specific_other")
    # }
  }else{
    result_vec =c(result_vec,rep(NA,3))
    names(result_vec)[2:4] = c("Region","Region_score","Region_other")
  }
  # result
  result_region_list[[cluster]] = result_vec
}

# combine to result
result_region = as.data.frame(do.call(rbind,result_region_list))

##########
### Add QC and power metrics of prediction
##########

## add qc
tmp_qc = hypoMap_region_annotation_full$gene_set_power
tmp_qc$cluster = rownames(tmp_qc)
result_region = dplyr::left_join(result_region,tmp_qc,by="cluster")

# filter out cluster prediction with low confidence:
result_region$Region_curated = result_region$Region
result_region$Region_curated[result_region$Region_score < parameter_list$min_score_region_summary] = NA
result_region$Region_curated[result_region$pct_of_genes < parameter_list$min_pct_of_genes | result_region$number_ids < parameter_list$min_number_of_ids] = NA
result_region$Region_other[is.na(result_region$Region_curated)] = paste0(result_region$Region[is.na(result_region$Region_curated)]," | ",result_region$Region_other[is.na(result_region$Region_curated)])

##########
### Save results
##########

message(Sys.time(),": Save region prediction .." )

# main results
data.table::fwrite(result_region,file = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_region_prediction.txt") ,sep="\t")

# all other information as rds
full_res_list = list(full_region_prediction_result = hypoMap_region_annotation_full,
                     dataset_pct_adjusted_per_region = dataset_info_wide,
                     regionweight_per_cluster = regionweight_per_cluster,
                     dataset_weight_enrichment_per_cluster = dataset_weight_enrichment_per_cluster)
saveRDS(full_res_list,paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_region_prediction_all_results.rds"))

message(Sys.time(),": Finalized .." )





