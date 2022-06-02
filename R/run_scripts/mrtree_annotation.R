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
# parameter_list = jsonlite::read_json("data/parameters_harmonization_v2_1_test.json")
# parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})
# parameter_list$marker_suffix = "pruned"
# parameter_list$new_name_suffix=paste0(parameter_list$new_name_suffix,"_curated")
# parameter_list$start_node = "C2-1"

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

hypoMap_test_curated_C21_markers_all_pruned = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_harmonization_test/hypoMap_test_curated_C2-1_markers_all_pruned.tsv")
hypoMap_test_curated_C21_markers_siblings_pruned = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_harmonization_test/hypoMap_test_curated_C2-1_markers_siblings_pruned.tsv")

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
### annotate_tree
##########

find_children = function(nodes,edges){
  current_children = edges$to[edges$from %in% nodes]
  #print(paste0(current_children,collapse = "|"))
  if(length(current_children)>0){
    all_children = c(current_children,find_children(current_children,edges))
  }else{
    all_children = current_children
  }
  return(all_children)
}

#' Add gene absed annotation to tree labels
#' @param seurat_object seurat_object to call FindMarkers
#' @param edgelist edgelist from mrtree corresponding to labelmat
#' @param labelmat dataframe of clusterlabels per cell with one column per level
#' @param markers_comparisons_all cluster markers for all clusters in labelmat. expects format from run_marker_detection.R
#' @param markers_comparisons_siblings cluster markers vs sibling for all clusters in labelmat. expects format from run_marker_detection.R
#' @param manual_names a named vector with names of clusters in edgelist --> will overwrite with manual names
#' @param manual_exclude_genes a set of genes that should not be included in anno
#' @param max_pval_adj max pavlue allowed for a marker to be considered
#' @param min_specificity min_specificity allowed for a marker to be considered
#' @param min_specificity_sibling_children maximum specificity allowed for a marker in sibling comparisons of children to be kept --> strong markers between children are probably not good for whole cluster
#' @param min_pct2_score when calculating the specificty score add this pseudo-count to pct.2 to avoid favouring of highly specific but low abundant genes
#' @param scale_preferred factor to scale scoe of preferred genes
#' @param preferred_genes a vector with genes that should be preferred
#' @param min_cells min_cells
#' @param limit_factor scale sibling or children_siblings scores down to be not larger than limit_factor*global_score
#' @param max_score_siblings_children todo: add explanation
#' @return dataframe with new labels for each cluster in edgelist$to

annotate_tree = function(edgelist,labelmat,markers_comparisons_all,markers_comparisons_siblings,preferred_genes=character(0),manual_names=c(),manual_exclude_genes=character(0),
                         max_pval_adj=0.0001, min_specificity = 0.5,min_specificity_sibling_children=10,scale_preferred=2,min_pct2_score=0.01,min_cells=20,limit_factor=5,max_score_siblings_children=20,
                         reverse_order = FALSE){
  # init edgelist and all_nodes
  edgelist = edgelist[,c("from","to","level")]
  cluster_levels = as.data.frame(labelmat) %>% tidyr::pivot_longer(everything(),names_to = "clusterlevel", values_to = "cluster")  %>% dplyr::group_by(cluster) %>%
    dplyr::add_count(name="ncells") %>% dplyr::distinct(clusterlevel,cluster,ncells)
  edgelist = dplyr::left_join(edgelist,cluster_levels,by=c("to"="cluster")) %>% dplyr::arrange(level)
  all_nodes = unique(edgelist$to)

  # annotation
  edgelist = edgelist %>% dplyr::arrange(level) # ensure top down order for tree traversal with for loop!+

  # list to store intermediate results
  annotation_list = list()
  descriptive_markers_list = list()

  message("Running annotation")

  # for each node in edgelist:
  for(n in 1:length(all_nodes)){

    # get information
    current_node = all_nodes[n]
    parent_node = edgelist$from[edgelist$to==current_node]
    sibling_nodes = edgelist$to[edgelist$from==parent_node & edgelist$to != current_node]
    children_nodes = find_children(current_node,edgelist)
    direct_children_nodes = edgelist$to[edgelist$from==current_node]
    current_level = edgelist$clusterlevel[edgelist$to==current_node]

    # get specific genes
    potential_descriptive_markers =markers_comparisons_all %>% dplyr::filter(cluster_id == current_node) %>% dplyr::arrange(desc(specificity)) %>%
      dplyr::filter(p_val_adj< max_pval_adj & specificity > min_specificity)

    # calculate score as adjusted specificity with a minimum value on pct.2 to put more emphasis on pct.1 (abundant markers)
    potential_descriptive_markers$pct.2_min = potential_descriptive_markers$pct.2
    potential_descriptive_markers$pct.2_min[potential_descriptive_markers$pct.2_min < min_pct2_score] = min_pct2_score
    potential_descriptive_markers$score = (potential_descriptive_markers$pct.1 / potential_descriptive_markers$pct.2_min) * potential_descriptive_markers$avg_log2FC

    # initiate vector with genes that should be excluded!
    exclude_genes=c(manual_exclude_genes)

    message("nrow(potential_descriptive_markers) 1: ",nrow(potential_descriptive_markers))

    # eliminate parent and sibling names
    reserved_genes=c()
    parent_name=""
    if(length(annotation_list)>0){
      parent_name = annotation_list[[parent_node]]
      if(parent_node %in% names(annotation_list)){
        reserved_genes = c(reserved_genes,stringr::str_split(parent_name,pattern = "\\.")[[1]])
      }
      if(length(sibling_nodes)>0){
        for(s in 1:length(sibling_nodes)){
          if(sibling_nodes[s] %in% names(annotation_list)){
            reserved_genes = c(reserved_genes,stringr::str_split(annotation_list[[sibling_nodes[s]]],pattern = "\\.")[[1]])
          }
        }
      }
    }else{
      reserved_genes=c()
    }
    # if a gene is in preferred_genes
    if(nrow(potential_descriptive_markers)>0){
      potential_descriptive_markers$score[potential_descriptive_markers$gene %in% preferred_genes] = potential_descriptive_markers$score[potential_descriptive_markers$gene %in% preferred_genes] * scale_preferred
    }
    ####
    # ranks in siblings and in children to compare with global rank in own markers
    ####

    if(nrow(potential_descriptive_markers)>0){
      # get sibling_markers of children and get inveretd ranks
      sibling_markers_children =markers_comparisons_siblings %>% dplyr::filter(cluster_id %in% children_nodes) %>% dplyr::arrange(desc(specificity))%>%
        dplyr::filter(p_val_adj<max_pval_adj & specificity > min_specificity_sibling_children)
      if(nrow(sibling_markers_children)>0){
        # calculate score and only keep highest children cluster
        sibling_markers_children$pct.2_min = sibling_markers_children$pct.2
        sibling_markers_children$pct.2_min[sibling_markers_children$pct.2_min < min_pct2_score] = min_pct2_score
        sibling_markers_children$score = (sibling_markers_children$pct.1 / sibling_markers_children$pct.2_min) * sibling_markers_children$avg_log2FC
        sibling_markers_children = sibling_markers_children %>% dplyr::group_by(gene) %>% dplyr::filter(score == max(score)) %>%
          dplyr::distinct(gene,.keep_all =TRUE)
        # filter by specificity
        sibling_markers_children =sibling_markers_children %>% dplyr::filter(! gene %in% exclude_genes) %>% dplyr::arrange(desc(score)) #  & score > min_specificity
        # join score with potential_descriptive_markers
        if(nrow(sibling_markers_children)>0){
          sibling_markers_children$score_siblings_children = sibling_markers_children$score
          potential_descriptive_markers = dplyr::left_join(potential_descriptive_markers,sibling_markers_children[,c("gene","score_siblings_children")],by="gene")
          potential_descriptive_markers$score_siblings_children[is.na(potential_descriptive_markers$score_siblings_children)] = 0 # if not a marker within children gene gets rank 1
        }else{
          potential_descriptive_markers$score_siblings_children = 0
        }
      }else{
        potential_descriptive_markers$score_siblings_children = 0
      }
      message("nrow(potential_descriptive_markers): ",nrow(potential_descriptive_markers))

      # check siblings and filter to genes that are also markers to siblings!
      if(length(sibling_nodes)>0){
        # repeat most steps for sibling markers
        sibling_markers =markers_comparisons_siblings %>% dplyr::filter(cluster_id == current_node) %>% dplyr::arrange(desc(specificity)) %>% dplyr::filter(p_val_adj<max_pval_adj & specificity > min_specificity)
        if(nrow(sibling_markers)>0){
          # calculate score, scale and filter
          sibling_markers$pct.2_min = sibling_markers$pct.2
          sibling_markers$pct.2_min[sibling_markers$pct.2_min < min_pct2_score] = min_pct2_score
          sibling_markers$score = (sibling_markers$pct.1 / sibling_markers$pct.2_min) * sibling_markers$avg_log2FC
          sibling_markers$score[sibling_markers$gene %in% preferred_genes] = sibling_markers$score[sibling_markers$gene %in% preferred_genes] * scale_preferred
          sibling_markers =sibling_markers %>% dplyr::filter(! gene %in% exclude_genes & score > min_specificity) %>% dplyr::arrange(desc(score))
          # join
          if(nrow(sibling_markers)>0){
            sibling_markers$score_siblings = sibling_markers$score
            potential_descriptive_markers = dplyr::left_join(potential_descriptive_markers,sibling_markers[,c("gene","score_siblings")],by="gene")
            potential_descriptive_markers$score_siblings[is.na(potential_descriptive_markers$score_siblings)] = 0
            # potential_descriptive_markers$score_siblings = potential_descriptive_markers$score_siblings * potential_descriptive_markers$children_fraction
          }else{
            potential_descriptive_markers$score_siblings = 0
          }
        }else{
          potential_descriptive_markers$score_siblings =0
        }
      }else{
        # fallback if nor sibling markers are available
        potential_descriptive_markers$score_siblings = 0
      }

      # get ranks by score
      if(nrow(potential_descriptive_markers)>0){
        # limit score_siblings_children and score_siblings to omit strong outliers
        potential_descriptive_markers$score_siblings_children[potential_descriptive_markers$score_siblings_children > limit_factor*potential_descriptive_markers$score] = limit_factor*potential_descriptive_markers$score[potential_descriptive_markers$score_siblings_children > limit_factor*potential_descriptive_markers$score]
        potential_descriptive_markers$score_siblings[potential_descriptive_markers$score_siblings > limit_factor*potential_descriptive_markers$score] = limit_factor*potential_descriptive_markers$score[potential_descriptive_markers$score_siblings > limit_factor*potential_descriptive_markers$score]

        # calculate average score as: score + score_siblings - score_siblings_children
        potential_descriptive_markers$avg_score = potential_descriptive_markers$score + potential_descriptive_markers$score_siblings  - potential_descriptive_markers$score_siblings_children
        potential_descriptive_markers$mult_score = potential_descriptive_markers$score * (1+potential_descriptive_markers$score_siblings) /  (1+potential_descriptive_markers$score_siblings_children)


        # prepare for selection: exclude all genes in exclude_genes vector, filter score again after all adjustements and arrange by score
        potential_descriptive_markers =potential_descriptive_markers %>%
          dplyr::filter(! gene %in% exclude_genes & avg_score > min_specificity  & score_siblings_children < max_score_siblings_children) %>% dplyr::arrange(desc(avg_score))
      }
    }
    # take top result as cluster descriptive gene
    if(nrow(potential_descriptive_markers)>0){
      # don't select sibling or parent gene
      potential_descriptive_markers_sel =potential_descriptive_markers %>% dplyr::filter(! gene %in% reserved_genes) %>% dplyr::arrange(desc(avg_score))
      #reserved_genes
      # select top gene
      name_gene = potential_descriptive_markers_sel$gene[1]
      # if node has no siblings: don't need annotation!
      if(length(sibling_nodes)==0){
        name_gene="d"
      }
    }else{
      # if merge take name from sibling!
      name_gene = "problematic"
    }
    # check manual overwrite vector
    if(current_node %in% names(manual_names)){
      name_gene = manual_names[current_node]
    }
    # update full cluster_name
    cluster_name = paste0(parent_name,".",name_gene)
    message(current_node,": name: ",cluster_name)

    # add to list
    annotation_list[[current_node]] = cluster_name
    if(nrow(potential_descriptive_markers)){
      descriptive_markers_list[[current_node]] = potential_descriptive_markers %>% dplyr::select(cluster_id,parent,gene,p_val_adj,avg_logFC,pct.1,pct.2,specificity,score,score_siblings,score_siblings_children,avg_score,mult_score)
    }
  }
  # also save a list of the good marker genes
  descriptive_markers_df = as.data.frame(do.call(rbind,descriptive_markers_list))

  # bind results and format
  annotation_df = as.data.frame(do.call(rbind,annotation_list))
  annotation_df = annotation_df %>% dplyr::rename(Map_CellType=V1) %>% dplyr::mutate(cluster_id = rownames(annotation_df))
  if(reverse_order){
    annotation_df$Map_CellType = sapply(annotation_df$Map_CellType,function(x){paste0(rev(strsplit(x,split = "\\.")[[1]]),collapse = ".")})
    annotation_df$clean_names = sub("\\.$","",gsub("d\\.","",annotation_df$Map_CellType))#sub(".","",gsub("\\.d","",annotation_df$Map_CellType))
  }else{
    annotation_df$clean_names = sub("\\.","",gsub("\\.d","",annotation_df$Map_CellType))
  }
  annotation_df = dplyr::left_join(annotation_df,cluster_levels,by=c("cluster_id"="cluster"))
  annotation_df$clean_names[is.na(annotation_df$clean_names)]="all"
  # named_edgelist
  # named_edgelist = dplyr::left_join(edgelist,annotation_df[,c("clean_names","cluster_id")],by=c("from"="cluster_id"))
  # named_edgelist = dplyr::left_join(named_edgelist,annotation_df[,c("clean_names","cluster_id")],by=c("to"="cluster_id"))
  # named_edgelist$from = named_edgelist$clean_names.x
  # named_edgelist$from[is.na(named_edgelist$from)] = "all"
  # named_edgelist$to = named_edgelist$clean_names.y
  # named_edgelist = named_edgelist[,c("from","to","level")]
  # return
  message("Returning results")
  return(list(annotation_df = annotation_df,descriptive_markers_df = descriptive_markers_df))
}

##########
### Run annotation
##########
parameter_list$manual_names_annotation = c("C2-1" = "Neurons","C2-2"="Non-Neurons")
parameter_list$min_specificity = 1.5 # ?
parameter_list$min_specificity_sibling_children = 2.5
parameter_list$limit_factor = 5
parameter_list$max_score_siblings_children = 20
parameter_list$reverse_order = TRUE

annotation_results = annotate_tree(edgelist = edgelist,
                                   labelmat = labelmat,
                                   markers_comparisons_all = markers_comparisons_all,
                                   markers_comparisons_siblings = markers_comparisons_siblings,
                                   preferred_genes=character(0),
                                   manual_names= parameter_list$manual_names_annotation,
                                   manual_exclude_genes=features_exclude_list,
                                   max_pval_adj= parameter_list$max_pvalue_prune,
                                   min_specificity = parameter_list$min_specificity,
                                   min_specificity_sibling_children= parameter_list$min_specificity_sibling_children,
                                   scale_preferred=1,
                                   limit_factor=parameter_list$limit_factor,
                                   max_score_siblings_children = parameter_list$max_score_siblings_children,
                                   reverse_order = parameter_list$reverse_order)

message("Formating annotation results")

# remove existing if running again on same object:
#seurat_object_harmonized@meta.data = seurat_object_harmonized@meta.data[,!grepl("K*\\_pruned",colnames(seurat_object_harmonized@meta.data))]
# add to seurat
seurat_object_harmonized@meta.data = cbind(seurat_object_harmonized@meta.data,temp)

## get annotation names
annotation_df = annotation_results$annotation_df
annotation_df$clean_names[annotation_df$clean_names==""] = "hypothalamus"
cell_cluster_map =seurat_object_harmonized@meta.data[,c("Cell_ID",new_pruned_names)] %>% tidyr::gather(-Cell_ID,key="clusterlevel",value="cluster_id")
cell_cluster_map$clusterlevel = gsub("_pruned","",cell_cluster_map$clusterlevel)
annotation_df_wide = annotation_df%>% dplyr::left_join(cell_cluster_map,by=c("clusterlevel"="clusterlevel","cluster_id"="cluster_id")) %>%
  dplyr::select(Cell_ID,clusterlevel,clean_names)  %>% tidyr::spread(key = clusterlevel,value = clean_names)
colnames(annotation_df_wide)[2:ncol(annotation_df_wide)] = paste0(colnames(annotation_df_wide)[2:ncol(annotation_df_wide)],"_named")

# remove existing if running again on same object:
seurat_object_harmonized@meta.data = seurat_object_harmonized@meta.data[,!grepl("K*\\_named",colnames(seurat_object_harmonized@meta.data))]
# add to seurat
seurat_object_harmonized@meta.data= dplyr::left_join(seurat_object_harmonized@meta.data,annotation_df_wide,by="Cell_ID")
rownames(seurat_object_harmonized@meta.data) =  seurat_object_harmonized@meta.data$Cell_ID

#data.table::fwrite(annotation_df,file = paste0(harmonization_file_path,project_name,"_annotation_labels.txt") ,sep="\t")

# add annotation_df to seurat misc
annotation_df_add = annotation_df %>% dplyr::select(cluster_id,cluster_name = clean_names, clusterlevel, ncells)
seurat_object_harmonized@misc$annotations = annotation_df_add

## clean markers
descriptive_markers_df = annotation_results$descriptive_markers_df
descriptive_markers_df = dplyr::left_join(descriptive_markers_df,annotation_df,by=c("cluster_id"="cluster_id"))
seurat_object_harmonized@misc$curated_markers = descriptive_markers_df
#data.table::fwrite(descriptive_markers_df,file = paste0(harmonization_file_path,project_name,"_curated_markers.txt") ,sep="\t")

