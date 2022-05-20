##########
### Load parameters and packages
##########

message(Sys.time(),": Starting fast marker detection .." )

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

# read features to excludes
features_exclude_list= jsonlite::read_json(parameter_list$genes_to_exclude_file)
features_exclude_list = lapply(features_exclude_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# load seurat
harmonized_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,".rds"))

# load mrtree
# TODO

# load markers
# TODO


##########
### Run mrtree pruning based on markers
##########

# TODO: need to update all the below code

# init edgelist and all_nodes
edgelist = edgelist[,c("from","to","level")]
cluster_levels = as.data.frame(labelmat) %>% tidyr::pivot_longer(everything(),names_to = "clusterlevel", values_to = "cluster")  %>% dplyr::group_by(cluster) %>%
  dplyr::add_count(name="ncells") %>% dplyr::distinct(clusterlevel,cluster,ncells)
edgelist = dplyr::left_join(edgelist,cluster_levels,by=c("to"="cluster")) %>% dplyr::arrange(level)
all_nodes = unique(edgelist$to)

##
message("Consolidate clusters before annotation")
# for each node in edgelist:
for(n in 1:length(all_nodes)){

  # get information
  current_node = all_nodes[n]
  #  message("At: ",current_node," with ",length(labelmat[labelmat[,current_level] == current_node,current_level])," cells")
  parent_node = edgelist$from[edgelist$to==current_node]
  sibling_nodes = edgelist$to[edgelist$from==parent_node & edgelist$to != current_node]
  current_level = edgelist$clusterlevel[edgelist$to==current_node]

  # ncells
  sibling_nodes = sibling_nodes[sibling_nodes %in% edgelist$to[edgelist$ncells>=min_cells]]
  # merge with sibling if too few cells
  merge=FALSE
  ncells=edgelist$ncells[edgelist$to==current_node]
  if(ncells<min_cells){
    merge=TRUE
  }
  # sibling markers if siblings exist:
  if(length(sibling_nodes)>0){
    sibling_markers =markers_comparisons_siblings %>% dplyr::filter(cluster_1 %in% sibling_nodes) %>% dplyr::arrange(desc(specificity))%>%
      dplyr::filter(p_val_adj<0.001 & specificity > min_specificity)
    if(nrow(sibling_markers)==0){
      #message("no sibling markers detected!")
      merge=TRUE
    }
  }
  # decide what to merge with:
  merge_node = ""
  if(length(sibling_nodes)==1 & merge){
    #message("  one sibling to merge with")
    merge_node =sibling_nodes[1]
  }else if(length(sibling_nodes)==0 & merge){
    # this should not happen
    merge_node = "cannot merge"
  }else if(length(sibling_nodes)>1 & merge){
    #message("  multiple siblings: determining which to use")
    # this is the difficult case:
    # I am using the intersection size of global markers as a criterion
    node_markers =markers_comparisons_all %>% dplyr::filter(cluster_1 == current_node) %>% dplyr::arrange(desc(specificity)) %>% dplyr::filter(p_val_adj<0.01 & specificity > min_specificity)
    intersection_lengths =c()
    for( s in sibling_nodes){
      sibling_markers = markers_comparisons_all %>% dplyr::filter(cluster_1 == s) %>% dplyr::arrange(desc(specificity)) %>% dplyr::filter(p_val_adj<0.01 & specificity > min_specificity)
      intersection_lengths[s] = length(base::intersect(node_markers$gene,sibling_markers$gene))
    }
    merge_node = names(intersection_lengths)[intersection_lengths==max(intersection_lengths)]
    merge_node = merge_node[1] # if still the same
  }
  # merge & update edgelist
  if(merge){
    message("  Merging: ",current_node," into ",merge_node)
    edgelist$from[edgelist$from == current_node] = merge_node
    edgelist$to[edgelist$to == current_node] = merge_node
    # group
    edgelist = edgelist %>% group_by(from,to,clusterlevel,level) %>% dplyr::summarise(ncells=sum(ncells)) # merge duplicate rows to one entry and updating the cell total
    # update labelmat
    labelmat[labelmat[,current_level] == current_node,current_level] = merge_node
  }
}

##########
### Creat pruned version
##########

# TODO

# TODO: save which clusters changed!

##########
### Save
##########



message("Finalized pruning")




