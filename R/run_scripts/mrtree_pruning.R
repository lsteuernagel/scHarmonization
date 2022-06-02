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
#parameter_list = jsonlite::read_json("data/parameters_harmonization_v2_1_test.json")
#parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})
#parameter_list$marker_suffix = "raw"
#parameter_list$new_name_suffix=paste0(parameter_list$new_name_suffix,"_curated")
#parameter_list$start_node = "K2-0"

# load seurat
harmonized_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,".rds"))

# load mrtree clustering
mrtree_result = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_mrtree_clustering_results",".rds"))


# read all available marker tables if marker detection was run on multiple subsets !!
markers_comparisons_all_list=list()
markers_comparisons_siblings_list =list()
for(current_start_node in parameter_list$start_nodes_pruning_markers){
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
message("Sibling markers for: ",length(unique(markers_comparisons_siblings$cluster_id))," clusters available")


##########
### prepare raw edgelists
##########

# get key objects
edgelist = mrtree_result$edgelist
labelmat = mrtree_result$labelmat

# add to seurat:
harmonized_seurat_object@meta.data = cbind(harmonized_seurat_object@meta.data,labelmat)

# find children in tree recusrively based on simple edgelist
# find_children = function(nodes,edges){
#   current_children = edges$to[edges$from %in% nodes]
#   #print(paste0(current_children,collapse = "|"))
#   if(length(current_children)>0){
#     all_children = c(current_children,find_children(current_children,edges))
#   }else{
#     all_children = current_children
#   }
#   return(all_children)
# }
#
# # walk through tree:
# all_children = find_children(nodes = parameter_list$start_node, edges = edgelist[,1:2])
#
# # subset =
# edgelist = edgelist[edgelist$to %in% all_children,]


##########
### Run mrtree pruning based on markers
##########

message(Sys.time(),": Prune clusters .." )

# init edgelist and all_nodes
edgelist = edgelist[,c("from","to","level")]
cluster_levels = as.data.frame(labelmat) %>% tidyr::pivot_longer(everything(),names_to = "clusterlevel", values_to = "cluster")  %>% dplyr::group_by(cluster) %>%
  dplyr::add_count(name="ncells") %>% dplyr::distinct(clusterlevel,cluster,ncells)
edgelist = dplyr::left_join(edgelist,cluster_levels,by=c("to"="cluster")) %>% dplyr::arrange(level)
all_nodes = unique(edgelist$to)

## iterate over all nodes
merge_list = list()
for(n in 1:length(all_nodes)){

  # get information
  current_node = all_nodes[n]
  #  message("At: ",current_node," with ",length(labelmat[labelmat[,current_level] == current_node,current_level])," cells")
  parent_node = edgelist$from[edgelist$to==current_node]
  sibling_nodes = edgelist$to[edgelist$from==parent_node & edgelist$to != current_node]
  current_level = edgelist$clusterlevel[edgelist$to==current_node]
  if(edgelist$level[edgelist$to==current_node] < parameter_list$min_prune_level){next}

  # filter siblings by ncells
  sibling_nodes = sibling_nodes[sibling_nodes %in% edgelist$to[edgelist$ncells >= parameter_list$min_cells]]
  if(length(sibling_nodes) == 0){next}
  #ncells of current
  ncells_current_node=edgelist$ncells[edgelist$to==current_node]

  # get number of sibling markers:
  sibling_markers = markers_comparisons_siblings %>% dplyr::filter(cluster_id == current_node) %>% dplyr::arrange(desc(specificity)) %>%
    dplyr::filter( p_val_adj < parameter_list$max_pvalue_prune & specificity > parameter_list$min_specificity)

  # get global markers:
  global_markers = markers_comparisons_all %>% dplyr::filter(cluster_id == current_node) %>% dplyr::arrange(desc(specificity)) %>%
    dplyr::filter( p_val_adj < parameter_list$max_pvalue_prune & specificity > parameter_list$min_specificity)

  # parent global markers:
  parent_global_markers = markers_comparisons_all %>% dplyr::filter(cluster_id == parent_node) %>% dplyr::arrange(desc(specificity)) %>%
    dplyr::filter( p_val_adj < parameter_list$max_pvalue_prune & specificity > parameter_list$min_specificity)

  # if any n_sibling_markers are < min_sibling_markers
  if(nrow(sibling_markers) < parameter_list$min_sibling_markers | ncells_current_node < parameter_list$min_cells){

    if(nrow(parent_global_markers) > 3){
      # get the expression of the parent markers in current cluster:
      current_gene_expression = FetchData(harmonized_seurat_object,vars = parent_global_markers$gene,cells = harmonized_seurat_object@meta.data$Cell_ID[harmonized_seurat_object@meta.data[,current_level] == current_node])
      current_gene_expression_mean = colMeans(current_gene_expression)

      # iterate over all sibling markers:
      intersection_lengths = c()
      coexpr = c()
      for(sib in sibling_nodes){
        # print(sib)
        # global_sibling_markers = markers_comparisons_all %>% dplyr::filter(cluster_id == sib) %>% dplyr::arrange(desc(specificity)) %>%
        #   dplyr::filter( p_val_adj < parameter_list$max_pvalue_prune & specificity > parameter_list$min_specificity)
        # print(head(global_sibling_markers))
        # intersection_lengths[sib] = length(base::intersect(global_sibling_markers$gene,global_markers$gene))
        #
        sibling_gene_expression = FetchData(harmonized_seurat_object,vars = parent_global_markers$gene,cells = harmonized_seurat_object@meta.data$Cell_ID[harmonized_seurat_object@meta.data[,current_level] == sib])
        sibling_gene_expression_mean = colMeans(sibling_gene_expression)
        coexpr[sib] = cor(current_gene_expression_mean,sibling_gene_expression_mean)
      }

      # add node with highest number of shared global markers
      merge_nodes = names(coexpr)[coexpr == max(coexpr)]
    }else{
      merge_nodes = sibling_nodes[1] # fall back
    }
    message("Merging node ",current_node," (",n,"/",length(all_nodes),") into node(s) ",paste0(merge_nodes,sep=", "))
    # print(coexpr)
    all_nodes_to_merge = c(merge_nodes,current_node)
    merge_list[[paste0(all_nodes_to_merge,collapse = "_")]] = all_nodes_to_merge
  }


}

##########
### Create pruned edgelist
##########

message(Sys.time(),": Update labels and save .." )

edgelist_updated = edgelist
labelmat_updated = labelmat
merge_list2 = merge_list

for(i in 1:length(merge_list2)){
  nodes_to_merge = merge_list2[[i]]
  print(nodes_to_merge)
  # only update edgelist and labelmat if there are truly 2 unique labels remaining in the current entry
  if(length(unique(nodes_to_merge)) > 1){
    # use the first node to label
    merge_node = as.character(sort((nodes_to_merge))[1])
    print(paste0(" >> merge to ",merge_node))
    # update edglist
    edgelist_updated$from[edgelist_updated$from %in% nodes_to_merge] = merge_node
    edgelist_updated$to[edgelist_updated$to %in% nodes_to_merge] = merge_node
    # remove repeated entries
    edgelist_updated = edgelist_updated %>% distinct(from,to,clusterlevel) # merge duplicate rows to one entry and updating the cell total
    # update labelmat
    current_level = edgelist$clusterlevel[edgelist$to==merge_node]
    labelmat_updated[labelmat_updated[,current_level] %in% nodes_to_merge,current_level] = merge_node
    # change all occurrences in merge_list
    merge_list2 = lapply(merge_list2,function(x,remove_nodes,new_node){x[x %in% remove_nodes] = new_node;return(x)},remove_nodes=nodes_to_merge,new_node=merge_node)

  }
}

##########
### Update labels in edgelist and labelmat
##########

old_prefix = parameter_list$old_prefix # "K"
new_prefix = parameter_list$new_prefix # "C"
all_cluster_levels_updated = edgelist_updated %>% dplyr::group_by(clusterlevel) %>% dplyr::count()
all_cluster_levels_updated$clusterlevel_new = paste0(all_cluster_levels_updated$clusterlevel %>% stringr::str_extract(old_prefix)  %>% stringr::str_replace(old_prefix, new_prefix),all_cluster_levels_updated$n)

# create new labels using the pruned numbers of labels per level
all_rows = list()
for(c in 1:nrow(all_cluster_levels_updated)){
  edgelist_updated_current_level = edgelist_updated[edgelist_updated$clusterlevel == all_cluster_levels_updated$clusterlevel[c],]
  for(i in 1:nrow(edgelist_updated_current_level)){
    new_node_label = data.frame(old_node = edgelist_updated_current_level$to[i],
                                new_node = paste0(all_cluster_levels_updated$clusterlevel_new[c],"-",i),
                                new_cluster_level = all_cluster_levels_updated$clusterlevel_new[c])
    all_rows[[new_node_label$new_node]] = new_node_label
  }
}
new_labels = as.data.frame(do.call(rbind,all_rows))

# make new edgelist
edgelist_updated_new_labels = dplyr::left_join(edgelist_updated,new_labels[,c("old_node","new_node")],by=c("from"="old_node"))
edgelist_updated_new_labels = dplyr::left_join(edgelist_updated_new_labels,new_labels,by=c("to"="old_node"))
edgelist_updated_new_labels = edgelist_updated_new_labels %>% dplyr::select(from = new_node.x, to = new_node.y,clusterlevel = new_cluster_level)
edgelist_updated_new_labels$from[is.na(edgelist_updated_new_labels$from)] = "all"

# make new labelmat
labelmat_updated_new_labels = apply(labelmat_updated,2,function(x,new_labels){
  new_labels$new_node[match(x,new_labels$old_node)]
},new_labels=new_labels)

colnames(labelmat_updated_new_labels) = stringr::str_extract(labelmat_updated_new_labels[1,],pattern = paste0(parameter_list$new_prefix,"[0-9]+"))

##########
### Creat pruned result version
##########

# save in data.tree format
require(data.tree)
df = as.data.frame(unique(labelmat_updated_new_labels), stringsAsFactors = F)
df$pathString = apply(df, 1, function(x) paste(c('all', x), collapse='/'))
tree.datatree = data.tree::as.Node(df)

# export edgelist
edges= data.tree::ToDataFrameNetwork(tree.datatree,"isLeaf","level","count","totalCount","height")
nodes = data.frame(id = c("all",as.character(unique(edges$to))),label=c("all",as.character(unique(edges$to))))
nodes = rbind(c("all",FALSE,1,5,max(edges$height)+1),edges[,2:ncol(edges)]) %>% dplyr::rename(id = to) %>% dplyr::mutate(label = id)

# make cluster object ? and add to misc
updated_cluster_object = list(labelmat = labelmat_updated_new_labels,
                              edgelist = edges ,
                              nodelist = nodes,
                              data_tree = tree.datatree)
##########
### Save
##########

saveRDS(updated_cluster_object,paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_pruned_mrtree_clustering_results",".rds"))

message("Finalized pruning")




