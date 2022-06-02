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

# load mrtree clustering
mrtree_result = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_mrtree_clustering_results",".rds"))

# important: start node
start_node = parameter_list$start_node

if(is.null(start_node)){
  start_node = "all"
}

# markers
n_cores_markers = parameter_list$n_cores_markers
assay_markers = parameter_list$assay_markers
assay_slot =parameter_list$assay_slot
logfc.threshold = parameter_list$logfc.threshold
min.pct =parameter_list$min.pct
min.diff.pct = parameter_list$min.diff.pct
max.cells.per.ident = parameter_list$max.cells.per.ident
min.cells.feature = parameter_list$min.cells.feature
min.cells.group =  parameter_list$min.cells.group
base = parameter_list$base
only.pos = parameter_list$only.pos
batch_var = parameter_list$batch_var
test.use = parameter_list$test.use
specificity_base = parameter_list$specificity_base
message("test.use ",test.use," ... ",parameter_list$test.use)
if(is.null(base)){base = 2}

##########
### findMarkers_tree2
##########

#' Walk through a mrtree / hierchical clustering tree using an edgelist and calculate Marker genes using Seurat's FindMarkers or a stratified version:
#' https://github.com/KChen-lab/stratified-tests-for-seurat
#' Two results: Vs sibling and Vs all . Runs in parallel to be faster
#' @param seurat_object seurat_object to call FindMarkers
#' @param edgelist minimum number of cells to keep cluster
#' @param labelmat
#' @param n_cores doParallel cores
#' @param use_stratified use stratfied
#' @param batch_var batch_var for stratfied mode
#' @param assay
#' @param slot
#' @param ... passed to detection functions
#' @return updated vector of labels

findMarkers_tree2 = function(seurat_object,edgelist,labelmat,n_cores=1,use_stratified=TRUE,test.use="wilcox-stratified",batch_var="Batch_ID",genes_to_include=NULL,assay="RNA",slot="data",...){

  require(doParallel)
  require(tidyselect)
  require(foreach)

  # info:
  message(n_cores)

  # set edglist
  edgelist = edgelist[,c("from","to")]

  # genes
  if(is.null(genes_to_include)){
    genes_to_include = rownames(seurat_object@assays[[assay]]$counts)
  }

  # add labelmat to seurat
  if(any(colnames(labelmat) %in% colnames(seurat_object@meta.data))){
    if(!all(colnames(labelmat) %in% colnames(seurat_object@meta.data))){
      stop("Cannot handle duplicated column names in seurat_object@meta.data and labelmat from mrtree. Either include all labelmat columns beforehand or none at all.")
    }
  }else{
    seurat_object@meta.data = cbind(seurat_object@meta.data,labelmat )
  }

  # init
  current_node="all"
  colnames(edgelist) = c("from","to")
  label_mat_long = as.data.frame(labelmat) %>% tidyr::pivot_longer(everything(),names_to = "clusterlevel", values_to = "cluster")
  edgelist = dplyr::left_join(edgelist,label_mat_long,by=c("to"="cluster")) %>% dplyr::distinct(from,to,clusterlevel)
  all_nodes = unique(edgelist[,2])
  #all_nodes = all_nodes[!all_nodes %in% c("all","root")]

  comparisons_siblings = NULL
  comparisons_all = NULL

  registerDoParallel(cores=n_cores)
  comparison_types = c("Sibling","All")
  return_list=list()
  for(comp in comparison_types){
    if(comp == "Sibling"){
      message("Sibling Comparisons")
    }
    if(comp == "All"){
      message("All Comparisons")
    }
    comparisons_main <- foreach(n = 1:length(all_nodes),.errorhandling = 'remove', .combine='rbind') %dopar% {
      current_node = all_nodes[n]
      parent_node = edgelist$from[edgelist$to==current_node]
      sibling_nodes = edgelist$to[edgelist$from==parent_node & edgelist$to != current_node]
      current_level = edgelist$clusterlevel[edgelist$to==current_node]
      # set ident to current level!
      Idents(seurat_object) = current_level

      #decide what to cal against:
      cluster_1 = current_node
      if(comp == "Sibling"){
        cluster_2 = sibling_nodes
      }
      if(comp == "All"){
        cluster_2 = all_nodes[! all_nodes %in% current_node]
        cluster_2 = cluster_2[cluster_2 %in% edgelist$to[edgelist$clusterlevel == current_level]] # need to ensure that only current level is used!
      }
      #calculate markers
      #  current_markers <- tryCatch({
      if(length(cluster_2)>0){
        if(test.use == "wilcox-stratified"){
          current_markers = FindMarkers2.Seurat(object = seurat_object,
                                                ident.1 = cluster_1,
                                                ident.2 = cluster_2,
                                                assay = assay_markers,
                                                features = genes_to_include,
                                                logfc.threshold = logfc.threshold,
                                                slot = assay_slot,
                                                test.use = "VE",
                                                min.pct = min.pct,
                                                min.diff.pct = min.diff.pct,
                                                max.cells.per.ident=max.cells.per.ident,
                                                min.cells.feature = min.cells.feature,
                                                min.cells.group = min.cells.group,
                                                return.thresh = 1,
                                                base = base,
                                                only.pos = only.pos,
                                                latent.vars = batch_var,
                                                genre = "locally-best")
          if(base==2){
            colnames(current_markers)[colnames(current_markers)=="avg_logFC"] = "avg_log2FC"
          }
        }else{
          current_markers=Seurat::FindMarkers(object = seurat_object,
                                              ident.1 = cluster_1,
                                              ident.2 = cluster_2,
                                              assay = assay_markers,
                                              logfc.threshold = logfc.threshold,
                                              features = genes_to_include,
                                              slot = assay_slot,
                                              test.use =test.use,
                                              min.pct = min.pct,
                                              min.diff.pct = min.diff.pct,
                                              max.cells.per.ident=max.cells.per.ident,
                                              min.cells.feature = min.cells.feature,
                                              min.cells.group = min.cells.group,
                                              base = base,
                                              only.pos = only.pos)
        }
        current_markers$gene = rownames(current_markers)
        current_markers$cluster_id = cluster_1
        current_markers$comparison = comp
        current_markers$parent = parent_node
        current_markers
      }else{
        NULL
      }
      # },error=function(cond) {
      #   message(cond)
      #   # Choose a return value in case of error
      #   return(NULL)
      # })

      # return
      #current_markers
    }
    return_list[[paste0("comparisons_",comp)]] = comparisons_main
  }
  # return
  return(return_list)
}

##########
### Prepare edgelist and labelmat
##########

# get key objects
edgelist = mrtree_result$edgelist[,1:2]
message("loaded edglist with ",nrow(edgelist)," rows.")
labelmat = mrtree_result$labelmat
message("loaded labelmat with ",nrow(labelmat)," rows.")

message("Object cells: ",ncol(harmonized_seurat_object),"  Labelmat cells: ",nrow(labelmat))

# if start node if specified:
# find subtree and subset edgelist:

# walk through tree:
all_children = scUtils::find_children(nodes = start_node, edges = edgelist[,1:2])
if(length(all_children) < 3){
  message("Cannot find enough children for starting node ",start_node," . Using all nodes in tree instead")
  all_children =  unique(edges[,2])
}

# subset
edgelist = edgelist[edgelist$to %in% all_children,]

# whichc ells to include when downsampleing for feature abundance check
if(start_node != "all"){
  cells_to_check = rownames(harmonized_seurat_object@meta.data)[labelmat[,which(apply(labelmat,2,function(x,target){target %in% x},target=start_node))] == start_node]
}else{
  cells_to_check = rownames(harmonized_seurat_object@meta.data)
}

##########
### Calculate markers between leaf-siblings ?  and merge ?
##########

## subset features to be tested;
if(ncol(harmonized_seurat_object@assays[['RNA']]@data) > 30000){
  set.seed(parameter_list$global_seed)
  subset = sample(colnames(harmonized_seurat_object@assays[['RNA']]@data[,cells_to_check]),size = 30000)
  gene_expr_dataset = harmonized_seurat_object@assays[['RNA']]@data[,subset]
}else{
  gene_expr_dataset = harmonized_seurat_object@assays[['RNA']]@data
}
gene_expr_dataset[gene_expr_dataset != 0] <- 1 # set to 1 for occ
gene_sums = Matrix::rowSums(gene_expr_dataset)
rm(gene_expr_dataset)
genes_to_include = names(gene_sums)[gene_sums > min.cells.feature ]
genes_to_include = genes_to_include[! genes_to_include %in% features_exclude_list]
message("Testing ",length(genes_to_include)," genes as cluster markers")

# run marker detection on tree using findMarkers_tree2
message(Sys.time(),": Start marker detection on mrtree" )
all_markers_mrtree_list = findMarkers_tree2(harmonized_seurat_object,
                                            edgelist=edgelist[,1:2],
                                            labelmat = labelmat,
                                            n_cores = n_cores_markers,
                                            assay=assay_markers,
                                            slot=assay_slot,
                                            batch_var=batch_var,
                                            genes_to_include = genes_to_include,
                                            logfc.threshold = logfc.threshold,
                                            min.pct = min.pct,
                                            min.diff.pct = min.diff.pct,
                                            max.cells.per.ident=max.cells.per.ident,
                                            min.cells.feature = min.cells.feature,
                                            min.cells.group = min.cells.group,
                                            base = base,
                                            only.pos = only.pos)

message("names(all_markers_mrtree_list): ",names(all_markers_mrtree_list))

#
comparisons_siblings = as.data.frame(all_markers_mrtree_list$comparisons_Sibling)
comparisons_all = as.data.frame(all_markers_mrtree_list$comparisons_All)

# save objects as txt#
data.table::fwrite(comparisons_siblings,paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$start_node,"_markers_siblings_",parameter_list$marker_suffix,".tsv"),sep="\t")
data.table::fwrite(comparisons_all,paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$start_node,"_markers_all_",parameter_list$marker_suffix,".tsv"),sep="\t")


## assumes base 2!
comparisons_siblings$specificity = ((comparisons_siblings$pct.1 + specificity_base) / (comparisons_siblings$pct.2 + specificity_base)) * comparisons_siblings$avg_log2FC
comparisons_all$specificity = ((comparisons_all$pct.1 + specificity_base) / (comparisons_all$pct.2 + specificity_base)) * comparisons_all$avg_log2FC


##########
### Save result
##########

# save objects as txt
data.table::fwrite(comparisons_siblings,paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$start_node,"_markers_siblings_",parameter_list$marker_suffix,".tsv"),sep="\t")
data.table::fwrite(comparisons_all,paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$start_node,"_markers_all_",parameter_list$marker_suffix,".tsv"),sep="\t")

message(Sys.time(),": Finalized marker detection." )
