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

# load clusters
hypoMap_test_initial_leiden_clustering = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_initial_leiden_clustering.txt"),data.table = F)
temp_meta = dplyr::left_join(harmonized_seurat_object@meta.data,hypoMap_test_initial_leiden_clustering[,c(1,ncol(hypoMap_test_initial_leiden_clustering))],by="Cell_ID")
rownames(temp_meta) = temp_meta$Cell_ID
harmonized_seurat_object@meta.data = temp_meta
cluster_column = colnames(hypoMap_test_initial_leiden_clustering)[ncol(hypoMap_test_initial_leiden_clustering)]

# define output file name
output_file_name = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,parameter_list$basic_marker_filename,".txt")

# markers
assay_markers = parameter_list$assay_markers
assay_slot =parameter_list$assay_slot
test.use = parameter_list$test.use
logfc.threshold = parameter_list$logfc.threshold
min.pct =parameter_list$min.pct
min.diff.pct = parameter_list$min.diff.pct
max.cells.per.ident = parameter_list$max.cells.per.ident
min.cells.feature = parameter_list$min.cells.feature
min.cells.group =  parameter_list$min.cells.group
base = parameter_list$base
only.pos = parameter_list$only.pos
use_stratified=parameter_list$use_stratified
batch_var = parameter_list$batch_var
test.use = parameter_list$test.use
message("test.use ",test.use," ... ",parameter_list$test.use)
if(is.null(base)){base = 2}



mrtree_object = readRDS(mrtree_object_filepath)

edgelist = mrtree_object$edgelist[,1:2]
labelmat = mrtree_object$labelmat

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
#' @param ...
#' @return updated vector of labels

findMarkers_tree2 = function(seurat_object,edgelist,labelmat,n_cores=1,use_stratified=TRUE,test.use="",batch_var="Batch_ID",assay="RNA",slot="data",...){
  require(doParallel)
  edgelist = edgelist[,c("from","to")]
  message(n_cores)

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
  message("Sibling Comparisons")
  comparisons_siblings <- foreach(n = 1:length(all_nodes),.errorhandling = 'remove', .combine='rbind') %dopar% {
    current_node = all_nodes[n]
    parent_node = edgelist$from[edgelist$to==current_node]
    sibling_nodes = edgelist$to[edgelist$from==parent_node & edgelist$to != current_node]
    current_level = edgelist$clusterlevel[edgelist$to==current_node]
    # set ident to current level!
    Idents(seurat_object) = current_level

    # calculate markers vs siblings
    current_markers <- tryCatch({
      cluster_1 = current_node
      cluster_2 = sibling_nodes
      if(test.use == "wilcox-stratified"){
        current_markers =FindMarkers2.Seurat(seurat_object,
                                             ident.1 = cluster_1,
                                             ident.2 = cluster_2,
                                             latent.vars = batch_var,
                                             test.use = "VE",
                                             genre = "locally-best",
                                             assay =assay,
                                             slot = slot,...)
        #current_markers =FindMarkers2.Seurat(seurat_object, ident.1 = cluster_1,ident.2 = cluster_2, latent.vars = batch_var, test.use = "VE", genre = "locally-best",only.pos=TRUE,logfc.threshold=0.2,min.pct=0.1,min.diff.pct=0.1)
      }else{
        current_markers=FindMarkers(seurat_object,
                                    ident.1 = cluster_1,
                                    ident.2 = cluster_2,
                                    assay =assay,
                                    slot = slot,
                                    test.use = test.use,...)
      }
      current_markers$gene = rownames(current_markers)
      current_markers$cluster_1 = cluster_1
      current_markers$comparison = "siblings"
      current_markers$parent = parent_node
      current_markers
    },error=function(cond) {
      message(cond)
      # Choose a return value in case of error
      return(NULL)
    })

    # return
    current_markers
  }
  message("All Comparisons")
  comparisons_all <- foreach(n = 1:length(all_nodes),.errorhandling = 'remove', .combine='rbind') %dopar% {
    #for(n in 1:length(all_nodes)){
    current_node = all_nodes[n]
    parent_node = edgelist$from[edgelist$to==current_node]
    sibling_nodes = edgelist$to[edgelist$from==parent_node & edgelist$to != current_node]
    current_level = edgelist$clusterlevel[edgelist$to==current_node]
    # set ident to current level!
    Idents(seurat_object) = current_level

    # calculate markers vs all
    current_markers <- tryCatch({
      cluster_1 = current_node
      cluster_2 = all_nodes[! all_nodes %in% current_node]
      #markers vs all
      if(test.use == "wilcox-stratified"){
        current_markers =FindMarkers2.Seurat(seurat_object,
                                             ident.1 = cluster_1,
                                             ident.2 = cluster_2,
                                             latent.vars = batch_var,
                                             test.use = "VE",
                                             genre = "locally-best",
                                             assay =assay,
                                             slot = slot,...)
      }else{
        current_markers=FindMarkers(seurat_object,
                                    ident.1 = cluster_1,
                                    ident.2 = cluster_2,
                                    assay =assay,
                                    slot = slot,
                                    test.use = test.use,...)
      }
      current_markers$gene = rownames(current_markers)
      current_markers$cluster_1 = cluster_1
      current_markers$comparison = "all"
      current_markers$parent = parent_node
      current_markers
    },error=function(cond) {
      message(cond)
      # Choose a return value in case of error
      return(NULL)
    })

    # return
    current_markers
  }
  # return
  return_list = list(comparisons_siblings=comparisons_siblings,comparisons_all=comparisons_all)
  return(return_list)
}



##########
### Calculate markers between leaf-siblings ?  and merge ?
##########


# TODO: needs to be able to only work with a subset of markers

message(Sys.time(),": Start marker detection on mrtree" )
all_markers_stratified = findMarkers_tree2(seurat_object_harmonized,edgelist=edgelist[,1:2],labelmat = labelmat,n_cores = n_cores_markers,
                                           assay=assay_markers,slot=assay_slot,use_stratified=use_stratified,batch_var=batch_var,max.cells.per.ident = max.cells.per.ident,
                                           only.pos=TRUE,logfc.threshold=logfc.threshold,min.pct=min.pct,min.diff.pct=min.diff.pct)

comparisons_siblings = all_markers_stratified$comparisons_siblings
comparisons_siblings$specificity = (comparisons_siblings$pct.1 / comparisons_siblings$pct.2) * comparisons_siblings$avg_logFC
comparisons_all = all_markers_stratified$comparisons_all
comparisons_all$specificity = (comparisons_all$pct.1 / comparisons_all$pct.2) * comparisons_all$avg_logFC


##########
### Calculate markers
##########


message(Sys.time(),": Start marker detection" )
Idents(harmonized_seurat_object) = cluster_column # set ident!!

# genes to test
#seurat_raw@assays[["RNA"]]@counts
gene_expr_dataset = harmonized_seurat_object@assays[['RNA']]@counts
gene_expr_dataset[gene_expr_dataset!=0] = 1
gene_sums = Matrix::rowSums(gene_expr_dataset)
genes_to_include = names(gene_sums)[gene_sums > min.cells.feature ]
genes_to_include = genes_to_include[! genes_to_include %in% features_exclude_list]
message("Testing ",length(genes_to_include)," genes as cluster markers")

if(test.use=="wilcox-stratified"){
  # check batch var
  if(! batch_var %in% colnames(harmonized_seurat_object@meta.data)){stop("Error: Cannot find ",batch_var," in object")}
  message("Using van elteren Test for batch-aware marker detection.")
  # simple for loop
  all_clusters = unique(harmonized_seurat_object@meta.data[,cluster_column])
  all_markers_list = list()
  for(current_cluster in all_clusters){
    message("Calculating cluster ",current_cluster)
    if(length(harmonized_seurat_object@meta.data[harmonized_seurat_object@meta.data[,cluster_column] == current_cluster,cluster_column]) >= min.cells.group){
      # FindAll stratified and use batch_var
      # see stratified_wilcoxon_functions.R
      all_markers_list[[current_cluster]] =FindMarkers2.Seurat(object = harmonized_seurat_object,
                                                               ident.1 = current_cluster,
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
      all_markers_list[[current_cluster]]$cluster = current_cluster
      all_markers_list[[current_cluster]]$gene = rownames(all_markers_list[[current_cluster]])
    }
  }
  all_markers = do.call(rbind, all_markers_list)
}else{
  message("Using Seurat::FindMarkers")
  all_clusters = unique(harmonized_seurat_object@meta.data[,cluster_column])
  all_markers_list = list()
  for(current_cluster in all_clusters){
    message("Calculating cluster ",current_cluster)
    # FindAll stratified and use batch_var
    # see stratified_wilcoxon_functions.R
    all_markers_list[[current_cluster]] <- tryCatch({
      temp_df = Seurat::FindMarkers(object = harmonized_seurat_object,
                                    ident.1 = current_cluster,
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
      temp_df$cluster = current_cluster
      temp_df$gene = rownames(temp_df)
      temp_df
    },
    error=function(cond) {
      message("Cannot calculate markers. Error:",cond)
      return(NULL)
    })
  }
  all_markers = do.call(rbind, all_markers_list)
}

#all_markers$specificity = (all_markers$pct.1 / all_markers$pct.2) * all_markers$avg_logFC

##########
### Save result
##########

# save objects as txt
data.table::fwrite(all_markers,file=output_file_name)

message(Sys.time(),": Finalized marker detection." )
