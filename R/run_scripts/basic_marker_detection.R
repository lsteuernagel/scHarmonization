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
harmonized_seurat_object = readRDS(paste0(parameter_list$merged_file))

# load clusters

### get all parameters:
# general
seurat_file_name =parameter_list$seurat_file_name
output_file_name =parameter_list$output_file_name
cluster_column = parameter_list$cluster_column
global_seed = parameter_list$global_seed

# markers
#n_cores_markers = parameter_list$n_cores_markers
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

if(is.null(base)){base = 2}

##########
### Calculate markers
##########

message(Sys.time(),": Start marker detection" )
Idents(seurat_object_harmonized) = cluster_column # set ident!!

if(use_stratified){
  # check batch var
  if(! batch_var %in% colnames(seurat_object_harmonized@meta.data)){stop("Error: Cannot find ",batch_var," in object")}
  message("Using van elteren Test for batch-aware marker detection.")
  # simple for loop
  all_clusters = unique(seurat_object_harmonized@meta.data[,cluster_column])
  all_markers_list = list()
  for(current_cluster in all_clusters){
    message("Calculating cluster ",current_cluster)
    if(length(seurat_object_harmonized@meta.data[seurat_object_harmonized@meta.data[,cluster_column] == current_cluster,cluster_column]) >= min.cells.group){
      # FindAll stratified and use batch_var
      # see stratified_wilcoxon_functions.R
      all_markers_list[[current_cluster]] =FindMarkers2.Seurat(object = seurat_object_harmonized,
                                                               ident.1 = current_cluster,
                                                               assay = assay_markers,
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
  all_clusters = unique(seurat_object_harmonized@meta.data[,cluster_column])
  all_markers_list = list()
  for(current_cluster in all_clusters){
    message("Calculating cluster ",current_cluster)
    # FindAll stratified and use batch_var
    # see stratified_wilcoxon_functions.R
    all_markers_list[[current_cluster]] <- tryCatch({
      temp_df = Seurat::FindMarkers(object = seurat_object_harmonized,
                                    ident.1 = current_cluster,
                                    assay = assay_markers,
                                    logfc.threshold = logfc.threshold,
                                    slot = assay_slot,
                                    test.use = "wilcox",
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
