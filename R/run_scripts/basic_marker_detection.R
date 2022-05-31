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
hypoMap_test_initial_leiden_clustering = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_leiden_clustering.txt"),data.table = F)
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

##########
### Calculate markers
##########

message(Sys.time(),": Start marker detection" )
Idents(harmonized_seurat_object) = cluster_column # set ident!!

# genes to test
if(ncol(harmonized_seurat_object@assays[['RNA']]@data) > 30000){
  set.seed(parameter_list$global_seed)
  subset = sample(colnames(harmonized_seurat_object@assays[['RNA']]@data),size = 30000)
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
counter=1
if(test.use=="wilcox-stratified"){
  # check batch var
  if(! batch_var %in% colnames(harmonized_seurat_object@meta.data)){stop("Error: Cannot find ",batch_var," in object")}
  message("Using van elteren Test for batch-aware marker detection.")
  # simple for loop
  all_clusters = unique(harmonized_seurat_object@meta.data[,cluster_column])
  all_markers_list = list()
  for(current_cluster in all_clusters){
    message("Calculating cluster ",current_cluster," (",counter,"/",length(all_clusters),")")
    counter = counter+1
    if(length(harmonized_seurat_object@meta.data[harmonized_seurat_object@meta.data[,cluster_column] == current_cluster,cluster_column]) >= min.cells.group){
      # FindAll stratified and use batch_var
      # see stratified_wilcoxon_functions.R
      all_markers_list[[paste0("c_",as.character(current_cluster))]] =FindMarkers2.Seurat(object = harmonized_seurat_object,
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
      all_markers_list[[paste0("c_",current_cluster)]]$cluster = current_cluster
      all_markers_list[[paste0("c_",current_cluster)]]$gene = rownames(all_markers_list[[paste0("c_",current_cluster)]])
    }
  }
  all_markers = do.call(rbind, all_markers_list)
}else{
  message("Using Seurat::FindMarkers")
  all_clusters = unique(harmonized_seurat_object@meta.data[,cluster_column])
  all_markers_list = list()
  for(current_cluster in all_clusters){
    message("Calculating cluster ",current_cluster," (",counter,"/",length(all_clusters),")")
    counter = counter+1
    all_markers_list[[paste0("c_",as.character(current_cluster))]] <- tryCatch({
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

all_markers$cluster = stringr::str_remove(all_markers$cluster,pattern = "c_")

##########
### Save result
##########

# save objects as txt
data.table::fwrite(all_markers,file=output_file_name)

message(Sys.time(),": Finalized marker detection." )
