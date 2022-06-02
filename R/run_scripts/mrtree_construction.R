##########
### Load parameters and packages
##########

message(Sys.time(),": Starting mtrtree construction .." )

message(" Load parameters and packages ")

require(tidyverse)
require(Seurat)
require(Matrix)

source("R/harmonization_functions.R")
source("R/mrtree_functions.R")

# get params-filename from commandline
command_args<-commandArgs(TRUE)
param_file = command_args[1]
# read all parameters and filepaths
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# load seurat
#harmonized_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_curated",".rds"))

# load clusters
cluster_matrix_for_mrtree = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$clusters_for_mrtree_file),data.table = F)
# ensure that cells are in right order!
# TODO


##########
### Run mrtree
##########

prefix=paste0("leiden_clusters_")

message(Sys.time(),": Prepare clustering" )

# ensure that we have good numerics
cluster_matrix_for_mrtree = as.matrix(cluster_matrix_for_mrtree)
cluster_matrix_for_mrtree <- apply(cluster_matrix_for_mrtree, 2, as.character)
cluster_matrix_for_mrtree <- apply(cluster_matrix_for_mrtree, 2, as.numeric)
#print(head(cluster_matrix_for_mrtree))

message(Sys.time(),": Build mrtree" )
# feed into mrtree
# function from 'mrtree_functions.R' which is copied from original repo
mrtree_res <- mrtree(cluster_matrix_for_mrtree,
                     prefix = prefix,
                     suffix = NULL,
                     max.k = Inf,
                     consensus = FALSE,
                     sample.weighted = FALSE, # parameter_list$mr_tree_weighted,
                     augment.path = FALSE,
                     verbose = FALSE,
                     n.cores = parameter_list$n_cores)

saveRDS(mrtree_res,paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_curated","_mrtree_res_raw",".rds"))

message(Sys.time(),": Create mrtree ouput list" )

# make labelmat with colnames included
labelmat=mrtree_res$labelmat.mrtree
n=nrow(labelmat)
backup = colnames(labelmat)
labelmat = matrix(paste(matrix(rep(colnames(labelmat),each=n), nrow = n), labelmat, sep='-'), nrow = n)
colnames(labelmat)=backup
df = as.data.frame(unique(labelmat), stringsAsFactors = F)

# add to seuratobject metadata
# seurat_object_harmonized@meta.data = cbind(seurat_object_harmonized@meta.data,labelmat)

# save in data.tree format
require(data.tree)
df$pathString = apply(df, 1, function(x) paste(c('all', x), collapse='/'))
tree.datatree = data.tree::as.Node(df)
# export edgelist
edges= data.tree::ToDataFrameNetwork(tree.datatree,"isLeaf","level","count","totalCount","height")
nodes = data.frame(id = c("all",as.character(unique(edges$to))),label=c("all",as.character(unique(edges$to))))
nodes = rbind(c("all",FALSE,1,5,max(edges$height)+1),edges[,2:ncol(edges)]) %>% dplyr::rename(id = to) %>% dplyr::mutate(label = id)

# make cluster object ? and add to misc
cluster_object = list(labelmat = labelmat,
                      edgelist = edges ,
                      nodelist = nodes,
                      data_tree = tree.datatree,
                      mrtree_output = mrtree_res)
#seurat_object_harmonized@misc$mrtree_clustering_results = cluster_object

##########
###Save
##########

# save object as rds
saveRDS(cluster_object,paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_mrtree_clustering_results",".rds"))

message("Finalized")
