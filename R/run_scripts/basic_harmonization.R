##########
### Load parameters and packages
##########
message(" Load parameters and packages ")

library(magrittr)
library(scUtils)
library(Seurat)

source("R/harmonization_functions")

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

#get metadata and save
seurat_metadata = harmonized_seurat_object@meta.data
data.table::fwrite(x = seurat_metadata,file = paste0(parameter_list$harmonization_folder_path,paste0(parameter_list$new_name_suffix,"_metadata.txt")),sep = "\t")


##########
### Raw data & new object
##########

# put scvi imputed into another assay
message("Adding counts as new assay from: ",parameter_list$corrected_counts_scvi)
corrected_counts_scvi = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"scVI_corrected.txt"),data.table = F)
rownames(corrected_counts_scvi) = corrected_counts_scvi[,1]
corrected_counts_scvi = t(as.matrix(corrected_counts_scvi[,2:ncol(corrected_counts_scvi)]))
scvi_assay = CreateAssayObject(data=corrected_counts_scvi, min.cells = 0, min.features = 0)
harmonized_seurat_object@assays[[paste0(parameter_list$integration_name,'_corrected')]] = scvi_assay
harmonized_seurat_object@assays[[paste0(parameter_list$integration_name,'scvi_corrected')]]@key =paste0(parameter_list$integration_name,"_corrected")

# # add hvgs
# message("Adding variable feature names to assays")
# harmonized_seurat_object@misc$var_features = as.character(unlist(parameter_list$var_features))
# # harmonized_seurat_object@assays[['RNA']]@var.features=parameter_list$var_features
# # harmonized_seurat_object@assays[[paste0(parameter_list$integration_name,'_corrected')]]@var.features=parameter_list$var_features


##########
### Reductions
##########

message("Adding integrated embedding and calculating default umap")
# add embedding manually
current_embedding = read_embedding(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_scVI_reduction.txt")),harmonized_seurat_object)

# make dim red
dimred <- Seurat::CreateDimReducObject(
  embeddings = as.matrix(current_embedding),
  stdev = as.numeric(apply(current_embedding, 2, stats::sd)),
  assay = "RNA",
  key = parameter_list$integration_name
)
# add  to object
harmonized_seurat_object@reductions[[parameter_list$integration_name]] = dimred


##########
### UMAP
##########

# run umap and save model
message(Sys.time(),": Build UMAP..." )
harmonized_seurat_object = RunUMAP(harmonized_seurat_object,
                                   reduction = parameter_list$integration_name,
                                   seed.use= parameter_list$global_seed,
                                   dims=1:ncol(harmonized_seurat_object@reductions[[parameter_list$integration_name]]@cell.embeddings),
                                   reduction.name=paste0("umap_",parameter_list$integration_name),
                                   reduction.key = paste0("umap_",parameter_list$integration_name),
                                   verbose=F,
                                   n.neighbors = parameter_list$k_param,
                                   return.model = TRUE)

##########
### basic SNN:
##########

# run seurat SNN annoy
message(Sys.time(),": Build SNN.." )
seurat_object_harmonized = FindNeighbors(seurat_object_harmonized,reduction=reduction,dims = 1:ncol(harmonized_seurat_object@reductions[[parameter_list$integration_name]]@cell.embeddings),
                                         k.param = parameter_list$k_param,
                                         nn.method="annoy",annoy.metric=parameter_list$dist_type,
                                         graph.name = paste0("SNN_",parameter_list$integration_name), verbose=TRUE)


##########
### Save object
##########

file_name_prefix = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix)

# save data to rds
saveRDS(harmonized_seurat_object,paste0(file_name_prefix,".rds"))

# save h5seurat
SeuratDisk::SaveH5Seurat(object = seurat_merged,filename = paste0(file_name_prefix,".h5seurat"), overwrite = TRUE, verbose = TRUE)

# save to anndata
SeuratDisk::Convert( paste0(file_name_prefix,".h5seurat"), dest =  paste0(file_name_prefix,".h5ad"),assay="RNA",verbose=TRUE,overwrite=TRUE)
system(paste0("rm ",paste0(file_name_prefix,".h5seurat")))

# metadata
#data.table::fwrite(harmonized_seurat_object@meta.data,file = paste0(parameter_list$harmonization_path,parameter_list$project_name,"_metadata.txt") ,sep="\t")

message("Finalized basic seurat object harmonization.")
