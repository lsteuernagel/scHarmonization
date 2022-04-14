##########
### Load parameters and packages
##########
message(" Load parameters and packages ")

library(magrittr)
library(scUtils)

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
seurat_merged = readRDS(paste0(parameter_list$merged_file))

#get metadata and save
seurat_metadata = seurat_merged@meta.data
data.table::fwrite(x = seurat_metadata,file = paste0(parameter_list$harmonization_folder_path,paste0(parameter_list$new_name_suffix,"_metadata.txt")),sep = "\t")


##########
### Raw data & new object
##########

# load raw data
message("Loading Seurat object from: ",parameter_list$raw_seurat_object_path)
harmonized_seurat_object <- SeuratDisk::LoadH5Seurat(parameter_list$raw_seurat_object_path,assays='RNA')
harmonized_seurat_object@project.name = parameter_list$project_name

# put scvi imputed into another assay
message("Adding counts as new assay from: ",parameter_list$corrected_counts_scvi)
corrected_counts_scvi = data.table::fread(parameter_list$corrected_counts_scvi,data.table = F)
rownames(corrected_counts_scvi) = corrected_counts_scvi[,1]
corrected_counts_scvi = t(as.matrix(corrected_counts_scvi[,2:ncol(corrected_counts_scvi)]))
scvi_assay = CreateAssayObject(data=corrected_counts_scvi, min.cells = 0, min.features = 0)
harmonized_seurat_object@assays[[paste0(parameter_list$integration_name,'_corrected')]] = scvi_assay
harmonized_seurat_object@assays[[paste0(parameter_list$integration_name,'_corrected')]]@key ="scvi_corrected"

# add hvgs
message("Adding variable feature names to assays")
harmonized_seurat_object@misc$var_features = as.character(unlist(parameter_list$var_features))
# harmonized_seurat_object@assays[['RNA']]@var.features=parameter_list$var_features
# harmonized_seurat_object@assays[[paste0(parameter_list$integration_name,'_corrected')]]@var.features=parameter_list$var_features


##########
### Reductions
##########

message("Adding integrated embedding and calculating default umap")
# add embedding manually
current_embedding = read_embedding(paste0(parameter_list$reduction_path),harmonized_seurat_object)

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
harmonized_seurat_object = RunUMAP(harmonized_seurat_object,
                                   reduction = parameter_list$integration_name,
                                   seed.use=global_seed,
                                   dims=1:ncol(harmonized_seurat_object@reductions[[parameter_list$integration_name]]@cell.embeddings),
                                   reduction.name=paste0("umap_",parameter_list$integration_name),
                                   reduction.key = paste0("umap_",parameter_list$integration_name),
                                   verbose=F,
                                   n.neighbors = parameter_list$k_param_umap,
                                   return.model = TRUE)
##########
### Add other results
##########

if(parameter_list$add_mappedCelltypes){
  message("Adding celltype matching results.")
  cellType_matching_res=readRDS(parameter_list$cellType_matching_path)
  tmp_meta=harmonized_seurat_object@meta.data
  #tmp_meta=tmp_meta[,which(!colnames(tmp_meta) %in% c("Curated_CellType","Curated_CellType_Confidence","Curated_CellType_HighConf"))]
  celltype_per_dataset = cellType_matching_res$celltype_per_dataset
  colnames(celltype_per_dataset) = paste0("predicted_",colnames(celltype_per_dataset))
  celltype_per_dataset$Cell_ID =rownames(celltype_per_dataset)
  tmp_meta = dplyr::left_join(tmp_meta,celltype_per_dataset,by="Cell_ID")
  rownames(tmp_meta) = tmp_meta$Cell_ID
  harmonized_seurat_object@meta.data = tmp_meta

}

##########
### Save object
##########

message("Saving harmonized Seurat object.")

h5Seurat_filename = save_h5Seurat(seurat_object=harmonized_seurat_object,
                                  filepath=parameter_list$harmonization_path,
                                  suffix=parameter_list$project_name,
                                  verbose=F,name_with_project = F)

# metadata
data.table::fwrite(harmonized_seurat_object@meta.data,file = paste0(parameter_list$harmonization_path,parameter_list$project_name,"_metadata.txt") ,sep="\t")

message("Finalized Seurat object harmonization.")
