

##########
### Load seurat object
##########

# load json file with all other information
parameter_list = jsonlite::read_json("data/parameters_harmonization_v2_4.json")
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# load seurat
harmonized_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,"hypoMap_harmonized",".rds"))

##########
### Removed
##########

## make plot
harmonized_seurat_object@meta.data$removed = "yes"
harmonized_seurat_object@meta.data$removed[harmonized_seurat_object@meta.data$Cell_ID %in% curated_seurat_object@meta.data$Cell_ID] = "no"

p1 = DimPlot(harmonized_seurat_object,group.by = "removed",raster = F,label=TRUE)+NoLegend()
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)


p1 = DimPlot(harmonized_seurat_object,group.by = "Author_Class",raster = F,label=TRUE)+NoLegend()
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)


# First round:
hypoMap_harmonized_leiden_clustering_initial = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_harmonization/hypoMap_harmonized_leiden_clustering.txt",data.table = F)
dim(hypoMap_harmonized_leiden_clustering_initial)
dim(harmonized_seurat_object)
harmonized_seurat_object$initial_curation_clustering = hypoMap_harmonized_leiden_clustering_initial$leiden_clusters_13


p1 = DimPlot(harmonized_seurat_object,group.by = "initial_curation_clustering",raster = F,label=TRUE)+NoLegend()
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)


### second round:
hypoMap_harmonized_leiden_clustering_initial = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_harmonization/hypoMap_harmonized_leiden_clustering.txt",data.table = F)
