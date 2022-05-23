##########
### Load params
##########

# load json file with all other information
parameter_list = jsonlite::read_json("data/parameters_harmonization_v2_1.json")
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})


##########
### Load seurat object + clustering + amrker detection results
##########

# load seurat
harmonized_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,".rds"))

# load clusters
hypoMap_test_initial_leiden_clustering = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_initial_leiden_clustering.txt"),data.table = F)
temp_meta = dplyr::left_join(harmonized_seurat_object@meta.data,hypoMap_test_initial_leiden_clustering[,c(1,ncol(hypoMap_test_initial_leiden_clustering))],by="Cell_ID")
rownames(temp_meta) = temp_meta$Cell_ID
harmonized_seurat_object@meta.data = temp_meta
cluster_column = colnames(hypoMap_test_initial_leiden_clustering)[ncol(hypoMap_test_initial_leiden_clustering)]

# markers
#data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,parameter_list$basic_marker_filename,".txt"))


##########
### Visualization
##########

library(Seurat)
harmonized_seurat_object@assays$scvi_corrected@counts = harmonized_seurat_object@assays$scvi_corrected@data
harmonized_seurat_object =NormalizeData(harmonized_seurat_object,assay = "scvi_corrected")

p1 = DimPlot(harmonized_seurat_object,group.by = "Author_Class",raster = F)
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

p1 = DimPlot(harmonized_seurat_object,group.by = "Dataset",raster = F,shuffle = TRUE)
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)


p1 = DimPlot(harmonized_seurat_object,group.by = cluster_column,raster = F,label=TRUE)+NoLegend()
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)
gene="Sim1"
DefaultAssay(harmonized_seurat_object) <- "RNA"
p1=FeaturePlot(harmonized_seurat_object,features = gene,raster = F,order=TRUE)
p1=scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

DefaultAssay(harmonized_seurat_object) <- "scvi_corrected"
p2=FeaturePlot(harmonized_seurat_object,features = gene,raster = F,order=TRUE)
p2=scUtils::rasterize_ggplot(p2,pixel_raster = 2048,pointsize = 1.8)

cowplot::plot_grid(p1,p2)

#####
DefaultAssay(harmonized_seurat_object) <- "RNA"

p1=FeaturePlot(harmonized_seurat_object,features = "Oxt",raster = F,order=TRUE)
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)
p1

p1 = DimPlot(harmonized_seurat_object,group.by = "Dataset",raster = F,shuffle = TRUE)
p1+theme_grey()+NoLegend()


###
harmonized_seurat_object = RunUMAP(harmonized_seurat_object,reduction = "scvi",dims = 1:85,n.neighbors = 20,seed.use = 123456,reduction.name = "umap20",reduction.key = "umap20")

p1 = DimPlot(harmonized_seurat_object,group.by = "Dataset",raster = F,shuffle = TRUE,reduction = "umap20")
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

p1 = DimPlot(harmonized_seurat_object,group.by = cluster_column,raster = F,label=TRUE,reduction = "umap20")+NoLegend()
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

p1 = DimPlot(harmonized_seurat_object,group.by = cluster_column,raster = F,label=TRUE,reduction = "umap_scvi")+NoLegend()
p2 = DimPlot(harmonized_seurat_object,group.by = cluster_column,raster = F,label=TRUE,reduction = "umap20")+NoLegend()
cowplot::plot_grid(p1,p2)

##########
### Curate clusters
##########

hypoMap_seurat_meta = harmonized_seurat_object_test@meta.data
hypoMap_seurat_meta[,cluster_column] = as.character(hypoMap_seurat_meta[,cluster_column])

hypoMap_seurat_meta$Author_Class[hypoMap_seurat_meta$Author_Class=="Vascular"] ="Endothelial"
hypoMap_seurat_meta$Author_Class[hypoMap_seurat_meta$Author_Class %in% c("Unassigned","Mixed")] = NA
class_per_cluster = hypoMap_seurat_meta %>% dplyr::group_by(!!sym(cluster_column)) %>% dplyr::add_count(name="cluster_total") %>%
  dplyr::group_by(!!sym(cluster_column),Author_Class,cluster_total)%>%
  dplyr::count(name="cluster_class_count") %>% dplyr::mutate(pct = round(cluster_class_count / cluster_total *100,4)) %>%
  dplyr::filter(!is.na(Author_Class)) %>% dplyr::group_by(!!sym(cluster_column)) %>% dplyr::mutate(annotated_pct = sum(pct)) %>%
  dplyr::ungroup() %>%  dplyr::filter(pct > 1) %>% as.data.frame()
class_per_cluster[,cluster_column] = as.factor(as.character(class_per_cluster[,cluster_column]))

class_per_cluster_top = class_per_cluster %>% dplyr::group_by(!!sym(cluster_column))  %>% dplyr::filter(pct == max(pct)) %>%
  dplyr::mutate(fraction_of_annotated = pct / annotated_pct) %>% dplyr::distinct(!!sym(cluster_column),.keep_all = TRUE)

class_per_cluster_top$Author_Class_Curated = class_per_cluster_top$Author_Class
class_per_cluster_top$Author_Class_Curated[class_per_cluster_top$fraction_of_annotated < 0.67] = "Unknown"
class_per_cluster_top$Author_Class_Curated[class_per_cluster_top$fraction_of_annotated < 0.8 & class_per_cluster_top$annotated_pct < 10] = "Unknown"

hypoMap_seurat_meta = hypoMap_seurat_meta %>% dplyr::rename(Author_Class_Curated_orig = Author_Class_Curated)
temp_meta = dplyr::left_join(hypoMap_seurat_meta,class_per_cluster_top[,c(cluster_column,"Author_Class_Curated")],by=cluster_column) %>% as.data.frame()
rownames(temp_meta) = temp_meta$Cell_ID
temp_meta$Author_Class_Curated[is.na(temp_meta$Author_Class_Curated)] ="Unknown"
harmonized_seurat_object_test@meta.data = temp_meta

##########
### Curate clusters manually
##########


# TODO

##########
### Subset object
##########

keep_cells = harmonized_seurat_object_test@meta.data$Cell_ID[!harmonized_seurat_object_test@meta.data$Author_Class_Curated %in% c("Doublet","Unknown") ]
curated_seurat_object = subset(harmonized_seurat_object_test,cells = keep_cells)

##########
### Updated curated object NN trees
##########

## also re-run

run umap and save model
message(Sys.time(),": Build UMAP with ",parameter_list$k_param," n.neighbors ..." )
curated_seurat_object = RunUMAP(curated_seurat_object,
                                   reduction = parameter_list$integration_name,
                                   seed.use= parameter_list$global_seed,
                                   dims=1:ncol(curated_seurat_object@reductions[[parameter_list$integration_name]]@cell.embeddings),
                                   reduction.name=paste0("umap_",parameter_list$integration_name),
                                   reduction.key = paste0("umap_",parameter_list$integration_name),
                                   verbose=F,
                                   n.neighbors = parameter_list$k_param,
                                   return.model = TRUE)

# run seurat SNN annoy
message(Sys.time(),": Build SNN with ",parameter_list$k_param," n.neighbors ..." )
curated_seurat_object = FindNeighbors(curated_seurat_object,
                                         reduction=parameter_list$integration_name,
                                         dims = 1:ncol(curated_seurat_object@reductions[[parameter_list$integration_name]]@cell.embeddings),
                                         k.param = parameter_list$k_param,
                                         nn.method="annoy",
                                         annoy.metric=parameter_list$dist_type,
                                         graph.name = paste0("SNN_",parameter_list$integration_name), verbose=TRUE)

##########
### Save
##########

# name
file_name_prefix = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_curated")

# save data to rds
saveRDS(curated_seurat_object,paste0(file_name_prefix,".rds"))


