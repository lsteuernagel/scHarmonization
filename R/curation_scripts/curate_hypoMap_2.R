##########
### Load params
##########

# load json file with all other information
parameter_list = jsonlite::read_json("data/parameters_harmonization_v2_7.json")
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})


##########
### Load seurat object + clustering + amrker detection results
##########

# load seurat
harmonized_seurat_object2 = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,".rds"))

# load clusters
hypoMap_test_initial_leiden_clustering = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_leiden_clustering.txt"),data.table = F)
temp_meta = dplyr::left_join(harmonized_seurat_object2@meta.data,hypoMap_test_initial_leiden_clustering[,c(1,ncol(hypoMap_test_initial_leiden_clustering))],by="Cell_ID")
rownames(temp_meta) = temp_meta$Cell_ID
harmonized_seurat_object2@meta.data = temp_meta
cluster_column = colnames(hypoMap_test_initial_leiden_clustering)[ncol(hypoMap_test_initial_leiden_clustering)]

# markers
basic_markers = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,parameter_list$basic_marker_filename,".txt"))
colnames(basic_markers)[colnames(basic_markers)=="avg_logFC"] = "avg_log2FC"
basic_markers$specificity = ((basic_markers$pct.1 + specificity_base) / (basic_markers$pct.2 + parameter_list$specificity_base)) * basic_markers$avg_log2FC
basic_markers = basic_markers %>% dplyr::arrange(desc(specificity))
# basic filtering
basic_markers = basic_markers[basic_markers$specificity > 0.5 & basic_markers$p_val_adj < 0.001,]
basic_markers = basic_markers[,c("gene","cluster","specificity","avg_log2FC","pct.1","pct.2","p_val_adj")]

##########
### Visualization
##########

library(Seurat)
#harmonized_seurat_object2@assays$scvi_corrected@counts = harmonized_seurat_object2@assays$scvi_corrected@data
#harmonized_seurat_object2 =NormalizeData(harmonized_seurat_object2,assay = "scvi_corrected")

p1 = DimPlot(harmonized_seurat_object2,group.by = "Author_Class",raster = F)
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

p1 = DimPlot(harmonized_seurat_object2,group.by = "Dataset",raster = F,shuffle = TRUE)
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)


p1 = DimPlot(harmonized_seurat_object2,group.by = cluster_column,raster = F,label=TRUE)+NoLegend()
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

gene="Lef1"
DefaultAssay(harmonized_seurat_object2) <- "RNA"
p1=FeaturePlot(harmonized_seurat_object2,features = gene,raster = F,order=TRUE)
p1=scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)
p1

# DefaultAssay(harmonized_seurat_object2) <- "scvi_corrected"
# p2=FeaturePlot(harmonized_seurat_object2,features = gene,raster = F,order=TRUE)
# p2=scUtils::rasterize_ggplot(p2,pixel_raster = 2048,pointsize = 1.8)
#
# cowplot::plot_grid(p1,p2)

#####
DefaultAssay(harmonized_seurat_object2) <- "RNA"

p1=FeaturePlot(harmonized_seurat_object2,features = "Sst",raster = F,order=TRUE)
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)
#p1

p1 = DimPlot(harmonized_seurat_object2,group.by = "Dataset",raster = F,shuffle = TRUE)
p1+theme_grey()+NoLegend()


###
harmonized_seurat_object2 = RunUMAP(harmonized_seurat_object2,reduction = "scvi",dims = 1:85,n.neighbors = 20,seed.use = 123456,reduction.name = "umap20",reduction.key = "umap20")

p1 = DimPlot(harmonized_seurat_object2,group.by = "Dataset",raster = F,shuffle = TRUE,reduction = "umap20")
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

p1 = DimPlot(harmonized_seurat_object2,group.by = cluster_column,raster = F,label=TRUE,reduction = "umap20")+NoLegend()
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

p1 = DimPlot(harmonized_seurat_object2,group.by = cluster_column,raster = F,label=TRUE,reduction = "umap_scvi")+NoLegend()
p2 = DimPlot(harmonized_seurat_object2,group.by = cluster_column,raster = F,label=TRUE,reduction = "umap20")+NoLegend()
cowplot::plot_grid(p1,p2)

##########
### Load other version
##########

hypoMap_v2_final = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_final/hypoMap_v2.rds")

# get relevant columns
hypoMap_v2_final_meta = hypoMap_v2_final@meta.data %>% dplyr::select(Cell_ID,C2,C7,C23,C62,C180,C280,C446,Author_Class_Curated)
temp = dplyr::left_join(harmonized_seurat_object2@meta.data,hypoMap_v2_final_meta,by="Cell_ID")
rownames(temp) = temp$Cell_ID
# left join to to new
harmonized_seurat_object2@meta.data = temp

# check result
p1 = DimPlot(harmonized_seurat_object2,group.by = "C23",raster = F,label=TRUE)+NoLegend()
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

# subset to all cell with Author Cll

p1 = DimPlot(hypoMap_v2_final,group.by = "C23",raster = F,label=TRUE)+NoLegend()
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

##########
### class_curated_per_cluster
##########
cluster_column="leiden_clusters_13.y"

class_curated_per_cluster = harmonized_seurat_object2@meta.data %>% dplyr::group_by(!!sym(cluster_column)) %>% dplyr::add_count(name="cluster_total") %>%
  dplyr::group_by(!!sym(cluster_column),Author_Class_Curated,cluster_total)%>%
  dplyr::count(name="cluster_class_count") %>% dplyr::mutate(pct = round(cluster_class_count / cluster_total *100,4)) %>%
  dplyr::filter(!is.na(Author_Class_Curated)) %>% dplyr::group_by(!!sym(cluster_column)) %>% dplyr::mutate(annotated_pct = sum(pct)) %>%
  dplyr::ungroup() %>%  dplyr::filter(pct > 1) %>% as.data.frame()
#class_curated_per_cluster[,cluster_column] = as.factor(as.character(class_curated_per_cluster[,cluster_column]))

updated_version = class_curated_per_cluster  %>% dplyr::group_by(!!sym(cluster_column)) %>% dplyr::filter(pct == max(pct)) %>%
  dplyr::select(cluster_number = leiden_clusters_13.y, Curated_Class_updated = Author_Class_Curated)

# add updated Curated_Class_updated version
temp = dplyr::left_join(harmonized_seurat_object2@meta.data,updated_version,by=c("leiden_clusters_13.y"="cluster_number"))
rownames(temp) = temp$Cell_ID
harmonized_seurat_object2@meta.data = temp

# remove three problematic clusterss:
remove_cells = harmonized_seurat_object2@meta.data$Cell_ID[harmonized_seurat_object2@meta.data$leiden_clusters_13.y %in% c("196","185","208")]

# remove cells:
remove_cells = c(remove_cells,harmonized_seurat_object2@meta.data$Cell_ID[is.na(harmonized_seurat_object2@meta.data$Author_Class_Curated)])

# delete neurons that are now different
remove_cells = c(remove_cells,harmonized_seurat_object2@meta.data$Cell_ID[harmonized_seurat_object2@meta.data$Author_Class_Curated == "Neurons" &
                                                                            harmonized_seurat_object2@meta.data$Curated_Class_updated %in% c("Astrocytes","Immune","NG/OPC","Oligodendrocytes")])

# SUBSET
remove_cells = unique(remove_cells)


##########
### Subset object
##########

keep_cells = harmonized_seurat_object2@meta.data$Cell_ID[!harmonized_seurat_object2@meta.data$Cell_ID %in% remove_cells]
nrow(harmonized_seurat_object2@meta.data) - length(keep_cells) # remove 677 cells
# subset
curated_seurat_object = subset(harmonized_seurat_object2,cells = keep_cells)

# update curated column
curated_seurat_object@meta.data$Author_Class_Curated = curated_seurat_object@meta.data$Curated_Class_updated
select_meta = curated_seurat_object@meta.data %>% dplyr::select(- leiden_clusters_13.x,-preliminary_clusters,-Curated_Class_updated) %>% dplyr::rename(preliminary_clusters = leiden_clusters_13.y)
rownames(select_meta) = select_meta$Cell_ID
curated_seurat_object@meta.data = select_meta

# manual adjst
curated_seurat_object@meta.data$Author_Class_Curated[curated_seurat_object@meta.data$preliminary_clusters == "200"] = "Differentiating"
curated_seurat_object@meta.data$Author_Class_Curated[curated_seurat_object@meta.data$Author_Class_Curated == "Differentiating" ] = "Dividing"

p1 = DimPlot(curated_seurat_object,group.by = "Author_Class_Curated",raster = F)
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)


##########
### Updated curated object NN trees
##########

## clean up metadata:
# colnames(curated_seurat_object@meta.data)
# temp = curated_seurat_object@meta.data  %>% dplyr::select(-seurat_clusters,-log_nCount_RNA,-leiden_clusters_0,-leiden_clusters_0.01,-leiden_clusters_0.05,-leiden_clusters_0.25,-leiden_clusters_1, -leiden_clusters_5, -leiden_clusters_27) %>%
#   dplyr::rename()
# rownames(temp) = temp$Cell_ID
# curated_seurat_object@meta.data = temp

## also re-run
parameter_list$k_param = 25

#run umap and save model
message(Sys.time(),": Build UMAP with ",parameter_list$k_param," n.neighbors ..." )
curated_seurat_object = RunUMAP(curated_seurat_object,
                                reduction = parameter_list$integration_name,
                                seed.use= parameter_list$global_seed,
                                dims=1:ncol(curated_seurat_object@reductions[[parameter_list$integration_name]]@cell.embeddings),
                                reduction.name=paste0("umap_",parameter_list$integration_name),
                                reduction.key = paste0("umap_",parameter_list$integration_name),
                                verbose=F,
                                n.neighbors = parameter_list$k_param, # parameter_list$k_param,
                                return.model = TRUE)

# run seurat SNN annoy
message(Sys.time(),": Build SNN with ",parameter_list$k_param," n.neighbors ..." )
curated_seurat_object = FindNeighbors(curated_seurat_object,
                                      reduction=parameter_list$integration_name,
                                      dims = 1:ncol(curated_seurat_object@reductions[[parameter_list$integration_name]]@cell.embeddings),
                                      k.param = parameter_list$k_param,
                                      nn.method="annoy",
                                      annoy.metric=parameter_list$dist_type,
                                      graph.name = c(paste0("NN_",parameter_list$integration_name),paste0("SNN_",parameter_list$integration_name)),
                                      verbose=TRUE)

p1 = DimPlot(curated_seurat_object,group.by = "Author_Class_Curated",raster = F)
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

p1=FeaturePlot(curated_seurat_object,features = "Crabp1",raster = F,order=TRUE)
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

##########
### Save
##########

select_meta = curated_seurat_object@meta.data %>% dplyr::select(- Final_Doublet, -Final_Exclude, -C2, -C7, -C23, -C62, -C180, -C280, -C446) #%>% dplyr::rename(preliminary_clusters = leiden_clusters_13.y)
rownames(select_meta) = select_meta$Cell_ID
curated_seurat_object@meta.data = select_meta

# name
file_name_prefix = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_curated")

# remove scvi corrected to make things smaller and faster
curated_seurat_object[['scvi_corrected']] <- NULL
# remove scale data
dummy=matrix(data = as.numeric())
curated_seurat_object@assays[["RNA"]]@scale.data <- dummy[,-1] # error is okay

# save data to rds
saveRDS(curated_seurat_object,paste0(file_name_prefix,".rds"))

# remove NN to save SNN when converting to anndata!
curated_seurat_object@graphs[[paste0("NN_",parameter_list$integration_name)]] =NULL

# save h5seurat
SeuratDisk::SaveH5Seurat(object = curated_seurat_object,filename = paste0(file_name_prefix,".h5seurat"), overwrite = TRUE, verbose = TRUE)

# save to anndata
SeuratDisk::Convert( paste0(file_name_prefix,".h5seurat"), dest =  paste0(file_name_prefix,".h5ad"),assay="RNA",verbose=TRUE,overwrite=TRUE)
system(paste0("rm ",paste0(file_name_prefix,".h5seurat")))


