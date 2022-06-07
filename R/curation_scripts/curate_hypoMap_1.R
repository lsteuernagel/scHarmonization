##########
### Load params
##########

# load json file with all other information
parameter_list = jsonlite::read_json("data/parameters_harmonization_v2_4.json")
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})


##########
### Load seurat object + clustering + amrker detection results
##########

# load seurat
harmonized_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,".rds"))

# load clusters
hypoMap_test_initial_leiden_clustering = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_leiden_clustering.txt"),data.table = F)
temp_meta = dplyr::left_join(harmonized_seurat_object@meta.data,hypoMap_test_initial_leiden_clustering[,c(1,ncol(hypoMap_test_initial_leiden_clustering))],by="Cell_ID")
rownames(temp_meta) = temp_meta$Cell_ID
harmonized_seurat_object@meta.data = temp_meta
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
#harmonized_seurat_object@assays$scvi_corrected@counts = harmonized_seurat_object@assays$scvi_corrected@data
#harmonized_seurat_object =NormalizeData(harmonized_seurat_object,assay = "scvi_corrected")

p1 = DimPlot(harmonized_seurat_object,group.by = "Author_Class",raster = F)
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

p1 = DimPlot(harmonized_seurat_object,group.by = "Dataset",raster = F,shuffle = TRUE)
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)


p1 = DimPlot(harmonized_seurat_object,group.by = cluster_column,raster = F,label=TRUE)+NoLegend()
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

gene="Pln"
DefaultAssay(harmonized_seurat_object) <- "RNA"
p1=FeaturePlot(harmonized_seurat_object,features = gene,raster = F,order=TRUE)
p1=scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)
p1

# DefaultAssay(harmonized_seurat_object) <- "scvi_corrected"
# p2=FeaturePlot(harmonized_seurat_object,features = gene,raster = F,order=TRUE)
# p2=scUtils::rasterize_ggplot(p2,pixel_raster = 2048,pointsize = 1.8)
#
# cowplot::plot_grid(p1,p2)

#####
DefaultAssay(harmonized_seurat_object) <- "RNA"

p1=FeaturePlot(harmonized_seurat_object,features = "Rab7b",raster = F,order=TRUE)
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)
#p1

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

hypoMap_seurat_meta = harmonized_seurat_object@meta.data
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

length(unique(as.character(class_per_cluster$leiden_clusters_13)))
length(unique(as.character(class_per_cluster_top$leiden_clusters_13)))

#hypoMap_seurat_meta = hypoMap_seurat_meta %>% dplyr::rename(Author_Class_Curated_orig = Author_Class_Curated)
temp_meta = dplyr::left_join(hypoMap_seurat_meta,class_per_cluster_top[,c(cluster_column,"Author_Class_Curated")],by=cluster_column) %>% as.data.frame()
rownames(temp_meta) = temp_meta$Cell_ID
temp_meta$Author_Class_Curated[is.na(temp_meta$Author_Class_Curated)] ="Unknown"
harmonized_seurat_object@meta.data = temp_meta


##########
### load preliminaryinfo
##########

# hypoMap_v2_preliminary_clusters_090522 = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/hypoMap_v2_preliminary_clusters_090522.tsv",data.table = F)
# hypoMap_v2_preliminary_clusters_090522 = dplyr::left_join(hypoMap_v2_preliminary_clusters_090522,hypoMap_seurat_meta[,c("Cell_ID",cluster_column)] ,by="Cell_ID")
#
# class_curated_per_cluster = hypoMap_v2_preliminary_clusters_090522 %>% dplyr::group_by(!!sym(cluster_column)) %>% dplyr::add_count(name="cluster_total") %>%
#   dplyr::group_by(!!sym(cluster_column),Author_Class_Curated,cluster_total)%>%
#   dplyr::count(name="cluster_class_count") %>% dplyr::mutate(pct = round(cluster_class_count / cluster_total *100,4)) %>%
#   dplyr::filter(!is.na(Author_Class_Curated)) %>% dplyr::group_by(!!sym(cluster_column)) %>% dplyr::mutate(annotated_pct = sum(pct)) %>%
#   dplyr::ungroup() %>%  dplyr::filter(pct > 1) %>% as.data.frame()
# class_curated_per_cluster[,cluster_column] = as.factor(as.character(class_curated_per_cluster[,cluster_column]))


##########
### Curate clusters manually
##########


# show problematic clusters
class_per_cluster_top[class_per_cluster_top$Author_Class_Curated=="Unknown",]
class_curated_per_cluster[class_curated_per_cluster$leiden_clusters_13 %in% class_per_cluster_top$leiden_clusters_13[class_per_cluster_top$Author_Class_Curated=="Unknown"],]

## helper code toc check on clusters:
cluster_target = "194"
class_per_cluster[class_per_cluster$leiden_clusters_13 == cluster_target,]
head(basic_markers[basic_markers$cluster == cluster_target,],10)
tmp_markers = basic_markers[basic_markers$cluster == cluster_target,]
table(harmonized_seurat_object@meta.data$Dataset[harmonized_seurat_object@meta.data$leiden_clusters_13 == cluster_target])
cellsh = harmonized_seurat_object@meta.data$Cell_ID[harmonized_seurat_object@meta.data$leiden_clusters_13 == cluster_target]
p1 = DimPlot(harmonized_seurat_object,group.by = "Dataset",raster = F,cells.highlight = cellsh,sizes.highlight = 0.1)+NoLegend()
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)


####### manually add annoation
harmonized_seurat_object@meta.data$Author_Class_Curated[harmonized_seurat_object@meta.data[,cluster_column]=="122"] = "Fibroblast"
harmonized_seurat_object@meta.data$Author_Class_Curated[harmonized_seurat_object@meta.data[,cluster_column]=="146"] = "Differentiating"
harmonized_seurat_object@meta.data$Author_Class_Curated[harmonized_seurat_object@meta.data[,cluster_column]=="198"] = "Neurons"
harmonized_seurat_object@meta.data$Author_Class_Curated[harmonized_seurat_object@meta.data[,cluster_column]=="201"] = "Doublet"
harmonized_seurat_object@meta.data$Author_Class_Curated[harmonized_seurat_object@meta.data[,cluster_column]=="202"] = "Differentiating"
harmonized_seurat_object@meta.data$Author_Class_Curated[harmonized_seurat_object@meta.data[,cluster_column]=="203"] = "Doublet"
harmonized_seurat_object@meta.data$Author_Class_Curated[harmonized_seurat_object@meta.data[,cluster_column]=="204"] = "Hypendymal"
harmonized_seurat_object@meta.data$Author_Class_Curated[harmonized_seurat_object@meta.data[,cluster_column]=="205"] = "Erythroid-like"


harmonized_seurat_object@meta.data$Author_Class_Curated[harmonized_seurat_object@meta.data[,cluster_column]=="23"] = "Astrocytes"
harmonized_seurat_object@meta.data$Author_Class_Curated[harmonized_seurat_object@meta.data[,cluster_column]=="47"] = "Tanycytes"
harmonized_seurat_object@meta.data$Author_Class_Curated[harmonized_seurat_object@meta.data[,cluster_column]=="84"] = "Mural"#
harmonized_seurat_object@meta.data$Author_Class_Curated[harmonized_seurat_object@meta.data[,cluster_column]=="195"] = "Immune"#
harmonized_seurat_object@meta.data$Author_Class_Curated[harmonized_seurat_object@meta.data[,cluster_column]=="193"] = "Immune"#
harmonized_seurat_object@meta.data$Author_Class_Curated[harmonized_seurat_object@meta.data[,cluster_column]=="191"] = "Krt19+ cells"#

# others:
harmonized_seurat_object@meta.data$Author_Class_Curated[harmonized_seurat_object@meta.data[,cluster_column]=="187"] = "Exclude"
harmonized_seurat_object@meta.data$Author_Class_Curated[harmonized_seurat_object@meta.data[,cluster_column]=="197"] = "Doublet"#

# show:
p1 = DimPlot(harmonized_seurat_object,group.by = "Author_Class_Curated",raster = F)
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)


##########
### Subset object
##########

keep_cells = harmonized_seurat_object@meta.data$Cell_ID[!harmonized_seurat_object@meta.data$Author_Class_Curated %in% c("Doublet","Unknown","Exclude") ]
nrow(harmonized_seurat_object@meta.data) - length(keep_cells) # remove 677 cells
curated_seurat_object = subset(harmonized_seurat_object,cells = keep_cells)


##########
### second curation round
##########

# Ideally I should ahve used a even higeh resolution for the inital clustering above
# When building the clusters for mrtree, I realized that I should take an additional loook at some small clusters:

## load other clusters:
hypoMap_test_additional_leiden_clustering = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_leiden_clustering_400.txt"),data.table = F)
clusters_to_add = hypoMap_test_additional_leiden_clustering[,c("Cell_ID","leiden_clusters_27")] %>% dplyr::rename(additional_clustering = leiden_clusters_27 )

temp_meta = dplyr::left_join(curated_seurat_object@meta.data,clusters_to_add,by="Cell_ID") %>% as.data.frame()
rownames(temp_meta) = temp_meta$Cell_ID
curated_seurat_object@meta.data = temp_meta

cluster_column = "additional_clustering"

#  c("360","330","386","364","377","321","325","314")

## helper code toc check on clusters:
cluster_target = "298"
#class_per_cluster[class_per_cluster[,cluster_column] == cluster_target,]
#head(basic_markers[basic_markers$cluster == cluster_target,],10)
#tmp_markers = basic_markers[basic_markers$cluster == cluster_target,]
table(curated_seurat_object@meta.data$Dataset[curated_seurat_object@meta.data[,cluster_column] == cluster_target])
table(curated_seurat_object@meta.data$Author_Class_Curated[curated_seurat_object@meta.data[,cluster_column] == cluster_target])
cellsh = curated_seurat_object@meta.data$Cell_ID[curated_seurat_object@meta.data[,cluster_column] == cluster_target]
p1 = DimPlot(curated_seurat_object,group.by = "Dataset",raster = F,cells.highlight = cellsh,sizes.highlight = 0.1)+NoLegend()
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

curated_seurat_object@meta.data$Author_Class_Curated[curated_seurat_object@meta.data[,cluster_column]=="360"] = "Doublet"
curated_seurat_object@meta.data$Author_Class_Curated[curated_seurat_object@meta.data[,cluster_column]=="330"] = "Doublet"
curated_seurat_object@meta.data$Author_Class_Curated[curated_seurat_object@meta.data[,cluster_column]=="386"] = "Differentiating"#
curated_seurat_object@meta.data$Author_Class_Curated[curated_seurat_object@meta.data[,cluster_column]=="364"] = "Differentiating"#
curated_seurat_object@meta.data$Author_Class_Curated[curated_seurat_object@meta.data[,cluster_column]=="377"] = "Differentiating"#
curated_seurat_object@meta.data$Author_Class_Curated[curated_seurat_object@meta.data[,cluster_column]=="321"] = "Astrocytes"# unknown !?
curated_seurat_object@meta.data$Author_Class_Curated[curated_seurat_object@meta.data[,cluster_column]=="325"] = "Unknown"#
curated_seurat_object@meta.data$Author_Class_Curated[curated_seurat_object@meta.data[,cluster_column]=="314"] = "Doublet"#

#curated_seurat_object@meta.data$log_nCount_RNA = log2(curated_seurat_object@meta.data$nCount_RNA)
p1=FeaturePlot(curated_seurat_object,features = "Slc32a1",raster = F,order=TRUE)
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)


##########
### Subset object again
##########

keep_cells = curated_seurat_object@meta.data$Cell_ID[!curated_seurat_object@meta.data$Author_Class_Curated %in% c("Doublet","Unknown","Exclude") ]
nrow(curated_seurat_object@meta.data) - length(keep_cells) # remove 677 cells
curated_seurat_object = subset(curated_seurat_object,cells = keep_cells)

p1 = DimPlot(curated_seurat_object,group.by = "Author_Class_Curated",raster = F)
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

# manually remove another 26 cells
# cells_selected = CellSelector(p1)
#cell_to_remove=c()
cell_to_remove = unique(c(cell_to_remove,cells_selected))
keep_cells = curated_seurat_object@meta.data$Cell_ID[!curated_seurat_object@meta.data$Cell_ID %in% cell_to_remove]
curated_seurat_object = subset(curated_seurat_object,cells = keep_cells)

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



##########
### Save
##########

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


