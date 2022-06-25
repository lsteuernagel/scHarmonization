##########
### Load params
##########

library(Seurat)
library(magrittr)
library(dplyr)
source("R/harmonization_functions.R")

# load json file with all other information
parameter_list = jsonlite::read_json("data/parameters_harmonization_v2_7.json")
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

file_name_prefix = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_curated")

# load clusters
hypoMap_clustering = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_curated","_leiden_clustering.txt"),data.table = F)
hypoMap_clustering = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2c_harmonization/hypoMap_harmonized_curated_selected_fast_leiden_clustering.txt",data.table = F)
#hypoMap_clustering = hypoMap_clustering[,c(1,3:ncol(hypoMap_clustering))]

# load seurat
#curated_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_curated",".rds"))

# get nn
neighbor_nn = as.Neighbor(curated_seurat_object@graphs$NN_scvi)
nn_idx = neighbor_nn@nn.idx

hypoMap_clustering_clear = apply(hypoMap_clustering[,2:ncol(hypoMap_clustering)],2,clear_clustering,min_cells = parameter_list$min_cells_valid,nn_idx = nn_idx) %>% as.data.frame()
#length(table(hypoMap_clustering_clear$leiden_clusters_0.1))
sort(apply(hypoMap_clustering_clear,2,function(x){length(table(x))}))

# I drop th first level because it contains only 1 cluster
#hypoMap_clustering_clear = hypoMap_clustering_clear[,2:ncol(hypoMap_clustering_clear)]

# add cell id
hypoMap_clustering_clear = dplyr::bind_cols(Cell_ID = hypoMap_clustering$Cell_ID,hypoMap_clustering_clear)
# drop
#hypoMap_clustering_clear = hypoMap_clustering_clear[,!grepl("_48|_19|_7",colnames(hypoMap_clustering_clear))]

curated_seurat_object@meta.data$temp_cluster = hypoMap_clustering_clear$leiden_clusters_45
p1 = DimPlot(curated_seurat_object,group.by = "temp_cluster",raster = F,label=TRUE,label.size = 2,repel = TRUE)+NoLegend()
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

p1 = DimPlot(curated_seurat_object,group.by = "Author_Class_Curated",raster = F,label=TRUE)+NoLegend()
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

p1=FeaturePlot(curated_seurat_object,features = "Pax6",raster = F,order=TRUE)
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

a1 = scUtils::gene_pct_cluster(curated_seurat_object,genes = c("Pax6","Trh","Gng8","Slc17a6"),col_name = "temp_cluster")
a1$cluster = rownames(a1)

##########
### define cluster for mrtree
##########

# make summary of clusters per level
hypoMap_clustering_long = hypoMap_clustering_clear %>% dplyr::select(-Cell_ID) %>% tidyr::gather(key="cluster_level",value="cluster_id")
hypoMap_clustering_long_summary = hypoMap_clustering_long %>% dplyr::group_by(cluster_level) %>%
  dplyr::summarise(n=length(unique(cluster_id))) %>% dplyr::arrange(n)
hypoMap_clustering_long_summary$cluster_level = factor(hypoMap_clustering_long_summary$cluster_level,levels = hypoMap_clustering_long_summary$cluster_level)

# take neuron and non-neuron
# decide which other cluster levels to use:
# take the lowest cluster
start_value = hypoMap_clustering_long_summary$n[hypoMap_clustering_long_summary$n == min(hypoMap_clustering_long_summary$n)][1]
start_value=8
# take the highest cluster
last_value = hypoMap_clustering_long_summary$n[hypoMap_clustering_long_summary$n == max(hypoMap_clustering_long_summary$n)]
# 4 other levels
x = start_value
target_clusters = x
while(x < last_value){
  x = x*2 + target_clusters[1] # this line controls how clusters are selected
  x = hypoMap_clustering_long_summary$n[abs(hypoMap_clustering_long_summary$n - x) == min(abs(hypoMap_clustering_long_summary$n - x))][1]
  if((last_value - x) > (x - target_clusters[length(target_clusters)])){target_clusters=c(target_clusters,x)}
}
target_clusters = c(target_clusters,last_value)
target_clusters

## make labels
cluster_cols = as.character(hypoMap_clustering_long_summary$cluster_level)[hypoMap_clustering_long_summary$n %in% target_clusters]
cluster_cols
# mrtree_input_labels = hypoMap_clustering_clear[,cluster_cols]

# add
cluster_cols[cluster_cols == "leiden_clusters_0.25"] = "leiden_clusters_0.5"
cluster_cols[cluster_cols == "leiden_clusters_0.75"] = "leiden_clusters_6"
cluster_cols = cluster_cols[cluster_cols!="leiden_clusters_3"]
cluster_cols[cluster_cols == "leiden_clusters_15"] = "leiden_clusters_16"
cluster_cols[cluster_cols == "leiden_clusters_50"] = "leiden_clusters_45"
#hypoMap_clustering_clear2 = hypoMap_clustering_clear[hypoMap_clustering_clear$Cell_ID %in% curated_seurat_object@meta.data$Cell_ID,]

mrtree_input_labels = hypoMap_clustering_clear[,cluster_cols]

##########
### manually add top level annotation
##########

neuron_label = curated_seurat_object@meta.data$Author_Class_Curated
neuron_label[neuron_label != "Neurons"] = 1
neuron_label[neuron_label == "Neurons"] = 0
mrtree_input_labels = cbind(neuron_label,mrtree_input_labels)
mrtree_input_labels <- apply(mrtree_input_labels, 2, as.character)
mrtree_input_labels <- apply(mrtree_input_labels, 2, as.numeric)
mrtree_input_labels = as.data.frame(mrtree_input_labels)
colnames(mrtree_input_labels)[1] = "leiden_clusters_0"

apply(mrtree_input_labels,2,function(x){length(table(x))})


##########
### manually adjust some clusters
##########



##########
### manually adjust some clusters
##########

## move tac2 to other vglut2+
mrtree_input_labels$leiden_clusters_0.0075[mrtree_input_labels$leiden_clusters_0.0075 == "6"] = "1"
# move Pomc to otehr vglu2+
mrtree_input_labels$leiden_clusters_0.0075[mrtree_input_labels$leiden_clusters_0.05 == "12"] = "1"

# Tbx19/pirt _--- > no fine
# mrtree_input_labels$leiden_clusters_0.5[mrtree_input_labels$leiden_clusters_6 == "79"] = "34"

# sst
mrtree_input_labels$leiden_clusters_0.5[mrtree_input_labels$leiden_clusters_6 == "169"] = "41"
mrtree_input_labels$leiden_clusters_0.5[mrtree_input_labels$leiden_clusters_6 == "40"] = "41"


# Grnh1+
mrtree_input_labels$leiden_clusters_0.05[mrtree_input_labels$leiden_clusters_0.5 == "66"] = "23" # make its own neuron cluster
mrtree_input_labels$leiden_clusters_0.0075[mrtree_input_labels$leiden_clusters_0.5 == "66"] = "1" # also move to neurons here
mrtree_input_labels$leiden_clusters_0[mrtree_input_labels$leiden_clusters_0.5 == "66"] = "0" # also move to neurons here

# pmch:
mrtree_input_labels$leiden_clusters_0.05[mrtree_input_labels$leiden_clusters_0.5 == "56"] = "24" # make its own neuron cluster
mrtree_input_labels$leiden_clusters_0.0075[mrtree_input_labels$leiden_clusters_0.5 == "56"] = "1" # also move to neurons here

# hypendymal Wfdc2+
mrtree_input_labels$leiden_clusters_0.5[mrtree_input_labels$leiden_clusters_6 == "199"] = "25" #
mrtree_input_labels$leiden_clusters_0.05[mrtree_input_labels$leiden_clusters_6 == "199"] = "8" #
mrtree_input_labels$leiden_clusters_0.0075[mrtree_input_labels$leiden_clusters_6 == "199"] = "3" #


# strange immune cells to immune:
#mrtree_input_labels$leiden_clusters_0.05[mrtree_input_labels$leiden_clusters_0.5 == "60"] = "7"

# # strange immune cells to immune:
# mrtree_input_labels$leiden_clusters_0.05[mrtree_input_labels$leiden_clusters_0.5 == "60"] = "7"

# check

#a2= as.data.frame(table(curated_seurat_object@meta.data$Author_Class_Curated,curated_seurat_object@meta.data$leiden_clusters_6))
#
# a3= as.data.frame(table(curated_seurat_object@meta.data$leiden_clusters_0.0075,curated_seurat_object@meta.data$leiden_clusters_0.05))

##########
### Plot
##########

# remove existing if running again on same object:
curated_seurat_object@meta.data = curated_seurat_object@meta.data[,!grepl("leiden_clusters_",colnames(curated_seurat_object@meta.data))]
#add
curated_seurat_object@meta.data = cbind(curated_seurat_object@meta.data,mrtree_input_labels)

p1 = DimPlot(curated_seurat_object,group.by = "leiden_clusters_0.5",raster = F,label=TRUE,label.size = 2.5)+NoLegend() #,label.size = 2
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)


p1=FeaturePlot(curated_seurat_object,features = "Wfdc2",raster = F,order=TRUE)
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

##########
### Save
##########

data.table::fwrite(mrtree_input_labels,paste0(parameter_list$harmonization_folder_path,parameter_list$clusters_for_mrtree_file))



