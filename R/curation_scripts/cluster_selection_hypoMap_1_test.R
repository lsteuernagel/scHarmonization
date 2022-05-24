##########
### Load params
##########

# load json file with all other information
parameter_list = jsonlite::read_json("data/parameters_harmonization_v2_1_test.json")
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

file_name_prefix = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_curated")

# load clusters
hypoMap_clustering = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_curated","_leiden_clustering.txt"),data.table = F)
#hypoMap_clustering = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_harmonization_test/hypoMap_test_curated_initial_leiden_clustering.txt",data.table = F)

##########
### define cluster for mrtree
##########

# make summary of clusters per level
hypoMap_clustering_long = hypoMap_clustering %>% dplyr::select(-Cell_ID) %>% tidyr::gather(key="cluster_level",value="cluster_id")
hypoMap_clustering_long_summary = hypoMap_clustering_long %>% dplyr::group_by(cluster_level) %>%
  dplyr::summarise(n=length(unique(cluster_id))) %>% dplyr::arrange(n)
hypoMap_clustering_long_summary$cluster_level = factor(hypoMap_clustering_long_summary$cluster_level,levels = hypoMap_clustering_long_summary$cluster_level)

# take neuron and non-neuron
# decide which other cluster levels to use:
# take the lowest cluster
start_value = hypoMap_clustering_long_summary$n[hypoMap_clustering_long_summary$n == min(hypoMap_clustering_long_summary$n)]
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
mrtree_input_labels = hypoMap_clustering[,cluster_cols]

##########
### manually add top level annotation
##########

# load seurat
curated_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_curated",".rds"))

neuron_label = curated_seurat_object@meta.data$Author_Class_Curated
neuron_label[neuron_label != "Neurons"] = 1
neuron_label[neuron_label == "Neurons"] = 0
mrtree_input_labels = cbind(neuron_label,mrtree_input_labels)
mrtree_input_labels <- apply(mrtree_input_labels, 2, as.character)
mrtree_input_labels <- apply(mrtree_input_labels, 2, as.numeric)
mrtree_input_labels = as.data.frame(mrtree_input_labels)
summary(mrtree_input_labels)
colnames(mrtree_input_labels)[1] = "leiden_clusters_0"

##########
### Save
##########

data.table::fwrite(mrtree_input_labels,paste0(parameter_list$harmonization_folder_path,parameter_list$clusters_for_mrtree_file))



