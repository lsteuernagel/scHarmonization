


parameter_list = jsonlite::read_json("data/parameters_harmonization_v2_8.json")
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})
parameter_list$marker_suffix = "pruned"
parameter_list$new_name_suffix = "hypoMap_harmonized_curated"
# load seurat
#harmonized_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,".rds"))

# load mrtree clustering
mrtree_result = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_mrtree_clustering_results",".rds"))

mrtree_output_raw = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_harmonization/hypoMap_harmonized_curated_curated_mrtree_res_raw.rds")

##
edgelist = mrtree_result$edgelist#[,1:2]
message("loaded edglist with ",nrow(edgelist)," rows.")
labelmat = mrtree_result$labelmat
message("loaded labelmat with ",nrow(labelmat)," rows.")

#a1 = mrtree_output_raw$labelmat.flat

apply(labelmat,2,function(x){length(unique(x))})

############

# make labelmat with colnames included
# consensus =FALSE
# labelmat.recon = mrtree_output_raw$labelmat.recon
# Ks.recon = apply(labelmat.recon, 2, function(y) length(unique(y)))
# # unique.idx = which(!duplicated(Ks.recon[-length(Ks.recon)], MARGIN = 2)) # remove the last column in labelmat.reco
# unique.idx = 1:length(Ks.recon) # don#t remove
# labelmat.tree = labelmat.recon[, unique.idx, drop=F]
# colnames(labelmat.tree) = paste0("K", Ks.recon[unique.idx])
#
# labelmat=labelmat.tree
# n=nrow(labelmat)
# backup = colnames(labelmat)
# labelmat = matrix(paste(matrix(rep(colnames(labelmat),each=n), nrow = n), labelmat, sep='-'), nrow = n)
# colnames(labelmat)=backup
df = as.data.frame(unique(labelmat), stringsAsFactors = F)
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
                      mrtree_output = mrtree_output_raw)

############ visualize

# find children in tree recusrively based on simple edgelist
find_ancestors = function(nodes,edges){
  current_ancestors = edges$from[edges$to %in% nodes]
  #print(paste0(current_children,collapse = "|"))
  if(length(current_ancestors)>0){
    all_ancestors = c(current_ancestors,find_ancestors(current_ancestors,edges))
  }else{
    all_ancestors = current_ancestors
  }
  return(all_ancestors)
}

library(igraph)
library(ggplot2)
library(scales)
library(treeio)
library(tidytree)
library(ggtree)
#library(ggnewscale)

edges = edges[edges$level <= 7,]

## convert to treedata
# only take
tree_data_igraph = base::suppressWarnings(igraph::graph_from_edgelist(as.matrix(edges[,1:2])))
tree_data_phylo = base::suppressWarnings(treeio::as.phylo(tree_data_igraph))
tree_data_tibble <- dplyr::as_tibble(tree_data_phylo)

# update additional columns
tree_data_tibble$nodesize = 1 # default node size
tree_data_tibble$n_children = sapply(tree_data_tibble$label,function(x,el){length(el$to[el$from==x])},el=edges[,1:2]) # count children number
tree_data_tibble$n_siblings = sapply(tree_data_tibble$label,function(x,el){ # count siblings
  parent = el$from[el$to==x]
  return(length(el$to[el$from==parent])-1)
},el=edges)

# convert back to treedata
tree_data = suppressWarnings(tidytree::as.treedata(tree_data_tibble))


# circular_tree +
#   geom_nodelab(aes(x=branch, label=first_cluster_name), size=label_size,vjust= -.5, color="darkred")+
#   geom_tiplab(ggplot2::aes(x=branch, label=first_cluster_name), size=label_size,vjust= -.5,color="darkred")
# circular_tree

label = "C180-140"#"K42-7"#"K160-112"
id_offset = 1.5
label_size =3.5
# add selection
tree_data_tibble$selected_edges = "not"
all_selected_nodes = c(label,find_ancestors(label,edges[,c("from","to")]))
tree_data_tibble$selected_edges[tree_data_tibble$label %in% all_selected_nodes] = "selected"
#plot circular tree
tree_data = suppressWarnings(tidytree::as.treedata(tree_data_tibble))
circular_tree =ggtree(tree_data,layout = 'circular', branch.length='none',aes(color=selected_edges))+
  geom_nodelab(ggplot2::aes(x=branch, label=label),color="grey20", size=label_size,hjust= -1*id_offset)+
  scale_color_manual(values = c("selected" = "blue","not" = 'grey60'))+guides(color= "none")+
  geom_nodepoint(aes(subset = n_children > 1))+geom_tippoint(aes(subset = selected_edges == "selected"))
# plot selection
# circular_tree = circular_tree +
#   geom_nodelab(ggplot2::aes(x=branch, label=label,color=selected_edges), size=label_size,hjust= -1*id_offset) + #hjust= .5
#   ggtree::geom_(ggplot2::aes(x=branch, label=label,color=selected_edges), size=label_size,hjust= -1*id_offset) #hjust= .5
circular_tree

##### plot on UMAP:
library(Seurat)
# load clusters
# hypoMap_clustering = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_curated","_leiden_clustering.txt"),data.table = F)

# load seurat
#curated_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_curated",".rds"))
curated_seurat_object@meta.data = curated_seurat_object@meta.data[,!grepl("C[0-9]+",colnames(curated_seurat_object@meta.data))]
curated_seurat_object@meta.data = cbind(curated_seurat_object@meta.data,labelmat)

apply(labelmat,2,function(x){length(unique(x))})

p1 = DimPlot(curated_seurat_object,group.by = "C185",raster = F,label=TRUE,label.size = 2.5,repel = TRUE)+NoLegend()
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

cluster_column = "C286"
set.seed(1234)
curated_seurat_object@meta.data[,cluster_column] = factor(curated_seurat_object@meta.data[,cluster_column], levels = sample(unique(curated_seurat_object@meta.data[,cluster_column]),length(unique(curated_seurat_object@meta.data[,cluster_column]))))
p1 = DimPlot(curated_seurat_object,group.by = cluster_column,raster = F,label=TRUE,label.size = 2.5)+NoLegend()
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)


p1=FeaturePlot(curated_seurat_object,features = "Sst",raster = F,order=TRUE)
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

a22 = CellSelector(p1)
a22_meta = curated_seurat_object@meta.data[curated_seurat_object@meta.data$Cell_ID %in% a22,]

