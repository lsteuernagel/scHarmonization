
##########
### Check on eliden input
##########

# TODO:
# need to load suerat
# need to load mretree res
# need to cbind labelmat!

# check Kim et all annoatations vs my clustering:
subset_seurat = subset(curated_seurat_object,subset = Dataset == "Kim10x")
# subset_seurat@meta.data$leiden_clusters_27 = as.character(subset_seurat@meta.data$leiden_clusters_27)
# subset_seurat@meta.data$leiden_clusters_5 = as.character(subset_seurat@meta.data$leiden_clusters_5)
# create cluster overview using mapscvi:
overview_clustering = mapscvi::compare_clustering(query_seura_object = subset_seurat,
                                                  clustering_1 = "Author_CellType",
                                                  clustering_2 = "C316",
                                                  min_cells = 5,
                                                  min_pct = 0.05,
                                                  return_data=TRUE)
# sankey
# i focus on the vmh part of the clusters
clustering_2_filter = unique(subset_seurat@meta.data$C316[subset_seurat@meta.data$C23 == "C23-2"]) # neurons only
#clustering_2_filter = unique(subset_seurat@meta.data$leiden_clusters_27[subset_seurat@meta.data$leiden_clusters_0.01 == "1"]) # neurons only
clustering_1_filter = NULL
sankey_clusters = mapscvi::plot_sankey_comparison(overview_clustering,
                                             clustering_1_filter = clustering_1_filter,
                                             clustering_2_filter = clustering_2_filter,
                                             text_size=20,
                                             use_and = FALSE,
                                             light_factor = 0.45)
sankey_clusters

# length(unique(overview_clustering$clustering_1))
# length(unique(overview_clustering$clustering_2))
#

# make a subset of the vmh subclusters
vmh_seurat = subset(curated_seurat_object,subset = leiden_clusters_0.05 == "3")
vmh_seurat = subset(vmh_seurat,subset = umapscvi_1 < -3 & umapscvi_2 > -3 & umapscvi_2 < 6)
DimPlot(vmh_seurat,group.by = "leiden_clusters_27",label=TRUE)+NoLegend()

vmh_seurat@meta.data$Author_CellType_Kim10x =NA
vmh_seurat@meta.data$Author_CellType_Kim10x[vmh_seurat@meta.data$Dataset=="Kim10x"] = vmh_seurat@meta.data$Author_CellType[vmh_seurat@meta.data$Dataset=="Kim10x"]

DimPlot(vmh_seurat,group.by = "Author_CellType_Kim10x",label=TRUE)+NoLegend()



