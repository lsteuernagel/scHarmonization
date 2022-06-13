
color_value_vector =unlist(jsonlite::read_json("data/region_prediction_mapping_colors.json"))
# leaf_lvel
leaf_level_column= "C180"

# make data for first heatmap with percentages per dataset
heatmap_data = curated_seurat_object@meta.data %>% dplyr::select(Cell_ID,Dataset,!!sym(leaf_level_column)) %>% dplyr::group_by(!!sym(leaf_level_column),Dataset) %>% #dplyr::filter(predicted_Campbell!="NA")
  dplyr::count(name = "presence")  %>% dplyr::group_by(!!sym(leaf_level_column)) %>% dplyr::mutate(presence = presence / sum(presence)*100) %>% dplyr::ungroup() %>% #%>%  dplyr::left_join(tree_data_tibble[,c("label","node")],by=c("K169"="label"))
  tidyr::spread(key = Dataset,value=presence) %>% as.data.frame()
heatmap_matrix = as.matrix(heatmap_data[,2:ncol(heatmap_data)])
rownames(heatmap_matrix) = heatmap_data[,leaf_level_column]
heatmap_matrix[is.na(heatmap_matrix)] = 0

# make data for second heatmap with regions
heatmap_data2 = curated_seurat_object@meta.data %>% dplyr::select(Cell_ID,Region_summarized,!!sym(leaf_level_column)) %>%
  dplyr::group_by(!!sym(leaf_level_column),Region_summarized) %>% dplyr::count() %>% dplyr::group_by(!!sym(leaf_level_column)) %>%
  dplyr::top_n(n = 1,wt = n) %>% ungroup() %>%
  dplyr::distinct(!!sym(leaf_level_column),Region_summarized,.keep_all=TRUE) %>% as.data.frame()
heatmap_matrix2 = as.matrix(heatmap_data2[,"Region_summarized",drop=F])
rownames(heatmap_matrix2) = heatmap_data2[,leaf_level_column]
colnames(heatmap_matrix2) = "Region"

## make anno_df outside of function:
#"cluster_id","clusterlevel","cluster_name","first_cluster_name"
anno_df = annotation_result %>% dplyr::select(cluster_id,clusterlevel,cluster_name = clean_names)
anno_df$first_cluster_name = sapply(anno_df$cluster_name,function(x){strsplit(x,"\\.")[[1]][1]})


require(RColorBrewer)
## expand palette size
colourCount <- length(unique(heatmap_data2$suggested_region_curated)) # number of levels
getPalette <- colorRampPalette(brewer.pal(12, "Set3"))
bg_col = "grey90"
leaf_level = 6 # needs to be the matching one to leaf_level_column

# plot tree with heatmaps
circular_tree_heat = plot_cluster_tree(edgelist = curated_seurat_object@misc$clustering_edgelist,
                                       heatmap_matrix=heatmap_matrix,
                                       heatmap_matrix2 = heatmap_matrix2,
                                       leaf_level=leaf_level,
                                       anno_df = anno_df ,
                                       metadata=curated_seurat_object@meta.data,
                                       label_size = 2, show_genes = TRUE, legend_title_1 = "Pct", legend_title_2 = "Region",
                                       matrix_offset = 0.2, matrix_width =0.4,matrix_width_2 = 0.1,heatmap_colnames = TRUE,
                                       manual_off_second = 2,legend_text_size = 8,heatmap_text_size = 2,colnames_angle=0,hjust_colnames=0.5,
                                       heatmap_colors =c(bg_col,"darkred")) +
  scale_fill_brewer(palette = "Paired",na.value = bg_col)
#circular_tree_heat
require(ggtree)
circular_tree_heat_rotated = rotate_tree(circular_tree_heat, -90)
circular_tree_heat_rotated + ggplot2::scale_fill_manual(values = color_value_vector,na.value = "grey80")


