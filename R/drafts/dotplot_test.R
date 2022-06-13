
target_col="C180_named"
panel_col = "C23_named"
cluster_to_include = unique(curated_seurat_object@meta.data[curated_seurat_object@meta.data$C2_named=="C2-1: Neurons",target_col])
#cluster_to_include = cluster_to_include[1:10]
cluster_to_include_wo = stringr::str_extract(cluster_to_include,"C180-[0-9]+")
top_features = curated_seurat_object@misc$marker_genes_all %>% dplyr::filter(specificity > 5 & cluster_id %in% cluster_to_include_wo) %>%
  dplyr::filter(!grepl("Rik|Gm",gene)) %>%
  dplyr::arrange(desc(specificity)) %>% dplyr::group_by(cluster_id) %>% dplyr::top_n(n = 1,wt = specificity)

# TODO: dont use top n features but instead the genes used for annotation !!!!!!!!

# need to reorder factor level in seurat
# order by number --> then it also is ordered by tree
unordered = unique(curated_seurat_object@meta.data[,target_col])
unordered = as.numeric(stringr::str_remove(stringr::str_extract(unordered,"-[0-9]+"),"-"))
names(unordered) = unique(curated_seurat_object@meta.data[,target_col])
curated_seurat_object@meta.data[,target_col] = factor(curated_seurat_object@meta.data[,target_col],levels = names(sort(unordered)))
# also reorder panel_col ?

Idents(curated_seurat_object) = target_col
dp1 = DotPlot(object = curated_seurat_object,
              features = unique(top_features$gene),
              #  group.by = target_col,
              idents= cluster_to_include,
              scale = FALSE,
              # split.by = panel_col,
              cluster.idents = F)#+facet_grid( ~ C23_named)
dp1 + theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust=0.75))

