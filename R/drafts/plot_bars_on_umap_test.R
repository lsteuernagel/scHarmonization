

example_barplot_data = curated_seurat_object@meta.data %>% dplyr::group_by(C180,Dataset) %>% dplyr::summarise(dataset_pct = n()) %>% dplyr::group_by(C180) %>%
  mutate(dataset_pct = dataset_pct / sum(dataset_pct)) %>% dplyr::filter(C180=="C180-1" & Dataset %in%  c("Rupp10x","ChenDropseq"))
example_barplot = ggplot(example_barplot_data,aes(x=Dataset,y=dataset_pct,fill=Dataset))+geom_col(show.legend = FALSE)+theme_minimal()+xlab("")+ylab("")
example_barplot


umap_plot = DimPlot(curated_seurat_object,group.by = "C180",raster = F)+NoLegend()
umap_plot + annotation_custom(ggplotGrob(example_barplot),xmin = -5,xmax = -3,ymin = 5,ymax = 7)


#

#simulated_data = data.frame(cluster_id = sample(unique(curated_seurat_object@meta.data$C180),15),
#                            enrichment = rnorm(n = 15,mean = 0.5,sd=0.25))
a1=scUtils::gene_pct_cluster(curated_seurat_object,col_name = "C180",genes = "Pomc")
a2 = a1 %>% dplyr::arrange(desc(Pomc)) %>% head(n=6)
simulated_data = data.frame(cluster_id= rownames(a2))
simulated_data$enrichment = rnorm(n = nrow(simulated_data),mean = 0.5,sd=0.25)
#simulated_data$enrichment[simulated_data$cluster_id == "C180-48"] = 0.6

get_coordinates = function(label_vector,label_column,seurat_object,reduction_name = "umap_scvi"){
  plotdata = bind_cols(seurat_object@reductions[[reduction_name]]@cell.embeddings,label=seurat_object@meta.data[,label_column]) %>%
    dplyr::filter(label %in% label_vector)
  label_centers = plotdata %>% dplyr::group_by(label) %>% dplyr::summarise_at(vars(matches("umap")), median)
  return(label_centers)
}

coordinate_centers = get_coordinates(simulated_data$cluster_id,label_column = "C180",seurat_object = curated_seurat_object)
simulated_data = dplyr::left_join(simulated_data,coordinate_centers,by=c("cluster_id"="label"))

# make base umap_plot
max_height = 2
below_center_y = 0.2
max_width = 0.5
color = "darkblue"
#umap_plot = DimPlot(curated_seurat_object,group.by = "C180",raster = F,label=TRUE,label.size = 2)+NoLegend()
umap_plot =FeaturePlot(curated_seurat_object,features = "Pomc",raster = F,order=TRUE)

for(i in 1:nrow(simulated_data)){
  temp_barplot = ggplot(simulated_data[i,],aes(x=cluster_id,y=enrichment))+geom_col(show.legend = FALSE,fill=color)+
    theme_bw()  + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),legend.position="none",
                        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),plot.background=element_blank())+
    xlab("")+ylab("")
  # add
  umap_plot = umap_plot +
    annotation_custom(ggplotGrob(temp_barplot),
                      xmin = simulated_data$umapscvi_1[i]-(max_width/2),
                      xmax = simulated_data$umapscvi_1[i]+(max_width/2),
                      ymin = simulated_data$umapscvi_2[i] - below_center_y,
                      ymax = simulated_data$umapscvi_2[i] + (max_height*simulated_data$enrichment[i]) - below_center_y)

}

scUtils::rasterize_ggplot(umap_plot,pixel_raster = 2048,pointsize = 1.8)

# TODO: maybe add some way to show max value (grey bar and only filled o x percent ?) # or small second bar with max vlue
# TODO: maybe add text labels for specific clusters
# TODO: maybe color text of cluster verfied by RNAscope

