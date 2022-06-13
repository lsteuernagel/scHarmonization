

#hypoMap_region_annotation_df
hypoMap_region_annotation_df$Region_curated = hypoMap_region_annotation_df$Region
hypoMap_region_annotation_df$Region_curated[hypoMap_region_annotation_df$Region_curated %in% c("Lateral mammillary nucleus","Medial mammillary nucleus")] = "Mammillary nucleus"
hypoMap_region_annotation_df$Region_curated[hypoMap_region_annotation_df$Region_curated %in% c("Median eminence")] = "Arcuate hypothalamic nucleus"
hypoMap_region_annotation_df$Region_curated[hypoMap_region_annotation_df$Region_curated %in% c("Periventricular hypothalamic nucleus, intermediate part","Anteroventral periventricular nucleus","Periventricular hypothalamic nucleus, preoptic part","Periventricular hypothalamic nucleus, posterior part")] = "Periventricular hypothalamic nucleus"
hypoMap_region_annotation_df$Region_curated[hypoMap_region_annotation_df$Region_confidence == 0] = NA

curated_seurat_object@meta.data = curated_seurat_object@meta.data[,1:43]
curated_seurat_object@meta.data = cbind(curated_seurat_object@meta.data ,mrtree_result$labelmat)


temp = dplyr::left_join(curated_seurat_object@meta.data,result_region[,c("cluster","Region_curated_Color")],by=c("C280"="cluster"))
# temp = dplyr::left_join(curated_seurat_object@meta.data,hypoMap_region_annotation_df[hypoMap_region_annotation_df$Region_confidence > 0.2,c("gene_set","Region_curated")],by=c("C196"="gene_set"))
#temp$Region[is.natemp$Region]
rownames(temp) = temp$Cell_ID
curated_seurat_object@meta.data = temp

sort(table(result_region$Region_curated_Color),decreasing = TRUE)


p1 = DimPlot(curated_seurat_object,group.by = "Region_curated_Color",raster = F)+ggplot2::scale_color_manual(values = color_value_vector,na.value = "grey60")
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)


p1 = DimPlot(curated_seurat_object,group.by = "C180",raster = F,label.size = 2,label = TRUE)+NoLegend()
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

p1 = DimPlot(curated_seurat_object,group.by = "C23_named",raster = F,label.size = 4,label = TRUE)+NoLegend()
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)


p1 = DimPlot(curated_seurat_object,group.by = "Dataset",raster = F,label.size = 2,label = F,shuffle = TRUE)+NoAxes()#+NoLegend()
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)


######### Figure 1 & 3

## for campbell
curated_seurat_object@meta.data$campbell_anno = NA
curated_seurat_object@meta.data$campbell_anno[curated_seurat_object@meta.data$Dataset=="CampbellDropseq"] = curated_seurat_object@meta.data$Author_CellType[curated_seurat_object@meta.data$Dataset=="CampbellDropseq"]
p1 = DimPlot(curated_seurat_object,group.by = "campbell_anno",raster = F,label.size = 4,label = TRUE,order = TRUE)+NoLegend()+NoAxes()+scale_color_discrete(na.value="grey90")+ggtitle("Campbell ARH celltypes")
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

## for dowsett
curated_seurat_object@meta.data$dowsett_anno = NA
curated_seurat_object@meta.data$dowsett_anno[curated_seurat_object@meta.data$Dataset=="Dowsett10xnuc"] = curated_seurat_object@meta.data$C180[curated_seurat_object@meta.data$Dataset=="Dowsett10xnuc"]
p1 = DimPlot(curated_seurat_object,group.by = "dowsett_anno",raster = F,label.size = 4,label = F,order = TRUE)+NoLegend()+NoAxes()+scale_color_discrete(na.value="grey90")+ggtitle("Nuc-seq data in HypoMap")
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

p1=FeaturePlot(curated_seurat_object,features = "Th",raster = F,order=TRUE)
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

# dotplot
# need to subset to
DotPlot()

pruned_annotation_result = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_harmonization/hypoMap_harmonized_curated_pruned_annotation_result.txt",data.table = F)

table(curated_seurat_object@meta.data$Dataset[curated_seurat_object@meta.data$C280 == "C280-163"])

cellsh = curated_seurat_object@meta.data$Cell_ID[ curated_seurat_object@meta.data$Dataset == "Moffit10x"]
cellsh = curated_seurat_object@meta.data$Cell_ID[ curated_seurat_object@meta.data$C180 == "C180-18"]
p1 = DimPlot(curated_seurat_object,group.by = "C280",raster = F,label.size = 2,label = TRUE,cells.highlight = cellsh,sizes.highlight = 0.1)+NoLegend()
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

######### Figure 1 UMAPs

## make feature plot with multiple genes:
p <- FeaturePlot(curated_seurat_object,features = c("Slc32a1","Slc17a6","Slc18a2","Trh","Vip","Crh","Npy","Cck"), combine = FALSE,raster = F,order=TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
  p[[i]] <- scUtils::rasterize_ggplot(p[[i]],pixel_raster = 2048,pointsize = 1.8)
}
cowplot::plot_grid(plotlist = p,ncol = 2)

######### agrp dowsett zoom plot:

dowsett_subset = curated_seurat_object@meta.data[curated_seurat_object@meta.data$Dataset=="Dowsett10xnuc","Cell_ID"]
dowsett_subset = subset(curated_seurat_object,cells=dowsett_subset)
dowsett_subset_agrp = subset(dowsett_subset,subset = C62 == "C62-25" & umapscvi_1 > 2 & umapscvi_2 > 2)

FeaturePlot(dowsett_subset_agrp,"1700016P03Rik",split.by = "Diet",keep.scale = "all")

DimPlot(dowsett_subset_agrp,group.by = "C180",label = TRUE)
DimPlot(dowsett_subset_agrp,group.by = "Diet",label = TRUE)

full_subset_agrp = subset(curated_seurat_object,subset = C62 == "C62-25" & umapscvi_1 > 2 & umapscvi_2 > 2)

FeaturePlot(full_subset_agrp,"1700016P03Rik")
DimPlot(full_subset_agrp,group.by = "Dataset",label = TRUE)

rupp_subset = subset(curated_seurat_object,subset = Dataset=="Rupp10x")
rupp_subset_agrp = subset(rupp_subset,subset = C62 == "C62-25" & umapscvi_1 > 2 & umapscvi_2 > 2)
DimPlot(rupp_subset_agrp,group.by = "Author_Condition",label = TRUE)

### for figure plot:
cellsh = curated_seurat_object@meta.data$Cell_ID[curated_seurat_object@meta.data$C62=="C62-25"]
DimPlot(dowsett_subset, cells.highlight = cellsh,sizes.highlight = 0.1,cols.highlight = "orange")
FeaturePlot(dowsett_subset_agrp,"Agrp")
FeaturePlot(dowsett_subset_agrp,"Fos",split.by = "Diet",keep.scale = "all")
FeaturePlot(dowsett_subset_agrp,"1700016P03Rik",split.by = "Diet",keep.scale = "all")


######## supplement

## make feature plot with multiple genes:
all_genes = c("Ghrh","Tbx19","Oxt","Pomc","Sst","Nkx2-4")
# p <- list()
# for(i in 1:length(all_genes)) {
#   p[[i]] <- FeaturePlot(curated_seurat_object,features = c("Glp1r",all_genes[i]), combine = TRUE,raster = F,order=TRUE,blend = TRUE,blend.threshold = 0.1)+NoLegend()
#   p[[i]] <- p[[i]][[3]] + NoLegend() + NoAxes()
#   p[[i]] <- scUtils::rasterize_ggplot(p[[i]],pixel_raster = 2048,pointsize = 1.8)
# }
# cowplot::plot_grid(plotlist = p,ncol = 2)
cols_for_feature_plot = c("grey90","#ad5c00") # "#0b3ebd"
p <- list()
for(i in 1:length(all_genes)) {
  print(i)
  curated_seurat_object@meta.data$temp_score = CalculateMultScore(curated_seurat_object,features = c("Glp1r",all_genes[i]))
  p[[i]] <- FeaturePlot(curated_seurat_object,"temp_score",order=TRUE,raster=FALSE,cols = cols_for_feature_plot)+ NoLegend() + NoAxes()+ggtitle(paste0("Glp1r + ",all_genes[i]))
  p[[i]] <- scUtils::rasterize_ggplot(p[[i]],pixel_raster = 2048,pointsize = 1.8)
}
cowplot::plot_grid(plotlist = p,ncol = 2)


CalculateMultScore = function(object,features,cells=NULL,allowed_zeros = 0){
  gene_data = FetchData(object,vars = features,cells = cells)
  if(allowed_zeros == 0){
    multScore = apply(gene_data,1,prod)
    multScore = multScore^(1/ncol(gene_data)) # normalize
  }else{
    multScore = apply(gene_data,1,function(x,subtr){
      zero_idx = which(x == 0) # get all zeros
      # decide which ones to drop (maybe not all)
      if(length(zero_idx)>subtr){
        drop_idx = zero_idx[(length(zero_idx)-subtr+1):length(zero_idx)]
      }else{
        drop_idx = zero_idx
      }
      y = x[!drop_idx] # exclude subtr zeros from vector
      z = prod(y)^(1/length(y)) # normalize
    },subtr=allowed_zeros)
  }
  return(multScore)
}
curated_seurat_object@meta.data$temp_score <- CalculateMultScore(object = curated_seurat_object, features = c("Lepr","Pomc"))
FeaturePlot(object = curated_seurat_object, features = "temp_score",order=TRUE,raster=FALSE,cols = cols_for_feature_plot)+ NoAxes()+ggtitle(paste0("Lepr + Pomc"))

