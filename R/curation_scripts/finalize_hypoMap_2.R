##########
### Scrip Aim
##########

# This script gathers all data from the hypoMap harmoinzation and curation, to save the most relevant files and a final seurat object (and anndata) into a sepratae folder.

##########
### Load parameters and packages
##########

message(" Load parameters and packages ")

require(tidyverse)
require(Seurat)
require(Matrix)
require(scCoco)

source("R/harmonization_functions.R")

#load parameters from annotation:
parameter_list = jsonlite::read_json("data/parameters_annotation_v2_2.json")
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# read features to excludes
features_exclude_list= unlist(jsonlite::read_json(parameter_list$genes_to_exclude_file))

# define new folder for output
parameter_list$final_output_folder = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2c_final/" # can also be subfolder of parameter_list$harmonization_folder_path
parameter_list$file_name_prefix = "hypoMap_v2"

min_specificity_markers = 0.5
min_specificity_markers_seurat = 1.5

##########
### Load objects
##########

# load seurat
curated_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,".rds"))

#curated_seurat_object@meta.data = curated_seurat_object@meta.data[,1:43]

# load mrtree clustering
mrtree_result = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_mrtree_clustering_results",".rds"))

# annotation results
annotation_result = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_annotation_result.txt"),data.table = F)
annotation_labelmat =  data.table::fread(file = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_annotation_labelmat.txt"),data.table = F)
annotation_markers_filtered =  data.table::fread(file =paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_annotation_markers_filtered.txt"),data.table = F)

# region prediction result
region_prediction = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_region_prediction.txt"),data.table = F)

##########
### Adjust annotation and add to seurat
##########

# sort annoation labelmat to metadata
annotation_labelmat = annotation_labelmat[match(curated_seurat_object@meta.data$Cell_ID,annotation_labelmat$Cell_ID),]

# add labels with annotation
label_df = cbind(mrtree_result$labelmat,annotation_labelmat %>% dplyr::select(-Cell_ID)) %>% as.data.frame()

## add labels to seurat:
curated_seurat_object@meta.data = cbind(curated_seurat_object@meta.data , label_df)

# check that named labels are ordered correctly:
p1 = DimPlot(curated_seurat_object,group.by = "C185_named",raster = F,label.size = 2,label = TRUE,repel = TRUE)+NoLegend()
scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)


##########
### Load and filter marker objects
##########

# read all available marker tables if marker detection was run on multiple subsets !!
markers_comparisons_all_list=list()
markers_comparisons_siblings_list =list()
for(current_start_node in parameter_list$start_nodes_annotation_markers){
  # load markers all
  filename=paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",current_start_node,"_markers_all_",parameter_list$marker_suffix,".tsv")
  if(file.exists(filename)){
    markers_comparisons_all_list[[current_start_node]] = data.table::fread(filename,data.table = F)
  }else{
    message("Cannot find markers stored in : ",filename)
  }
  # siblings
  filename_sib=paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",current_start_node,"_markers_siblings_",parameter_list$marker_suffix,".tsv")
  if(file.exists(filename_sib)){
    markers_comparisons_siblings_list[[current_start_node]] = data.table::fread(filename_sib,data.table = F)
  }else{
    message("Cannot find markers stored in : ",filename_sib)
  }
}
markers_comparisons_all = as.data.frame(do.call(rbind,markers_comparisons_all_list))
markers_comparisons_siblings  = as.data.frame(do.call(rbind,markers_comparisons_siblings_list))
message("All markers for: ",length(unique(markers_comparisons_all$cluster_id))," clusters available")
message("Sibling markers for: ",length(unique(markers_comparisons_siblings$cluster_id))," clusters available")

# load markers for exclusion ad further filter
# load additional remove
additional_remove_genes = unlist(jsonlite::read_json(unlist(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"raw","_additionally_removed_markers.json"))))
# gather some other genes that are not informative during annoation:
other_genes_remove = rownames(curated_seurat_object@assays$RNA@counts)[grepl("RP",rownames(curated_seurat_object@assays$RNA@counts))]
# make list of all genes that should be removed:
all_exclusion_genes = unique(c(features_exclude_list,additional_remove_genes,other_genes_remove))

# filter marker genes:
markers_comparisons_all = markers_comparisons_all %>% dplyr::filter(! gene %in% all_exclusion_genes &
                                                                      specificity >= min_specificity_markers &
                                                                      p_val_adj < parameter_list$max_pval_adj &
                                                                      pct.2 < 0.5
)

# filter marker genes:
markers_comparisons_siblings = markers_comparisons_siblings %>% dplyr::filter(! gene %in% all_exclusion_genes &
                                                                      specificity >= min_specificity_markers &
                                                                      p_val_adj < parameter_list$max_pval_adj &
                                                                      pct.2 < 0.5
)


# join labels with annotation to marker tables
markers_comparisons_all_final = dplyr::left_join(markers_comparisons_all,annotation_result %>% dplyr::select(cluster_id,cluster_name = clean_names_withID),by="cluster_id") %>%
  dplyr::select(gene,cluster_id,cluster_name,parent_id = parent,specificity,avg_log2FC,pct.1,pct.2,p_val_adj) # %>% dplyr::arrange(cluster_id,desc(specificity))
markers_comparisons_siblings_final = dplyr::left_join(markers_comparisons_siblings,annotation_result %>% dplyr::select(cluster_id,cluster_name = clean_names_withID),by="cluster_id") %>%
  dplyr::select(gene,cluster_id,cluster_name,parent_id = parent,specificity,avg_log2FC,pct.1,pct.2,p_val_adj) # %>% dplyr::arrange(cluster_id,desc(specificity))

# save marker table is @misc slot
curated_seurat_object@misc$marker_genes_all = markers_comparisons_all_final[markers_comparisons_all_final$specificity > min_specificity_markers_seurat,]
curated_seurat_object@misc$marker_genes_siblings = markers_comparisons_siblings_final[markers_comparisons_siblings_final$specificity > min_specificity_markers_seurat,]

# store edgelist in @misc
curated_seurat_object@misc$annotation_result = annotation_result
curated_seurat_object@misc$clustering_edgelist = mrtree_result$edgelist


##########
### Adjust predicted regions and summarize
##########

# hypoMap_old = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_final/hypoMap_v2.rds")

anno_comparison = dplyr::left_join(curated_seurat_object@meta.data[,c("Cell_ID","C286")],
                                   hypoMap_old@meta.data[,c("Cell_ID","C280","Region_predicted")],by="Cell_ID") %>%
  dplyr::group_by(C286,C280) %>% dplyr::add_count() %>% dplyr::group_by(C286) %>% dplyr::filter(n==max(n)) %>% dplyr::distinct(C286,C280,Region_predicted,n)# %>% dplyr::select(-Cell_ID)

# anno_comparison = dplyr::left_join(curated_seurat_object@misc$annotation_result, hypoMap_old@misc$annotation_result,by="cluster_id",suffix=c("new","old") )
# anno_comparison = dplyr::left_join(anno_comparison,hypoMap_old@meta.data[,c("")],by=c(""=""))

# read full prediction result:
full_region_prediction = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_region_prediction_all_results.rds"))
scores_per_target_level_region_all=full_region_prediction$full_region_prediction_result$scores_per_target_level_region_all %>% as.data.frame()

## check in there:
target = "C286-20"
scores_per_target_level_region_all[,c("topname",target)] %>% dplyr::arrange(desc(!!sym(target))) %>% top_n(n = 10,wt=!!sym(target))
markers_comparisons_all[markers_comparisons_all$cluster_id == target,] %>% dplyr::arrange(desc(specificity)) %>% dplyr::top_n(10,wt = specificity) %>% dplyr::select(gene,specificity,pct.1,pct.2)
annotation_result_with_region[annotation_result_with_region$cluster_id == target,] %>% t()
sort(table(curated_seurat_object@meta.data$Dataset[curated_seurat_object@meta.data$C286 == target]),decreasing = TRUE)

# read color scheme --> important !!
#color_value_vector =unlist(jsonlite::read_json("data/region_prediction_mapping_colors.json"))
annotation_result_with_region = dplyr::left_join(annotation_result,region_prediction,by=c("cluster_id"="cluster")) %>% dplyr::filter(clusterlevel=="C286")


#### Manually adjust:

# chat:
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-197"] = NA#"Arcuate hypothalamic nucleus"

# gnrh1
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-92"] = NA

# overwrite Ghrh C286-121
#"C286-121" = "Arcuate hypothalamic nucleus"
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-175"] = "Arcuate hypothalamic nucleus"
#"C286-12" ="Lateral hypothalamic area" # keep as lateral!

# overwrite C286-21 Prdm8
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5782816/
#"C286-21" = "Periventricular hypothalamic nucleus, preoptic part"
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-19"] = "Periventricular hypothalamic nucleus, preoptic part"

# Qrfp: Lateral hypothalamic area
# https://www.researchgate.net/figure/QRFP-producing-neurons-are-distributed-in-and-around-the-LHA-and-innervated-by-fibers_fig3_309957432
#"C286-68" = "Lateral hypothalamic area"
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-70"] = "Lateral hypothalamic area"

# npvf
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4364833/
# seems to be mostly in DMH and DMH is also highly enriched in annotaion process
#"C286-67" = "Dorsomedial nucleus of the hypothalamus"
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-69"] = "Dorsomedial nucleus of the hypothalamus"

# overwrite -14 ?

# "C286-3" = NA # NSCs
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-3"] = NA
# "C286-10" = "Lateral hypothalamic area" --> Ebf3/Sim1 . Not clear from prediction could be PVH or LH
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-10"] = NA #"Lateral hypothalamic area"
# "C286-8" = "Lateral hypothalamic area" NPW neurons
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-11"] = "Lateral hypothalamic area"
#"C286-28" = NA # not clear mostly Flyn/Dowsett/rupp, Medial mammillary nucleus seems likely
# see: https://pubmed.ncbi.nlm.nih.gov/19111912/ Lmx1a (--> and Ngfr)
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-7"] = "Medial mammillary nucleus"
#"C286-29" = "Ventromedial hypothalamic nucleus" --> Fezf1+ but too few markers for reliable prediction
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-30"] = "Ventromedial hypothalamic nucleus"
#"C286-42" = "Ventromedial hypothalamic nucleus" # clearly VMH with Esr1 /Cd40+ but ARC prediction was slightly higher
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-44"] = "Ventromedial hypothalamic nucleus"
# C286-52" = NA # not clear
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-55"] = NA

# We found neuropeptide VF precursor (Npvf),PR domain containing 13 (Prdm13), and SK1 family transcriptional corepressor (Skor1) as DMH-enriched genes. Particularly, Prdm13,
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4364833/
#"C286-53" = "Dorsomedial nucleus of the hypothalamus" # Prdm13
#"C286-54" = "Dorsomedial nucleus of the hypothalamus" # Prdm13
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-59"] = "Dorsomedial nucleus of the hypothalamus"
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-60"] = "Dorsomedial nucleus of the hypothalamus"
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-61"] = "Dorsomedial nucleus of the hypothalamus"

# others:
# NA # not clear Pthlh neurons but nt sure from where exatcly
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-134"] = NA
#" NA # not clear, maybe lateral or
# annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-109"] = NA # TODO
#= NA # not clear
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-144"] = NA
#NA # not clear, probably posterior Satb1+
# https://www.sciencedirect.com/science/article/pii/S0168010211001519?via%3Dihub
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-131"] = NA
# "Posterior hypothalamic nucleus"unclear!
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-119"] = NA#"Posterior hypothalamic nucleus"
#"Zona incerta" # just based on ABA prediction
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-116"] = "Zona incerta"
# NA # not clear (ARC ME predicted but many SCN cells ?)
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-129"] = NA
#""Lateral hypothalamic area" Crhbp+ and many LH dataset cells
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-98"] =  "Lateral hypothalamic area"
#"C286-38" = NA --> prediction too low
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-39"] =  NA
# not really clear from where (prediction more anterior, but datasets more posterior)
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-150"] =  NA
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-151"] =  NA
# annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-152"] =  NA # leave as is
#"C286-41" = "Medial mammillary nucleus" # overwrite --> don't seems to be VMH cluster!!
# annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-41"] =  "Medial mammillary nucleus"

# hdc nueorns in uerbomamillary
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5618271/
# 12 are in fact not hdc+!!
#annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-126"] =  "Tuberomamillary nucleus"
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-157"] =  "Tuberomamillary nucleus"
# Medial mammillary nucleus seems best, also gets the highest prediction !
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-129"] =  "Medial mammillary nucleus"

# some PVH clusters with Sim1+:
# Trh/Crh+
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-20"] =  "Paraventricular hypothalamic nucleus"
# PVH is second highest score and Trh+
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-13"] =  "Paraventricular hypothalamic nucleus"
# not so sure about this one !
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-15"] =  "Paraventricular hypothalamic nucleus"
# C286-18 is indeed Lateral hypothalamic area

# set these to NA
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-10"] =  NA
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-12"] =  NA

# some SCN clusters:
# Suprachiasmatic nucleus" # nms neurons SCN
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-123"] =  "Suprachiasmatic nucleus"
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-111"] =  "Suprachiasmatic nucleus"
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-115"] =  "Suprachiasmatic nucleus"
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-114"] =  "Suprachiasmatic nucleus"
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-112"] =  "Suprachiasmatic nucleus"
# not really clear as prediction and origin are so different
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-122"] =  NA
#"Suprachiasmatic nucleus" # vip neurons SCN
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-124"] =  "Suprachiasmatic nucleus"
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-125"] =  "Suprachiasmatic nucleus"
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-127"] =  "Suprachiasmatic nucleus"
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-129"] =  NA
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-186"] =  NA
annotation_result_with_region$Region_curated[annotation_result_with_region$cluster_id == "C286-187"] =  NA

#all
annotation_result_with_region$Region_curated[annotation_result_with_region$Region_curated==""] =NA

# add region for coloring:
hypothalamus_regions_mapping = data.table::fread("data/region_prediction_mapping.tsv",data.table = FALSE)
hypothalamus_regions_mapping$name[hypothalamus_regions_mapping$name == "Tuberomammillary nucleus"] = "Tuberomamillary nucleus"
annotation_result_with_region = dplyr::left_join(annotation_result_with_region,hypothalamus_regions_mapping %>% dplyr::select(name,Region_curated_Color=name_summarized),by=c("Region_curated"="name"))

## add to map
temp = dplyr::left_join(curated_seurat_object@meta.data,annotation_result_with_region %>% dplyr::select(cluster_id,Region_predicted = Region_curated , Region_summarized = Region_curated_Color),by=c("C286"="cluster_id"))
rownames(temp) = temp$Cell_ID
curated_seurat_object@meta.data = temp

# check:
# color_value_vector =unlist(jsonlite::read_json("data/region_prediction_mapping_colors.json"))
# p1 = DimPlot(curated_seurat_object,group.by = "Region_summarized",raster = F)+ggplot2::scale_color_manual(values = color_value_vector,na.value = "grey80")
# scUtils::rasterize_ggplot(p1,pixel_raster = 2048,pointsize = 1.8)

##########
### Clean up metadata
##########

metadata_clean_up = curated_seurat_object@meta.data[,!grepl("leiden_",colnames(curated_seurat_object@meta.data))]
colnames(metadata_clean_up)
#metadata_clean_up = metadata_clean_up %>% dplyr::select(-preliminary_clusters,-Processing_clusters,-Final_Doublet,-Final_Exclude,-additional_clustering,-temp_cluster,-Doublet)
metadata_clean_up = metadata_clean_up %>% dplyr::select(-preliminary_clusters,-seurat_clusters,-Processing_clusters,-Doublet)
rownames(metadata_clean_up) = metadata_clean_up$Cell_ID

curated_seurat_object@meta.data = metadata_clean_up

##########
### save files
##########

# check that folder exists
system(paste0("mkdir -p ",parameter_list$final_output_folder))
print(parameter_list$final_output_folder)

# save other relevant files
# markers:
data.table::fwrite(markers_comparisons_all_final,paste0(parameter_list$final_output_folder,parameter_list$file_name_prefix,"_markers_all.tsv"),sep="\t")
data.table::fwrite(markers_comparisons_siblings_final,paste0(parameter_list$final_output_folder,parameter_list$file_name_prefix,"_markers_siblings.tsv"),sep="\t")
# clustering tree:
data.table::fwrite(as.data.frame(mrtree_result$labelmat),paste0(parameter_list$final_output_folder,parameter_list$file_name_prefix,"_mrtree_labelmat.tsv"),sep="\t")
data.table::fwrite(mrtree_result$edgelist,paste0(parameter_list$final_output_folder,parameter_list$file_name_prefix,"_mrtree_edgelist.tsv"),sep="\t")
data.table::fwrite(annotation_result,paste0(parameter_list$final_output_folder,parameter_list$file_name_prefix,"_annotation_result.tsv"),sep="\t")
data.table::fwrite(annotation_result_with_region,paste0(parameter_list$final_output_folder,parameter_list$file_name_prefix,"_annotation_result_with_region.tsv"),sep="\t")

# copy model
system(paste0("cp -r ",parameter_list$harmonization_folder_path,"hypoMap_harmonized_scVI_model","/  ",parameter_list$final_output_folder))

# save data to rds
saveRDS(curated_seurat_object,paste0(parameter_list$final_output_folder,parameter_list$file_name_prefix,".rds"))

# save h5seurat
SeuratDisk::SaveH5Seurat(object = curated_seurat_object,filename = paste0(parameter_list$final_output_folder,parameter_list$file_name_prefix,".h5seurat"), overwrite = TRUE, verbose = TRUE)

# save to anndata
SeuratDisk::Convert( paste0(parameter_list$final_output_folder,parameter_list$file_name_prefix,".h5seurat"), dest =  paste0(parameter_list$final_output_folder,parameter_list$file_name_prefix,".h5ad"),assay="RNA",verbose=TRUE,overwrite=TRUE)
system(paste0("rm ",paste0(parameter_list$final_output_folder,parameter_list$file_name_prefix,".h5seurat")))

