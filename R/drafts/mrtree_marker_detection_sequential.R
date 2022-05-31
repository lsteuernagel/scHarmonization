##########
### Load parameters and packages
##########

message(Sys.time(),": Starting fast marker detection .." )

message(" Load parameters and packages ")

require(tidyverse)
require(Seurat)
require(Matrix)

source("R/harmonization_functions.R")
source("R/stratified_wilcoxon_functions.R")

# get params-filename from commandline
param_file = "data/parameters_harmonization_v2_1_test.json"
# read all parameters and filepaths
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# read features to excludes
features_exclude_list= jsonlite::read_json(parameter_list$genes_to_exclude_file)
features_exclude_list = lapply(features_exclude_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# load seurat
harmonized_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,".rds"))

# load mrtree clustering
mrtree_result = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_mrtree_clustering_results",".rds"))

# define output file name
output_file_name = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,parameter_list$basic_marker_filename,".txt")

# important: start node
start_node = parameter_list$start_node

if(is.null(start_node)){
  start_node = "all"
}

# markers
n_cores_markers = parameter_list$n_cores_markers
assay_markers = parameter_list$assay_markers
assay_slot =parameter_list$assay_slot
logfc.threshold = parameter_list$logfc.threshold
min.pct =parameter_list$min.pct
min.diff.pct = parameter_list$min.diff.pct
max.cells.per.ident = parameter_list$max.cells.per.ident
min.cells.feature = parameter_list$min.cells.feature
min.cells.group =  parameter_list$min.cells.group
base = parameter_list$base
only.pos = parameter_list$only.pos
batch_var = parameter_list$batch_var
test.use = parameter_list$test.use
specificity_base = parameter_list$specificity_base
message("test.use ",test.use," ... ",parameter_list$test.use)
if(is.null(base)){base = 2}

##########
### Prepare edgelist and labelmat
##########

# get key objects
edgelist = mrtree_result$edgelist[,1:2]
labelmat = mrtree_result$labelmat

# if start node if specified:
# find subtree and subset edgelist:

# find children in tree recusrively based on simple edgelist
find_children = function(nodes,edges){
  current_children = edges$to[edges$from %in% nodes]
  #print(paste0(current_children,collapse = "|"))
  if(length(current_children)>0){
    all_children = c(current_children,find_children(current_children,edges))
  }else{
    all_children = current_children
  }
  return(all_children)
}

# walk through tree:
all_children = find_children(nodes = start_node, edges = edgelist[,1:2])

# subset =
edgelist = edgelist[edgelist$to %in% all_children,]

# whichc ells to include when downsampleing for feature abundance check
cells_to_check = rownames(harmonized_seurat_object@meta.data)[labelmat[,which(apply(labelmat,2,function(x,target){target %in% x},target=start_node))] == start_node]

##########
### Calculate markers between leaf-siblings ?  and merge ?
##########

## subset features to be tested;
if(ncol(harmonized_seurat_object@assays[['RNA']]@data) > 30000){
  set.seed(parameter_list$global_seed)
  subset = sample(colnames(harmonized_seurat_object@assays[['RNA']]@data[,cells_to_check]),size = 30000)
  gene_expr_dataset = harmonized_seurat_object@assays[['RNA']]@data[,subset]
}else{
  gene_expr_dataset = harmonized_seurat_object@assays[['RNA']]@data
}
gene_expr_dataset[gene_expr_dataset != 0] <- 1 # set to 1 for occ
gene_sums = Matrix::rowSums(gene_expr_dataset)
rm(gene_expr_dataset)
genes_to_include = names(gene_sums)[gene_sums > min.cells.feature ]
genes_to_include = genes_to_include[! genes_to_include %in% features_exclude_list]
message("Testing ",length(genes_to_include)," genes as cluster markers")

##########
### findMarkers_tree2
##########

require(tidyselect)

seurat_object =harmonized_seurat_object

# set edglist
edgelist = edgelist[,c("from","to")]

# genes
if(is.null(genes_to_include)){
  genes_to_include = rownames(seurat_object@assays[[assay]]$counts)
}

# add labelmat to seurat
if(any(colnames(labelmat) %in% colnames(seurat_object@meta.data))){
  if(!all(colnames(labelmat) %in% colnames(seurat_object@meta.data))){
    stop("Cannot handle duplicated column names in seurat_object@meta.data and labelmat from mrtree. Either include all labelmat columns beforehand or none at all.")
  }
}else{
  seurat_object@meta.data = cbind(seurat_object@meta.data,labelmat )
}

# init
current_node=start_node
colnames(edgelist) = c("from","to")
label_mat_long = as.data.frame(labelmat) %>% tidyr::pivot_longer(everything(),names_to = "clusterlevel", values_to = "cluster")
edgelist = dplyr::left_join(edgelist,label_mat_long,by=c("to"="cluster")) %>% dplyr::distinct(from,to,clusterlevel)
all_nodes = unique(edgelist[,2])
#all_nodes = all_nodes[!all_nodes %in% c("all","root")]

comparisons_siblings = NULL
comparisons_all = NULL

comparison_types = c("Sibling","All")
return_list=list()
for(comp in comparison_types){
  if(comp == "Sibling"){
    message("Sibling Comparisons")
  }
  if(comp == "All"){
    message("All Comparisons")
  }
  all_res = list()
  for(n in 1:length(all_nodes)){
    message("Running current_node ","(",n,"/",length(all_nodes),")")
    current_node = all_nodes[n]
    parent_node = edgelist$from[edgelist$to==current_node]
    sibling_nodes = edgelist$to[edgelist$from==parent_node & edgelist$to != current_node]
    current_level = edgelist$clusterlevel[edgelist$to==current_node]
    # set ident to current level!
    Idents(seurat_object) = current_level

    #decide what to cal against:
    cluster_1 = current_node
    if(comp == "Sibling"){
      cluster_2 = sibling_nodes
    }
    if(comp == "All"){
      cluster_2 = all_nodes[! all_nodes %in% current_node]
    }
    #calculate markers
    #  current_markers <- tryCatch({
    if(length(cluster_2)>0){
      if(test.use == "wilcox-stratified"){
        current_markers = FindMarkers2.Seurat(object = seurat_object,
                                              ident.1 = cluster_1,
                                              ident.2 = cluster_2,
                                              assay = assay_markers,
                                              features = genes_to_include,
                                              logfc.threshold = logfc.threshold,
                                              slot = assay_slot,
                                              test.use = "VE",
                                              min.pct = min.pct,
                                              min.diff.pct = min.diff.pct,
                                              max.cells.per.ident=max.cells.per.ident,
                                              min.cells.feature = min.cells.feature,
                                              min.cells.group = min.cells.group,
                                              return.thresh = 1,
                                              base = base,
                                              only.pos = only.pos,
                                              latent.vars = batch_var,
                                              genre = "locally-best")

      }else{
        current_markers=Seurat::FindMarkers(object = seurat_object,
                                            ident.1 = cluster_1,
                                            ident.2 = cluster_2,
                                            assay = assay_markers,
                                            logfc.threshold = logfc.threshold,
                                            features = genes_to_include,
                                            slot = assay_slot,
                                            test.use =test.use,
                                            min.pct = min.pct,
                                            min.diff.pct = min.diff.pct,
                                            max.cells.per.ident=max.cells.per.ident,
                                            min.cells.feature = min.cells.feature,
                                            min.cells.group = min.cells.group,
                                            base = base,
                                            only.pos = only.pos)
      }
      current_markers$gene = rownames(current_markers)
      current_markers$cluster_id = cluster_1
      current_markers$comparison = comp
      current_markers$parent = parent_node
      # },error=function(cond) {
      #   message(cond)
      #   # Choose a return value in case of error
      #   return(NULL)
      # })
    }
    # return
    all_res[[n]] = current_markers
  }
  comparisons_main = do.call(rbind,all_res)
  return_list[[paste0("comparisons_",comp)]] = comparisons_main
}



# run marker detection on tree using findMarkers_tree2
message(Sys.time(),": Start marker detection on mrtree" )
all_markers_mrtree_list = findMarkers_tree2(harmonized_seurat_object,
                                            edgelist=edgelist[,1:2],
                                            labelmat = labelmat,
                                            n_cores = n_cores_markers,
                                            assay=assay_markers,
                                            slot=assay_slot,
                                            batch_var=batch_var,
                                            genes_to_include = genes_to_include,
                                            logfc.threshold = logfc.threshold,
                                            min.pct = min.pct,
                                            min.diff.pct = min.diff.pct,
                                            max.cells.per.ident=max.cells.per.ident,
                                            min.cells.feature = min.cells.feature,
                                            min.cells.group = min.cells.group,
                                            base = base,
                                            only.pos = only.pos)

## assumes base 2!
comparisons_siblings = all_markers_mrtree_list$comparisons_Siblings
comparisons_siblings$specificity = ((comparisons_siblings$pct.1 + specificity_base) / (comparisons_siblings$pct.2 + specificity_base)) * comparisons_siblings$avg_log2FC
comparisons_all = all_markers_mrtree_list$comparisons_All
comparisons_all$specificity = ((comparisons_all$pct.1 + specificity_base) / (comparisons_all$pct.2 + specificity_base)) * comparisons_all$avg_log2FC


##########
### Save result
##########

# save objects as txt
data.table::fwrite(comparisons_all,paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_curated",parameter_list$start_node,"_markers_all_",parameter_list$marker_suffix,".tsv"),sep="\t")
data.table::fwrite(comparisons_siblings,paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_curated",parameter_list$start_node,"_markers_siblings_",parameter_list$marker_suffix,".tsv"),sep="\t")

message(Sys.time(),": Finalized marker detection." )
