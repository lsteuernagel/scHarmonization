##########
### read_embedding
##########

#' Load an emebedding with cells x lowDims from flatfile, ensuring consistency with a Seurat object (or metadata only for faster usage)
#' @param filename_withpath filepath
#' @param seurat_object seuratobject associated with current embedding. If specified metadata does not have to be set explicitly.
#' @param seurat_object_metadata metadata only of seuratobject associated with current embedding
#' @return

read_embedding = function(filename_withpath,seurat_object=NULL,seurat_object_metadata=NULL){

  #get metadata
  if(!is.null(seurat_object_metadata)){
    metadata = seurat_object_metadata
  }else{
    if(!is.null(seurat_object)){
      metadata = seurat_object@meta.data
    }else{
      stop("Please provide either a dataframe with metadata or a seurat object with metadata that can be exctracted!")
    }
  }
  # load
  current_embedding = data.table::fread(filename_withpath,data.table = F)
  # use first col as rownames
  if(is.character(current_embedding[,1])){
    rnames = current_embedding[,1]
    current_embedding = current_embedding[,2:ncol(current_embedding)]
    rownames(current_embedding)=rnames
    # reorder to align with rest of object
    if(any(is.na(match(rownames(metadata),rownames(current_embedding))))){
      message("Found ",length(rnames)," rows in embedding and ",length(rownames(metadata))," rows in metadata.")
      stop("Cell names from loaded reduction and new object are not matching exactly. Stopping import.")
    }
    current_embedding = current_embedding[match(rownames(metadata),rownames(current_embedding)),]
  }else{
    warning("First column of loaded file is not of type character, using rownames of metadata as rownames of added reduction. This can induce bugs if the order changed due to split/merge of the Seurat object!")
    rownames(current_embedding) = rownames(metadata)
  }
  return(current_embedding)
}

##########
### clear_clustering
##########

#' Eliminate small clusters from a vector of labels substituting with the labels of NN
#' @param x vector of labels
#' @param min_cells minimum number of cells to keep cluster
#' @param nn_idx matrix of cells x k NN idx --> e.g. output of annoy or rann
#' @return updated vector of labels

clear_clustering = function(x,min_cells,nn_idx){
  x = as.character(x)
  new_x=x
  # which clusters are too small ?
  small_clusters = names(table(x)[table(x) < min_cells])
  # go through small clusters and move cells to laregr clusters based on neighbors
  if(length(small_clusters)>0){
    #message("Removing ",length(small_clusters)," clusters")
    for(i in 1:length(small_clusters)){
      current_cluster = small_clusters[i]
      which_idx = which(x == current_cluster)
      # get idx for k NN
      neighbor_idx = nn_idx[which_idx,]
      # substitute with cluster labels
      if(length(which_idx)>1){
        neighbor_clusters = apply(neighbor_idx,1,function(z,cluster_vector){return(cluster_vector[z])},cluster_vector=x)
        # extract that most common label in neighbors
        clusters_vote = apply(neighbor_clusters,2,function(z,exclude_clusters){
          return(names(sort(table(z),decreasing = TRUE))[! names(sort(table(z),decreasing = TRUE)) %in% exclude_clusters][1])
        },exclude_clusters=small_clusters)
      }else{
        neighbor_clusters = x[neighbor_idx]
        clusters_vote = names(sort(table(neighbor_clusters),decreasing = TRUE))[! names(sort(table(neighbor_clusters),decreasing = TRUE)) %in% small_clusters][1]
      }
      # overwrite cluster label
      new_x[which_idx]=clusters_vote
    }
  }

  return(new_x)
}


##########
### annotate_tree
##########


#' Add gene absed annotation to tree labels
#' @param seurat_object seurat_object to call FindMarkers
#' @param edgelist edgelist from mrtree corresponding to labelmat
#' @param labelmat dataframe of clusterlabels per cell with one column per level
#' @param markers_comparisons_all cluster markers for all clusters in labelmat. expects format from run_marker_detection.R
#' @param markers_comparisons_siblings cluster markers vs sibling for all clusters in labelmat. expects format from run_marker_detection.R
#' @param manual_names a named vector with names of clusters in edgelist --> will overwrite with manual names
#' @param overwrite_with_manual whether to overwrite the full tree until this point with the manual name or just append similar to the best gene. defaults to FALSE
#' @param manual_exclude_genes a set of genes that should not be included in anno
#' @param max_pval_adj max pavlue allowed for a marker to be considered
#' @param min_specificity min_specificity allowed for a marker to be considered
#' @param min_specificity_sibling_children maximum specificity allowed for a marker in sibling comparisons of children to be kept --> strong markers between children are probably not good for whole cluster
#' @param min_pct2_score when calculating the specificty score add this pseudo-count to pct.2 to avoid favouring of highly specific but low abundant genes
#' @param scale_preferred factor to scale scoe of preferred genes
#' @param preferred_genes a vector with genes that should be preferred
#' @param min_cells min_cells
#' @param limit_factor scale sibling or children_siblings scores down to be not larger than limit_factor*global_score
#' @param max_score_siblings_children todo: add explanation
#' @return dataframe with new labels for each cluster in edgelist$to

annotate_tree = function(edgelist,labelmat,markers_comparisons_all,markers_comparisons_siblings,preferred_genes=character(0),manual_names=c(),overwrite_with_manual=FALSE,manual_exclude_genes=character(0),
                         max_pval_adj=0.0001, min_specificity = 0.5,min_specificity_sibling_children=10,scale_preferred=2,min_pct2_score=0.01,min_cells=20,limit_factor=5,max_score_siblings_children=20,
                         reverse_order = FALSE){
  # init edgelist and all_nodes
  edgelist = edgelist[,c("from","to","level")]
  cluster_levels = as.data.frame(labelmat) %>% tidyr::pivot_longer(everything(),names_to = "clusterlevel", values_to = "cluster")  %>% dplyr::group_by(cluster) %>%
    dplyr::add_count(name="ncells") %>% dplyr::distinct(clusterlevel,cluster,ncells)
  edgelist = dplyr::left_join(edgelist,cluster_levels,by=c("to"="cluster")) %>% dplyr::arrange(level)
  all_nodes = unique(edgelist$to)

  # annotation
  edgelist = edgelist %>% dplyr::arrange(level) # ensure top down order for tree traversal with for loop!+

  # list to store intermediate results
  annotation_list = list()
  descriptive_markers_list = list()

  message("Running annotation")
  message("Using ",length(manual_names)," manually provided names to overwrite annotation for specific clusters.")

  # for each node in edgelist:
  for(n in 1:length(all_nodes)){

    # get information
    current_node = all_nodes[n]
    parent_node = edgelist$from[edgelist$to==current_node]
    sibling_nodes = edgelist$to[edgelist$from==parent_node & edgelist$to != current_node]
    children_nodes = scUtils::find_children(current_node,edgelist)
    direct_children_nodes = edgelist$to[edgelist$from==current_node]
    current_level = edgelist$clusterlevel[edgelist$to==current_node]

    message("n: ",n," ",current_node)
    # get specific genes
    potential_descriptive_markers =markers_comparisons_all %>% dplyr::filter(cluster_id == current_node) %>% dplyr::arrange(desc(specificity)) %>%
      dplyr::filter(p_val_adj< max_pval_adj & specificity > min_specificity)
    # calculate score as adjusted specificity with a minimum value on pct.2 to put more emphasis on pct.1 (abundant markers)
    potential_descriptive_markers$pct.2_min = potential_descriptive_markers$pct.2
    potential_descriptive_markers$pct.2_min[potential_descriptive_markers$pct.2_min < min_pct2_score] = min_pct2_score
    potential_descriptive_markers$score = (potential_descriptive_markers$pct.1 / potential_descriptive_markers$pct.2_min) * potential_descriptive_markers$avg_log2FC

    # initiate vector with genes that should be excluded!
    exclude_genes=c(manual_exclude_genes)
    # message("nrow(potential_descriptive_markers) 1: ",nrow(potential_descriptive_markers))

    # eliminate parent and sibling names
    reserved_genes=c()
    parent_name=""
    if(length(annotation_list)>0){
      parent_name = annotation_list[[parent_node]]
      if(parent_node %in% names(annotation_list)){
        reserved_genes = c(reserved_genes,stringr::str_split(parent_name,pattern = "\\.")[[1]])
      }
      if(length(sibling_nodes)>0){
        for(s in 1:length(sibling_nodes)){
          if(sibling_nodes[s] %in% names(annotation_list)){
            reserved_genes = c(reserved_genes,stringr::str_split(annotation_list[[sibling_nodes[s]]],pattern = "\\.")[[1]])
          }
        }
      }
    }else{
      reserved_genes=c()
    }
    # if a gene is in preferred_genes
    if(nrow(potential_descriptive_markers)>0){
      potential_descriptive_markers$score[potential_descriptive_markers$gene %in% preferred_genes] = potential_descriptive_markers$score[potential_descriptive_markers$gene %in% preferred_genes] * scale_preferred
    }
    ####
    # ranks in siblings and in children to compare with global rank in own markers
    ####
    if(nrow(potential_descriptive_markers)>0){
      # get sibling_markers of children and get inveretd ranks
      sibling_markers_children =markers_comparisons_siblings %>% dplyr::filter(cluster_id %in% children_nodes) %>% dplyr::arrange(desc(specificity))%>%
        dplyr::filter(p_val_adj<max_pval_adj & specificity > min_specificity_sibling_children)
      if(nrow(sibling_markers_children)>0){
        # calculate score and only keep highest children cluster
        sibling_markers_children$pct.2_min = sibling_markers_children$pct.2
        sibling_markers_children$pct.2_min[sibling_markers_children$pct.2_min < min_pct2_score] = min_pct2_score
        sibling_markers_children$score = (sibling_markers_children$pct.1 / sibling_markers_children$pct.2_min) * sibling_markers_children$avg_log2FC
        sibling_markers_children = sibling_markers_children %>% dplyr::group_by(gene) %>% dplyr::filter(score == max(score)) %>%
          dplyr::distinct(gene,.keep_all =TRUE)
        # filter by specificity
        sibling_markers_children =sibling_markers_children %>% dplyr::filter(! gene %in% exclude_genes) %>% dplyr::arrange(desc(score)) #  & score > min_specificity
        # join score with potential_descriptive_markers
        if(nrow(sibling_markers_children)>0){
          sibling_markers_children$score_siblings_children = sibling_markers_children$score
          potential_descriptive_markers = dplyr::left_join(potential_descriptive_markers,sibling_markers_children[,c("gene","score_siblings_children")],by="gene")
          potential_descriptive_markers$score_siblings_children[is.na(potential_descriptive_markers$score_siblings_children)] = 0 # if not a marker within children gene gets rank 1
        }else{
          potential_descriptive_markers$score_siblings_children = 0
        }
      }else{
        potential_descriptive_markers$score_siblings_children = 0
      }
      #message("nrow(potential_descriptive_markers): ",nrow(potential_descriptive_markers))
      # check siblings and filter to genes that are also markers to siblings!
      if(length(sibling_nodes)>0){
        # repeat most steps for sibling markers
        sibling_markers =markers_comparisons_siblings %>% dplyr::filter(cluster_id == current_node) %>% dplyr::arrange(desc(specificity)) %>% dplyr::filter(p_val_adj<max_pval_adj & specificity > min_specificity)
        if(nrow(sibling_markers)>0){
          # calculate score, scale and filter
          sibling_markers$pct.2_min = sibling_markers$pct.2
          sibling_markers$pct.2_min[sibling_markers$pct.2_min < min_pct2_score] = min_pct2_score
          sibling_markers$score = (sibling_markers$pct.1 / sibling_markers$pct.2_min) * sibling_markers$avg_log2FC
          sibling_markers$score[sibling_markers$gene %in% preferred_genes] = sibling_markers$score[sibling_markers$gene %in% preferred_genes] * scale_preferred
          sibling_markers =sibling_markers %>% dplyr::filter(! gene %in% exclude_genes & score > min_specificity) %>% dplyr::arrange(desc(score))
          # join
          if(nrow(sibling_markers)>0){
            sibling_markers$score_siblings = sibling_markers$score
            potential_descriptive_markers = dplyr::left_join(potential_descriptive_markers,sibling_markers[,c("gene","score_siblings")],by="gene")
            potential_descriptive_markers$score_siblings[is.na(potential_descriptive_markers$score_siblings)] = 0
            # potential_descriptive_markers$score_siblings = potential_descriptive_markers$score_siblings * potential_descriptive_markers$children_fraction
          }else{
            potential_descriptive_markers$score_siblings = 0
          }
        }else{
          potential_descriptive_markers$score_siblings =0
        }
      }else{
        # fallback if nor sibling markers are available
        potential_descriptive_markers$score_siblings = 0
      }

      # get ranks by score
      if(nrow(potential_descriptive_markers)>0){
        # limit score_siblings_children and score_siblings to omit strong outliers
        potential_descriptive_markers$score_siblings_children[potential_descriptive_markers$score_siblings_children > limit_factor*potential_descriptive_markers$score] = limit_factor*potential_descriptive_markers$score[potential_descriptive_markers$score_siblings_children > limit_factor*potential_descriptive_markers$score]
        potential_descriptive_markers$score_siblings[potential_descriptive_markers$score_siblings > limit_factor*potential_descriptive_markers$score] = limit_factor*potential_descriptive_markers$score[potential_descriptive_markers$score_siblings > limit_factor*potential_descriptive_markers$score]

        # calculate average score as: score + score_siblings - score_siblings_children
        potential_descriptive_markers$avg_score = potential_descriptive_markers$score + potential_descriptive_markers$score_siblings  - potential_descriptive_markers$score_siblings_children
        potential_descriptive_markers$mult_score = potential_descriptive_markers$score * (1+potential_descriptive_markers$score_siblings) /  (1+potential_descriptive_markers$score_siblings_children)


        # prepare for selection: exclude all genes in exclude_genes vector, filter score again after all adjustements and arrange by score
        potential_descriptive_markers =potential_descriptive_markers %>%
          dplyr::filter(! gene %in% exclude_genes & avg_score > min_specificity  & score_siblings_children < max_score_siblings_children) %>% dplyr::arrange(desc(avg_score))
      }
    }
    # take top result as cluster descriptive gene
    if(nrow(potential_descriptive_markers)>0){
      # don't select sibling or parent gene
      potential_descriptive_markers_sel =potential_descriptive_markers %>% dplyr::filter(! gene %in% reserved_genes) %>% dplyr::arrange(desc(avg_score))
      #reserved_genes
      # select top gene
      name_gene = potential_descriptive_markers_sel$gene[1]
      # if node has no siblings: don't need annotation!
      if(length(sibling_nodes)==0){
        name_gene="m"
      }
    }else{
      # if merge take name from sibling!
      name_gene = "problematic"
    }
    # check manual overwrite vector
    if(current_node %in% names(manual_names)){
      name_gene = manual_names[current_node]
    }
    # update full cluster_name
    if(overwrite_with_manual & current_node %in% names(manual_names)){
      cluster_name = name_gene
      message(current_node,": name: ",cluster_name)
    }else{
      cluster_name = paste0(parent_name,".",name_gene)
      message(current_node,": name: ",cluster_name)
    }

    # add to list
    annotation_list[[current_node]] = cluster_name
    if(nrow(potential_descriptive_markers)){
      descriptive_markers_list[[current_node]] = potential_descriptive_markers %>% dplyr::select(cluster_id,parent,gene,p_val_adj,avg_log2FC,pct.1,pct.2,specificity,score,score_siblings,score_siblings_children,avg_score,mult_score)
    }
  }
  # also save a list of the good marker genes
  descriptive_markers_df = as.data.frame(do.call(rbind,descriptive_markers_list))

  # bind results and format
  annotation_df = as.data.frame(do.call(rbind,annotation_list))
  colnames(annotation_df) = "V1"
  #print(head(annotation_df))
  annotation_df = annotation_df %>% dplyr::rename(Map_CellType=V1) %>% dplyr::mutate(cluster_id = rownames(annotation_df))
  if(reverse_order){
    annotation_df$Map_CellType = sapply(annotation_df$Map_CellType,function(x){paste0(rev(strsplit(x,split = "\\.")[[1]]),collapse = ".")})
    annotation_df$clean_names = sub("\\.$","",gsub("m\\.","",annotation_df$Map_CellType))#sub(".","",gsub("\\.d","",annotation_df$Map_CellType))
  }else{
    annotation_df$clean_names = sub("\\.","",gsub("\\.m","",annotation_df$Map_CellType))
  }
  annotation_df = dplyr::left_join(annotation_df,cluster_levels,by=c("cluster_id"="cluster"))
  annotation_df$clean_names[is.na(annotation_df$clean_names)]="all"
  # named_edgelist
  # named_edgelist = dplyr::left_join(edgelist,annotation_df[,c("clean_names","cluster_id")],by=c("from"="cluster_id"))
  # named_edgelist = dplyr::left_join(named_edgelist,annotation_df[,c("clean_names","cluster_id")],by=c("to"="cluster_id"))
  # named_edgelist$from = named_edgelist$clean_names.x
  # named_edgelist$from[is.na(named_edgelist$from)] = "all"
  # named_edgelist$to = named_edgelist$clean_names.y
  # named_edgelist = named_edgelist[,c("from","to","level")]
  # return
  message("Returning results")
  return(list(annotation_df = annotation_df,descriptive_markers_df = descriptive_markers_df))
}


##########
### findMarkers_tree2
##########

#' Walk through a mrtree / hierchical clustering tree using an edgelist and calculate Marker genes using Seurat's FindMarkers or a stratified version:
#' https://github.com/KChen-lab/stratified-tests-for-seurat
#' Two results: Vs sibling and Vs all . Runs in parallel to be faster
#' @param seurat_object seurat_object to call FindMarkers
#' @param edgelist minimum number of cells to keep cluster
#' @param labelmat
#' @param n_cores doParallel cores
#' @param use_stratified use stratfied
#' @param batch_var batch_var for stratfied mode
#' @param assay
#' @param slot
#' @param ... passed to detection functions
#' @return updated vector of labels

findMarkers_tree2 = function(seurat_object,edgelist,labelmat,n_cores=1,use_stratified=TRUE,test.use="wilcox-stratified",batch_var="Batch_ID",genes_to_include=NULL,assay="RNA",slot="data",...){

  require(doParallel)
  require(tidyselect)
  require(foreach)

  # info:
  message(n_cores)

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
  current_node="all"
  colnames(edgelist) = c("from","to")
  label_mat_long = as.data.frame(labelmat) %>% tidyr::pivot_longer(everything(),names_to = "clusterlevel", values_to = "cluster")
  edgelist = dplyr::left_join(edgelist,label_mat_long,by=c("to"="cluster")) %>% dplyr::distinct(from,to,clusterlevel)
  all_nodes = unique(edgelist[,2])
  #all_nodes = all_nodes[!all_nodes %in% c("all","root")]

  comparisons_siblings = NULL
  comparisons_all = NULL

  registerDoParallel(cores=n_cores)
  comparison_types = c("Sibling","All")
  return_list=list()
  for(comp in comparison_types){
    if(comp == "Sibling"){
      message("Sibling Comparisons")
    }
    if(comp == "All"){
      message("All Comparisons")
    }
    comparisons_main <- foreach(n = 1:length(all_nodes),.errorhandling = 'remove', .combine='rbind') %dopar% {
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
        cluster_2 = cluster_2[cluster_2 %in% edgelist$to[edgelist$clusterlevel == current_level]] # need to ensure that only current level is used!
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
          if(base==2){
            colnames(current_markers)[colnames(current_markers)=="avg_logFC"] = "avg_log2FC"
          }
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
        current_markers
      }else{
        NULL
      }
      # },error=function(cond) {
      #   message(cond)
      #   # Choose a return value in case of error
      #   return(NULL)
      # })

      # return
      #current_markers
    }
    return_list[[paste0("comparisons_",comp)]] = comparisons_main
  }
  # return
  return(return_list)
}
