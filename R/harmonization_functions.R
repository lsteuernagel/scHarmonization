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
