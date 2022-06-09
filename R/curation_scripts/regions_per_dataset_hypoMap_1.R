
# This script is a helper script to create a json with ABA ish reion names per cluster which can be used as an additional information when predicting regions of origin per dataset


# helper function
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

# what regions exist on level 8
mba_ontology_flatten = data.table::fread("data/mba_ontology_flatten.tsv",data.table = FALSE)
mba_ontology_flatten_edges = data.frame(from=as.character(mba_ontology_flatten$parent_structure_id), to = as.character(mba_ontology_flatten$id),stringsAsFactors = F)

# need to give parent structure id (1097 for hypothalamus)
target_structure_id = "1097"
# need to give target level
target_level = 8

# find all children of target structure
structures_selected = find_children(nodes=target_structure_id,edges=mba_ontology_flatten_edges) # todo: move find_children to scUtils !?
# what regions exist above target level
structures_above_target_level = structures_selected[structures_selected %in% mba_ontology_flatten$id[mba_ontology_flatten$st_level < as.numeric(target_level)]]
# add all children of these
children_structures_above_target_level =  mba_ontology_flatten_edges$to[mba_ontology_flatten_edges$from %in%structures_above_target_level]
# select all structures on target level + potentially some above
structures_on_target_level = children_structures_above_target_level[children_structures_above_target_level %in% mba_ontology_flatten$id[mba_ontology_flatten$st_level >= as.numeric(target_level)]]
# sort by name to print
relevant_structures = sort(mba_ontology_flatten$name[mba_ontology_flatten$id %in% structures_on_target_level])
relevant_structures

table(curated_seurat_object@meta.data$Dataset)

# shortcuts for HypoMap v2:
mediobasal_hypo = c("Ventromedial hypothalamic nucleus","Arcuate hypothalamic nucleus","Median eminence","Dorsomedial nucleus of the hypothalamus","Periventricular hypothalamic nucleus, posterior part",
                    "Ventral premammillary nucleus","Paraventricular hypothalamic nucleus","Periventricular hypothalamic nucleus, intermediate part","Tuberal nucleus","Anterior hypothalamic nucleus") # added peri and paraventricular hypo because they are typicall also included
lateral_hypo = c("Lateral hypothalamic area","Zona incerta","Tuberal nucleus","Subthalamic nucleus","Parasubthalamic nucleus","Preparasubthalamic nucleus","Retrochiasmatic area")
scn_hypo = c("Subparaventricular zone","Suprachiasmatic nucleus","Anterior hypothalamic nucleus","Suprachiasmatic preoptic nucleus","Retrochiasmatic area","Paraventricular hypothalamic nucleus","Periventricular hypothalamic nucleus, preoptic part","Medial preoptic area","Medial preoptic nucleus","Periventricular hypothalamic nucleus, anterior part")

suggested_region_per_dataset = list(Affinati10x = mediobasal_hypo,
                                    Anderson10x = c(mediobasal_hypo,"Lateral hypothalamic area"),
                                    CampbellDropseq = c(mediobasal_hypo,"Suprachiasmatic nucleus"),
                                    ChenDropseq = relevant_structures,
                                    Dowsett10xnuc = relevant_structures,
                                    Flynn10x = c("Posterior hypothalamic nucleus","Tuberomammillary nucleus","Supramammillary nucleus","Medial mammillary nucleus","Lateral mammillary nucleus","Ventral premammillary nucleus", "Dorsal premammillary nucleus","Periventricular hypothalamic nucleus, posterior part"),
                                    Kim10x = mediobasal_hypo,
                                    KimDev10x =  relevant_structures,
                                    LeeDropseq = relevant_structures,
                                    Mickelsen10x = lateral_hypo,
                                    Moffit10x = c("Lateral preoptic area","Medial preoptic area", "Median preoptic nucleus","Medial preoptic nucleus","Periventricular hypothalamic nucleus preoptic part",
                                                  "Periventricular hypothalamic nucleus, anterior part","Ventrolateral preoptic nucleus","Posterodorsal preoptic nucleus","Anteroventral periventricular nucleus",
                                                  "Supraoptic nucleus","Parastrial nucleus","Anterodorsal preoptic nucleus" ,"Paraventricular hypothalamic nucleus","Suprachiasmatic preoptic nucleus" ),
                                    Morris10x = scn_hypo,
                                    Mousebrainorg10x = relevant_structures,
                                    RomanovDev10x = relevant_structures,
                                    RossiDropseq = lateral_hypo,
                                    Rupp10x = c(mediobasal_hypo,"Lateral hypothalamic area"),
                                    Wen10x = scn_hypo,
                                    wenDropseq = scn_hypo)

## save to file
scUtils::writeList_to_JSON(suggested_region_per_dataset,filename = "data/hypoMap_suggested_region_per_dataset.json")






