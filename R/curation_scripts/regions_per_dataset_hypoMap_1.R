

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




target_structure_id = "1097"
target_level = 8
# what regions exist on level 8
mba_ontology_flatten = data.table::fread("data/mba_ontology_flatten.tsv",data.table = FALSE)
mba_ontology_flatten_edges = data.frame(from=as.character(mba_ontology_flatten$parent_structure_id), to = as.character(mba_ontology_flatten$id),stringsAsFactors = F)


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

suggested_region_per_dataset = list("Kim" = c("Ventromedial hypothalamic nucleus","Arcuate hypothalamic nucleus","Median eminence","Dorsomedial nucleus of the hypothalamus","Periventricular hypothalamic nucleus, posterior part","Tuberomammillary nucleus","Ventral premammillary nucleus","Paraventricular hypothalamic nucleus","Periventricular hypothalamic nucleus, intermediate part"), #,"Paraventricular hypothalamic nucleus"
                                    "Moffit" = c("Lateral preoptic area","Medial preoptic area", "Median preoptic nucleus","Medial preoptic nucleus","Periventricular hypothalamic nucleus preoptic part","Periventricular hypothalamic nucleus, anterior part","Ventrolateral preoptic nucleus","Posterodorsal preoptic nucleus","Anteroventral periventricular nucleus","Supraoptic nucleus","Parastrial nucleus","Anterodorsal preoptic nucleus" ,"Paraventricular hypothalamic nucleus","Suprachiasmatic preoptic nucleus" ),
                                    "Campbell" = c("Ventromedial hypothalamic nucleus","Arcuate hypothalamic nucleus","Median eminence","Dorsomedial nucleus of the hypothalamus","Periventricular hypothalamic nucleus, posterior part","Tuberomammillary nucleus","Ventral premammillary nucleus","Periventricular hypothalamic nucleus, intermediate part"), #,"Paraventricular hypothalamic nucleus"
                                    "wenDropSeq" = c("Subparaventricular zone","Suprachiasmatic nucleus","Anterior hypothalamic nucleus","Suprachiasmatic preoptic nucleus","Retrochiasmatic area","Paraventricular hypothalamic nucleus","Periventricular hypothalamic nucleus, preoptic part","Medial preoptic area","Medial preoptic nucleus","Periventricular hypothalamic nucleus, anterior part"),
                                    "wen10x" =  c("Subparaventricular zone","Suprachiasmatic nucleus","Anterior hypothalamic nucleus","Suprachiasmatic preoptic nucleus","Retrochiasmatic area","Paraventricular hypothalamic nucleus","Periventricular hypothalamic nucleus, preoptic part","Medial preoptic area","Medial preoptic nucleus","Periventricular hypothalamic nucleus, anterior part"),
                                    "Rossi" = c("Lateral hypothalamic area","Zona incerta","Tuberal nucleus","Subthalamic nucleus","Parasubthalamic nucleus" ),
                                    "Mickelsen" = c("Lateral hypothalamic area","Zona incerta","Tuberal nucleus","Subthalamic nucleus","Parasubthalamic nucleus" ),
                                    "Flynn" =c("Posterior hypothalamic nucleus","Tuberomammillary nucleus","Supramammillary nucleus","Medial mammillary nucleus","Lateral mammillary nucleus","Ventral premammillary nucleus", "Dorsal premammillary nucleus","Periventricular hypothalamic nucleus, posterior part"),
                                    "Chen" = relevant_structures,
                                    "kimDev" = relevant_structures,
                                    "RomanovDev" = relevant_structures,
                                    "Lee_Idol" = relevant_structures,
                                    "Mousebrainorg" = relevant_structures)




