##########
### Load
##########

# source("R/harmonization_functions.R")
singularity_path = "~/Documents/r_scvi_015.simg"

# direct output and logs to some files on the local filesystem:
# where to store temporary json's with params for jobs:
param_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/hypoMap_v2_params/"
# where to save log files --> use this path in the slurm.sh files!
log_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/hypoMap_v2_slurmlogs/"

# load json file with all other information
params_harmonization = jsonlite::read_json("data/parameters_harmonization_v2_1_test.json")
# if some fields are lists --> unlist
params_harmonization = lapply(params_harmonization,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

### try to creat dir if necessary:
system(paste0("mkdir -p ",paste0(param_path)))
system(paste0("mkdir -p ",paste0(log_path)))
system(paste0("mkdir -p ",paste0(params_harmonization$harmonization_folder_path)))


# define helper function
writeList_to_JSON = function (list_with_rows, filename){
  jsonfile = jsonlite::toJSON(list_with_rows, pretty = TRUE, auto_unbox = TRUE, digits = NA)
  writeLines(jsonfile, con = paste0(filename))
}


##########
### [6] Updated advanced object
##########

## also re-run

# run umap and save model
# message(Sys.time(),": Build UMAP with ",parameter_list$k_param," n.neighbors ..." )
# harmonized_seurat_object = RunUMAP(harmonized_seurat_object,
#                                    reduction = parameter_list$integration_name,
#                                    seed.use= parameter_list$global_seed,
#                                    dims=1:ncol(harmonized_seurat_object@reductions[[parameter_list$integration_name]]@cell.embeddings),
#                                    reduction.name=paste0("umap_",parameter_list$integration_name),
#                                    reduction.key = paste0("umap_",parameter_list$integration_name),
#                                    verbose=F,
#                                    n.neighbors = parameter_list$k_param,
#                                    return.model = TRUE)
#
# # run seurat SNN annoy
# message(Sys.time(),": Build SNN with ",parameter_list$k_param," n.neighbors ..." )
# harmonized_seurat_object = FindNeighbors(harmonized_seurat_object,
#                                          reduction=parameter_list$integration_name,
#                                          dims = 1:ncol(harmonized_seurat_object@reductions[[parameter_list$integration_name]]@cell.embeddings),
#                                          k.param = parameter_list$k_param,
#                                          nn.method="annoy",
#                                          annoy.metric=parameter_list$dist_type,
#                                          graph.name = paste0("SNN_",parameter_list$integration_name), verbose=TRUE)


##########
### [7] Clustering for hierachical tree after curation
##########

# set additional parameters for scvi
param_set = params_harmonization
# make unique id:
job_id=digest::digest(param_set)
param_set$job_id = job_id
# write to JSON as transfer file
param_file = paste0(param_path,"leiden_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# execute job
script_path = "python/basic_leiden_clustering.py"
# set sbatch params:
jobname = paste0("leiden_scanpy_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = c(slurm_id_2)
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes python/run_Python_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_3 = stringr::str_remove(output_message,pattern = "Submitted batch job ")

##########
### [8] Define cluster levels for mr_tree
##########


## Run selection script and save file in

# parameter_list$clusters_for_mrtree_file,".rds"


##########
### [9] Hierachical tree
##########

# set params
param_set = params_harmonization
# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"mrtree_construction_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "R/run_scripts/mrtree_construction.R"
# set sbatch params:
jobname = paste0("mrtree_construction_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = c(slurm_id_3)
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes R/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_9 = stringr::str_remove(output_message,pattern = "Submitted batch job ")

##########
### [10] Hierachical tree cluster markers + pruning
##########

# set params
param_set = params_harmonization
# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"mrtree_markers_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "R/run_scripts/mrtree_marker_detection.R"
# set sbatch params:
jobname = paste0("mrtree_markers_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = c(slurm_id_3)
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes R/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_9 = stringr::str_remove(output_message,pattern = "Submitted batch job ")

##########
### [11] Prune tree
##########

# set params
param_set = params_harmonization
# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"mrtree_construction_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "R/run_scripts/mrtree_construction.R"
# set sbatch params:
jobname = paste0("mrtree_construction_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = c(slurm_id_3)
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes R/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_9 = stringr::str_remove(output_message,pattern = "Submitted batch job ")


##########
### [12] Hierachical tree cluster markers re-run
##########

# TODO: only re run on clusters that actually changed!

# set params
param_set = params_harmonization
# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"mrtree_markers_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "R/run_scripts/mrtree_marker_detection.R"
# set sbatch params:
jobname = paste0("mrtree_markers_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = c(slurm_id_3)
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes R/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_9 = stringr::str_remove(output_message,pattern = "Submitted batch job ")

##########
### [13] scMRMR
##########

# TODO: still need to decide whether to include this part


##########
### [14] propagate author cell types
##########

# TODO
# maybe use approach from mapscvi to make a column per annotated dataset

##########
### [15] regionPrediction with scCoco
##########

#TODO: add region_prediction_new and scCoco !


##########
### [16] final curation
##########

# TODO: need to gather everything
# TODO: clean up final object
# todo: export final objects


