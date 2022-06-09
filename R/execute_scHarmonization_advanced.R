##########
### Load
##########

# Important: add your curated name as file !

# source("R/harmonization_functions.R")
singularity_path = "~/Documents/r_scvi_015.simg"

# direct output and logs to some files on the local filesystem:
# where to store temporary json's with params for jobs:
param_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/hypoMap_v2_params/"
# where to save log files --> use this path in the slurm.sh files!
log_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/hypoMap_v2_slurmlogs/"

# load json file with all other information
params_harmonization = jsonlite::read_json("data/parameters_harmonization_v2_4.json")
# if some fields are lists --> unlist
params_harmonization = lapply(params_harmonization,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# update suffix name:
params_harmonization$new_name_suffix = paste0(params_harmonization$new_name_suffix,"_curated")

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
### [1] Clustering for hierachical tree after curation
##########

# set additional parameters for scvi
param_set = params_harmonization
## full clustering extra
# param_set$target_clusterN = 600
# param_set$start_res = 28
# param_set$end_res = 50
# param_set$step_size = 1
# param_set$include_low_res = FALSE
# param_set$min_cells_valid = 5
# param_set$additional_clustering_suffix = "_extracluster28"
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
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --kill-on-invalid-dep=yes python/run_Python_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_1 = stringr::str_remove(output_message,pattern = "Submitted batch job ")

##########
### [2] Define cluster levels for mr_tree
##########


## Run selection script and save file in

# parameter_list$clusters_for_mrtree_file


##########
### [3] Hierachical tree
##########

# set params
param_set = params_harmonization
param_set$marker_suffix = "raw"
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
dependency_ids = c(slurm_id_1)
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes R/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_3 = stringr::str_remove(output_message,pattern = "Submitted batch job ")

##########
### [4] Hierachical tree cluster markers
##########

# set params
param_set = params_harmonization
param_set$n_cores_markers = 4
param_set$marker_suffix = "raw"
param_set$start_node = "K2-0" # "all" for everything --> then 4b is not necessary
#param_set$marker_suffix = "nonneuron"
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
slurm_id_4 = stringr::str_remove(output_message,pattern = "Submitted batch job ")

##########
### [4b] Hierachical tree cluster markers ---> for HypoMap where we split into two subsets (neuron and non neuron)
##########

# set params
param_set = params_harmonization
param_set$n_cores_markers = 4
param_set$marker_suffix = "raw"
param_set$start_node = "K2-1" # "all" for everything
#param_set$marker_suffix = "nonneuron"
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
slurm_id_4b = stringr::str_remove(output_message,pattern = "Submitted batch job ")

##########
### [5] Prune tree
##########

# set params
param_set = params_harmonization
param_set$marker_suffix = "raw"
# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"mrtree_pruning_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "R/run_scripts/mrtree_pruning.R"
# set sbatch params:
jobname = paste0("mrtree_pruning_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = c(slurm_id_4,slurm_id_4b) ## Might need to adjust this !
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes R/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_5 = stringr::str_remove(output_message,pattern = "Submitted batch job ")


##########
### [6] Hierachical tree cluster markers re-run
##########

# TODO: only re run on clusters that actually changed ?

# set params
param_set = params_harmonization
param_set$n_cores_markers = 4
param_set$marker_suffix = "pruned"
param_set$start_node = "C2-1" # "all" for everything

#param_set$marker_suffix = "nonneuron"
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
dependency_ids = c(slurm_id_5)
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes R/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_6 = stringr::str_remove(output_message,pattern = "Submitted batch job ")

##########
### [6b] Hierachical tree cluster markers re-run
##########

# TODO: only re run on clusters that actually changed ?

# set params
param_set = params_harmonization
param_set$n_cores_markers = 4
param_set$marker_suffix = "pruned"
param_set$start_node = "C2-2" # "all" for everything

#param_set$marker_suffix = "nonneuron"
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
dependency_ids = c(slurm_id_5)
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes R/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_6b = stringr::str_remove(output_message,pattern = "Submitted batch job ")


##########
### Read parameters fro annotation
##########

# load json file with all other information
params_annotation= jsonlite::read_json("data/parameters_annotation_v2_1.json")
# if some fields are lists --> unlist
params_annotation = lapply(params_annotation,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

##########
### [7] cluster annotation
##########

# set params
param_set = params_annotation
param_set$manual_names_annotation = as.list(param_set$manual_names_annotation)
# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"mrtree_annotation_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "R/run_scripts/mrtree_annotation.R"
# set sbatch params:
jobname = paste0("mrtree_annotation_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = c(slurm_id_6,slurm_id_6b) ## Might need to adjust this !
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes R/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_7 = stringr::str_remove(output_message,pattern = "Submitted batch job ")

##########
### [8] scMRMR
##########

# TODO: still need to decide whether to include this part
# probably not for now !

##########
### [9] propagate author cell types
##########

# TODO
# maybe use approach from mapscvi to make a column per annotated dataset

##########
### [10] regionPrediction with scCoco
##########

# set params
param_set = params_annotation
# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"region_prediction_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "R/run_scripts/region_prediction.R"
# set sbatch params:
jobname = paste0("region_prediction_annotation_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = c(slurm_id_6,slurm_id_6b) ## Might need to adjust this !
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes R/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_7 = stringr::str_remove(output_message,pattern = "Submitted batch job ")

##########
### [11] final curation
##########

# TODO: need to gather everything
# TODO: clean up final object
# todo: export final objects


