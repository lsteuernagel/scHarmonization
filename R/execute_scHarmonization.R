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
### [0] prepare
##########

# set params
param_set = params_harmonization
# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"prepare_harmonization_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "R/run_scripts/prepare_harmonization.R"
# set sbatch params:
jobname = paste0("prepare_harmonization_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = ""
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --kill-on-invalid-dep=yes R/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_0 = stringr::str_remove(output_message,pattern = "Submitted batch job ")

##########
### [1] Integrate with scVI
##########

# set additional parameters for scvi
param_set = params_harmonization
# param_set$output_folder = harmonization_scvi_folder
# param_set$data_filepath_full = paste0(params_harmonization$harmonization_folder_path,param_set$new_name_suffix,".h5ad")
# param_set$feature_set_file = paste0(params_harmonization$harmonization_folder_path,"features/feature_sets.json")
# param_set$hvgs_set_name = paste0(params_harmonization$assay_name,".log.", "vst", ".split_", params_harmonization$batch_var, ".features.",params_harmonization$feature_set_size)

# make unique id:
job_id=digest::digest(param_set)
param_set$job_id = job_id
# write to JSON as transfer file
param_file = paste0(param_path,"harmonization_scvi_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# execute job
script_path = "python/integrate_scVI_v015.py"
# set sbatch params:
jobname = paste0("harmonization_scvi_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = c(slurm_id_0)
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes python/run_Python_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_1 = stringr::str_remove(output_message,pattern = "Submitted batch job ")

##########
### [2] Basic harmonization
##########

# set params
param_set = params_harmonization
# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"basic_harmonization_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "R/run_scripts/basic_harmonization.R"
# set sbatch params:
jobname = paste0("basic_harmonization_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = c(slurm_id_1)
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes R/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_2 = stringr::str_remove(output_message,pattern = "Submitted batch job ")


##########
### [3] Initial clustering using python leiden
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
### [4] Initial marker detection
##########


# set params
param_set = params_harmonization
# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"basic_markers_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "R/run_scripts/basic_marker_detection.R"
# set sbatch params:
jobname = paste0("basic_harmonization_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = c(slurm_id_3)
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes R/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_4 = stringr::str_remove(output_message,pattern = "Submitted batch job ")


##########
### [5] Manual curation
##########



##########
### [6] Hierachical tree
##########


##########
### [7] Hierachical tree cluster markers + pruning
##########


##########
### [8] Hierachical tree cluster markers re-run
##########


