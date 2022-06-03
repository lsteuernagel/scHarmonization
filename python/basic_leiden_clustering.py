
# Import relevant modules
import pandas as pd
import numpy as np
import scanpy as sc
import louvain
import igraph
import re
import os
import sys
import json
import gc

print(" Read parameters")

# Parameters
# get path to parameters and path json
json_file = open(sys.argv[1])
# read json into dictionary
json_str = json_file.read()
parameter_dict = json.loads(json_str)

# general params from dict (could also save these lines, but I think this way it is easier to digest)
batch=parameter_dict["batch_var"]
data_filepath_full = parameter_dict["harmonization_folder_path"]+parameter_dict["new_name_suffix"]+".h5ad"
results_path = parameter_dict["harmonization_folder_path"]
new_name_suffix = parameter_dict["new_name_suffix"]
global_seed = parameter_dict["global_seed"]
target_clusterN = parameter_dict["target_clusterN"]
job_id = parameter_dict["job_id"]
start_res = parameter_dict["start_res"]
end_res = parameter_dict["end_res"]
step_size = parameter_dict["step_size"]
include_low_res = parameter_dict["include_low_res"]
n_neighbors = parameter_dict["k_param"]
snn_name = "SNN_"+parameter_dict["integration_name"]
min_cells_valid = parameter_dict["min_cells_valid"]
additional_clustering_suffix = parameter_dict["additional_clustering_suffix"]

#define resolution range, hardcoded atm
resolutions = [round(x*step_size,3) for x in range(int(1/step_size)*start_res,int(1/step_size)*end_res+1)]
if(include_low_res):
  low_res_list = [0.001,0.005,0.01,0.05,0.1,0.175,0.25,0.5,0.75]
  resolutions = low_res_list + resolutions
  
  
print(" Read anndata")

# read adata
adata = sc.read_h5ad(data_filepath_full)

#ensure there are no bytestrings 
str_df = adata.obs
str_df = str_df.applymap(lambda x: x.decode() if isinstance(x, bytes) else x)
str_df = str_df.set_index('Cell_ID',drop=False)
adata.obs = str_df
# for features:
str_df = adata.var
str_df = str_df.applymap(lambda x: x.decode() if isinstance(x, bytes) else x)
str_df = str_df.set_index('features',drop=False)
adata.var = str_df

# when exporting from h5seurat: 'Adding SNN_scvi as neighbors'
# scanpy then throws this warbing: 'Moving element from .uns['neighbors']['distances'] to .obsp['distances']'
# --> so as long as we don't overwrite it neighbors contains the seurat graph .obsp['distances'] !
# however sc.tl.leiden does not like this so I also copy the same adjacency matrix into the connectivities slot:
adata.obsp['connectivities'] = adata.obsp['distances']

print(" Run clustering")

# clustering
for res in resolutions:
        key_name = "leiden_clusters_"+str(res)
        sc.tl.leiden(adata,resolution=res,key_added=key_name,random_state=global_seed,neighbors_key='neighbors') # could use neighbors_key
        value_counts = adata.obs[key_name].value_counts()
        print(" Ran leiden with resolution "+str(res)+" and found "+str(len(set(adata.obs[key_name])))+" total clusters with "+str(len(value_counts[value_counts > min_cells_valid]))+" valid clusters")
        if(len(value_counts[value_counts > min_cells_valid]) >= target_clusterN):
            print(" Reached "+str(len(value_counts[value_counts > min_cells_valid]))+" valid clusters")
            break    

#save
meta_subset = adata.obs.filter(regex=("leiden_clusters_"))
meta_subset.to_csv(results_path+new_name_suffix+additional_clustering_suffix+"_leiden_clustering.txt", sep='\t',index=True)

print(" Finalized clustering")





