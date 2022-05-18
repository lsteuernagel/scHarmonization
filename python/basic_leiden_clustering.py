
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

# find file helper
def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if pattern in name:
                result.append(os.path.join(root, name))
    return result

# Parameters
# get path to parameters and path json
json_file = open(sys.argv[1])
# read json into dictionary
json_str = json_file.read()
parameter_dict = json.loads(json_str)

# general params from dict (could also save these lines, but I think this way it is easier to digest)
batch=parameter_dict["batch_var"]
data_filepath_full = parameter_dict["data_filepath_full"]
results_path = parameter_dict["harmonization_folder_path"]
new_name_suffix = parameter_dict["new_name_suffix"]
global_seed = parameter_dict["global_seed"]
target_clusterN = parameter_dict["target_clusterN"]
job_id = parameter_dict["job_id"]
start_res = parameter_dict["start_res_asw"]
end_res = parameter_dict["end_res_asw"]
n_neighbors = parameter_dict["k_param"]
snn_name = "SNN_"+parameter_dict["integration_name"]

#define resolution range, hardcoded atm
resolutions = [x*0.5 for x in range(2*start_res, 2*end_res+1)]

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

# clustering
for res in resolutions:
        key_name = "leiden_clusters_"+str(res)
        sc.tl.leiden(adata,resolution=res,key_added=key_name,random_state=global_seed,neighbors_key=snn_name) # could use neighbors_key
        print(" Ran leiden with resolution "+str(res)+" and found "+str(len(set(adata.obs[key_name])))+" clusters")
        if(len(set(adata.obs[key_name])) >= target_clusterN):
            print(" Reached "+str(len(set(adata.obs[key_name])))+" clusters")
            break    

#save
meta_subset = adata.obs.filter(regex=("Cell_ID|leiden_clusters_"))
meta_subset.to_csv(results_path+new_name_suffix+"_initial_leiden_clustering.txt", sep='\t',index=True)

print(" Finalized clustering")





