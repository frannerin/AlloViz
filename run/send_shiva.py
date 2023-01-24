# RUN IN SHIVA, saved in an accessible, shared folder of /gpcr/users
# Previously, create an environment in the same folder called "env" (/gpcr/users/whatever/env), e.g. with dev/dev_env_ubuntu.yml

import os
from contextlib import redirect_stdout, redirect_stderr
from GPCRmd_dictionaries_shiva import d



# Save the path from where this script is stored and launched to a variable
running_dir = os.path.dirname(os.path.realpath(__file__)) # /gpcr/users//mlopez/allosteric/AlloViz_to_GPCRmd/AlloViz/run
# Create a folder into the running_dir to store the slurm logs an other files
os.makedirs(f"{running_dir}/slurm_files", exist_ok=True)


# Create a file to store the calculation times for each system if it doesn't exist yet
completion_times = f"{running_dir}/completion.times"
if not os.path.isfile(completion_times):
    with open(completion_times, "a") as f:
        f.write("dynid,node,cores,taskcpus,num_trajs,num_frames,num_residues,time\n") #,pkgs_list


# Space-separated list, as a single string, of the dynids to run
dyns = "9".split()
dyns = " ".join([f"dyn{dynid}" for dynid in dyns]).split()
# print('Computing dynID:', dyns)

### HYDRA CLUSTER
# # Put in "use" the nodes of gpcr_gpu queue (unique queue) in which calculations can be sent to
# nodes = set("gimli kili legolas thorin fili dwalin arwen balin bifur bombur aragorn".split())
# use = set(["aragorn"]) #set("kili legolas fili aragorn".split())

# # Cores and taskcpus for the jobs
# cores = 12 # Aragorn number of CPU cores
# taskcpus = 4 # cores per replica
# # 4taskCPUS * 3replicates = 12cores


### SHIVA CLUSTER
# Put in "use" the nodes of wanted queue ('short', 'normal',...) in which calculations can be sent to:

## 'short' queue
partition = "short"
nodes = set("node01 node02 node08 node09 node16".split())
use = set(["node09"]) 

## 'normal' queue
# partition = "normal"
# nodes = set("node10 node11 node12 node18 node19 node20 node25 node26 node27 node28 node29".split())
# use = set(["node11"])

## 'long' queue
# partition = "long"
# nodes = set("node13 node14 node15 node22 node23".split())
# use = set(["node15"])


# Cores and taskcpus for the jobs
cores = 12
taskcpus = 4 # cores per replica
# 8taskCPUS * 3replicates = 24cores



# String with the .sh file that will be sent to the slurm queue (one string per dyn)
# CHANGE THE --partition (name of the Shiva queue selected)
send = lambda dyn: f"""#!/bin/bash
#SBATCH --job-name={dyn}
#SBATCH --partition={partition}
#SBATCH --exclude={','.join(nodes - use)}
#SBATCH --nodes=1-1
#SBATCH --ntasks={cores}
#SBATCH --ntasks-per-core=1
#SBATCH --mem-per-cpu=2500
#SBATCH --output={running_dir}/slurm_files/{dyn}/%x_%j.%N.out # output to shared gpcr folder
#SBATCH --error={running_dir}/slurm_files/{dyn}/%x_%j.%N.err

module load Miniconda3
eval "$(conda shell.bash hook)"
conda activate {running_dir}/env

mkdir -p {dyn}
cd {dyn}

python3 {running_dir}/work_shiva.py {dyn} {cores} {taskcpus}

exit"""



for dyn in dyns:
    # Retrieve the information for this dynid from the GPCRmd compl_dict
    dynd = d[dyn]
    # print('Build sh for dynID:', dyn)
    
    # Create a dynid folder to store slurm files in the shared folder
    os.makedirs(f"{running_dir}/slurm_files/{dyn}", exist_ok=True)
    # Lambda function to create running-related files for this dynid
    exe = lambda ext: f"{running_dir}/slurm_files/{dyn}/{dyn}.{ext}"
    
    # Create the sending .sh file
    with open(exe("sh"), "w") as f:
        f.write(send(dyn))
    
    # Retrieve the names of IDs and files and create a .files file
    # it includes the path of the files (in ori) that need to be rsync'ed to the computer in which the 
    # calculations are running to use them: 
    # pdb (struc_f), trajectories (traj_f) and GetContacts files (one per each traj_f id)
    trajsl = dynd["traj_f"] # names are like 10180_trj_10.xtc
    trajids = [traj.split("/")[-1].split("_")[0] for traj in trajsl]
    with open(exe("files"), "w") as f:
        for file in trajsl + [dynd["struc_f"]]:
            f.write(file + "\n")
        for ix in trajids:
            file = f"Precomputed/get_contacts_files/dynamic_symlinks/{dyn}/{dyn}-{ix}_dynamic.tsv"
            f.write(file + "\n")
    
    # Send the job to the slurm queue
    os.system(f"sbatch -D /users/gpcr/{os.environ.get('USER')} {exe('sh')}")