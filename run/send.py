# RUN IN HYDRA, saved in an accessible, shared folder of /gpcr/users

import os
from contextlib import redirect_stdout, redirect_stderr
from GPCRmd_dictionaries import d

running_dir = os.path.dirname(os.path.realpath(__file__)) # /gpcr/users/frann/AlloViz_GPCRmd
os.makedirs(f"{running_dir}/slurm_files", exist_ok=True)

with open(f"{running_dir}/combinations.csv", "a") as f:
    pass

completion_times = f"{running_dir}/completion.times"
if not os.path.isfile(completion_times):
    with open(completion_times, "a") as f:
        f.write("dynid,node,cores,taskcpus,num_trajs,num_frames,num_residues,time\n") #,pkgs_list



dyns = "9".split() # add dyn9
dyns = " ".join([f"dyn{dynid}" for dynid in dyns]).split()

nodes = set("gimli kili legolas thorin fili dwalin arwen balin bifur bombur aragorn".split())
use = set(["legolas"]) #set("kili legolas fili aragorn".split())

cores = 1
taskcpus = 12



send = lambda dyn: f"""#!/bin/bash
#SBATCH --job-name={dyn}
#SBATCH --partition=gpcr_gpu
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

python3 {running_dir}/work.py {dyn} {cores} {taskcpus}

exit"""





for dyn, dynd in d.items():
    if dyn in dyns:
        os.makedirs(f"{running_dir}/slurm_files/{dyn}", exist_ok=True)
        exe = lambda ext: f"{running_dir}/slurm_files/{dyn}/{dyn}.{ext}"
        with open(exe("sh"), "w") as f:
            f.write(send(dyn))

        trajsl = dynd["traj_f"]
        trajids = [traj.split("/")[-1].split("_")[0] for traj in trajsl]

        with open(exe("files"), "w") as f:
            for file in trajsl + [dynd["struc_f"]]:
                f.write(file + "\n")
            for ix in trajids:
                file = f"Precomputed/get_contacts_files/dynamic_symlinks/{dyn}/{dyn}-{ix}_dynamic.tsv"
                f.write(file + "\n")

        os.system(f"sbatch -D /home/{os.environ.get('USER')} {exe('sh')}")
    
    
    #copy results to shared folder and maybe also delete the originals
    #future: copy results to ori and delete the originals