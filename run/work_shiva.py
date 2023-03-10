import os, sys, random, datetime, shutil

from GPCRmd_dictionaries_shiva import d, abrvs
from rsync_shiva import rsync

print('\n##########################')
print('Running work.py')
print('##########################\n')

running_dir = os.path.dirname(os.path.realpath(__file__))

_, dynid, cores, taskcpus = sys.argv
dynd = d[dynid]


# The generic numbering saved in the GPCRmd dictionary is a dictionary itself that looks like: { 2x47:"77-P-A" }
# The mapper translates the residue names from AlloViz 
# (for the GPCRmd implementation, they look like: "A:PRO:77"; the normal AlloViz implementation doesn't include the chain)
from Bio.SeqUtils import seq3
print('### Residue name mapping ###\n')
gen_num = dynd["gpcr_pdb"]
mapper = {f"{chain}:{seq3(code).upper()}:{num}": key for key, val in gen_num.items()
                                             for num, chain, code in [val.split("-")]}



# rsync the files from the .files file from ori to the folder where the job is being run
print('### Calling rsync function ###\n')
print('Calling rsync function...')
print(f'rsync -z --no-R --inplace --files-from={running_dir}/slurm_files/{dynid}/{dynid}.files ori.prib.upf.edu:/protwis/sites/files .')
rsync(f'rsync -z --no-R --inplace --files-from={running_dir}/slurm_files/{dynid}/{dynid}.files ori.prib.upf.edu:/protwis/sites/files .')
if not any(["_dynamic.tsv" in f for f in os.listdir()]):
    raise Exception(f"ERROR: {dynid} doesn't have contacts data available yet")

# Change the names and locations (symlink) of the GetContacts' .tsvs so that AlloViz can recognize them when GetContacts calculations are launched
print('### Changing names and locations of GetContacts ###\n')
tsvs = sorted([f for f in os.listdir() if "_dynamic.tsv" in f], key=lambda trajid: trajid.split("_")[0].split("-")[-1]) # filenames look like: {dynid}-{trajid}_dynamic.tsv
for num, f in enumerate(tsvs, 1):
    if os.path.isfile(f"{num}.tsv"):
        os.remove(f"{num}.tsv")
    os.symlink(f, f"{num}.tsv")


# Extract and save the information about the ligand resname/chain from the GPCRmd dictionary so that only the protein is 
# processed later by AlloViz (and not e.g. a peptidic ligand)
print('### Extracting ligand information ###\n')
extra_sel = ""
lig = dynd["lig_sname"]
if lig and len(lig) > 0:
    extra_sel = f" and (not chainid {lig.replace(':', '')})" if ":" in lig else f" and (not resname {lig})"
protein_sel = f"(same segid as protein) and (not segid LIG) and (not chainid L){extra_sel}"
print("protein_sel", protein_sel)
    
    

import AlloViz

# pkgs = ['MDTASK', 'pytraj_CA', 'pytraj_CB', 'dynetan', #'dynetan_COM',
# 'correlationplus_CA_Pear', 'correlationplus_COM_Pear', 'correlationplus_CA_LMI', 'correlationplus_COM_LMI',
# 'correlationplus_Phi', 'correlationplus_Psi', 'correlationplus_Omega', 
# #'correlationplus_Dihs',
# 'AlloViz_Phi', 'AlloViz_Psi', 'AlloViz_Omega', 
# #'AlloViz_Dihs',
# 'MDEntropy_Dihs', 'MDEntropy_AlphaAngle', #'MDEntropy_Contacts',
# 'GetContacts', 'PyInteraph2_Contacts', 'PyInteraph2_Energy',
# 'g_correlation_CA_MI', 'g_correlation_COM_MI', 'g_correlation_CA_LMI', 'g_correlation_COM_LMI']

# csvize = lambda pkgs, end="\n": ",".join(pkgs) + end
# with open(f"{running_dir}/combinations.csv", "r+") as tested:
#     random.shuffle(pkgs)
#     while csvize(pkgs) in tested.readlines():
#         random.shuffle(pkgs)
#     tested.write(csvize(pkgs))
    

t = datetime.datetime.now()
print('### Running Alloviz ###')
print('# Preparing the system...')
dyn = AlloViz.Protein(pdb = dynd['struc_fname'],
                      trajs = sorted(dynd['traj_fnames'], key=lambda trajid: trajid.split("_")[0]), # filenames look like: {trajid}_trj_{dynid}.xtc
                      name = dynid,
                      #path = out + dynid,
                      protein_sel = protein_sel,
                      special_res = {"HSP": "H"}, # for dyn10
                      GPCR = mapper, # In the GPCRmd AlloViz implementation, the GPCR parameter is the generic numbering mapper dictionary
                     )
print('# Done.\n')

print('# Computing networks...\n')
# dyn.calculate(pkgs=["all"], cores=int(cores), taskcpus=int(taskcpus)) #"all"
# mypkgs=["MDTASK","pytraj_CA","pytraj_CB","dynetan","correlationplus_CA_Pear","correlationplus_COM_Pear","correlationplus_CA_LMI","correlationplus_COM_LMI","correlationplus_Phi","correlationplus_Psi","correlationplus_Backbone_Dihs_Avg","correlationplus_Backbone_Dihs_Max","AlloViz_Phi","AlloViz_Psi","AlloViz_Backbone_Dihs_Avg","AlloViz_Backbone_Dihs_Max","AlloViz_Chi1","AlloViz_Chi2","AlloViz_Chi3","AlloViz_Chi4","AlloViz_Sidechain_Dihs_Avg","AlloViz_Sidechain_Dihs_Max","AlloViz_Dihs_Avg","MDEntropy_Phi","MDEntropy_Psi","MDEntropy_Dihs","MDEntropy_AlphaAngle","PyInteraph2_Atomic_Contacts_Occurrence","PyInteraph2_Atomic_Contacts_Strength","PyInteraph2_COM_Contacts","PyInteraph2_COM_Contacts_Corrected","PyInteraph2_Energy", "GetContacts"]
mypkgs = [pkg for pkg in AlloViz.AlloViz.info.wrappers if \
 "g_correlation" not in pkg and \
 "CARDS" not in pkg and \
 "_Max" not in pkg and \
 "MDEntropy_Contacts" not in pkg]

dyn.calculate(pkgs=mypkgs, cores=int(cores), taskcpus=int(taskcpus)) 
print('# Done.\n')


print('# Computing filters...\n')
dyn.filter(pkgs="all", filterings=["All", "GetContacts_edges", "GPCR_Interhelix", "Spatially_distant"], GetContacts_threshold=0.5) #
# dyn.filter(pkgs="all", filterings=["All", "GPCR_Interhelix", "Spatially_distant"]) #
print('# Done.\n')

print('# Analysing networks...\n')
dyn.analyze(pkgs="all", filterings="all", elements="edges", metrics="all", normalize=True, cores=int(cores)) 
# dyn.analyze(pkgs="pytraj_CA", filterings="Spatially_distant", elements="edges", metrics="all", normalize=True, cores=int(cores))
print('# Done.\n\n')


print('### Preparations to change AlloViz standard 3-letter residue codes to the original PDB ones ###\n')
# Added to be able to map AlloViz-standardized 3-letter residue names back to the original GPCRmd PDB resnames
from MDAnalysis import Universe
from Bio.SeqUtils import seq1, seq3
process = lambda name: seq3(seq1(name, custom_map=dyn._standard_resdict)).upper()
# df.index = df.index.map(lambda idx: tuple(":".join([chain, process(name), num]) for res in idx for chain, name, num in [res.split(":")]))
# Mapper to transform AlloViz-standardized 3-letter residue names back to the original GPCRmd PDB resnames for saving
original_resmap = {f"{aa.atoms[0].chainID}:{process(aa.resname)}:{aa.resid}": f"{aa.atoms[0].chainID}:{aa.resname}:{aa.resid}" \
                         for aa in Universe(dynd['struc_fname']).select_atoms(dyn._protein_sel).residues}



print('### Reorganizing outfiles ###\n')
shutil.move(f"{dynid}.times", f"{running_dir}/slurm_files/{dynid}/{dynid}.times")
for file in [f for f in os.listdir() if os.path.isfile(f)]:
    os.remove(file)
    

    
print('### Building csv files ###\n') 
# abrvs['GetContacts'] = abrvs['getcontacts']
abrvs['MI'] = abrvs['Mutual Information (MI)']
abrvs['Linear MI (LMI)'] = abrvs['Linear Mutual Information (LMI)']
abrvs['LMI'] = abrvs['Linear Mutual Information (LMI)']
abrvs["Pearson's"] = abrvs['Pearson']
abrvs["Carbon \u03B1"] = abrvs['Alpha-carbon']
abrvs["Carbon \u03B2"] = abrvs['Beta-carbon']
# abrvs["Residue COM contacts"] = abrvs['Residue COM']
# abrvs["Backbone dihedrals"] = abrvs['All dihedrals']
# abrvs["Alpha angle"] = abrvs['Alpha']
abrvs["Atoms' displacements"] = abrvs["Atoms' movement correlation"]
abrvs["Dihedral angles"] = abrvs["Dihedrals' movement correlation"]
abrvs["All backbone dihedrals (Phi and psi)"] = abrvs["All backbone dihedrals"]
# abrvs["All backbone dihedrals (Phi and psi) (max. value)"] = abrvs["All backbone dihedrals (max. value)"]

w = AlloViz.AlloViz.info.wrappers


def get_csv_name(pkg, filterby):
    # Order of strings in w is: Info, Pkg, Metric, Atom/angle
    value = list(w[pkg])
    
    # Order for csv name is: (info)_[(subinfo for dihs)]_(atom/angle)_(metric)_(pkg)_(filterby).csv
    # Except for int. ene. and cont. freq: (info)_(pkg)_(filterby).csv
    if value[0] == "Interaction energy":
        order = [0,1]
            
    elif value[0] == "Contacts":
        order = [0,3,1]
            
    elif value[0] == "Atoms' displacements":
        order = [0,3,2,1]
            
    elif value[0] == "Dihedral angles":
        if any([i in pkg for i in ["Phi", "Psi", "AlphaAngle", "Backbone_Dihs"]]): #"Omega",
            value.insert(4, "Backbone dihedrals")
            
        elif "Chi" in pkg or "Sidechain" in pkg:
            value.insert(4, "Side-chain dihedrals")
            
        elif "AlloViz_Dihs" in pkg or ("CARDS" in pkg and "_MI_Dihs" in pkg):
            value.insert(4, "All dihedrals")
            
        order = [0,4,3,2,1]
        
    strings = [value[ix] for ix in order]
    
    name = abrvs[strings.pop(0)]
    for string in strings:
        name = name + "_" + abrvs[string]
        
    return name + "_" + abrvs[filterby]




get_gen_num = lambda res: mapper[res] if res in mapper else "-"

for pkg in [pkg for pkg in dyn.__dict__ if pkg in w]:
    for filterby in [fb for fb in getattr(dyn, pkg).__dict__ if fb in AlloViz.AlloViz.utils.filteringsl]:
        df = AlloViz.AlloViz.utils.rgetattr(dyn, pkg, filterby, "edges")
        if not isinstance(df, bool):
            df.rename(columns={"weight": "rwd_weight", "weight_avg": "rwd_weight",
                               "cfb": "cfb_weight", "cfb_avg": "cfb_weight",
                               "btw": "btw_weight", "btw_avg": "btw_weight"},
                     inplace=True)
            cols = [col for col in df.columns if ("weight" in col or "std" in col)]
            df = df[cols]
            
            # Get the generic numbers for the residues and save them; the gen num mapper was built with standard residue names
            gen_num_ix = df.index.map( lambda ix: (get_gen_num(ix[0]), get_gen_num(ix[1])) )
            
            # Change the standardized 3-letter residue names back to the original GPCRmd PDB resnames
            df.index = df.index.map(lambda idx: tuple(original_resmap[ix] for ix in idx))
            
            # Append the saved generic numbers to the de-standardized index
            df.set_index(gen_num_ix, append=True, inplace=True)

            # Save
            df.to_csv(f"{get_csv_name(pkg, filterby)}.csv", index_label=["resid1", "resid2", "resid1_gennum", "resid2_gennum"])


# with open(f"{running_dir}/slurm_files/{dynid}/save_{dynid}.log", "w") as f:

    
print('### Removing files from cluster ###\n')   
rsync(f'rsync -z ./* ori.prib.upf.edu:/protwis/sites/files/Precomputed/allosteric_com/{dynid}/') # --remove-source-files
rsync(f'rsync -zr data/ ori.prib.upf.edu:/protwis/sites/files/Precomputed/allosteric_com/data/{dynid}/') # --remove-source-files

print('### Writting completion.times ###\n')
with open(f"{running_dir}/completion.times", "a") as f:
    tf = datetime.datetime.now()
    trajs = len(dyn.u.trajectory.readers) if hasattr(dyn.u.trajectory, "readers") else 1
    frames = dyn.u.trajectory.n_frames
    residues = dyn.u.atoms.residues.n_residues
    f.write(f"{dynid},{cores},{taskcpus},{trajs},{frames},{residues},{tf - t}\n") #{os.environ['HOSTNAME']} removed, not working in shiva

print('### WORK FINISHED ###\n')
