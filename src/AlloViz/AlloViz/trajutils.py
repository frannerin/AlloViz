import os

import requests, parmed
import MDAnalysis as mda
import numpy as np
from multiprocess import Pool

from . import utils




def download_GPCRmd_files(Protein):
    import sys, re
    import tarfile, fileinput
    from bs4 import BeautifulSoup
    import urllib.request as pwget

    web = "https://submission.gpcrmd.org"

    html = requests.get(f"{web}/dynadb/dynamics/id/{Protein.GPCR}/")
    soup = BeautifulSoup(html.text, features="html.parser").find_all('a')
    links = [link.get('href') for link in soup if re.search("(xtc|pdb|psf|prm)", link.get('href'))]
    get_name = lambda link: f"{Protein._path}/{link.rsplit('/')[-1]}"


    mypool = Pool()
    utils.pool = mypool        

    for link in links:
        fname = get_name(link)
        print(f"downloading {fname}")
        utils.get_pool().apply_async(pwget.urlretrieve, args=(f"{web}{link}", fname))

    mypool.close()
    mypool.join()
    utils.pool = utils.dummypool()


    fname = next(get_name(link) for link in links if "prm" in get_name(link))
    with tarfile.open(fname) as tar:
        tar.extractall(Protein._path)
    os.remove(fname)

    for line in fileinput.input(f"{Protein._path}/parameters", inplace=True):
        if line.strip().startswith('HBOND'):
            line = 'HBOND CUTHB 0.5\n'
        elif line.strip().startswith('END'):
            line = 'END\n'
        sys.stdout.write(line)

    return







def standardize_resnames(protein, **kwargs):
    from Bio.SeqUtils import seq1, seq3
    from Bio.SeqUtils import IUPACData
    from MDAnalysis.lib.util import inverse_aa_codes

    res_d = {}
    res_d.update(IUPACData.protein_letters_3to1_extended)
    res_d.update(inverse_aa_codes)
    if "special_res" in kwargs and isinstance(kwargs["special_res"], dict):
        res_d.update(kwargs["special_res"])

    try:
        for res in protein.residues:
            res.resname = seq3(seq1(res.resname, custom_map = res_d)).upper()
    except:
        print(f"""
resname(s) {[res for res in set(protein.residues.resnames) if res not in res_d]} could not be translated to a standard name.
Please provide a mapping from 3/4-letter code to 1-letter code as a dictionary with keyword argument 'special_res' upon initialization.
              """)
        raise
        
        
def get_GPCRdb_numbering(protein):
    # print "GPCRdb script labeling generic numbers on the annotated pdb structure\nKey bindings for labels:\nF1 - show generic numbers\nF2 - show Ballesteros-Weinstein numbers\nF3 - clear labels"
    # label n. CA or n. N
    # cmd.set_key('F1', 'label n. CA & (b >-8.1 and  b <8.1), str("%1.2f" %b).replace(".","x") if b > 0 else str("%1.3f" %(-b + 0.001)).replace(".", "x")')
    # cmd.set_key('F2', 'label n. N & (b > 0 and  b <8.1), "%1.2f" %b')
    # cmd.set_key('F3', 'label n. CA or n. N')
    from io import StringIO as IO
    
    ### USING CONTEXT MANAGERS
    def get_GPCRdb_result(protein):
        with IO("") as protf, mda.lib.util.NamedStream(protf, "prot.pdb") as f:
            protein.write(f, file_format="pdb")
            GPCRdb_output = requests.post('https://gpcrdb.org/services/structure/assign_generic_numbers',
                                          files = {'pdb_file': protf.getvalue()})
            return GPCRdb_output.text    
    
    out = get_GPCRdb_result(protein)
    with mda.lib.util.NamedStream(IO(out), "gen_numb.pdb") as gen_numb:
        nums = mda.Universe(gen_numb).select_atoms("name CA").tempfactors
        GPCRdb_nums = np.array([round(bfactor, 2) if bfactor>1 else round(-bfactor+0.001, 3) if bfactor<0 else 0 for bfactor in nums])
        protein.select_atoms("name CA").tempfactors = GPCRdb_nums
        
    return

    ### WITHOUT CONTEXT MANAGERS
#     protf = IO()
#     f = mda.lib.util.NamedStream(protf, "prot.pdb")
    
#     protein.write(f, file_format="pdb")
#     GPCRdb_output = requests.post('https://gpcrdb.org/services/structure/assign_generic_numbers',
#                                   files = {'pdb_file': protf.getvalue()})
#     out = GPCRdb_output.text
    
#     gen_numb = mda.lib.util.NamedStream(IO(out), "gen_numb.pdb")
#     nums = mda.Universe(gen_numb).select_atoms("name CA").tempfactors
#     GPCRdb_nums = np.array([round(bfactor, 2) if bfactor>1 else round(-bfactor+0.001, 3) if bfactor<0 else 0 for bfactor in nums])
#     protein.select_atoms("name CA").tempfactors = GPCRdb_nums
        
        
        
def write_protein_trajs(whole, protein, trajs):
    def write_protein_traj(protein, reader, trajf):
        with mda.Writer(trajf, protein.n_atoms) as W:
                for ts in reader:
                    W.write(protein.atoms)
    
    
    mypool = Pool()
    utils.pool = mypool

    for num, trajf in trajs.items():
        if not os.path.isfile(trajf):
            reader = whole.trajectory.readers[num-1] if hasattr(whole.trajectory, "readers") else whole.trajectory
            utils.get_pool().apply_async(write_protein_traj, args=(protein, reader, trajf))

    mypool.close()
    mypool.join()
    utils.pool = utils.dummypool()
    
        


def process_input(Protein, **kwargs):
    whole = mda.Universe(Protein.pdb)
    protein = whole.select_atoms(Protein._protein_sel)

    # Rename all residues in protein_sel to standard names
    standardize_resnames(protein, **kwargs)

    # Retrieve GPCRdb residue generic numbering if it's a GPCR
    if Protein.GPCR:
        get_GPCRdb_numbering(protein)

    # Write protein pdb file
    if not os.path.isfile(Protein._pdbf):
        protein.write(Protein._pdbf)

    # Write protein psf file if it exists
    psff = Protein._pdbf.replace("pdb", "psf")
    if "psf" in kwargs and not os.path.isfile(psff):
        psf = parmed.load_file(kwargs["psf"])[protein.indices]
        psf.title = psff
        psf.write_psf(psff)
        setattr(Protein, "_psff", psff)

    # Write protein trajectory(ies) file(s)
    if any([not os.path.isfile(f) for f in Protein._trajs.values()]):
        whole.load_new(Protein.trajs)
        write_protein_trajs(whole, protein, Protein._trajs)
    
    
    
    
    
    
def get_bonded_cys(Protein):
    "Identify disulfide bond-forming cysteines' sulphur atoms"
    bonded_cys = []
    pdb = parmed.read_pdb(Protein._pdbf)
    mask = parmed.amber.mask.AmberMask(pdb, ':CY?@SG')
    for sel in mask.Selected():
        atom = pdb.atoms[sel]
        bonded = [a.name for a in atom.bond_partners]
        if "SG" in bonded:
            bonded_cys.append(atom.residue.idx)

    return bonded_cys
                    
                    
                    
                    
           
        
def get_COM_trajs(Protein):
    compath = f"{Protein._datadir}/COM_trajs"
    if not os.path.isdir(compath): os.makedirs(compath, exist_ok=True)

    compdb = f"{compath}/ca.pdb"
    if not os.path.isfile(compdb):
        print(f"Making trajectories of residue COM for {Protein.name}")
        prot = Protein.protein.select_atoms("name CA")
        prot.write(compdb)
    setattr(Protein, "_compdbf", compdb)

    comtrajs = [f"{compath}/{xtc}.xtc" for xtc in Protein._trajs]
    setattr(Protein, "_comtrajs", {num: traj for num, traj in enumerate(comtrajs, 1)})

    for xtc, comtraj in Protein._comtrajs.items():
        if not os.path.isfile(comtraj):
            prot = Protein.protein
            traj = Protein.U.trajectory.readers[xtc-1] if hasattr(Protein.U.trajectory, "readers") else Protein.U.trajectory
            # traj = next(traj for traj in Protein.mdau.trajectory.readers if traj.filename == Protein._trajs[xtc])
            arr = np.empty((prot.n_residues, traj.n_frames, 3))
            for ts in traj:
                arr[:, ts.frame] = prot.center_of_mass(compound='residues')

            cau = mda.Universe(compdb, arr, format=mda.coordinates.memory.MemoryReader, order='afc')
            with mda.Writer(comtraj, cau.atoms.n_atoms) as W:
                for ts in cau.trajectory:
                    W.write(cau.atoms)



    




# def get_mdau(state, protein_sel):
#     pdb = mda.Universe(state._pdbf)#, *state._trajs.values())  
#     prot = pdb.select_atoms(protein_sel)
#     mdau = mda.Universe(state._pdbf, *state._trajs.values())
    
#     protcopy = mda.core.universe.Merge(prot.copy()).atoms
#     # protfile = f"{state._datadir}/prot.pdb"
#     protfile = f"{state._datadir}/prot_{state._pdbf.rsplit('/', 1)[-1]}"
#     setattr(state, "_protpdb", protfile)
    
#     try:
#         from Bio.SeqUtils import seq1, seq3
#         for res in protcopy.residues:
#             res.resname = seq3(seq1(res.resname, custom_map = state._res_dict)).upper()
#     except:
#         print(f"""
# resname(s) {[res for res in set(protcopy.residues.resnames) if res not in state._res_dict]} could not be translated to a standard name.
# Please provide a mapping from 3/4-letter code to 1-letter code as a dictionary with keyword argument 'special_res' upon initialization.
#               """)
#         raise 

    
#     if state.GPCR:#hasattr(state, "_gpcrmdid"):
#         # prot_numsf = f"{state._datadir}/{state._pdbf.replace(f'{state._path}/', 'gpcrdb_')}"
#         prot_numsf = state._protpdb

#         if not os.path.isfile(prot_numsf):
#             print(f"retrieving {prot_numsf}")
#             with io.StringIO() as protf:
#                 with mda.lib.util.NamedStream(protf, "prot.pdb") as f:
#                     protcopy.write(f, file_format="pdb")
#                 response = requests.post('https://gpcrdb.org/services/structure/assign_generic_numbers', 
#                                          files = {'pdb_file': protf.getvalue()})
#                 with open(prot_numsf, "w") as prot_nums:
#                     prot_nums.write(response.text)
        
#         # print "GPCRdb script labeling generic numbers on the annotated pdb structure\nKey bindings for labels:\nF1 - show generic numbers\nF2 - show Ballesteros-Weinstein numbers\nF3 - clear labels"
#         # label n. CA or n. N
#         # cmd.set_key('F1', 'label n. CA & (b >-8.1 and  b <8.1), str("%1.2f" %b).replace(".","x") if b > 0 else str("%1.3f" %(-b + 0.001)).replace(".", "x")')
#         # cmd.set_key('F2', 'label n. N & (b > 0 and  b <8.1), "%1.2f" %b')
#         # cmd.set_key('F3', 'label n. CA or n. N')
#         nums = mda.Universe(prot_numsf).select_atoms(f"{protein_sel} and name CA").tempfactors
#         GPCRdb_nums = np.array([round(bfactor, 2) if bfactor>1 else round(-bfactor+0.001, 3) if bfactor<0 else 0 for bfactor in nums])
        
#         mdau.select_atoms(f"{protein_sel} and name CA").tempfactors = GPCRdb_nums#.round(2)
#         pdb.select_atoms(f"{protein_sel} and name CA").tempfactors = GPCRdb_nums
#         prot.select_atoms(f"{protein_sel} and name CA").tempfactors = GPCRdb_nums
    
#     else:
#         protcopy.write(state._protpdb)
    
    
#     return pdb, prot, mdau




# def get_dihedral_residx(prot):
#     res_arrays = np.split(prot.residues.resindices, np.where(np.diff(prot.residues.resnums) != 1)[0]+1)
#     dihedral_residx = lambda end=-1: [elem for arr in res_arrays for elem in arr[1:end]]
    
#     return dihedral_residx



# def get_res_dict(state, **kwargs):
#     from Bio.SeqUtils import IUPACData
#     from MDAnalysis.lib.util import inverse_aa_codes
    
#     res_d = {}
#     res_d.update(IUPACData.protein_letters_3to1_extended)
#     res_d.update(inverse_aa_codes)
#     if "special_res" in kwargs and isinstance(kwargs["special_res"], dict):
#         res_d.update(kwargs["special_res"])
        
#     return res_d



#######3 change state to Protein and use Protein's "u"; don't use protein sel







# def get_bonded_cys(state, pdb):
#     import parmed
    
#     indices = []

#     pdb = parmed.read_pdb(pdb)
#     mask = parmed.amber.mask.AmberMask(pdb, ':CY?@SG')
#     for sel in mask.Selected():
#         atom = pdb.atoms[sel]
#         bonded = [a.name for a in atom.bond_partners]
#         if "SG" in bonded:
#             indices.append(atom.residue.idx)

#     return indices