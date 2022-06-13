import sys, os, re, io, requests

import MDAnalysis as mda
import numpy as np

from . import utils



def download_files(state):
    from multiprocess import Pool
    from bs4 import BeautifulSoup
    import urllib.request as pwget
    import tarfile, fileinput

    web = "https://submission.gpcrmd.org"

    html = requests.get(f"{web}/dynadb/dynamics/id/{state._gpcrmdid}/")
    soup = BeautifulSoup(html.text, features="html.parser").find_all('a')
    links = [link.get('href') for link in soup if re.search("(xtc|pdb|psf|prm)", link.get('href'))]
    get_name = lambda link: f"{state._path}/{link.rsplit('/')[-1]}"


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
        tar.extractall(state._path)
    os.remove(fname)

    for line in fileinput.input(f"{state._path}/parameters", inplace=True):
        if line.strip().startswith('HBOND'):
            line = 'HBOND CUTHB 0.5\n'
        elif line.strip().startswith('END'):
            line = 'END\n'
        sys.stdout.write(line)

    return




def get_mdau(state, protein_sel):
    pdb = mda.Universe(state._pdbf)#, *state._trajs.values())  
    prot = pdb.select_atoms(protein_sel)
    mdau = mda.Universe(state._pdbf, *state._trajs.values())
    
    protcopy = mda.core.universe.Merge(prot.copy()).atoms
    # protfile = f"{state._datadir}/prot.pdb"
    protfile = f"{state._datadir}/prot_{state._pdbf.rsplit('/', 1)[-1]}"
    setattr(state, "_protpdb", protfile)
    
    try:
        from Bio.SeqUtils import seq1, seq3
        for res in protcopy.residues:
            res.resname = seq3(seq1(res.resname, custom_map = state._res_dict)).upper()
    except:
        print(f"""
resname(s) {[res for res in set(protcopy.residues.resnames) if res not in state._res_dict]} could not be translated to a standard name.
Please provide a mapping from 3/4-letter code to 1-letter code as a dictionary with keyword argument 'special_res' upon initialization.
              """)
        raise 

    
    if state.GPCR:#hasattr(state, "_gpcrmdid"):
        # prot_numsf = f"{state._datadir}/{state._pdbf.replace(f'{state._path}/', 'gpcrdb_')}"
        prot_numsf = state._protpdb

        if not os.path.isfile(prot_numsf):
            print(f"retrieving {prot_numsf}")
            with io.StringIO() as protf:
                with mda.lib.util.NamedStream(protf, "prot.pdb") as f:
                    protcopy.write(f, file_format="pdb")
                response = requests.post('https://gpcrdb.org/services/structure/assign_generic_numbers', 
                                         files = {'pdb_file': protf.getvalue()})
                with open(prot_numsf, "w") as prot_nums:
                    prot_nums.write(response.text)

        nums = mda.Universe(prot_numsf).select_atoms(f"{protein_sel} and name N").tempfactors
        mdau.select_atoms(f"{protein_sel} and name N").tempfactors = nums.round(2)
    
    else:
        protcopy.write(state._protpdb)
    
    
    return pdb, prot, mdau




def get_dihedral_residx(prot):
    res_arrays = np.split(prot.residues.resindices, np.where(np.diff(prot.residues.resnums) != 1)[0]+1)
    dihedral_residx = lambda end=-1: [elem for arr in res_arrays for elem in arr[1:end]]
    
    return dihedral_residx



def get_res_dict(state, **kwargs):
    from Bio.SeqUtils import IUPACData
    from MDAnalysis.lib.util import inverse_aa_codes
    
    res_d = {}
    res_d.update(IUPACData.protein_letters_3to1_extended)
    res_d.update(inverse_aa_codes)
    if "special_res" in kwargs and isinstance(kwargs["special_res"], dict):
        res_d.update(kwargs["special_res"])
        
    return res_d




def add_comtrajs(state, protein_sel):
    compath = f"{state._datadir}/COMtrajs"
    if not os.path.isdir(compath): os.makedirs(compath, exist_ok=True)

    compdb = f"{compath}/ca.pdb"
    if not os.path.isfile(compdb):
        print(f"Making trajectories of residue COM for {state._pdbf}")
        prot = state.mdau.select_atoms(f"({protein_sel}) and name CA")
        prot.write(compdb)
    setattr(state, "_compdbf", compdb)

    comtrajs = [f"{compath}/{xtc}.xtc" for xtc in state._trajs]
    setattr(state, "_comtrajs", {num: traj for num, traj in enumerate(comtrajs, 1)})

    for xtc, comtraj in state._comtrajs.items():
        if not os.path.isfile(comtraj):
            prot = state.mdau.select_atoms(protein_sel)
            traj = state.mdau.trajectory.readers[xtc-1] if hasattr(state.mdau.trajectory, "readers") else state.mdau.trajectory
            # traj = next(traj for traj in state.mdau.trajectory.readers if traj.filename == state._trajs[xtc])
            arr = np.empty((prot.n_residues, traj.n_frames, 3))
            for ts in traj:
                arr[:, ts.frame] = prot.center_of_mass(compound='residues')

            cau = mda.Universe(compdb, arr, format=mda.coordinates.memory.MemoryReader, order='afc')
            with mda.Writer(comtraj, cau.atoms.n_atoms) as W:
                for ts in cau.trajectory:
                    W.write(cau.atoms)
    return




def make_dcds(state):        
    dcdpath = f"{state._datadir}/dcds"
    if not os.path.isdir(dcdpath): os.makedirs(dcdpath, exist_ok=True)
    # setattr(state, "_protf", lambda ext: f"{dcdpath}/prot.{ext}")
    #
    # atoms = state.mdau.select_atoms("same segid as protein")
    #
    # dcdpdb = state._protf("pdb")
    # if not os.path.isfile(dcdpdb): atoms.write(dcdpdb) # could it be prot_nums? IT HAS THE LIGAND THOUGH
    
    setattr(state, "_protpsf", f"{dcdpath}/prot.psf")
    dcdpsf = state._protpsf
    if not os.path.isfile(dcdpsf):
        import parmed
        print(f"Making dcd trajectories for {state._pdbf}")
        psf = parmed.load_file(state._psff)[atoms.indices]
        psf.title = state._psff
        psf.write_psf(dcdpsf)


    dcds = [f"{dcdpath}/{xtc}.dcd" for xtc in state._trajs]
    setattr(state, "_dcds", {num: traj for num, traj in enumerate(dcds, 1)})

    import mdtraj
    for xtc, dcd in state._dcds.items():
        if not os.path.isfile(dcd):
            traj = mdtraj.load(state._trajs[xtc], top=state._pdbf, atom_indices=atoms.indices)
            traj.save_dcd(dcd)
    return






def get_bonded_cys(state, pdb):
    import parmed
    
    indices = []

    pdb = parmed.read_pdb(pdb)
    mask = parmed.amber.mask.AmberMask(pdb, ':CY?@SG')
    for sel in mask.Selected():
        atom = pdb.atoms[sel]
        bonded = [a.name for a in atom.bond_partners]
        if "SG" in bonded:
            indices.append(atom.residue.idx)

    return indices