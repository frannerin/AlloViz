import sys, os, re, io, requests

import MDAnalysis as mda
import numpy as np
from lazyasd import LazyObject

from . import utils



def _download_files(self):
    from multiprocess import Pool
    from bs4 import BeautifulSoup
    import urllib.request as pwget
    import tarfile, fileinput

    web = "https://submission.gpcrmd.org"

    html = requests.get(f"{web}/dynadb/dynamics/id/{self._gpcrmdid}/")
    soup = BeautifulSoup(html.text, features="html.parser").find_all('a')
    links = [link.get('href') for link in soup if re.search("(xtc|pdb|psf|prm)", link.get('href'))]
    get_name = lambda link: f"{self._path}/{link.rsplit('/')[-1]}"


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
        tar.extractall(self._path)
    os.remove(fname)

    for line in fileinput.input(f"{self._path}/parameters", inplace=True):
        if line.strip().startswith('HBOND'):
            line = 'HBOND CUTHB 0.5\n'
        elif line.strip().startswith('END'):
            line = 'END\n'
        sys.stdout.write(line)

    return




def _get_mdau(self):
    mdau = mda.Universe(self._pdbf, *self._trajs.values())
    prot = mdau.select_atoms("protein")

    if self.GPCR:#hasattr(self, "_gpcrmdid"):
        prot_numsf = f"{self._datadir}/gpcrdb_gennums.pdb"

        if not os.path.isfile(prot_numsf):
            print(f"retrieving {prot_numsf}")
            with io.StringIO() as protf:
                with mda.lib.util.NamedStream(protf, "prot.pdb") as f:
                    prot.write(f, file_format="pdb")
                response = requests.post('https://gpcrdb.org/services/structure/assign_generic_numbers', 
                                         files = {'pdb_file': protf.getvalue()})
                with open(prot_numsf, "w") as prot_nums:
                    prot_nums.write(response.text)

        nums = mda.Universe(prot_numsf).select_atoms("protein").tempfactors
        prot.tempfactors = nums.round(2)
    
    
    res_arrays = np.split(prot.residues.resindices, np.where(np.diff(prot.residues.resnums) != 1)[0]+1)
    dihedral_residx = lambda end=-1: [elem for arr in res_arrays for elem in arr[1:end]]

    return mdau, dihedral_residx




def _add_comtrajs(self):
    compath = f"{self._datadir}/COMtrajs"
    if not os.path.isdir(compath): os.makedirs(compath, exist_ok=True)

    compdb = f"{compath}/ca.pdb"
    if not os.path.isfile(compdb):
        print(f"Making trajectories of residue COM for {self._pdbf}")
        prot = self.mdau.select_atoms("protein and name CA")
        prot.write(compdb)
    setattr(self, "_compdbf", compdb)

    comtrajs = [f"{compath}/{xtc}.xtc" for xtc in self._trajs]
    setattr(self, "_comtrajs", {num: traj for num, traj in enumerate(comtrajs, 1)})

    for xtc, comtraj in self._comtrajs.items():
        if not os.path.isfile(comtraj):
            prot = self.mdau.select_atoms("protein")
            traj = next(traj for traj in self.mdau.trajectory.readers if traj.filename == self._trajs[xtc])
            arr = np.empty((prot.n_residues, traj.n_frames, 3))
            for ts in traj:
                arr[:, ts.frame] = prot.center_of_mass(compound='residues')

            cau = mda.Universe(compdb, arr, format=mda.coordinates.memory.MemoryReader, order='afc')
            with mda.Writer(comtraj, cau.atoms.n_atoms) as W:
                for ts in cau.trajectory:
                    W.write(cau.atoms)
    return




def _make_dcds(self):        
    dcdpath = f"{self._datadir}/dcds"
    if not os.path.isdir(dcdpath): os.makedirs(dcdpath, exist_ok=True)
    setattr(self, "_protf", lambda ext: f"{dcdpath}/prot.{ext}")

    atoms = self.mdau.select_atoms("protein")

    dcdpdb = self._protf("pdb")
    if not os.path.isfile(dcdpdb): atoms.write(dcdpdb) # could it be prot_nums? IT HAS THE LIGAND THOUGH

    dcdpsf = self._protf("psf")
    if not os.path.isfile(dcdpsf):
        import parmed
        print(f"Making dcd trajectories for {self._pdbf}")
        psf = parmed.load_file(self._psff)[atoms.indices]
        psf.title = self._psff
        psf.write_psf(dcdpsf)


    dcds = [f"{dcdpath}/{xtc}.dcd" for xtc in self._trajs]
    setattr(self, "_dcds", {num: traj for num, traj in enumerate(dcds, 1)})

    import mdtraj
    for xtc, dcd in self._dcds.items():
        if not os.path.isfile(dcd):
            traj = mdtraj.load(self._trajs[xtc], top=self._pdbf, atom_indices=atoms.indices)
            traj.save_dcd(dcd)
    return