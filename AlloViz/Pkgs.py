import sys, os, pandas, time#, importlib
import numpy as np
from .utils import *
from contextlib import redirect_stdout, redirect_stderr

from importlib import import_module
from lazyasd import LazyObject

imports = {
"_getcontacts_contacts": ".Packages.getcontacts.get_dynamic_contacts",
"_getcontacts_freqs": ".Packages.getcontacts.get_contact_frequencies",
"_dynetan": ".Packages.dynetan.dynetan.proctraj",
"_corrplus": ".Packages.correlationplus.correlationplus.calculate",
"_mdtask": ".Packages.MD-TASK.mdtask.calc_correlation",
"_pytraj": "pytraj",
"_pyinteraph": "pyinteraph.main",
"_grinn_args": ".Packages.gRINN_Bitbucket.source.grinn",
"_grinn_calc": ".Packages.gRINN_Bitbucket.source.calc",
"_grinn_corr": ".Packages.gRINN_Bitbucket.source.corr",
"_npeet_lnc": ".Packages.NPEET_LNC.lnc",
"_mda_dihedrals": "MDAnalysis.analysis.dihedrals",
"_mdtraj": "mdtraj",
"_mdentropy": ".Packages.mdentropy.mdentropy.metrics",
}

_extra_arg = lambda val: "package='AlloViz'" if 'Packages' in val else ''

for key, val in imports.items():
    exec(f"{key} = LazyObject(lambda: import_module('{val}', {_extra_arg(val)}), globals(), '{key}')")
    
    
    



class Pkg:
    def __new__(cls, state, d):#**kwargs):
        #print("new", os.getpid(), state, dir(state), dir())
        new = super().__new__(cls)
        new._name = new.__class__.__name__
        new.state = state
        new._d = d
        
        new._pdbf = d["_pdbf"]
        new._trajs = d["_trajs"]
        
        new._path = f"{d['_datadir']}/{new._name}/raw"
        os.makedirs(new._path, exist_ok=True)
        new._rawpq = lambda xtc: f"{new._path}/{xtc if isinstance(xtc, int) else xtc.rsplit('.', 1)[0]}.pq"
        
        return new
    
    
    def __getnewargs__(self):
        return self.state, self._d
        
        
    def __init__(self, *args):
        pqs = [self._rawpq(xtc) for xtc in self._trajs]
        no_exist = lambda pqs: [not os.path.isfile(pq) for pq in pqs]
        
        if any(no_exist(pqs)):
            for xtc in (xtc for xtc in self._trajs if no_exist(pqs)[xtc-1]):
                self._calculate(xtc)
                # time.sleep(5)
                
        
        def wait_calculate(pqs):
            while any(no_exist(pqs)):
                time.sleep(5)
            return pqs
                
        def get_raw(pqs):
            print(f"adding raw data of {self._name} for {self._pdbf}: ", pqs)
            flist = map(lambda pq: pandas.read_parquet(pq), pqs)
            df = pandas.concat(flist, axis=1)
            cols = [f"{num}" for num in self._trajs]
            df["weight_avg"] = df[cols].fillna(0).mean(axis=1)
            df["weight_std"] = df[cols].fillna(0).std(axis=1)
            return df
        
        add_raw = lambda pqs: setattr(self, "raw", get_raw(pqs))
        get_pool().apply_async(wait_calculate,
                               args=(pqs,),
                               callback=add_raw)
        
        
    
    def _calculate(self, xtc, *args):
        def send_and_log(xtc, *args):
            #print(f"sending {xtc}", os.getpid())
            with open(f"{self._path}/{self._name}.log", "a+") as f:
                with redirect_stdout(f), redirect_stderr(f):
                    return self._computation(xtc, *args)
                
        get_pool().apply_async(send_and_log,
                         args=(xtc, *args),
                         callback=self._save_pq)
        
    # def _computation
    # def _save_pq
    
    
    


class Multicorepkg(Pkg):
    def __new__(cls, state, d):
        new = super().__new__(cls, state, d)
        
        if "taskcpus" not in new._d:
            new.taskcpus = int(np.ceil(os.cpu_count()/2))
        else:
            new.taskcpus = new._d["taskcpus"]
        
        new._empties = new.taskcpus-1
            
        return new
    
    
    
    def _calculate(self, xtc, *args):
        super()._calculate(xtc, *args)
        
        
        def calculate_empty(pq):
            # print("sleeping", pq, os.getpid())
            while not os.path.isfile(pq):
                time.sleep(5)
            return
        
        for _ in range(self._empties): get_pool().apply_async(calculate_empty, args=(self._rawpq(xtc),))
    
    
    






class Matrixoutput(Pkg):
    def __new__(cls, state, d):
        new = super().__new__(cls, state, d)
        new._selection = "protein"
        return new
    
    
    def _save_pq(self, args):
        corr, xtc, *resl = args
        
        if len(resl) != 0:
            resl = resl[0]
            if corr.shape != (len(resl), len(resl)):
                corr = corr[np.ix_(resl, resl)]
        elif len(resl) == 0:
            resl = slice(0, corr.shape[0])
            
        resnames = [f"{aa.resname}:{aa.resid}" for aa in self._d["mdau"].select_atoms(self._selection).residues[resl]]
        
        df = pandas.DataFrame(corr, columns=resnames, index=resnames)
        df = df.where( np.triu(np.ones(df.shape), k=1).astype(bool) )
        df = pandas.DataFrame({f"{xtc}": df.stack()})
        # if not len(df[f"{xtc}"].unique()) == 1:
        df.to_parquet(self._rawpq(xtc))
        
                
        
        
class CombinedDihs(Matrixoutput):
    
    def _calculate(self, xtc):
        pkg = self._name.replace("Dihs", "")
        Dihl = ["Phi", "Psi", "Omega"]
        get_rawpq = lambda Dih: rgetattr(self, "state", f"{pkg}{Dih}", "_rawpq")(xtc)
        no_exist = lambda Dihl: [not rhasattr(self, "state", f"{pkg}{Dih}") for Dih in Dihl]
        
        if any(no_exist(Dihl)):
            for Dih in (Dih for Dih in Dihl if no_exist(Dihl)[Dihl.index(Dih)]):
                pkgclass = eval(capitalize(f"{pkg}{Dih}")) if isinstance(pkg, str) else pkg
                setattr(self.state, pkgclass.__name__, pkgclass(self.state, self._d))

        
        def wait_calculate(Dihl):
            not_finished = lambda Dihl: [not os.path.isfile(get_rawpq(Dih)) for Dih in Dihl]
            while any(not_finished(Dihl)):
                time.sleep(5)
            return Dihl
                
        def save_pq(Dihl):
            dfs = [pandas.read_parquet(get_rawpq(Dih))[f"{xtc}"].abs() for Dih in Dihl]
            
            final = None
            for df in dfs:
                if final is None:
                    final = df
                else:
                    final = final + df
            df = final / len(Dihl) # average of the absolute number

            # df = (final - final.min()) / (final.max() - final.min()) # This is done column-wise # This would be needed for absolute number sum; we are doing averaging
            pandas.DataFrame(df).to_parquet(self._rawpq(xtc))
            return
            
            
        get_pool().apply_async(wait_calculate,
                         args=(Dihl,),
                         callback=save_pq)
        
        




class COMpkg(Pkg):
    def __new__(cls, state, d):
        new = super().__new__(cls, state, d)
        
        new._pdbf = new._d["_compdbf"]
        new._trajs = new._d["_comtrajs"]
        
        return new
    
    

    
    
class dcdpkg(Matrixoutput):
    def __new__(cls, state, d):
        new = super().__new__(cls, state, d)
        
        new._pdbf = new._d["_protf"]("pdb")
        new._trajs = new._d["dcds"]
        
        return new




class Getcontacts(Multicorepkg):
            
    def _computation(self, xtc):#pdb, traj, xtc, pq, ctcs, freqs, taskcpus):
        path = self._path
        ctcs = f"{path}/{xtc}.tsv"
        freqs = f"{path}/{xtc}_freqs.tsv"
        
        if not os.path.isfile(freqs):# or ow:
            _getcontacts_contacts.main(f"--topology {self._pdbf} --trajectory {self._trajs[xtc]} --output {ctcs} --itypes all --cores {self.taskcpus}".split())
            _getcontacts_freqs.main(f"--input_files {ctcs} --output_file {freqs}".split())
        return freqs, xtc
        
        
    def _save_pq(self, args):
        freqs, xtc = args
        
        df = pandas.read_csv(freqs, sep="\t", skiprows=2,
                             index_col = (0, 1), names = [f"{xtc}"])
        df.index = df.index.map(lambda idx: tuple(sorted([res.split(":", 1)[-1] for res in idx], key = lambda res: int(res.split(":")[-1]))))
        df.to_parquet(self._rawpq(xtc))
        
    
    
    




class Dynetan(Matrixoutput, Multicorepkg):    
    
    def _computation(self, xtc):# pdb, traj, xtc, pq, taskcpus):
        obj = _dynetan.DNAproc()
        obj.loadSystem(self._pdbf, self._trajs[xtc]) # pdb.replace("pdb", "psf")
        prot = obj.getU().select_atoms("protein")

        protseg = list(prot.segments.segids)
        obj.setSegIDs(protseg)
        obj.selectSystem(withSolvent=False, userSelStr=f"protein")

        obj.setCustomResNodes({})
        obj.setUsrNodeGroups({})

        obj.setNumWinds(1)
        obj.alignTraj(inMemory=False)
        obj.prepareNetwork()
        obj.contactMatAll = np.triu(np.ones([obj.numWinds, obj.numNodes, obj.numNodes], dtype=int), k=1)

        obj.calcCor(ncores=self.taskcpus)
        
        return obj.corrMatAll[0], xtc



class DynetanCOM(COMpkg, Dynetan):
    pass






class Corrplus(Matrixoutput):    
    
    def _computation(self, xtc):
        corr = _corrplus.calcMDnDCC(self._pdbf, self._trajs[xtc], saveMatrix = False)
        return corr, xtc

        
        
        
class CorrplusLMI(Corrplus):
    
    def _computation(self, xtc):#pdb, traj, xtc, pq):
        corr = _corrplus.calcMD_LMI(self._pdbf, self._trajs[xtc], saveMatrix = False)
        return corr, xtc
        
        
        
class CorrplusCOM(COMpkg, Corrplus):
    pass
        

        
        
class CorrplusCOMLMI(COMpkg, CorrplusLMI):
    pass

        
        
        
class CorrplusPsi(Corrplus):
    def __new__(cls, state, d):
        new = super().__new__(cls, state, d)
        new._dih = "psi"
        return new
    
    def _computation(self, xtc):#pdb, traj, xtc, pq):
        corr = _corrplus.calcMDsingleDihedralCC(self._pdbf, self._trajs[xtc], dihedralType = self._dih, saveMatrix = False) # outputs a n_res x n_res matrix nevertheless
        return corr, xtc, self._d["_dihedral_residx"]()#[1, -1]
    
    
    
class CorrplusPhi(CorrplusPsi, Corrplus):
    def __new__(cls, state, d):
        new = super().__new__(cls, state, d)
        new._dih = "phi"
        return new

        

class CorrplusOmega(CorrplusPsi, Corrplus):
    def __new__(cls, state, d):
        new = super().__new__(cls, state, d)
        new._dih = "omega"
        return new

        
        
class CorrplusDihs(CombinedDihs, Corrplus):
    pass

        
        
        
        
        
        
class AlloVizPsi(Matrixoutput):
    def __new__(cls, state, d):
        new = super().__new__(cls, state, d)
        new._dih = "psi"
        return new
            
    
    def _computation(self, xtc):#pdb, traj, xtc, pq):
        prot = self._d["mdau"].select_atoms("protein")
        selected_res = self._d["_dihedral_residx"]()
        
        select_dih = lambda res: eval(f"res.{self._dih.lower()}_selection()")
        selected = [select_dih(res) for res in prot.residues[selected_res]]

        offset = 0
        get_frames = lambda numreader: self._d["mdau"].trajectory.readers[numreader].n_frames
        for numreader in range(xtc-1):
            offset += get_frames(numreader)

        values = _mda_dihedrals.Dihedral(selected).run(start=offset, stop=offset+get_frames(xtc-1)).results.angles.transpose()
        
        corr = np.zeros((len(selected_res), len(selected_res)))
        print("prot.n_residues:", prot.n_residues, "len(selected_res)", len(selected_res), "dihedrals array shape:", values.shape, "empty corr matrix shape:", corr.shape)

        iterator = set(range(len(selected_res)))
        for res1 in iterator:
            for res2 in iterator - set(range(res1)) - {res1}:
                corr[res1, res2] = _npeet_lnc.MI.mi_LNC([values[res1], values[res2]])
    
        return corr, xtc, selected_res
    
    
    
class AlloVizPhi(AlloVizPsi):
    def __new__(cls, state, d):
        new = super().__new__(cls, state, d)
        new._dih = "phi"
        return new

        

class AlloVizOmega(AlloVizPsi):
    def __new__(cls, state, d):
        new = super().__new__(cls, state, d)
        new._dih = "omega"
        return new

        
        
        
    
class AlloVizDihs(CombinedDihs, AlloVizPsi):
    pass
        

        

        
        
        
class MDEntropyContacts(Matrixoutput, Multicorepkg):
    def __new__(cls, state, d):
        new = super().__new__(cls, state, d)
        
        # # Opción si las necesidades de memoria incrementan con taskcpus
        # extra_taskcpus = int((self.taskcpus/4 - 1) * 4) if self.taskcpus>=4 else 0
        # taskcpus = 1 + extra_taskcpus # Minimum 4*2000 of memory (this taskcpus is like that to use 4 cpus-2000 mem in .sh files)
        # empties = 3
        
        # # Opción si las necesidades altas de memoria sólo son con el inicio del cálculo
        # extra_taskcpus = int((self.taskcpus/4 - 1) * 4) if self.taskcpus>=4 else 0
        # taskcpus = 1 + extra_taskcpus # Minimum 4*2000 of memory (this taskcpus is like that to use 4 cpus-2000 mem in .sh files)
        # empties = 3 - extra_taskcpus if extra_taskcpus<=3 else 0

        # Opción de 1 empty por tasckpu, pasando las extras a empties; 1 taskcpu requiere 3,2G; hay una necesidad mayor de memoria al principio pero i pretend i do not see
        half = int(np.floor(new.taskcpus/2))
        new.taskcpus = half if new.taskcpus > 1 else 1
        new._empties = half
        
        new._function = _mdentropy.ContactMutualInformation
        new._types = {}
        new._resl = []
        return new
        
        
        
    def _computation(self, xtc):#pdb, traj, xtc, pq, taskcpus):
        mytraj = _mdtraj.load(self._trajs[xtc], top=self._pdbf) # hopefully mdtraj is loaded from the Classes module
        mi = self._function(threads=self.taskcpus, **self._types) # n_bins=3, method='knn', normed=True
        corr = mi.partial_transform(traj=mytraj, shuffle=0, verbose=True)
        return corr, xtc, self._resl
    
        
        
        
class MDEntropyDihs(MDEntropyContacts):
    def __new__(cls, state, d):
        new = super().__new__(cls, state, d)        
        new._function = _mdentropy.DihedralMutualInformation
        new._types = {"types": ["phi", "psi", "omega"]}
        new._resl = new._d["_dihedral_residx"]()
        return new
    
    

class MDEntropyAlphaAngle(MDEntropyContacts):
    """
    The alpha angle of residue `i` is the dihedral formed by the four CA atoms
    of residues `i-1`, `i`, `i+1` and `i+2`.
    """
    def __new__(cls, state, d):
        new = super().__new__(cls, state, d)        
        new._function = _mdentropy.AlphaAngleMutualInformation
        new._types = {"types": ["alpha"]}
        new._resl = new._d["_dihedral_residx"](-2)
        return new
    
    
    def _computation(self, xtc):
        corr, xtc, dihedral_resl = super()._computation(xtc)
        
        length = corr.shape[0]
        to_insert = {0, length+1, length+2}
        new_length = length + len(to_insert)
        new_indices = set(range(new_length))
        
        new_corr = np.zeros((new_length, new_length), dtype=corr.dtype)
        corr_indices = np.array(list(new_indices - to_insert))
        new_corr[corr_indices.reshape(-1,1), corr_indices] = corr
        
        return new_corr, xtc, dihedral_resl
        
        
        
        
    
    
    
    

class MDTASK(Matrixoutput, Multicorepkg):
    def __new__(cls, state, d):
        new = super().__new__(cls, state, d)        
        new._empties = 2
        return new
    
    
    def _computation(self, xtc):#pdb, traj, xtc, pq):
        corr = _mdtask.correlate(_mdtask.parse_traj(traj = self._trajs[xtc], topology = self._pdbf))
        return corr, xtc









class PytrajCA(Matrixoutput):
    def __new__(cls, state, d):
        new = super().__new__(cls, state, d)
        new._mask = new._name[-2:]
        return new
    
    
    def _computation(self, xtc):#pdb, traj, mask, xtc, pq):
        top = _pytraj.load_topology(self._pdbf)
        traj = _pytraj.load(self._trajs[xtc], top, mask = f'@{self._mask}')
        corr = _pytraj.matrix.correl(traj, f'@{self._mask}')
        return corr, xtc

    
    
class PytrajCB(PytrajCA):
    def __new__(cls, state, d):
        new = super().__new__(cls, state, d)
        new._selection = "protein and not resname GLY"
        return new

        
        
        
        
        
        
class PyInteraph(Matrixoutput):
    def __new__(cls, state, d):
        new = super().__new__(cls, state, d)
        new._CLIargs = "-m --cmpsn-graph dummy"
        return new
                
        
    def _computation(self, xtc):#pdb, traj, xtc, pq, CLIargs):
        corr = _pyinteraph.main(f"-s {self._pdbf} -t {self._trajs[xtc]} {self._CLIargs}".split()) / 100
        return corr, xtc

    
    
class PyInteraphEne(PyInteraph):
    def __new__(cls, state, d):
        new = super().__new__(cls, state, d)
        new._CLIargs = "-p --kbp-graph dummy"
        return new

        
        


# only for local
class G_corrCAMI(Matrixoutput):
    def __new__(cls, state, d):
        new = super().__new__(cls, state, d)
        new._CLIargs = ""
        return new
                
        
    def _computation(self, xtc):#pdb, traj, xtc, pq):
        # Send g_correlation
        pq = self._rawpq(xtc)
        if not os.path.isfile(f"{pq}.dat"):
            os.system(f"""
module load g_correlation
export GMXLIB=/soft/EB_repo/bio/sequence/programs/noarch/gromacs/3.3.1/share/gromacs/top/
g_correlation -f {self._trajs[xtc]} -s {self._pdbf} -o {pq}.dat {self._CLIargs} &> {self._path}/{xtc}.log <<EOF
1
3
EOF
""")
        # Read output.dat
        size = self._d["mdau"].select_atoms("protein and name CA").n_atoms
        corr = np.empty([size, size])
        rown = 0
        row = []
        
        with open(f"{pq}.dat", "r") as f:
            for num, line in enumerate(f):
                if num == 0:
                    pass
                else:
                    row.extend( line.strip().split() )

                if len(row) >= size:
                    corr[rown, :] = row[:size]
                    rown += 1
                    del row[:size]
        
        
        return corr, xtc

    
    
    
class G_corrCOMMI(COMpkg, G_corrCAMI):
    def __new__(cls, state, d):
        new = super().__new__(cls, state, d)
        new._CLIargs = ""
        return new
        
        
        
        
class G_corrCALMI(G_corrCAMI):
    def __new__(cls, state, d):
        new = super().__new__(cls, state, d)
        new._CLIargs = "-linear"
        return new
    

    
class G_corrCOMLMI(COMpkg, G_corrCALMI):
    pass
        

        

        
# only for local
class GSAtools(Matrixoutput):
        
    def _computation(self, xtc):# pdb, traj, xtc, pq, out):
        # Send gsatools
        out = f"{self._path}/{xtc}"
        if not os.path.isfile(f"{out}/lf_nMImat.out"):
            os.system(f"""
module load GSAtools
g_sa_encode -s {self._pdbf} -f {self._trajs[xtc]} -rmsdlf {out}/lf_rmsd.xvg -strlf {out}/lf_str.out -log {out}/log.log 
g_sa_analyze -sa {out}/lf_str.out -MImatrix -MImat {out}/lf_MImat.out -eeMImat {out}/lf_eeMImat.out -jHmat {out}/lf_jHmat.out -nMImat {out}/lf_nMImat.out >>{out}/{xtc}.log 2>&1 
""")
        # Read output
        corr = np.loadtxt(f"{out}/lf_nMImat.out")        
        
        return corr, xtc      
        
        
        
        

        
        
class GRINN(dcdpkg, Multicorepkg):
    def __new__(cls, state, d):
        new = super().__new__(cls, state, d)
        
        if "namd" in new._d:
            new.namd = new._d["namd"]
        else:
            from distutils.spawn import find_executable
            new.namd = find_executable('namd2')
            if new.namd is None:
                raise Exception("namd executable for gRINN computation not found")
                
        if "auto_send" not in new._d:
            d = new._d.copy()
            d.update({"namd": new.namd, "auto_send": True})
            new.state._set_pkgclass(self.state, "gRINNcorr", d)
        
        return new
        
        
    def _computation(self, xtc):# pdb, traj, out, xtc, pq, psf, params, taskcpus):
        psf = self._d["_protf"]("psf")
        params = self._d["_paramf"]
        out = f"{self._path}/{xtc}"
        outf = f"{out}/energies_intEnMeanTotal.dat"
        
        if not os.path.isfile(outf):
            if os.path.isdir(out):
                from shutil import rmtree
                rmtree(out)
            
            _grinn_calc.getResIntEn(_grinn_args.arg_parser(f"-calc --pdb {self._pdbf} --top {psf} --traj {self._trajs[xtc]} --exe {self.namd} --outfolder {out} --numcores {self.taskcpus} --parameterfile {params}".split()))
            
        corr = np.loadtxt(outf)
        return corr, xtc
    

    
    
    
    
class GRINNcorr(GRINN):
    def __new__(cls, state, d):
        new = super().__new__(cls, state, d)
        
        if "auto_send" in new._d:
            new._auto_send = new._d["auto_send"]
        else:
            new._auto_send = False
        
        return new
    
    
    def __init__(self):
        if not rhasattr(self, "state", "GRINN") and not self._auto_send:
            raise Exception("Make sure to send GRINN calculation before GRINNcorr")
        
        no_exist = lambda files: [not os.path.isfile(file) for file in files]
        
        pqs = [self._rawpq(xtc) for xtc in self._trajs]
        
        outf = lambda xtc: f"{self._path.replace('GRINNcorr', 'GRINN')}/{xtc}/energies_intEnTotal.csv"
        outfs = [outf(xtc) for xtc in self._trajs]
        
        
        def wait_calculate(pqs):
            while any(no_exist(outfs)):
                time.sleep(30)
                
            self._initialize_real()
            
            while any(no_exist(pqs)):
                time.sleep(5)
            return pqs
                
            
        def get_raw(pqs):
            print(f"adding raw data of {self._name} for {self._pdbf}: ", pqs)
            flist = map(lambda pq: pandas.read_parquet(pq), pqs)
            df = pandas.concat(flist, axis=1)
            cols = [f"{num}" for num in self._trajs]
            df["weight_avg"] = df[cols].fillna(0).mean(axis=1)
            df["weight_std"] = df[cols].fillna(0).std(axis=1)
            return df
        
        add_raw = lambda pqs: setattr(self, "raw", get_raw(pqs))
        get_pool().apply_async(wait_calculate,
                               args=(pqs,),
                               callback=add_raw)
            
        
    def _initialize_real(self):
        pqs = [self._rawpq(xtc) for xtc in self._trajs]
        no_exist = lambda pqs: [not os.path.isfile(pq) for pq in pqs]
        
        if any(no_exist(pqs)):
            for xtc in (xtc for xtc in self._trajs if no_exist(pqs)[xtc-1]):
                self._calculate(xtc)
        
    
        
    def _computation(self, xtc):#pdb, out, xtc, pq, taskcpus):
        out = f"{self._path}/{xtc}"
        logFile = f"{out}/grinncorr.log"
        os.system(f"mkdir -p {out}; touch {logFile}")
        outf = f"{out}/energies_resCorr.dat"
        
        if not os.path.isfile(outf):
            _grinn_corr.getResIntCorr(_grinn_args.arg_parser(f"-corr --pdb {self._pdbf} --corrinfile {out.replace('GRINNcorr', 'GRINN')}/energies_intEnTotal.csv --corrprefix {out}/energies --numcores {self.taskcpus}".split()), logFile=logFile)
            
        corr = np.loadtxt(outf) 
        return corr, xtc
    


# class Carma(): # needs the dcd
# class Bio3D(): # R package; needs the dcd
# class wordom(): # doesn't have python bindings for croscorr and lmi but if it had it would've been great because it looks fast?
