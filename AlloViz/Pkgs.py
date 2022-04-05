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
"_pyinteraph": ".Packages.pyinteraph2.pyinteraph.main",
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
    def __init__(self, state, **kwargs):
        self.state = state        
        self._name = self.__class__.__name__
        
        self._path = f"{self.state._datadir}/{self._name}/raw"
        os.makedirs(self._path, exist_ok=True)
        self._rawpq = lambda xtc: f"{self._path}/{xtc if isinstance(xtc, int) else xtc.rsplit('.', 1)[0]}.pq"
        
        self._initialize()
        
        
    
    def _initialize(self):
        pqs = [self._rawpq(xtc) for xtc in self.state._trajs]
        no_exist = lambda pqs: [not os.path.isfile(pq) for pq in pqs]
        
        if any(no_exist(pqs)):
            for xtc in (xtc for xtc in self.state._trajs if no_exist(pqs)[xtc-1]):
                self._calculate(xtc)
                
        
        def wait_calculate(pqs):
            while any(no_exist(pqs)):
                time.sleep(5)
            return pqs
                
        def get_raw(pqs):
            print(f"adding raw data of {self._name} for {self.state._pdbf}: ", pqs)
            flist = map(lambda pq: pandas.read_parquet(pq), pqs)
            df = pandas.concat(flist, axis=1)
            cols = [f"{num}" for num in self.state._trajs]
            df["weight_avg"] = df[cols].fillna(0).mean(axis=1)
            df["weight_std"] = df[cols].fillna(0).std(axis=1)
            return df
        
        add_raw = lambda pqs: setattr(self, "raw", get_raw(pqs))
        get_pool().apply_async(wait_calculate,
                               args=(pqs,),
                               callback=add_raw)
        return
        
    
    def _calculate(self, xtc):
        pool = get_pool()
        pdb = self.state._pdbf
        traj = self.state._trajs[xtc]
        pq = self._rawpq(xtc)
        print(f"making {pq}")
        return pool, pdb, traj, pq
    
    
    def _send_and_log(self, *args):
        with open(f"{self._path}/{self._name}.log", "a+") as f:
            with redirect_stdout(f), redirect_stderr(f):
                return self._computation(*args)
        
    
    # abstractmethod
    def _computation(self, xtc):
        return info, xtc, pq
    
    # abstractmethod
    def _save_pq(self, args):
        info, xtc, pq = args
        # df.to_parquet(pq)
    
    
    


class Multicorepkg(Pkg):
    def __init__(self, state, **kwargs):
        if "taskcpus" not in kwargs:
            self.taskcpus = 3
        else:
            self.taskcpus = kwargs["taskcpus"]
            
        super().__init__(state, **kwargs)
    
    def _calculate_empty(self, pqf):
        print("sleeping", pqf, os.getpid())
        while not os.path.isfile(pqf):
            time.sleep(5)
        return






class Matrixoutput(Pkg):
    def __init__(self, state, **kwargs):
        self._selection = "protein"
        super().__init__(state, **kwargs)
    
    
    def _save_pq(self, args):
        if len(args) == 3:
            corr, xtc, pq = args
            ids = [0, corr.shape[0]]
        else:
            corr, xtc, pq, ids = args
        slicing = slice(*ids)
        
        length = len(list(range(*ids)))
        if corr.shape != (length, length):
            corr = corr[slicing, slicing]
            
        resl = [f"{aa.resname}:{aa.resid}" for aa in self.state.mdau.select_atoms(self._selection).residues[slicing]]
        
        df = pandas.DataFrame(corr, columns=resl, index=resl)
        df = df.where( np.triu(np.ones(df.shape), k=1).astype(bool) )
        df = pandas.DataFrame({f"{xtc}": df.stack()})
        # if not len(df[f"{xtc}"].unique()) == 1:
        df.to_parquet(pq)






class COMpkg(Pkg):
    def __init__(self, state, **kwargs):
        if not hasattr(state, "_comtrajs"):
            self.state._add_comtrajs()
            
        super().__init__(state, **kwargs)
    
    
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        pdb = self.state._compdbf
        traj = self.state._comtrajs[xtc]
        
        return pool, pdb, traj, pq






class Getcontacts(Multicorepkg):
    def __init__(self, state, **kwargs):        
        super().__init__(state, **kwargs)
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        
        path = self._path
        ctcs = f"{path}/{xtc}.tsv"
        freqs = f"{path}/{xtc}_freqs.tsv"
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq, ctcs, freqs, self.taskcpus),
                         callback=self._save_pq)
        
        for _ in range(self.taskcpus-1): pool.apply_async(self._calculate_empty, args=(pq,))
        
        
    def _computation(self, pdb, traj, xtc, pq, ctcs, freqs, taskcpus):
        if not os.path.isfile(freqs):# or ow:
            _getcontacts_contacts.main(f"--topology {pdb} --trajectory {traj} --output {ctcs} --itypes all --cores {taskcpus}".split())
            _getcontacts_freqs.main(f"--input_files {ctcs} --output_file {freqs}".split())
        return freqs, xtc, pq
        
        
    def _save_pq(self, args):
        freqs, xtc, pq = args
        
        df = pandas.read_csv(freqs, sep="\t", skiprows=2,
                             index_col = (0, 1), names = [f"{xtc}"])
        df.index = df.index.map(lambda idx: tuple(sorted([res.split(":", 1)[-1] for res in idx], key = lambda res: int(res.split(":")[-1]))))
        df.to_parquet(pq)
        
    
    
    










class Dynetan(Matrixoutput, Multicorepkg):
    def __init__(self, state, **kwargs):
        super().__init__(state, **kwargs)
        
    
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq, self.taskcpus),
                         callback=self._save_pq)
        
        for _ in range(self.taskcpus-1): pool.apply_async(self._calculate_empty, args=(pq,))
    
    
    def _computation(self, pdb, traj, xtc, pq, taskcpus):
        obj = _dynetan.DNAproc()
        obj.loadSystem(pdb, traj) # pdb.replace("pdb", "psf")
        prot = obj.getU().select_atoms("protein")

        protseg = list(prot.segments.segids)
        obj.setSegIDs(protseg)
        obj.selectSystem(withSolvent=False, userSelStr=f"segid {protseg[0]}")

        obj.setCustomResNodes({})
        obj.setUsrNodeGroups({})

        obj.setNumWinds(1)
        obj.alignTraj(inMemory=False)
        obj.prepareNetwork()
        obj.contactMatAll = np.triu(np.ones([obj.numWinds, obj.numNodes, obj.numNodes], dtype=int), k=1)

        obj.calcCor(ncores=taskcpus)
        
        return obj.corrMatAll[0], xtc, pq



class DynetanCOM(Dynetan, COMpkg):
    def __init__(self, state, **kwargs):
        super().__init__(state, **kwargs)
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = COMpkg._calculate(self, xtc)
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq, self.taskcpus),
                         callback=self._save_pq)
        
        for _ in range(self.taskcpus-1): pool.apply_async(self._calculate_empty, args=(pq,))










class Corrplus(Matrixoutput):
    def __init__(self, state, **kwargs):
        super().__init__(state, **kwargs)
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq),
                         callback=self._save_pq)
    
    
    def _computation(self, pdb, traj, xtc, pq):
        corr = _corrplus.calcMDnDCC(pdb, traj, saveMatrix = False)
        return corr, xtc, pq

        
        
        
class CorrplusLMI(Corrplus):
    def __init__(self, state, **kwargs):
        super().__init__(state, **kwargs)
    
    def _computation(self, pdb, traj, xtc, pq):
        corr = _corrplus.calcMD_LMI(pdb, traj, saveMatrix = False)
        return corr, xtc, pq
        
        
        
class CorrplusCOM(Corrplus, COMpkg):
    def __init__(self, state, **kwargs):
        super().__init__(state, **kwargs)
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = COMpkg._calculate(self, xtc)
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq),
                         callback=self._save_pq)
        

        
        
class CorrplusCOMLMI(CorrplusLMI, COMpkg):
    def __init__(self, state, **kwargs):
        super().__init__(state, **kwargs)
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = COMpkg._calculate(self, xtc)
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq),
                         callback=self._save_pq)

        
        
        
class CorrplusPsi(Corrplus):        
    def __init__(self, state, **kwargs):
        self._dih = "psi"
        super().__init__(state, **kwargs)
    
    def _computation(self, pdb, traj, xtc, pq):
        # dih = "psi"
        corr = _corrplus.calcMDsingleDihedralCC(pdb, traj, dihedralType = self._dih, saveMatrix = False)
        return corr, xtc, pq, [1, -1]
    
    
#     def _save_pq(self, args):
#         corr, xtc, pq = args
        
#         corr = corr[1:-1, 1:-1]
#         resl = [f"A:{aa.resname}:{aa.resid}" for aa in self.state.mdau.select_atoms(self._selection).residues[1:-1]]

#         df = pandas.DataFrame(corr, columns=resl, index=resl)
#         df = df.where( np.triu(np.ones(df.shape), k=1).astype(bool) )
#         df = pandas.DataFrame({f"{xtc}": df.stack()})
#         # if not len(df[f"{xtc}"].unique()) == 1:
#         df.to_parquet(pq)

    
    
class CorrplusPhi(CorrplusPsi, Corrplus):
    def __init__(self, state, **kwargs):
        self._dih = "phi"
        Corrplus.__init__(self, state, **kwargs)

        

class CorrplusOmega(CorrplusPsi, Corrplus):        
    def __init__(self, state, **kwargs):
        self._dih = "omega"
        Corrplus.__init__(self, state, **kwargs)

        
        
        
    
class CorrplusDihs(Corrplus):        
    def __init__(self, state, **kwargs):
        super().__init__(state, **kwargs)
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super(Corrplus, self)._calculate(xtc)
        
        Dihl = ["Phi", "Psi", "Omega"]
        get_rawpq = lambda Dih: rgetattr(self, "state", f"Corrplus{Dih}", "_rawpq")(xtc)
        no_exist = lambda Dihl: [not rhasattr(self, "state", f"Corrplus{Dih}") for Dih in Dihl]
        
        if any(no_exist(Dihl)):
            for Dih in (Dih for Dih in Dihl if no_exist(Dihl)[Dihl.index(Dih)]):
                self.state._set_pkgclass(self.state, f"Corrplus{Dih}")

        
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
            pandas.DataFrame(df).to_parquet(pq)
            return
            
            
        pool.apply_async(wait_calculate,
                         args=(Dihl,),
                         callback=save_pq)

        
        
        
        
class AlloVizPsi(Matrixoutput):        
    def __init__(self, state, **kwargs):
        self._dih = "psi"
        super().__init__(state, **kwargs)
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq),
                         callback=self._save_pq)
        
    
    def _computation(self, pdb, traj, xtc, pq):
        select_dih = lambda res: eval(f"res.{self._dih.lower()}_selection()")
    
        prot = self.state.mdau.select_atoms("protein")
        atomgroups = [select_dih(res) for res in prot.residues]
        no_dihs = {idx for idx, group in enumerate(atomgroups) if group is None}

        is_terminal = lambda idx: True if idx == 0 or idx == prot.n_residues-1 else False
        if all([is_terminal(idx) for idx in no_dihs]):
            print("All missing dihedrals belong to terminal residues")
        if any([not is_terminal(idx) for idx in no_dihs]):
            raise Exception("There is a missing dihedral that does not belong to either of the terminal residues.")

        no_dihs = {0, prot.n_residues-1}
        selected = atomgroups[1:-1]

        offset = 0
        get_frames = lambda numreader: self.state.mdau.trajectory.readers[numreader].n_frames
        for numreader in range(xtc-1):
            offset += get_frames(numreader)

        values = _mda_dihedrals.Dihedral(selected).run(start=offset, stop=offset+get_frames(xtc-1)).results.angles.transpose()
        for idx in no_dihs:
            values = np.insert(values, idx, np.array([1]*values.shape[1]), axis=0)

        corr = np.zeros((prot.n_residues, prot.n_residues))
        print("prot.n_residues:", prot.n_residues, "dihedrals array shape:", values.shape, "empty corr matrix shape:", corr.shape)

        iterator = set(range(prot.n_residues)) - no_dihs
        for res1 in iterator:
            for res2 in iterator - set(range(res1)) - {res1}:
                corr[res1, res2] = _npeet_lnc.MI.mi_LNC([values[res1], values[res2]])
    
        return corr, xtc, pq, [1,-1]
    
    

    
    
class AlloVizPhi(AlloVizPsi):
    def __init__(self, state, **kwargs):
        self._dih = "phi"
        super().__init__(state, **kwargs)

        

class AlloVizOmega(AlloVizPsi):        
    def __init__(self, state, **kwargs):
        self._dih = "omega"
        super().__init__(state, **kwargs)

        
        
        
    
class AlloVizDihs(AlloVizPsi):        
    def __init__(self, state, **kwargs):
        super().__init__(state, **kwargs)
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super(AlloVizPsi, self)._calculate(xtc)
        
        Dihl = ["Phi", "Psi", "Omega"]
        get_rawpq = lambda Dih: rgetattr(self, "state", f"AlloViz{Dih}", "_rawpq")(xtc)
        no_exist = lambda Dihl: [not rhasattr(self, "state", f"AlloViz{Dih}") for Dih in Dihl]
        
        if any(no_exist(Dihl)):
            for Dih in (Dih for Dih in Dihl if no_exist(Dihl)[Dihl.index(Dih)]):
                self.state._set_pkgclass(self.state, f"AlloViz{Dih}")

        
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
            pandas.DataFrame(df).to_parquet(pq)
            return
            
            
        pool.apply_async(wait_calculate,
                         args=(Dihl,),
                         callback=save_pq)
        

        

        
        
        
class MDEntropyContacts(Matrixoutput, Multicorepkg):        
    def __init__(self, state, **kwargs):
        super().__init__(state, **kwargs)
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        
        # # Opción si las necesidades de memoria incrementan con taskcpus
        # extra_taskcpus = int((self.taskcpus/4 - 1) * 4) if self.taskcpus>=4 else 0
        # taskcpus = 1 + extra_taskcpus # Minimum 4*2000 of memory (this taskcpus is like that to use 4 cpus-2000 mem in .sh files)
        # empties = 3
        
        # # Opción si las necesidades altas de memoria sólo son con el inicio del cálculo
        # extra_taskcpus = int((self.taskcpus/4 - 1) * 4) if self.taskcpus>=4 else 0
        # taskcpus = 1 + extra_taskcpus # Minimum 4*2000 of memory (this taskcpus is like that to use 4 cpus-2000 mem in .sh files)
        # empties = 3 - extra_taskcpus if extra_taskcpus<=3 else 0

        # Opción de 1 empty por tasckpu, pasando las extras a empties; 1 taskcpu requiere 3,2G; hay una necesidad mayor de memoria al principio pero i pretend i do not see
        half = int(np.floor(self.taskcpus/2))
        taskcpus = half if self.taskcpus > 1 else 1
        empties = half
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq, taskcpus),
                         callback=self._save_pq)
        
        for _ in range(empties): pool.apply_async(self._calculate_empty, args=(pq,))
        
    
    def _computation(self, pdb, traj, xtc, pq, taskcpus):
        mytraj = _mdtraj.load(traj, top=pdb) # hopefully mdtraj is loaded from the Classes module
        mi = _mdentropy.ContactMutualInformation(threads=taskcpus) # n_bins=3, method='knn', normed=True
        corr = mi.partial_transform(traj=mytraj, shuffle=0, verbose=True)
        return corr, xtc, pq
    
        
        
        
class MDEntropyDihs(MDEntropyContacts):        
    def __init__(self, state, **kwargs):
        super().__init__(state, **kwargs)
        
    
    def _computation(self, pdb, traj, xtc, pq, taskcpus):
        mytraj = _mdtraj.load(traj, top=pdb) # hopefully mdtraj is loaded from the Classes module
        mi = _mdentropy.DihedralMutualInformation(types=['phi', 'psi', 'omega'], threads=taskcpus) # n_bins=3, method='knn', normed=True
        corr = mi.partial_transform(traj=mytraj, shuffle=0, verbose=True)
        return corr, xtc, pq, [1, -1]
    
    

        
    
    

class MDEntropyAlphaAngle(MDEntropyContacts):        
    def __init__(self, state, **kwargs):
        """
        The alpha angle of residue `i` is the dihedral formed by the four CA atoms
        of residues `i-1`, `i`, `i+1` and `i+2`.
        """
        super().__init__(state, **kwargs)
        
    
    def _computation(self, pdb, traj, xtc, pq, taskcpus):
        mytraj = _mdtraj.load(traj, top=pdb) # hopefully mdtraj is loaded from the Classes module
        mi = _mdentropy.AlphaAngleMutualInformation(threads=taskcpus) # n_bins=3, method='knn', normed=True
        corr = mi.partial_transform(traj=mytraj, shuffle=0, verbose=True)
        return corr, xtc, pq, [mi.labels[0], mi.labels[-1]+1]
        
        
    
    
    
    

class MDTASK(Matrixoutput, Multicorepkg):
    def __init__(self, state, **kwargs):
        super().__init__(state, **kwargs)
    
    
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq),
                         callback=self._save_pq)
        
        for _ in range(2): pool.apply_async(self._calculate_empty, args=(pq,)) # Just 1 empty job to use its memory
    
    
    def _computation(self, pdb, traj, xtc, pq):
        corr = _mdtask.correlate(_mdtask.parse_traj(traj = traj, topology = pdb))
        return corr, xtc, pq









class PytrajCA(Matrixoutput):
    def __init__(self, state, **kwargs):
        super().__init__(state, **kwargs)
    
    
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        mask = self.__class__.__name__[-2:]
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, mask, xtc, pq),
                         callback=self._save_pq)
    
    
    def _computation(self, pdb, traj, mask, xtc, pq):
        top = _pytraj.load_topology(pdb)
        traj = _pytraj.load(traj, top, mask = f'@{mask}')
        corr = _pytraj.matrix.correl(traj, f'@{mask}')
        return corr, xtc, pq

    
    
class PytrajCB(PytrajCA):
    def __init__(self, state, **kwargs):
        super().__init__(state, **kwargs)
    
    
    def _calculate(self, xtc):
        pool, pdb, traj, pq = Matrixoutput._calculate(self, xtc)
        mask = self.__class__.__name__[-2:]
        self._selection = "protein and not resname GLY"
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, mask, xtc, pq),
                         callback=self._save_pq)


        
        
        
        
        
        
class PyInteraph(Matrixoutput):
    def __init__(self, state, **kwargs):        
        super().__init__(state, **kwargs)
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        
        CLIargs = "-m --cmpsn-graph dummy"
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq, CLIargs),
                         callback=self._save_pq)
        
        
    def _computation(self, pdb, traj, xtc, pq, CLIargs):
        corr = _pyinteraph.main(f"-s {pdb} -t {traj} {CLIargs}".split()) / 100
        return corr, xtc, pq

    
    
class PyInteraphEne(PyInteraph):
    def __init__(self, state, **kwargs):        
        super().__init__(state, **kwargs)
        
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = Matrixoutput._calculate(self, xtc)
        
        CLIargs = "-p --kbp-graph dummy"
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq, CLIargs),
                         callback=self._save_pq)

        
        


# only for local
class G_corrCAMI(Matrixoutput):
    def __init__(self, state, **kwargs):        
        super().__init__(state, **kwargs)
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        CLIargs = f"-f {traj} -s {pdb} -o {pq}.dat"
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq, CLIargs),
                         callback=self._save_pq)
        
        
    def _computation(self, pdb, traj, xtc, pq, CLIargs):
        # Send g_correlation
        if not os.path.isfile(f"{pq}.dat"):
            os.system(f"""
module load g_correlation
export GMXLIB=/soft/EB_repo/bio/sequence/programs/noarch/gromacs/3.3.1/share/gromacs/top/
g_correlation {CLIargs} &> {self._path}/{xtc}.log <<EOF
1
3
EOF
""")
        # Read output.dat
        size = self.state.mdau.select_atoms("protein and name CA").n_atoms
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
        
        
        return corr, xtc, pq

    
    
    
class G_corrCOMMI(G_corrCAMI, COMpkg):
    def __init__(self, state, **kwargs):        
        super().__init__(state, **kwargs)
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = COMpkg._calculate(self, xtc)
        CLIargs = f"-f {traj} -s {pdb} -o {pq}.dat"
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq, CLIargs),
                         callback=self._save_pq)
        
        
        
        
class G_corrCALMI(G_corrCAMI):
    def __init__(self, state, **kwargs):        
        super().__init__(state, **kwargs)
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super(G_corrCAMI, self)._calculate(xtc)
        CLIargs = f"-f {traj} -s {pdb} -o {pq}.dat -linear"
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq, CLIargs),
                         callback=self._save_pq)
    

    
class G_corrCOMLMI(G_corrCAMI, COMpkg):
    def __init__(self, state, **kwargs):        
        super().__init__(state, **kwargs)
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = COMpkg._calculate(self, xtc)
        CLIargs = f"-f {traj} -s {pdb} -o {pq}.dat -linear"
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq, CLIargs),
                         callback=self._save_pq)
        

        

        
# only for local
class GSAtools(Matrixoutput):
    def __init__(self, state, **kwargs):        
        super().__init__(state, **kwargs)
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        out = f"{self._path}/{xtc}"
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq, out),
                         callback=self._save_pq)
        
        
    def _computation(self, pdb, traj, xtc, pq, out):
        # Send gsatools
        if not os.path.isfile(f"{out}/lf_nMImat.out"):
            os.system(f"""
module load GSAtools
g_sa_encode -s {pdb} -f {traj} -rmsdlf {out}/lf_rmsd.xvg -strlf {out}/lf_str.out -log {out}/log.log 
g_sa_analyze -sa {out}/lf_str.out -MImatrix -MImat {out}/lf_MImat.out -eeMImat {out}/lf_eeMImat.out -jHmat {out}/lf_jHmat.out -nMImat {out}/lf_nMImat.out >>{out}/{xtc}.log 2>&1 
""")
        # Read output
        corr = np.loadtxt(f"{out}/lf_nMImat.out")        
        
        return corr, xtc, pq        
        
        
        
        
        
        
        
        
        
        

        
        
class dcdpkg(Matrixoutput):
    def __init__(self, state, **kwargs):
        if not hasattr(state, "_dcds"):
            self.state._make_dcds()
            
        super().__init__(state, **kwargs)
    
    
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        pdb = self.state._protf("pdb")
        traj = self.state._dcds[xtc]
        
        return pool, pdb, traj, pq # psf? cores! multicorepkg
        
        
        
        
        
class GRINN(dcdpkg, Multicorepkg):
    def __init__(self, state, **kwargs):
        if "namd" in kwargs:
            self.namd = kwargs["namd"]
        else:
            from distutils.spawn import find_executable
            self.namd = find_executable('namd2')
            if self.namd is None:
                raise Exception("namd executable for gRINN computation not found")
        
        super().__init__(state, **kwargs)
        
        if "auto_send" not in kwargs:
            self.state._set_pkgclass(self.state, "gRINNcorr", namd = self.namd, auto_send=True)
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        out = f"{self._path}/{xtc}"
        psf = self.state._protf("psf")
        params = self.state._paramf
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, out, xtc, pq, psf, params, self.taskcpus),
                         callback=self._save_pq)
        
        for _ in range(self.taskcpus-1): pool.apply_async(self._calculate_empty, args=(pq,))
        
        
    def _computation(self, pdb, traj, out, xtc, pq, psf, params, cores):
        outf = f"{out}/energies_intEnMeanTotal.dat"
        
        if not os.path.isfile(outf):
            if os.path.isdir(out):
                from shutil import rmtree
                rmtree(out)
            
            _grinn_calc.getResIntEn(_grinn_args.arg_parser(f"-calc --pdb {pdb} --top {psf} --traj {traj} --exe {self.namd} --outfolder {out} --numcores {cores} --parameterfile {params}".split()))
            
        corr = np.loadtxt(outf)
        return corr, xtc, pq
    

    
    
    
    
class GRINNcorr(GRINN):        
    def __init__(self, state, **kwargs):
        if "auto_send" in kwargs:
            self._auto_send = kwargs["auto_send"]
        else:
            self._auto_send = False
        
        super().__init__(state, **kwargs)
        
    
    def _initialize(self):
        if not rhasattr(self, "state", "GRINN") and not self._auto_send:
            raise Exception("Make sure to send GRINN calculation before GRINNcorr")
        
        no_exist = lambda files: [not os.path.isfile(file) for file in files]
        
        pqs = [self._rawpq(xtc) for xtc in self.state._trajs]
        
        outf = lambda xtc: f"{self._path.replace('GRINNcorr', 'GRINN')}/{xtc}/energies_intEnTotal.csv"
        outfs = [outf(xtc) for xtc in self.state._trajs]
        
        
        def wait_calculate(pqs):
            while any(no_exist(outfs)):
                time.sleep(30)
                
            self._initialize_real()
            
            while any(no_exist(pqs)):
                time.sleep(5)
            return pqs
                
            
        def get_raw(pqs):
            print(f"adding raw data of {self._name} for {self.state._pdbf}: ", pqs)
            flist = map(lambda pq: pandas.read_parquet(pq), pqs)
            df = pandas.concat(flist, axis=1)
            cols = [f"{num}" for num in self.state._trajs]
            df["weight_avg"] = df[cols].fillna(0).mean(axis=1)
            df["weight_std"] = df[cols].fillna(0).std(axis=1)
            return df
        
        add_raw = lambda pqs: setattr(self, "raw", get_raw(pqs))
        get_pool().apply_async(wait_calculate,
                               args=(pqs,),
                               callback=add_raw)
            
        
    def _initialize_real(self):
        pqs = [self._rawpq(xtc) for xtc in self.state._trajs]
        no_exist = lambda pqs: [not os.path.isfile(pq) for pq in pqs]
        
        if any(no_exist(pqs)):
            for xtc in (xtc for xtc in self.state._trajs if no_exist(pqs)[xtc-1]):
                self._calculate(xtc)
                
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super(GRINN, self)._calculate(xtc)
        out = f"{self._path}/{xtc}"
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, out, xtc, pq, self.taskcpus),
                         callback=self._save_pq)
        
        for _ in range(self.taskcpus-1): pool.apply_async(self._calculate_empty, args=(pq,))
        
    
        
    def _computation(self, pdb, out, xtc, pq, cores):
        logFile = f"{out}/grinncorr.log"
        os.system(f"mkdir -p {out}; touch {logFile}")
        outf = f"{out}/energies_resCorr.dat"
        
        if not os.path.isfile(outf):
            _grinn_corr.getResIntCorr(_grinn_args.arg_parser(f"-corr --pdb {pdb} --corrinfile {out.replace('GRINNcorr', 'GRINN')}/energies_intEnTotal.csv --corrprefix {out}/energies --numcores {cores}".split()), logFile=logFile)
            
        corr = np.loadtxt(outf) 
        return corr, xtc, pq
    


# class Carma(): # needs the dcd
# class Bio3D(): # R package; needs the dcd
# class gRINN(): # needs the dcd; could be used with the dcd + pdb and psf + NAMD binaries (in the docs it's recommended to remove non-protein)
# class wordom(): # doesn't have python bindings for croscorr and lmi but if it had it would've been great because it looks fast?
