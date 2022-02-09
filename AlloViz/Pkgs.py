import sys, os, pandas, time
import numpy as np
from .utils import *



def _cannot_import(pkgname):
        return f"{pkgname} can't be imported"




class Pkg: #abc.ABC
    def __init__(self, state):#, ¿pkg?): #metrics="all", filterby="whole", normalize=True, cores=None, ow=False
        args = locals()
        del args["self"]
        if "args" in args: del args["args"]
        self.__dict__.update(args)
        # self.pool = get_pool()
        
        self._name = self.__class__.__name__
        self._path = lambda filterby: f"{self.state.name}/data/{self._name}/{filterby}"
        self._rawpq = lambda xtc: f"{self._path('raw')}/{xtc if isinstance(xtc, int) else xtc.rsplit('.', 1)[0]}.pq"
        
        self.raw = self._initialize()
        
        
#     def _get_norm_str(self, normalize):
#         return "norm" if normalize else "no_norm"
    
#     def _datapn(self, filterby, normalize, pkg):
#         norm = self._get_norm_str(normalize)
#         return f"{self.name}/data/{filterby}/{norm}/{pkg}"
    
#     def _datafn(self, filterby, normalize, pkg, metric):
#         return f"{self._datapn(filterby, normalize, pkg)}/{metric}.pq"
    
    
        
        
        
    
    def _initialize(self):
        pqs = [self._rawpq(xtc) for xtc in self.state._trajs]
        no_exist = [not os.path.isfile(pq) for pq in pqs]
        
        if any(no_exist):
            os.makedirs(self._path("raw"), exist_ok=True)
            for xtc in (xtc for xtc in self.state._trajs if no_exist[xtc-1]):
                self._calculate(xtc)
        
        
        
        def wait_calculate(pqs):
            while any([not os.path.isfile(pq) for pq in pqs]):
                time.sleep(5)
                
        def get_raw(*args):
            print(f"adding raw data of {self._name} for {self.state.name}")
            flist = map(lambda pq: pandas.read_parquet(pq), pqs)
            df = pandas.concat(flist, axis=1)
            cols = [f"{num}" for num in self.state._trajs]
            df["weight_avg"] = df[cols].fillna(0).mean(axis=1)
            df["weight_std"] = df[cols].fillna(0).std(axis=1)
            return df
            
        add_raw = lambda _: setattr(self, "raw", get_raw())
        get_pool().apply_async(wait_calculate,
                               args=(pqs,),
                               callback=add_raw)
        return
        

        
    
    def _calculate(self, xtc):
        pool = get_pool()
        pdb = self.state._pdbf
        traj = f"{self.state.name}/{self.state._trajs[xtc]}"
        pq = self._rawpq(xtc)
        print(f"making {pq}")
        return pool, pdb, traj, pq
    
    
    
    # abstractmethod
    def _computation(self, xtc):
        return freqs, xtc, pq
    
    # abstractmethod
    def _save_pq(self, args):
        corr, xtc, pq = args
        # df.to_parquet(pq)





class Multicorepkg(Pkg):
    def __init__(self, state, taskcpus = 3):
        self.taskcpus = taskcpus
        super().__init__(state)
    
    def _calculate_empty(self, pqf):
        while not os.path.isfile(pqf):
            time.sleep(5)
        return






class Matrixoutput(Pkg):
    def __init__(self, state):
        self.selection = "protein"
        super().__init__(state)
    
    
    def _save_pq(self, args):
        corr, xtc, pq = args
        
        resl = [f"A:{aa.resname}:{aa.resid}" for aa in self.state.mdau.select_atoms(self.selection).residues]

        df = pandas.DataFrame(corr, columns=resl, index=resl)
        df = df.where( np.triu(np.ones(df.shape), k=1).astype(bool) )
        df = pandas.DataFrame({f"{xtc}": df.stack()})
        df.to_parquet(pq)






class COMpkg(Pkg):
    def __init__(self, state):
        if not hasattr(state, "_comtrajs"):
            self.state._add_comtrajs()
            
        super().__init__(state)
    
    
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        pdb = self.state._compdbf
        traj = self.state._comtrajs[xtc]
        
        return pool, pdb, traj, pq








# sys.path.append('getcontacts')
# from getcontacts import get_dynamic_contacts, get_contact_frequencies


class Getcontacts(Multicorepkg):
    try:
        sys.path.append('getcontacts')
        from getcontacts import get_dynamic_contacts, get_contact_frequencies
    except:
        print(_cannot_import(__qualname__))
        
        
    def __init__(self, state):        
        super().__init__(state)
        
        
        
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        
        path = self._path("raw")
        ctcs = f"{path}/{xtc}.tsv"
        freqs = f"{path}/{xtc}_freqs.tsv"
        
        pool.apply_async(self._computation,
                         args=(pdb, traj, xtc, pq, ctcs, freqs, self.taskcpus),
                         callback=self._save_pq)
        
        for _ in range(self.taskcpus-1): pool.apply_async(self._calculate_empty, args=(pq,))
        
        
    def _computation(self, pdb, traj, xtc, pq, ctcs, freqs, taskcpus):
        if not os.path.isfile(freqs):# or ow:
            self.get_dynamic_contacts.main(f"--topology {pdb} --trajectory {traj} --output {ctcs} --itypes all --cores {taskcpus}".split())
            self.get_contact_frequencies.main(f"--input_files {ctcs} --output_file {freqs}".split())
        return freqs, xtc, pq
        
        
    def _save_pq(self, args):
        freqs, xtc, pq = args
        
        df = pandas.read_csv(freqs, sep="\t", skiprows=2,
                             index_col = (0, 1), names = [f"{xtc}"])
        df.index = df.index.map(lambda idx: tuple(sorted(idx, key = lambda res: int(res.split(":")[-1]))))
        df.to_parquet(pq)










class Dynetan(Matrixoutput, Multicorepkg):
    try: 
        from dynetan.proctraj import DNAproc as dynetan
    except:
        print(_cannot_import(__qualname__))
            
    def __init__(self, state):
        super().__init__(state)
    
    
    
    
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        
        pool.apply_async(self._computation,
                         args=(pdb, traj, xtc, pq, self.taskcpus),
                         callback=self._save_pq)
        for _ in range(self.taskcpus-1): pool.apply_async(self._calculate_empty, args=(pq,))
    
    
    
    def _computation(self, pdb, traj, xtc, pq, taskcpus):
        obj = self.dynetan()
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
    def __init__(self, state):
        super().__init__(state)
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = COMpkg._calculate(self, xtc)
        
        pool.apply_async(self._computation,
                         args=(pdb, traj, xtc, pq, self.taskcpus),
                         callback=self._save_pq)
        for _ in range(self.taskcpus-1): pool.apply_async(self._calculate_empty, args=(pq,))










class Corrplus(Matrixoutput):
    try:
        import correlationplus.calculate as corrplus
    except:
        print(_cannot_import(__qualname__))
        
    def __init__(self, state):
        super().__init__(state)
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        
        pool.apply_async(self._computation,
                         args=(pdb, traj, xtc, pq),
                         callback=self._save_pq)
    
    
    def _computation(self, pdb, traj, xtc, pq):
        corr = self.corrplus.calcMDnDCC(pdb, traj, saveMatrix = False)
        return corr, xtc, pq

        
        
        
class CorrplusLMI(Corrplus):
    def __init__(self, state):
        super().__init__(state)
    
    def _computation(self, pdb, traj, xtc, pq):
        corr = self.corrplus.calcMD_LMI(pdb, traj, saveMatrix = False)
        return corr, xtc, pq
        
        
        
class CorrplusCOM(Corrplus, COMpkg):
    def __init__(self, state):
        super().__init__(state)
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = COMpkg._calculate(self, xtc)
        
        pool.apply_async(self._computation,
                         args=(pdb, traj, xtc, pq),
                         callback=self._save_pq)
        

        
        
class CorrplusCOMLMI(CorrplusLMI, COMpkg):
    def __init__(self, state):
        super().__init__(state)
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = COMpkg._calculate(self, xtc)
        
        pool.apply_async(self._computation,
                         args=(pdb, traj, xtc, pq),
                         callback=self._save_pq)









class MDTASK(Matrixoutput):
    try:
        import MDTASK.calc_correlation as mdtask
    except:
        print(_cannot_import(__qualname__))
    def __init__(self, state):
        super().__init__(state)
    
    
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        
        pool.apply_async(self._computation,
                         args=(pdb, traj, xtc, pq),
                         callback=self._save_pq)
    
    def _computation(self, pdb, traj, xtc, pq):
        corr = self.mdtask.correlate(mdtask.parse_traj(traj = traj, topology = pdb))
        return corr, xtc, pq









class PytrajCA(Matrixoutput):
    try:
        import pytraj
    except:
        print(_cannot_import(__qualname__))
            
    def __init__(self, state):
        super().__init__(state)
    
    
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        mask = self.__class__.__name__[-2:]
        
        pool.apply_async(self._computation,
                         args=(pdb, traj, mask, xtc, pq),
                         callback=self._save_pq)
    
    
    
    def _computation(self, pdb, traj, xtc, pq):
        top = self.pytraj.load_topology(pdb)
        traj = self.pytraj.load(traj, top, mask = f'@{mask}')
        corr = self.pytraj.matrix.correl(traj, f'@{mask}')
        return corr, xtc, pq

    
    
class PytrajCB(PytrajCA):
    def __init__(self, state):
        super().__init__(state)
    
    
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        mask = self.__class__.__name__[-2:]
        self.selection = "protein and not resname GLY"
        
        pool.apply_async(self._computation,
                         args=(pdb, traj, mask, xtc, pq),
                         callback=self._save_pq)


        
        
        
        
        
        
class Pyinteraph(Matrixoutput):
    try:
        print(dir(), __module__, __qualname__)
        from pyinteraph.main import main as pyinteraph
    except:
        print(_cannot_import(__qualname__))
        
        
    def __init__(self, state):        
        super().__init__(state)
        
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        
        CLIargs = "-m --cmpsn-graph dummy"
        
        pool.apply_async(self._computation,
                         args=(pdb, traj, xtc, pq, CLIargs),
                         callback=self._save_pq)
        
        
    def _computation(self, pdb, traj, xtc, pq, CLIargs):
        self.pyinteraph(f"-s {pdb} -t {traj} {CLIargs}".split())
        return freqs, xtc, pq

    
    
class PyinteraphEne(Pyinteraph):
    def __init__(self, state):        
        super().__init__(state)
        
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        
        CLIargs = "-p --kbp-graph dummy"
        
        pool.apply_async(self._computation,
                         args=(pdb, traj, xtc, pq, CLIargs),
                         callback=self._save_pq)
