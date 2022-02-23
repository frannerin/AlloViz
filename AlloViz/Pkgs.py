import sys, os, pandas, time
import numpy as np
from .utils import *
from contextlib import redirect_stdout, redirect_stderr



def _cannot_import(pkgname):
        return f"{pkgname} can't be imported"




class Pkg: #abc.ABC
    def __init__(self, state):#, ¿pkg?): #metrics="all", filterby="whole", normalize=True, cores=None, ow=False
        args = locals()
        del args["self"]
        if "args" in args: del args["args"]
        self.__dict__.update(args)
        
        self._name = self.__class__.__name__
        # self._path = lambda filterby: f"{self.state.name}/data/{self._name}/{filterby}" # lambda FILTERBY MIGHT BE NOT NEEDED 
        self._path = f"{self.state.name}/data/{self._name}/raw"
        os.makedirs(self._path, exist_ok=True)
        self._rawpq = lambda xtc: f"{self._path}/{xtc if isinstance(xtc, int) else xtc.rsplit('.', 1)[0]}.pq"
        
        if self.state.__class__.__name__ == "State":
            with open(f"{self._path}/{self._name}.log", "a+") as f:
                with redirect_stdout(f), redirect_stderr(f):
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
                
        def get_raw(pqs): # *args
            print(f"adding raw data of {self._name} for {self.state.name}: ", pqs)
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
        traj = f"{self.state.name}/{self.state._trajs[xtc]}"
        pq = self._rawpq(xtc)
        print(f"making {pq}")
        return pool, pdb, traj, pq
    
    # abstractmethod
    def _computation(self, xtc):
        return info, xtc, pq
    
    # abstractmethod
    def _save_pq(self, args):
        info, xtc, pq = args
        # df.to_parquet(pq)
    
    
    
    
    
#     def _analyze(self, ):
#         norm = self._state._get_norm_str(normalize)
        
        
#         os.makedirs(self._state._datapn(filterby, normalize, pkg), exist_ok=True)
#         raw = getattr(self.raw, pkg)
#         if filterby == "incontact" or filterby == "intercontact":
#             raw = raw.filter(self.raw.getcontacts.index, axis=0)
#             if filterby == "intercontact":
#                 get_resnum = lambda res: int(res.rsplit(":")[-1])
#                 raw = raw.filter(get_intercontacts(raw.index), axis=0)
#         if filterby == "whole":
#             metrics = [metric for metric in metrics if "subset" not in metric]
            
        
#         data = rgetattr(self, filterby, norm, pkg).df
        
#         if any([f"{metric}_avg" not in data.columns for metric in metrics]) or ow:
#             #print(f"analyzing {metrics} from {self._state.name} {pkg} data for {filterby} residues")
#             raw_data = raw[[col for col in raw.columns if col != "weight_std"]]
#             self._send_anal(raw_data, pkg, metrics, normalize, ow, filterby)
#             # setattr(self.incontact, pkg, Edges(self._state, data.join(analyzed)))
            
#         return
        
        
#     def _send_anal(self, rawdata, pkg, metrics, normalize, ow, filterby):
#         data = None
        
#         if filterby == "whole":
#             metrics = [metric for metric in metrics if "subset" not in metric]
        
#         for metric in metrics:
#             pqf = self._state._datafn(filterby, normalize, pkg, metric)
            
#             if not os.path.isfile(pqf) or ow:
#                 print(f"\tmaking {pqf}")
#                 newcolnames = {name: f"{metric}_{name}" for name in rawdata.columns}
#                 df = rawdata.apply(self._networkx_analysis, args=(metric, normalize)).rename(columns=newcolnames)
#                 df.to_parquet(pqf)
            
#         return   
        





class Multicorepkg(Pkg):
    def __init__(self, state, taskcpus = 3):
        self.taskcpus = taskcpus
        super().__init__(state)
    
    def _calculate_empty(self, pqf):
        print("sleeping", pqf, os.getpid())
        while not os.path.isfile(pqf):
            time.sleep(5)
        return






class Matrixoutput(Pkg):
    def __init__(self, state):
        self._selection = "protein"
        super().__init__(state)
    
    
    def _save_pq(self, args):
        corr, xtc, pq = args
        
        resl = [f"A:{aa.resname}:{aa.resid}" for aa in self.state.mdau.select_atoms(self._selection).residues]

        df = pandas.DataFrame(corr, columns=resl, index=resl)
        df = df.where( np.triu(np.ones(df.shape), k=1).astype(bool) )
        df = pandas.DataFrame({f"{xtc}": df.stack()})
        if not len(df[f"{xtc}"].unique()) == 1:
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
        
        path = self._path
        ctcs = f"{path}/{xtc}.tsv"
        freqs = f"{path}/{xtc}_freqs.tsv"
        
        pool.apply_async(self._computation,
                         args=(pdb, traj, xtc, pq, ctcs, freqs, self.taskcpus),
                         callback=self._save_pq)
        # update_pdict((self.name, self._name, xtc), p)
        
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
        self.taskcpus = 4
        pool.apply_async(self._computation,
                         args=(pdb, traj, xtc, pq, self.taskcpus),
                         callback=self._save_pq)
        #for _ in range(self.taskcpus-1): pool.apply_async(self._calculate_empty, args=(pq,))










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
        corr = self.mdtask.correlate(self.mdtask.parse_traj(traj = traj, topology = pdb))
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
    
    
    
    def _computation(self, pdb, traj, mask, xtc, pq):
        top = self.pytraj.load_topology(pdb)
        traj = self.pytraj.load(traj, top, mask = f'@{mask}')
        corr = self.pytraj.matrix.correl(traj, f'@{mask}')
        return corr, xtc, pq

    
    
class PytrajCB(PytrajCA):
    def __init__(self, state):
        super().__init__(state)
    
    
    def _calculate(self, xtc):
        pool, pdb, traj, pq = Matrixoutput._calculate(self, xtc)
        mask = self.__class__.__name__[-2:]
        self._selection = "protein and not resname GLY"
        
        pool.apply_async(self._computation,
                         args=(pdb, traj, mask, xtc, pq),
                         callback=self._save_pq)


        
        
        
        
        
        
class Pyinteraph(Matrixoutput):
    try:
        #print(dir(), __module__, __qualname__)
        # from pyinteraph.main import main as pyinteraph
        import pyinteraph.main as pyinteraph
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
        corr = self.pyinteraph.main(f"-s {pdb} -t {traj} {CLIargs}".split()) / 100
        return corr, xtc, pq

    
    
class PyinteraphEne(Pyinteraph):
    def __init__(self, state):        
        super().__init__(state)
        
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = Matrixoutput._calculate(self, xtc)
        
        CLIargs = "-p --kbp-graph dummy"
        
        pool.apply_async(self._computation,
                         args=(pdb, traj, xtc, pq, CLIargs),
                         callback=self._save_pq)

        
        
        
        
# class Carma(): # needs the dcd
# class Bio3D(): # R package; needs the dcd
# class wordom(): # doesn't have python bindings for croscorr and lmi but if it had it would've been great because it looks fast?

