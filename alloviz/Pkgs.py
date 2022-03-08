import sys, os, pandas, time, importlib
import numpy as np
from .utils import *
from contextlib import redirect_stdout, redirect_stderr
from lazyasd import LazyObject






class Pkg: #abc.ABC
    def __init__(self, state):#, ¿pkg?): #metrics="all", filterby="whole", normalize=True, cores=None, ow=False
        self.state = state        
        self._name = self.__class__.__name__
        
        # self._path = lambda filterby: f"{self.state.name}/data/{self._name}/{filterby}" # lambda FILTERBY MIGHT BE NOT NEEDED 
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
                
        def get_raw(pqs): # *args
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
        # pool = get_pool()
        pdb = self.state._pdbf
        traj = self.state._trajs[xtc]
        pq = self._rawpq(xtc)
        print(f"making {pq}")
        return pool, pdb, traj, pq
    
    
    def _send_and_log(self, *args):#func=None, args:tuple=(), callback=print):
        # if func is None:
        #     func = self._computation
        with open(f"{self._path}/{self._name}.log", "a+") as f:
            with redirect_stdout(f), redirect_stderr(f):
                return self._computation(*args)
                # get_pool().apply_async(func = func,
                #                        args = args,
                #                        callback = callback)
        
    
    # abstractmethod
    def _computation(self, xtc):
        return info, xtc, pq
    
    # abstractmethod
    def _save_pq(self, args):
        info, xtc, pq = args
        # df.to_parquet(pq)
    
    
    


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
        # if not len(df[f"{xtc}"].unique()) == 1:
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






class Getcontacts(Multicorepkg):
    # sys.path.append(f"{__file__.rsplit('/', 1)[0]}/Forks/getcontacts") #/getcontacts
    # get_dynamic_contacts = lazy_import.lazy_module("get_dynamic_contacts")
    # get_contact_frequencies = lazy_import.lazy_module("get_contact_frequencies")
    # from getcontacts import get_dynamic_contacts, get_contact_frequencies
    
    # from .Forks.getcontacts import get_dynamic_contacts, get_contact_frequencies
    get_dynamic_contacts = LazyObject(lambda: importlib.import_module(".Forks.getcontacts.get_dynamic_contacts", package="AlloViz"), globals(), "get_dynamic_contacts")
    get_contact_frequencies = LazyObject(lambda: importlib.import_module(".Forks.getcontacts.get_contact_frequencies", package="AlloViz"), globals(), "get_contact_frequencies")
    
    
    def __init__(self, state):        
        super().__init__(state)
        
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        
        path = self._path
        ctcs = f"{path}/{xtc}.tsv"
        freqs = f"{path}/{xtc}_freqs.tsv"
        
        pool.apply_async(self._send_and_log,
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
    # try: 
    #     from dynetan.proctraj import DNAproc as dynetan
    # except:
    #     print(_cannot_import(__qualname__))
    #sys.path.append(f"{__file__.rsplit('/', 1)[0]}/Forks/dynetan/dynetan")
    #dynetan = lazy_import.lazy_callable(".Forks.dynetan.dynetan.proctraj.DNAproc")
    #from proctraj import DNAproc as dynetan
    #from .Forks.dynetan.dynetan.proctraj import DNAproc as dynetan
    #from .Forks.dynetan.dynetan import proctraj
    
    # from .Forks.dynetan.dynetan.proctraj import DNAproc as dynetan
    dynetanf = LazyObject(lambda: importlib.import_module(".Forks.dynetan.dynetan.proctraj", package="AlloViz").DNAproc, globals(), "dynetanf")
            
    def __init__(self, state):
        super().__init__(state)
    
    
    
    
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq, self.taskcpus),
                         callback=self._save_pq)
        for _ in range(self.taskcpus-1): pool.apply_async(self._calculate_empty, args=(pq,))
    
    
    
    def _computation(self, pdb, traj, xtc, pq, taskcpus):
        obj = self.dynetanf()
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
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq, self.taskcpus),
                         callback=self._save_pq)
        #for _ in range(self.taskcpus-1): pool.apply_async(self._calculate_empty, args=(pq,))










class Corrplus(Matrixoutput):
    # try:
    #     import correlationplus.calculate as corrplus
    # except:
    #     print(_cannot_import(__qualname__))
    #sys.path.append(f"{__file__.rsplit('/', 1)[0]}/Forks/correlationplus/correlationplus")
    #corrplus = lazy_import.lazy_module(".Forks.correlationplus.correlationplus.calculate")
    #print(dir(), dir(corrpluscalc))
    # import correlationplus.calculate as corrplus
    
    # from .Forks.correlationplus.correlationplus import calculate as corrplus
    corrplusf = LazyObject(lambda: importlib.import_module(".Forks.correlationplus.correlationplus.calculate", package="AlloViz"), globals(), "corrplusf")
        
    def __init__(self, state):
        super().__init__(state)
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq),
                         callback=self._save_pq)
    
    
    def _computation(self, pdb, traj, xtc, pq):
        corr = self.corrplusf.calcMDnDCC(pdb, traj, saveMatrix = False)
        return corr, xtc, pq

        
        
        
class CorrplusLMI(Corrplus):
    def __init__(self, state):
        super().__init__(state)
    
    def _computation(self, pdb, traj, xtc, pq):
        corr = self.corrplusf.calcMD_LMI(pdb, traj, saveMatrix = False)
        return corr, xtc, pq
        
        
        
class CorrplusCOM(Corrplus, COMpkg):
    def __init__(self, state):
        super().__init__(state)
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = COMpkg._calculate(self, xtc)
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq),
                         callback=self._save_pq)
        

        
        
class CorrplusCOMLMI(CorrplusLMI, COMpkg):
    def __init__(self, state):
        super().__init__(state)
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = COMpkg._calculate(self, xtc)
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq),
                         callback=self._save_pq)

        
        
        
class CorrplusPsi(Corrplus):        
    def __init__(self, state):
        super().__init__(state)
    
    def _computation(self, pdb, traj, xtc, pq):
        dih = "psi"
        corr = self.corrplusf.calcMDsingleDihedralCC(pdb, traj, dihedralType = dih, saveMatrix = False)
        return corr, xtc, pq

    
class CorrplusPhi(Corrplus):        
    def __init__(self, state):
        super().__init__(state)
    
    def _computation(self, pdb, traj, xtc, pq):
        dih = "phi"
        corr = self.corrplusf.calcMDsingleDihedralCC(pdb, traj, dihedralType = dih, saveMatrix = False)
        return corr, xtc, pq


class CorrplusOmega(Corrplus):        
    def __init__(self, state):
        super().__init__(state)
    
    def _computation(self, pdb, traj, xtc, pq):
        dih = "omega"
        corr = self.corrplusf.calcMDsingleDihedralCC(pdb, traj, dihedralType = dih, saveMatrix = False)
        return corr, xtc, pq



    
    

class MDTASK(Matrixoutput):
    # try:
    #     import MDTASK.calc_correlation as mdtask
    # except:
    #     print(_cannot_import(__qualname__))
    # sys.path.append(f"{__file__.rsplit('/', 1)[0]}/Forks/MD-TASK/mdtask")
    # mdtask = lazy_import.lazy_module("calc_correlation") # REVERT mdtask FOLDER CREATION!
    # import mdtask.calc_correlation as mdtask
    # from .Forks.correlationplus.correlationplus import calculate as corrplus ADAPT TO MD-TASK WHEN FOLDER CREATION IS REVERSED
    # from .Forks.correlationplus.correlationplus import calculate as corrplus
    
    mdtaskf = LazyObject(lambda: importlib.import_module(".Forks.MD-TASK.mdtask.calc_correlation", package="AlloViz"), globals(), "mdtaskf")
    
    
    def __init__(self, state):
        super().__init__(state)
    
    
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq),
                         callback=self._save_pq)
    
    def _computation(self, pdb, traj, xtc, pq):
        corr = self.mdtaskf.correlate(self.mdtaskf.parse_traj(traj = traj, topology = pdb))
        return corr, xtc, pq









class PytrajCA(Matrixoutput):
    # try:
    #     import pytraj
    # except:
    #     print(_cannot_import(__qualname__))
    # corrplus = lazy_import.lazy_module("pytraj")
    
    #import pytraj
    pytrajf = LazyObject(lambda: importlib.import_module("pytraj"), globals(), "pytrajf")
            
    def __init__(self, state):
        super().__init__(state)
    
    
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        mask = self.__class__.__name__[-2:]
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, mask, xtc, pq),
                         callback=self._save_pq)
    
    
    
    def _computation(self, pdb, traj, mask, xtc, pq):
        top = self.pytrajf.load_topology(pdb)
        traj = self.pytrajf.load(traj, top, mask = f'@{mask}')
        corr = self.pytrajf.matrix.correl(traj, f'@{mask}')
        return corr, xtc, pq

    
    
class PytrajCB(PytrajCA):
    def __init__(self, state):
        super().__init__(state)
    
    
    def _calculate(self, xtc):
        pool, pdb, traj, pq = Matrixoutput._calculate(self, xtc)
        mask = self.__class__.__name__[-2:]
        self._selection = "protein and not resname GLY"
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, mask, xtc, pq),
                         callback=self._save_pq)


        
        
        
        
        
        
class Pyinteraph(Matrixoutput):
    # try:
    #     #print(dir(), __module__, __qualname__)
    #     # from pyinteraph.main import main as pyinteraph
    #     import pyinteraph.main as pyinteraph
    # except:
    #     print(_cannot_import(__qualname__))
    # sys.path.append(f"{__file__.rsplit('/', 1)[0]}/Forks/pyinteraph2/pyinteraph")
    # pyinteraph = lazy_import.lazy_module("main")
    #import pyinteraph.main as pyinteraph
    #from .Forks.pyinteraph2.pyinteraph import main as pyinteraph
    pyinteraphf = LazyObject(lambda: importlib.import_module("pyinteraph.main"), globals(), "pyinteraphf")
    
                
    def __init__(self, state):        
        super().__init__(state)
        
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        
        CLIargs = "-m --cmpsn-graph dummy"
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq, CLIargs),
                         callback=self._save_pq)
        
        
    def _computation(self, pdb, traj, xtc, pq, CLIargs):
        corr = self.pyinteraphf.main(f"-s {pdb} -t {traj} {CLIargs}".split()) / 100
        return corr, xtc, pq

    
    
class PyinteraphEne(Pyinteraph):
    def __init__(self, state):        
        super().__init__(state)
        
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = Matrixoutput._calculate(self, xtc)
        
        CLIargs = "-p --kbp-graph dummy"
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq, CLIargs),
                         callback=self._save_pq)

        
        


# only for local
class G_correlationCA(Matrixoutput):
    # os.environ["GMXDIR"] = "/home/frann/gromacs-3.3.1/"
    # os.environ["GMXLIB"] = "/home/frann/gromacs-3.3.1/gromacs/share/gromacs/top/"
                
    def __init__(self, state):        
        super().__init__(state)
        
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq),
                         callback=self._save_pq)
        
        
    def _computation(self, pdb, traj, xtc, pq):
        # Send g_correlation
        if not os.path.isfile(f"{pq}.dat"):
            # os.system("")
            os.system(f"""
module load g_correlation
export GMXLIB=/soft/EB_repo/bio/sequence/programs/noarch/gromacs/3.3.1/share/gromacs/top/
g_correlation -f {traj} -s {pdb} -o {pq}.dat <<EOF
1
3
EOF
""")
            # proc = pexpect.spawnu(f"/home/frann/g_correlation -f {traj} -s {pdb} -o {pq}.dat")
            # proc.logfile = sys.stdout
            # proc.expect('Select a group: ')
            # # proc.send('1')
            # proc.sendline('1')
            # proc.expect('Select a group: ')
            # # proc.send('3')
            # proc.sendline('3')
            # proc.sendline()
            # proc.sendline()
            # proc.wait()
        
        
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

    
    
    
class G_correlationCOM(G_correlationCA, COMpkg):
    def __init__(self, state):        
        super().__init__(state)
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = COMpkg._calculate(self, xtc)
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq),
                         callback=self._save_pq)
        
        
        
        

        
        
class dcdpkg(Matrixoutput):
    def __init__(self, state):
        if not hasattr(state, "_dcds"):
            self.state._make_dcds()
            
        super().__init__(state)
    
    
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        pdb = self.state._protf("pdb")
        traj = self.state._dcds[xtc]
        
        return pool, pdb, traj, pq # psf? cores! multicorepkg
        
        
        
        
        
class GRINN(dcdpkg, Multicorepkg):
    #from .Forks.gRINN_Bitbucket.source import grinn, calc
    grinnf = LazyObject(lambda: importlib.import_module(".Forks.gRINN_Bitbucket.source.grinn", package="AlloViz"), globals(), "grinnf")
    calcf = LazyObject(lambda: importlib.import_module(".Forks.gRINN_Bitbucket.source.calc", package="AlloViz"), globals(), "calcf")

    from .Forks.gRINN_Bitbucket import source
    from distutils.spawn import find_executable
    namd = find_executable('namd2')
    if namd is None:
        namd = f"{source.__file__.rsplit('/', 1)[0]}/NAMD_2.14_Linux-x86_64-multicore/namd2"
                
    def __init__(self, state):
        super().__init__(state)
        
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        psf = self.state._protf("psf")
        params = self.state._paramf
        
        pool.apply_async(self._send_and_log,
                         args=(pdb, traj, xtc, pq, psf, params, self.taskcpus), # psf, self.taskcpus
                         callback=self._save_pq)
        
        if self._name == "GRINN":
            for _ in range(self.taskcpus-1): pool.apply_async(self._calculate_empty, args=(pq,))
        
        
    def _computation(self, pdb, traj, xtc, pq, psf, params, cores):
        out = f"{self._path}/{xtc}"
        if os.path.isdir(out):
            from shutil import rmtree
            rmtree(out)
            
        self.calcf.getResIntEn(self.grinnf.arg_parser(f"-calc --pdb {pdb} --top {psf} --traj {traj} --exe {self.namd} --outfolder {out} --numcores {cores} --parameterfile {params}".split())) # calc.getResIntCorr(grinn.arg_parser(f"".split()), logfile=None)
        corr = np.loadtxt(f"{out}/energies_intEnMeanTotal.dat") # energies_resCorr.dat
        return corr, xtc, pq

    
    

class GRINNcorr(GRINN):
    def __init__(self, state):
        super().__init__(state)
        
    def _computation(self, pdb, traj, xtc, pq, psf, params, cores):
        out = f"{self._path}/{xtc}"
        # if os.path.isdir(out):
        #     import shutil.rmtree
        #     shutil.rmtree(out)
            
        self.calcf.getResIntCorr(self.grinnf.arg_parser(f"-corr --corrinfile {out}/energies_intEnTotal.csv".split()), logfile=None)
        corr = np.loadtxt(f"{out}/energies_resCorr.dat.dat") # 
        return corr, xtc, pq


# class Carma(): # needs the dcd
# class Bio3D(): # R package; needs the dcd
# class gRINN(): # needs the dcd; could be used with the dcd + pdb and psf + NAMD binaries (in the docs it's recommended to remove non-protein)
# class wordom(): # doesn't have python bindings for croscorr and lmi but if it had it would've been great because it looks fast?