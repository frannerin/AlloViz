import os

import pandas

import numpy as np

from .Base import Multicore, Use_dcd

from ..AlloViz.utils import lazy_import

imports = {
"_grinn_args": ".Packages.gRINN_Bitbucket.source.grinn",
"_grinn_calc": ".Packages.gRINN_Bitbucket.source.calc",
"_grinn_corr": ".Packages.gRINN_Bitbucket.source.corr",
}

for key, val in imports.items():
    exec(f"{key} = lazy_import(*{key, val})")
    
    
    


class gRINN(Use_dcd, Multicore):
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        
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
            new.protein._set_pkgclass(self.protein, "gRINN_corr", d)
        
        return new
        
        
    def _computation(self, xtc):# pdb, traj, out, xtc, pq, psf, params, taskcpus):
        psf = self._d["_protpsf"]
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
    

    
    
    
    
class gRINN_corr(gRINN):
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        
        if "auto_send" in new._d:
            new._auto_send = new._d["auto_send"]
        else:
            new._auto_send = False
        
        return new
    
    
    def __init__(self):
        if not rhasattr(self, "protein", "GRINN") and not self._auto_send:
            raise Exception("Make sure to send GRINN calculation before GRINN_corr")
        
        no_exist = lambda files: [not os.path.isfile(file) for file in files]
        
        pqs = [self._rawpq(xtc) for xtc in self._trajs]
        
        outf = lambda xtc: f"{self._path.replace('GRINN_corr', 'GRINN')}/{xtc}/energies_intEnTotal.csv"
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
            _grinn_corr.getResIntCorr(_grinn_args.arg_parser(f"-corr --pdb {self._pdbf} --corrinfile {out.replace('GRINN_corr', 'GRINN')}/energies_intEnTotal.csv --corrprefix {out}/energies --numcores {self.taskcpus}".split()), logFile=logFile)
            
        corr = np.loadtxt(outf) 
        return corr, xtc