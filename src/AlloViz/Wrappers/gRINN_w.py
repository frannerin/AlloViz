"""gRINN wrapper

It calculates interaction energies and the correlation between them.

"""

import os, time

import pandas
import numpy as np

from .Base import lazy_import, Multicore

from ..AlloViz.utils import get_pool

imports = {
"_grinn_args": ".Packages.gRINN_Bitbucket.source.grinn",
"_grinn_calc": ".Packages.gRINN_Bitbucket.source.calc",
"_grinn_corr": ".Packages.gRINN_Bitbucket.source.corr",
"_mdtraj": "mdtraj"
}

for key, val in imports.items():
    exec(f"{key} = lazy_import(*{key, val})")
    

    


class gRINN(Multicore):
    """gRINN's interaction energies
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        
        assert new._d["_psff"], "Provide a valid psf file on Protein object initialization to use gRINN"
        assert new._d["_paramf"], "Provide a valid force-field parameters file on Protein object initialization to use gRINN"
        
        new._psf = new._d["_psff"]
        new._params = new._d["_paramf"]
        
        
        if "namd" in new._d:
            new._namd = new._d["namd"]
        else:
            from distutils.spawn import find_executable
            new._namd = find_executable('namd2')
            if new._namd is None:
                raise Exception("namd2 executable for gRINN computation not found")
        
        return new
    
    
    def __init__(self, *args):
        super().__init__(*args)
        
        pqs = [self._rawpq(xtc) for xtc in self._trajs]
        no_exist = lambda pqs: [not os.path.isfile(pq) for pq in pqs]

        def wait_calculate(pqs):
            while any(no_exist(pqs)):
                time.sleep(5)
            return
        
        corrclass = eval("gRINN_corr")
        add_corr = lambda _: setattr(self.protein, "gRINN_corr", corrclass(self.protein, self._d))
        get_pool().apply_async(wait_calculate,
                               args=(pqs,),
                               callback=add_corr)
        
        
    def _computation(self, xtc):
        out = f"{self._path}/{xtc}"
        outf = f"{out}/energies_intEnMeanTotal.dat"
        dcd = f"{self._path}/{xtc}.dcd"
        traj = _mdtraj.load(self._trajs[xtc], top=self._pdbf)#, atom_indices=Protein.protein.indices)
        traj.save_dcd(dcd)
        
        if not os.path.isfile(outf):
            if os.path.isdir(out):
                from shutil import rmtree
                rmtree(out)
            
            _grinn_calc.getResIntEn(_grinn_args.arg_parser(f"-calc --pdb {self._pdbf} --top {self._psf} --traj {dcd} --exe {self._namd} --outfolder {out} --numcores {self.taskcpus} --parameterfile {self._params}".split()))
            
        corr = np.loadtxt(outf)
        return corr, xtc
    

    
    
    
    
class gRINN_corr(gRINN):
    """gRINN's interaction energies' correlation
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)        
        return new
    
    def _computation(self, xtc):
        out = f"{self._path}/{xtc}"
        logFile = f"{out}/grinncorr.log"
        os.makedirs(out, exist_ok=True)
        open(logFile, 'w').close()
        outf = f"{out}/energies_resCorr.dat"
        
        if not os.path.isfile(outf):
            _grinn_corr.getResIntCorr(_grinn_args.arg_parser(f"-corr --pdb {self._pdbf} --corrinfile {out.replace('GRINN_corr', 'GRINN')}/energies_intEnTotal.csv --corrprefix {out}/energies --numcores {self.taskcpus}".split()), logFile=logFile)
            
        corr = np.loadtxt(outf) 
        return corr, xtc
