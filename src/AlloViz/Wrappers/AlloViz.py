import numpy as np

from .Base import lazy_import, Multicore, Combined_Dihs

imports = {
"_npeet_lnc": ".Packages.NPEET_LNC.lnc",
"_mda_dihedrals": "MDAnalysis.analysis.dihedrals",
}

for key, val in imports.items():
    exec(f"{key} = lazy_import(*{key, val})")
    
    
    
def _calculate_row(in_data, values, nres):
    res, others = in_data
    row = np.zeros((nres,))

    for other in others:
        row[other] = _npeet_lnc.MI.mi_LNC([values[res], values[other]])
        
    return res, row




    
class AlloViz(Multicore):
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._empties = 0
        return new
            
    
    def _computation(self, xtc):#pdb, traj, xtc, pq):
        prot = self._d["u"].atoms
        selected_res = self._d["_dihedral_residx"]()
        
        select_dih = lambda res: eval(f"res.{self._dih.lower()}_selection()")
        selected = [select_dih(res) for res in prot.residues[selected_res]]
        
        offset = 0
        
        if hasattr(self._d["u"].trajectory, "readers"):
            get_frames = lambda numreader: self._d["u"].trajectory.readers[numreader].n_frames
            for numreader in range(xtc-1):
                offset += get_frames(numreader)
        else:
            get_frames = lambda _: self._d["u"].trajectory.n_frames

        values = _mda_dihedrals.Dihedral(selected).run(start=offset, stop=offset+get_frames(xtc-1)).results.angles.transpose()
        
        corr = np.zeros((len(selected_res), len(selected_res)))
        print("prot.n_residues:", prot.n_residues, "len(selected_res)", len(selected_res), "dihedrals array shape:", values.shape, "empty corr array shape:", corr.shape)#
        
        in_data = []
        iterator = set(range(len(selected_res)))
        for res1 in iterator:
            others = list(iterator - set(range(res1)) - {res1})
            if len(others) > 0:
                in_data.append((res1, others))
                
        from concurrent.futures import ProcessPoolExecutor as Pool
        from functools import partial
        with Pool(self.taskcpus) as p:
            results = list(p.map(partial(_calculate_row, values=values, nres=len(selected_res)), in_data))#, chunksize=len(in_data)/(self.taskcpus*2)))
            p.shutdown()
            
        for result in results:
            res, row = result
            corr[res, :] = row
            
            
        return corr, xtc, selected_res
    
    
    
    
    
class AlloViz_Psi(AlloViz):
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dih = "psi"
        return new

    
    
    
class AlloViz_Phi(AlloViz):
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dih = "phi"
        return new

        

class AlloViz_Omega(AlloViz):
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dih = "omega"
        return new

        
        
        
    
class AlloViz_Dihs(Combined_Dihs, AlloViz):
    _Phi = AlloViz_Phi
    _Psi = AlloViz_Psi
    _Omega = AlloViz_Omega
    pass