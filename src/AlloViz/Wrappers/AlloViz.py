import numpy as np

from .Base import Multicore, Use_COM, Combined_Dihs

from ..AlloViz.utils import lazy_import

imports = {
"_npeet_lnc": "..Packages.NPEET_LNC.lnc",
"_mda_dihedrals": "MDAnalysis.analysis.dihedrals",
}

for key, val in imports.items():
    exec(f"{key} = lazy_import(*{key, val})")
    
    
    


    
class AlloViz(Multicore):
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        return new
            
    
    def _computation(self, xtc):#pdb, traj, xtc, pq):
        prot = self._d["u"].select_atoms(self._d["_protein_sel"])
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
        
        

        from multiprocess import Queue, Process
    
        in_data = Queue()
        out_data = Queue()
        
        
        iterator = set(range(len(selected_res)))
        for res1 in iterator:
            others = []
            for res2 in iterator - set(range(res1)) - {res1}:
                others.append(res2)
                
            if len(others) > 0: 
                in_data.put((res1, others))
        
        
        def calculate_row(values, nres, in_data, out_data):
             while not in_data.empty():
                res, others = in_data.get()
                row = np.zeros((nres,))
                
                for other in others:
                    row[other] = _npeet_lnc.MI.mi_LNC([values[res], values[other]])
                
                out_data.put((res, row))
                
        
        ps = []
        for _ in range(self.taskcpus):
            p = Process(target=calculate_row,
                        args=(values, len(selected_res), in_data, out_data)) 
            p.start()
            ps.append(p)
        
        for _ in range(len(selected_res)-1):
            res, row = out_data.get()
            corr[res, :] = row
        
        for p in ps:
            p.join()
            
            
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
    pass