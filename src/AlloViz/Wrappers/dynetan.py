import numpy as np

from .Base import Multicore, Use_COM

from ..AlloViz.utils import lazy_import

imports = {
"_dynetan": ".Packages.dynetan.dynetan.proctraj",
}

for key, val in imports.items():
    exec(f"{key} = lazy_import(*{key, val})")
    
    


class dynetan(Multicore):    
    
    def _computation(self, xtc):
        obj = _dynetan.DNAproc()
        obj.loadSystem(self._pdbf, self._trajs[xtc])
        
        prot = obj.getU().select_atoms(self._d["_protein_sel"])
        from Bio.SeqUtils import seq1, seq3
        for res in prot.residues:
            res.resname = seq3(seq1(res.resname, custom_map = self._d["_res_dict"])).upper()
        
        protseg = list(prot.segments.segids)
        obj.setSegIDs(protseg)
        obj.selectSystem(withSolvent=False, userSelStr=self._d["_protein_sel"])

        obj.setCustomResNodes({})
        obj.setUsrNodeGroups({})

        obj.setNumWinds(1)
        obj.alignTraj(inMemory=False)
        obj.prepareNetwork()
        obj.contactMatAll = np.triu(np.ones([obj.numWinds, obj.numNodes, obj.numNodes], dtype=int), k=1)

        obj.calcCor(ncores=self.taskcpus)
        
        return obj.corrMatAll[0], xtc



class dynetan_COM(Use_COM, Dynetan):
    pass