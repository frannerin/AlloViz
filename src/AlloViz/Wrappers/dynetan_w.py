"""dynetan wrapper

It calculates the Mutual Information of the residues' movement and of their COMs.

"""

import numpy as np

from .Base import lazy_import, Multicore, Use_COM

imports = {
"_dynetan": ".Packages.dynetan.dynetan.proctraj",
}

for key, val in imports.items():
    exec(f"{key} = lazy_import(*{key, val})")
    
    


class dynetan(Multicore):
    """dyentan's MI of the residues
    """    
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        if "stride" in d:
            new.stride = d["stride"]
        return new
    
    def _computation(self, xtc):
        """"""
        step = {"in_memory": True, "in_memory_step": self.stride} if hasattr(self, "stride") else {}
        
        obj = _dynetan.DNAproc(notebookMode=False)
        obj.workU = _dynetan.mda.Universe(self._pdbf, self._trajs[xtc], **step)
        
        # prot = obj.getU().select_atoms(self._d["_protein_sel"])
        # from Bio.SeqUtils import seq1, seq3
        # for res in prot.residues:
        #     res.resname = seq3(seq1(res.resname, custom_map = self._d["_res_dict"])).upper()
        
        protseg = list(obj.getU().segments.segids)
        obj.setSegIDs(protseg)
        obj.selectSystem(withSolvent=False)#, userSelStr=self._d["_protein_sel"])

        obj.setCustomResNodes({})
        obj.setUsrNodeGroups({})

        obj.setNumWinds(1)
        obj.alignTraj(inMemory=False)
        obj.prepareNetwork()
        obj.contactMatAll = np.triu(np.ones([obj.numWinds, obj.numNodes, obj.numNodes], dtype=int), k=1)

        obj.calcCor(ncores=self.taskcpus)
        
        return obj.corrMatAll[0], xtc



class dynetan_COM(Use_COM, dynetan):
    """dyentan's MI of the residues' COM
    """    
    pass