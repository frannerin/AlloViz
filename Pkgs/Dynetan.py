from Pkgs import *
from dynetan.proctraj import DNAproc as dynetan




class Dynetan(Correlationpkg, Multicorepkg):
    def __init__(self, state):
        super().__init__(state)
    
    
    
    
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        
        pool.apply_async(self._computation,
                         args=(pdb, traj, xtc, pq, self.taskcpus),
                         callback=self._save_pq)
        for _ in range(self.taskcpus-1): pool.apply_async(self._calculate_empty, args=(pq,))
    
    
    
    def _computation(self, pdb, traj, xtc, pq, taskcpus):
        obj = dynetan()
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


# ### COMpkgs

# In[32]:





# In[33]:


class DynetanCOM(Dynetan, COMpkg):
    def __init__(self, state):
        super().__init__(state)
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = COMpkg._calculate(self, xtc)
        
        pool.apply_async(self._computation,
                         args=(pdb, traj, xtc, pq, self.taskcpus),
                         callback=self._save_pq)
        for _ in range(self.taskcpus-1): pool.apply_async(self._calculate_empty, args=(pq,))
