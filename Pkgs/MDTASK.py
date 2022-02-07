from Pkgs import *
import MDTASK.calc_correlation as mdtask




class MDTASK(Correlationpkg):
    def __init__(self, state):
        super().__init__(state)
    
    
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        
        pool.apply_async(self._computation,
                         args=(pdb, traj, xtc, pq),
                         callback=self._save_pq)
    
    def _computation(self, pdb, traj, xtc, pq):
        corr = mdtask.correlate(mdtask.parse_traj(traj = traj, topology = pdb))
        return corr, xtc, pq
