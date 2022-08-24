"""GSAtools wrapper

It calculates Mutual Information (MI).

"""

import os
import numpy as np

from .Base import Base




# only for local
class GSAtools(Base):
    """GSAtools' MI
    """        
    def _computation(self, xtc):# pdb, traj, xtc, pq, out):
        # Send gsatools
        out = f"{self._path}/{xtc}"
        if not os.path.isfile(f"{out}/lf_nMImat.out"):
            os.system(f"""
module load GSAtools
g_sa_encode -s {self._pdbf} -f {self._trajs[xtc]} -rmsdlf {out}/lf_rmsd.xvg -strlf {out}/lf_str.out -log {out}/log.log 
g_sa_analyze -sa {out}/lf_str.out -MImatrix -MImat {out}/lf_MImat.out -eeMImat {out}/lf_eeMImat.out -jHmat {out}/lf_jHmat.out -nMImat {out}/lf_nMImat.out >>{out}/{xtc}.log 2>&1 
""")
        # Read output
        corr = np.loadtxt(f"{out}/lf_nMImat.out")        
        
        return corr, xtc      