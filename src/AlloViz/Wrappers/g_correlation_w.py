"""g_correlation wrapper

It calculates the Mutual Information (MI) and Linear MI (LMI) of the residues' CA atoms. Only for local installations.

"""

import os
import numpy as np

from .Base import Base
    
    
    
    
# only for local
class g_correlation_CA_MI(Base):
    """g_correlation's MI of CA atoms
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._CLIargs = ""
        return new
                
        
    def _computation(self, xtc):
        # Send g_correlation
        pq = self._rawpq(xtc)
        if not os.path.isfile(f"{pq}.dat"):
            os.system(f"""
module load g_correlation
export GMXLIB=/soft/EB_repo/bio/sequence/programs/noarch/gromacs/3.3.1/share/gromacs/top/
g_correlation -f {self._trajs[xtc]} -s {self._pdbf} -o {pq}.dat {self._CLIargs} &> {self._path}/{xtc}.log <<EOF
1
3
EOF
""")
        # Read output.dat
        size = self._d["protein"].select_atoms("name CA").n_atoms
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
        
        
        return corr, xtc

    
    
    
# class g_correlation_COM_MI(Use_COM, g_correlation_CA_MI):
#     """g_correlation's MI of residues' COM
#     """
#     def __new__(cls, protein, d):
#         new = super().__new__(cls, protein, d)
#         new._CLIargs = ""
#         return new
        
        
        
        
class g_correlation_CA_LMI(g_correlation_CA_MI):
    """g_correlation's LMI of CA atoms
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._CLIargs = "-linear"
        return new