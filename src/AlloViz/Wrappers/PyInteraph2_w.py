"""PyInteraph2 wrapper

It calculates interaction energies and contact frequencies.

"""

import numpy as np

from .Base import lazy_import, Base


imports = {
"_pyinteraph": "pyinteraph.main",
}

for key, val in imports.items():
    exec(f"{key} = lazy_import(*{key, val})")
    
    
cmpsn_reslist = \
        ["ALA", "CYS", "ASP", "GLU", "PHE", "HIS", "ILE", "LYS", "LEU",
         "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP",
         "TYR"] + ["GLY"]


    

class PyInteraph2_Base(Base):
    """PyInteraph2 base class
    """    
    def _computation(self, xtc):
        """"""
        _pyinteraph.main(f"-s {self._pdbf} -t {self._trajs[xtc]} {self._CLIargs(xtc)}".split(" "))
        corr = np.loadtxt(f"{self._rawpq(xtc)}{'_all' if 'Energy' not in self._name else ''}.dat")
        return corr, xtc
    
    
    
class PyInteraph2_Contacts(PyInteraph2_Base):
    """PyInteraph2's side-chain contact frequencies
    
    Contact frequencies based on the fulfillment of a distance threshold (5 angstroms)
    requirement of each residue pair's COMs' distance
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._CLIargs = lambda xtc: f"-m --cmpsn-csv {new._rawpq(xtc)}.csv --cmpsn-graph {new._rawpq(xtc)}.dat --cmpsn-residues {','.join(cmpsn_reslist)}"
        return new
                
        
    def _computation(self, xtc):
        """"""
        corr, xtc = super()._computation(xtc)
        corr = corr / 100
        
        bonded_cys_indices = self._d["_bonded_cys"]
        
        if len(bonded_cys_indices) > 1:
            for res1 in bonded_cys_indices:
                for res2 in bonded_cys_indices:
                    if res1 != res2:
                        corr[res1, res2] = 0
        
        return corr, xtc
    
    
    
# class PyInteraph2_Contacts_Corrected(PyInteraph2_Contacts):
#     """PyInteraph2's corrected side-chain contact frequencies
    
#     Contact frequencies based on the fulfillment of a distance threshold requirement
#     of each residue pair's COMs' distance. A correction based on the radius of gyration
#     of the involved residues is applied to the COMs' distance and thus the distance
#     threshold is lowered to 2.5 angstroms, based on the authors' calculations.
#     """
#     def __new__(cls, protein, d):
#         new = super().__new__(cls, protein, d)
#         new._CLIargs = lambda xtc: f"-m --cmpsn-csv {new._rawpq(xtc)}.csv --cmpsn-graph {new._rawpq(xtc)}.dat --cmpsn-residues {','.join(cmpsn_reslist)} --cmpsn-correction rg --cmpsn-co 2.5"
#         return new
    
    
    
    
class PyInteraph2_Atomic_Contacts_Strength(PyInteraph2_Base):
    """PyInteraph2's atomic contacts strength
    
    Atom pair distance-based residue pairs' contacts strength, based on the work by
    Brinda and Vishveshwara, 2005 [1]_.
    
    .. [1] Brinda,K. V and Vishveshwara,S. (2005) A network representation of protein
       structures: implications for protein stability. Biophys. J., 89, 4159–70. 
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._CLIargs = lambda xtc: f"-a --acpsn-csv {new._rawpq(xtc)}.csv --acpsn-graph {new._rawpq(xtc)}.dat --acpsn-proxco 0 --acpsn-imin 0 --acpsn-nf-permissive --acpsn-ew strength"
        return new
    
    
class PyInteraph2_Atomic_Contacts_Occurrence(PyInteraph2_Contacts):
    """PyInteraph2's atomic contacts occurrence
    
     Atom pair distance-based residue pairs' contacts occurrence/"persistence"/frequency.
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._CLIargs = lambda xtc: f"-a --acpsn-csv {new._rawpq(xtc)}.csv --acpsn-graph {new._rawpq(xtc)}.dat --acpsn-proxco 0 --acpsn-imin 0 --acpsn-nf-permissive --acpsn-ew persistence"
        return new

    
    
class PyInteraph2_Energy(PyInteraph2_Base):
    """PyInteraph2's interaction energies
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._CLIargs = lambda xtc: f"-p --kbp-csv {new._rawpq(xtc)}.csv --kbp-graph {new._rawpq(xtc)}.dat"
        return new