"""pytraj wrapper

It calculates the Pearson's correlation of the residues' CA and CB atoms.

"""

from .Base import lazy_import, Base

imports = {
"_pytraj": "pytraj",
}

for key, val in imports.items():
    exec(f"{key} = lazy_import(*{key, val})")




class pytraj_CA(Base):
    """pytraj's Pearson's correlation of CA atoms
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._mask = new._name[-2:]
        return new
    
    
    def _computation(self, xtc):
        """"""
        top = _pytraj.load_topology(self._pdbf)
        traj = _pytraj.load(self._trajs[xtc], top, mask = f'@{self._mask}')
        corr = _pytraj.matrix.correl(traj, f'@{self._mask}')
        return corr, xtc

    
    
class pytraj_CB(pytraj_CA):
    """pytraj's Pearson's correlation of CB atoms
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._selection = "not resname GLY"
        return new
