import AlloViz
import pytest
import MDAnalysis as mda
import numpy as np



@pytest.fixture(scope="module")
def test_path(request):
    return request.fspath.dirname

@pytest.fixture(scope="module")
def noninit_protein(test_path):
    return AlloViz.Protein.__new__(AlloViz.Protein, pdb = f"{test_path}/data/protein.pdb",
                                   trajs = f"{test_path}/data/traj_1.xtc",
                                   path = test_path)




def test_pathname(test_path):
    assert "test" in test_path

@pytest.mark.filterwarnings("ignore::UserWarning:MDAnalysis")
def test_standardize_resnames():
    testu = mda.Universe.empty(4, n_residues=4, atom_resindex=range(4))
    testu.add_TopologyAttr('resnames')
    testu.residues.resnames = ["HSD", "HSE", "CYSP", "PROX"]
    AlloViz.AlloViz.trajutils.standardize_resnames(testu)
    
    assert (np.array(['HIS', 'HIS', 'CYS', 'PRO'], dtype=object) == testu.residues.resnames).all(), "trajutils.standardize_resnames didn't work as expected"