import os
import pytest
import AlloViz
import MDAnalysis as mda
import numpy as np

pytestmark = pytest.mark.filterwarnings("ignore:::MDAnalysis")


@pytest.fixture(scope="module")
def test_path(request):
    return request.fspath.dirname

@pytest.fixture(scope="module")
def noninit_protein(test_path):
    prot = AlloViz.Protein.__new__(AlloViz.Protein, pdb = f"{test_path}/data/protein.pdb",
                                   trajs = f"{test_path}/data/traj_1.xtc",
                                   path = test_path)
    yield prot
    del prot




def test_pathname(test_path):
    assert "test" in test_path

    
# @pytest.mark.filterwarnings("ignore::UserWarning:MDAnalysis")
def test_standardize_resnames():
    testu = mda.Universe.empty(4, n_residues=4, atom_resindex=range(4))
    testu.add_TopologyAttr('resnames')
    testu.residues.resnames = ["HSD", "HSE", "CYSP", "PROX"]
    AlloViz.AlloViz.trajutils.standardize_resnames(testu)
    
    assert (np.array(['HIS', 'HIS', 'CYS', 'PRO'], dtype=object) == testu.residues.resnames).all(), "trajutils.standardize_resnames didn't work as expected"
    

# @pytest.mark.filterwarnings("ignore::UserWarning:MDAnalysis")
def test_traj_processing(noninit_protein):
    p = noninit_protein
    pdb, trajs, psff = AlloViz.AlloViz.trajutils.process_input(p.GPCR, p.pdb, p._protein_sel, p._pdbf, p._compdbf, False, None, p.trajs, p._trajs, p._comtrajs)
    assert p.pdb == pdb, "PDB file wasn't correctly recognised as already processed"
    assert p.trajs == list(trajs.values()), "Trajectory file wasn't correctly recognised as already processed"
    assert os.path.isfile(p._compdbf), "COM topology file wasn't correctly generated"
    