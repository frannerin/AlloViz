import os, shutil
import pytest
import AlloViz
import numpy as np

pytestmark = pytest.mark.filterwarnings("ignore:::MDAnalysis")

@pytest.fixture(scope="module")
def path(tmp_path_factory):
    path = tmp_path_factory.mktemp("downloads")
    yield str(path)
    shutil.rmtree(str(path))

    
    
def test_download(path):
    AlloViz.AlloViz.trajutils.download_GPCRmd_files(10, path)
    files = os.listdir(path)
    assert ("10182_dyn_10.pdb" in files and "10180_trj_10.dcd" in files), "download_GPCRmd_files of dynid 10 weren't correctly downloaded"
    
    
    
# @pytest.mark.filterwarnings("ignore::UserWarning:MDAnalysis")
@pytest.fixture(scope="module", autouse=True)
def prot(path):
    prot = AlloViz.Protein(GPCR=10, path=path, special_res={"HSP": "H"})
    yield prot
    del prot


    
def test_pdbname(prot):
    assert "10182_dyn_10.pdb" in prot.pdb
    
    
def test_processing(path):
    files = os.listdir(path + "/data")
    assert ("protein.pdb" in files and "traj_1.xtc" in files), "Structure and trajectory files weren't correctly processed"
    

    
def test_gpcrdb_numbering(prot):
    np.testing.assert_allclose(
        prot.protein.residues[100:150].atoms.select_atoms("name CA").tempfactors,
        np.array([  3.49000001,  3.5       ,  3.50999999,  3.51999998,  3.52999997,
                    3.53999996,  3.54999995,  3.55999994, 34.5       , 34.50999832,
                    34.52000046, 34.52999878, 34.54000092, 34.54999924, 34.56000137,
                    34.56999969,  4.38000011,  4.38999987,  4.4000001 ,  4.40999985,
                    4.42000008,  4.42999983,  4.44000006,  4.44999981,  4.46000004,
                    4.46999979,  4.48000002,  4.48999977,  4.5       ,  4.51000023,
                    4.51999998,  4.53000021,  4.53999996,  4.55000019,  4.55999994,
                    4.57000017,  4.57999992,  4.59000015,  4.5999999 ,  4.61000013,
                    4.61999989,  4.63000011,  4.63999987,  0.        ,  0.        ,
                    0.        ,  0.        ,  0.        ,  0.        ,  0.        ]),
        err_msg = "There was a problem retrieving GPCR generic numbering from GPCRdb"
    )
    
def test_comfile(path):
    assert os.path.isfile(f"{path}/data/COM_trajs/ca.pdb"), "COM topology file wasn't generated"
    
def test_protein_sel(prot, path):
    with open(f"{path}/data/protein.pdb", "r") as processed, open(f"{path}/10182_dyn_10.pdb", "r") as original:
        assert len(processed.readlines()) < len(original.readlines()), "protein.pdb file hasn't reduced the atom number of the original topology"