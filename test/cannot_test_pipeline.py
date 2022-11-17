import os, shutil
import pytest
import AlloViz
import multiprocess


pytestmark = pytest.mark.filterwarnings("ignore:::MDAnalysis", "ignore::DeprecationWarning", "ignore::PendingDeprecationWarning")


@pytest.fixture(scope="module")
def path(tmp_path_factory):
    path = tmp_path_factory.mktemp("data")
    yield str(path)
    # shutil.rmtree(str(path))

@pytest.fixture(scope="module")
def test_path(request):
    return request.fspath.dirname

@pytest.fixture(scope="module", autouse=True)
def prot(test_path, path):
    prot = AlloViz.Protein(pdb = f"{test_path}/data/protein.pdb",
                           trajs = [f"{test_path}/data/traj_2.xtc",
                                    f"{test_path}/data/traj_3.xtc"],
                           path = f"{test_path}/pipetest")
    yield prot
    del prot
    
def test_calculation(prot, test_path, path):
    # prot.calculate(["correlationplus_CA_Pear", "GetContacts"], cores=2)
    assert prot.pdb == f"{test_path}/data/protein.pdb"
    assert os.path.isfile(prot.pdb)
    assert f"{test_path}/data/traj_2.xtc" in prot._trajs.values()
    assert os.path.isfile(f"{test_path}/data/traj_2.xtc")
    assert f"{test_path}/data/traj_3.xtc" in prot._trajs.values()
    assert os.path.isfile(f"{test_path}/data/traj_3.xtc")
    assert os.path.isfile(prot._comtrajs[1])
    
def test_pool(prot, test_path, path):
    def func(prot):
        assert os.path.isfile(prot._trajs[1])
    with multiprocess.Pool(2) as p:
        p.apply_async(func, (prot,))
        p.close()
        p.join()
    assert True
    
    # assert True
    
# def test_calculation_results(prot, test_path):
    
    
# def test_filterings(prot):
#     prot.filter(filterings=["GetContacts_edges", "All", ["Spatially_distant", "No_Sequence_Neighbors"]], Sequence_Neighbor_distance=6)
#     assert True

# def test_analysis(prot):
#     prot.correlationplus_CA_Pear.GetContacts_edges.analyze(["edges", "nodes"], cores=2)
#     assert True
    
# def test_viz(prot):
#     nv = prot.correlationplus_CA_Pear.GetContacts_edges.edges.view("btw")
#     assert nv.n_components == 21, "Visualization doesn't seem to be rendering correctly"
    