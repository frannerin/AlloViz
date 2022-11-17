import os, shutil
import pytest
import AlloViz
import numpy as np
import pandas

pytestmark = pytest.mark.filterwarnings("ignore:::MDAnalysis", "ignore::DeprecationWarning", "ignore::PendingDeprecationWarning")


@pytest.fixture(scope="module")
def path(tmp_path_factory):
    path = tmp_path_factory.mktemp("data")
    yield str(path)
    shutil.rmtree(str(path))

@pytest.fixture(scope="module")
def test_path(request):
    return request.fspath.dirname

@pytest.fixture(scope="module", autouse=True)
def single_traj(test_path, path):
    prot = AlloViz.Protein(pdb = f"{test_path}/data/protein.pdb",
                           trajs = f"{test_path}/data/traj_1.xtc",
                           path = path)
    yield prot
    del prot
    
@pytest.fixture(scope="module", autouse=True)
def single_traj(test_path, path):
    prot = AlloViz.Protein(pdb = f"{test_path}/data/protein.pdb",
                           trajs = f"{test_path}/data/traj_1.xtc",
                           path = path)
    yield prot
    del prot
    

pkgs = [ "MDTASK", "pytraj_CB", "dynetan",
         "correlationplus_COM_LMI", 
         "correlationplus_Phi", "correlationplus_Phi", 
         "correlationplus_Backbone_Dihs_Avg", "correlationplus_Backbone_Dihs_Max",
         "CARDS_MI_Phi", #"AlloViz_Chi4",
         #"MDEntropy_Phi", "MDEntropy_AlphaAngle",
         "GetContacts", "PyInteraph2_Energy",
         "PyInteraph2_Atomic_Contacts_Occurrence", "PyInteraph2_COM_Contacts_Corrected"]
######## MDEntropy-related and AlloViz_Chi4 can't be included in the regression tests because the data in the raw attribute is normalized in comparison to the saved .pq

    
# @pytest.mark.filterwarnings("ignore:::MDAnalysis")
# list(set(AlloViz.AlloViz.utils.pkgsl) - {"gRINN", "MDEntropy_Contacts"})
@pytest.mark.parametrize("pkg", pkgs)
def test_submodules(single_traj, pkg):
    single_traj.calculate(pkg, taskcpus=os.cpu_count())
    
    
@pytest.mark.parametrize("pkg", pkgs)
def test_calculation_reg(single_traj, test_path, pkg):
    reg = pandas.read_parquet(f"{test_path}/data/{pkg}/raw/1.pq")
    if pkg in ["GetContacts"]:
        getattr(single_traj, pkg).raw.sort_index(inplace=True)
        reg.sort_index(inplace=True)
    np.testing.assert_array_equal(getattr(single_traj, pkg).raw.index.values, reg.index.values)
    np.testing.assert_allclose(getattr(single_traj, pkg).raw.values, reg.values)