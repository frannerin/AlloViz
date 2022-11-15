import os, shutil
import pytest
import AlloViz
import MDAnalysis as mda
import numpy as np

pytestmark = pytest.mark.filterwarnings("ignore:::MDAnalysis", "ignore::DeprecationWarning")


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
    

# @pytest.mark.filterwarnings("ignore:::MDAnalysis")
# list(set(AlloViz.AlloViz.utils.pkgsl) - {"gRINN", "MDEntropy_Contacts"})
@pytest.mark.parametrize("pkg", ["MDTASK", "pytraj_CB", "dynetan",
                                 "correlationplus_COM_LMI", 
                                 "correlationplus_Phi", "correlationplus_Phi", 
                                 "correlationplus_Backbone_Dihs_Avg", "correlationplus_Backbone_Dihs_Max",
                                 "AlloViz_Chi4", "CARDS_MI_Phi",
                                 "MDEntropy_Phi", "MDEntropy_AlphaAngle",
                                 "GetContacts", "PyInteraph2_Energy",
                                 "PyInteraph2_Atomic_Contacts_Occurrence", "PyInteraph2_COM_Contacts_Corrected"])
def test_submodules(single_traj, pkg):
    single_traj.calculate(pkg, taskcpus=os.cpu_count())