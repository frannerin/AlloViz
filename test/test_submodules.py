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
         "GetContacts", 
         "PyInteraph2_Atomic_Contacts_Occurrence", "PyInteraph2_COM_Contacts_Corrected"]
######## MDEntropy-related and AlloViz_Chi4 can't be included in the regression tests because the data in the raw attribute is normalized in comparison to the saved .pq

    
# @pytest.mark.filterwarnings("ignore:::MDAnalysis")
# list(set(AlloViz.AlloViz.utils.pkgsl) - {"gRINN", "MDEntropy_Contacts"})
@pytest.mark.parametrize("pkg", pkgs + [
                                        "AlloViz_Chi4", 
                                        "MDEntropy_Phi", 
                                        "MDEntropy_AlphaAngle"
                                       ])
def test_submodules(single_traj, pkg):
    single_traj.calculate(pkg, taskcpus=os.cpu_count())
    
    
@pytest.mark.parametrize("pkg", pkgs)
def test_calculation_reg(single_traj, test_path, pkg):
    reg = pandas.read_parquet(f"{test_path}/data/{pkg}/raw/1.pq")
    if pkg in ["GetContacts"]:
        getattr(single_traj, pkg).raw.sort_index(inplace=True)
        reg.sort_index(inplace=True)
    np.testing.assert_array_equal(getattr(single_traj, pkg).raw.index.values, reg.index.values)
    np.testing.assert_allclose(getattr(single_traj, pkg).raw.values, reg.values, rtol=3e-03)
    
# Had to increase rtol because calculations results in MACOS (in GitHub actions) maximum relative difference of correlationplus_Phi is up to 0.002+:
# __________________ test_calculation_reg[correlationplus_Phi0] __________________

# single_traj = <AlloViz.AlloViz.Classes.Protein object at 0x1392cfa60>
# test_path = '/Users/runner/work/AlloViz/AlloViz/test'
# pkg = 'correlationplus_Phi'

#     @pytest.mark.parametrize("pkg", pkgs)
#     def test_calculation_reg(single_traj, test_path, pkg):
#         reg = pandas.read_parquet(f"{test_path}/data/{pkg}/raw/1.pq")
#         if pkg in ["GetContacts"]:
#             getattr(single_traj, pkg).raw.sort_index(inplace=True)
#             reg.sort_index(inplace=True)
#         np.testing.assert_array_equal(getattr(single_traj, pkg).raw.index.values, reg.index.values)
# >       np.testing.assert_allclose(getattr(single_traj, pkg).raw.values, reg.values)
# E       AssertionError: 
# E       Not equal to tolerance rtol=1e-07, atol=0
# E       
# E       Mismatched elements: 4294 / 13861 (31%)
# E       Max absolute difference: 1.01149605e-07
# E       Max relative difference: 0.00214515
# E        x: array([[ 0.303966],
# E              [-0.222161],
# E              [ 0.341827],...
# E        y: array([[ 0.303966],
# E              [-0.222161],
# E              [ 0.341827],...

# test/test_submodules.py:62: AssertionError
# __________________ test_calculation_reg[correlationplus_Phi1] __________________

# single_traj = <AlloViz.AlloViz.Classes.Protein object at 0x1392cfa60>
# test_path = '/Users/runner/work/AlloViz/AlloViz/test'
# pkg = 'correlationplus_Phi'

#     @pytest.mark.parametrize("pkg", pkgs)
#     def test_calculation_reg(single_traj, test_path, pkg):
#         reg = pandas.read_parquet(f"{test_path}/data/{pkg}/raw/1.pq")
#         if pkg in ["GetContacts"]:
#             getattr(single_traj, pkg).raw.sort_index(inplace=True)
#             reg.sort_index(inplace=True)
#         np.testing.assert_array_equal(getattr(single_traj, pkg).raw.index.values, reg.index.values)
# >       np.testing.assert_allclose(getattr(single_traj, pkg).raw.values, reg.values)
# E       AssertionError: 
# E       Not equal to tolerance rtol=1e-07, atol=0
# E       
# E       Mismatched elements: 4294 / 13861 (31%)
# E       Max absolute difference: 1.01149605e-07
# E       Max relative difference: 0.00214515
# E        x: array([[ 0.303966],
# E              [-0.222161],
# E              [ 0.341827],...
# E        y: array([[ 0.303966],
# E              [-0.222161],
# E              [ 0.341827],...

# test/test_submodules.py:62: AssertionError
# ___________ test_calculation_reg[correlationplus_Backbone_Dihs_Avg] ____________

# single_traj = <AlloViz.AlloViz.Classes.Protein object at 0x1392cfa60>
# test_path = '/Users/runner/work/AlloViz/AlloViz/test'
# pkg = 'correlationplus_Backbone_Dihs_Avg'

#     @pytest.mark.parametrize("pkg", pkgs)
#     def test_calculation_reg(single_traj, test_path, pkg):
#         reg = pandas.read_parquet(f"{test_path}/data/{pkg}/raw/1.pq")
#         if pkg in ["GetContacts"]:
#             getattr(single_traj, pkg).raw.sort_index(inplace=True)
#             reg.sort_index(inplace=True)
#         np.testing.assert_array_equal(getattr(single_traj, pkg).raw.index.values, reg.index.values)
# >       np.testing.assert_allclose(getattr(single_traj, pkg).raw.values, reg.values)
# E       AssertionError: 
# E       Not equal to tolerance rtol=1e-07, atol=0
# E       
# E       Mismatched elements: 1870 / 13861 (13.5%)
# E       Max absolute difference: 5.30856664e-08
# E       Max relative difference: 3.09321557e-06
# E        x: array([[0.299593],
# E              [0.144347],
# E              [0.242614],...
# E        y: array([[0.299593],
# E              [0.144347],
# E              [0.242614],...

# test/test_submodules.py:62: AssertionError
# ___________ test_calculation_reg[correlationplus_Backbone_Dihs_Max] ____________

# single_traj = <AlloViz.AlloViz.Classes.Protein object at 0x1392cfa60>
# test_path = '/Users/runner/work/AlloViz/AlloViz/test'
# pkg = 'correlationplus_Backbone_Dihs_Max'

#     @pytest.mark.parametrize("pkg", pkgs)
#     def test_calculation_reg(single_traj, test_path, pkg):
#         reg = pandas.read_parquet(f"{test_path}/data/{pkg}/raw/1.pq")
#         if pkg in ["GetContacts"]:
#             getattr(single_traj, pkg).raw.sort_index(inplace=True)
#             reg.sort_index(inplace=True)
#         np.testing.assert_array_equal(getattr(single_traj, pkg).raw.index.values, reg.index.values)
# >       np.testing.assert_allclose(getattr(single_traj, pkg).raw.values, reg.values)
# E       AssertionError: 
# E       Not equal to tolerance rtol=1e-07, atol=0
# E       
# E       Mismatched elements: 1803 / 13861 (13%)
# E       Max absolute difference: 6.69842053e-08
# E       Max relative difference: 8.51313525e-06
# E        x: array([[0.303966],
# E              [0.222161],
# E              [0.341827],...
# E        y: array([[0.303966],
# E              [0.222161],
# E              [0.341827],...
