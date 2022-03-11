import sys, os, io, re, pandas, time, requests
import MDAnalysis as mda
import numpy as np
#from multiprocess import Pool, set_start_method
#set_start_method("fork", force=True)
import multiprocess
import matplotlib, nglview#, ipywidgets, matplotlib.cm
from matplotlib import pyplot as pl

from networkx import from_pandas_edgelist as networkx_from_pandas
from networkx.algorithms.centrality import edge_betweenness_centrality, edge_betweenness_centrality_subset # edge_betweenness
from networkx.algorithms.centrality import edge_current_flow_betweenness_centrality, edge_current_flow_betweenness_centrality_subset

from . import Pkgs

from . import utils
rgetattr = utils.rgetattr
rhasattr = utils.rhasattr
capitalize = utils.capitalize
norm = utils.norm


# pkgsl = ["MDTASK", "getcontacts", "pyinteraph", "pyinteraphEne", "dynetan", "pytrajCA", "pytrajCB",
#          "corrplus", "corrplusLMI", "corrplusCOM", "corrplusCOMLMI", "corrplusPsi", "corrplusPhi", "corrplusOmega"] #"dynetanCOM", 
pkgsl = ["MDTASK", "getcontacts", "pyinteraph", "pyinteraphEne", "dynetan", #"pytrajCA", "pytrajCB",
         "corrplus", "corrplusLMI", "corrplusCOM", "corrplusCOMLMI", "corrplusPsi", "corrplusPhi", "corrplusOmega",
        "gRINN", "gRINNcorr", "g_correlationCA", "g_correlationCOM"] #"dynetanCOM", 
# metricsl = ["cfb", "cfb_subset", "btw", "btw_subset"]
metricsl = ["cfb", "btw"]
#filterbyl = ["whole", "incontact", "intercontact"]
filterbyl = ["whole"]





class Pair:
    def __init__(self, refstate, state2):
        self.state1, self.state2 = refstate, state2
        self.states = [self.state1, self.state2]
        for state in self.states:
            state._pair = self
        
        nodes_subset = self._add_nodes_subset()
        self.sources_subset, self.targets_subset = nodes_subset.values()
    
    
    
    def get_delta(self):
        delta = self.state1 - self.state2
        setattr(self, "delta", delta)
        setattr(self.delta, "pair", self)
        return delta
        
     
    
    
    def _add_nodes_subset(self):
        bfac = lambda atom: f"{atom.tempfactor:.2f}"
        is_wb = lambda atom: atom.name == "N" and -8.1 <= atom.tempfactor <= 8.1 and (bfac(atom) != "1.00" and bfac(atom) != "0.00")
        
        
        targets = [3.50, 6.30, 7.49, 7.50, 7.51, 7.52, 7.53]
        # In Ballesteros-Weinstein: Ionic lock 3.50 and 6.30; NPxxY 7.49-7.53
    
    
        def get_sources(states):
            resnum_to_wb = lambda nums, state: (bfac(atom)
                                                for residue in nums
                                                for atom in state.mdau.select_atoms(f"resnum {residue}")
                                                if is_wb(atom))
            
            has_lig = [state if ("LIG" in (seg.segid for seg in state.mdau.segments)) else False for state in self.states]
            if any(has_lig):
                state = next(item for item in has_lig if item != False)
                aas = state.mdau.select_atoms("(same residue as around 4 segid LIG) and protein").residues
                return list(resnum_to_wb(aas.resnums, state))
            else:
                return [3.32] # Conserved position
      
    
        nodes_subset = {"sources": [val for val in get_sources(self.states) if val not in targets],
                        "targets": list(map(lambda x: f"{x:.2f}", targets))}
        
        
        
        format_res = lambda aa: f"A:{aa.resname}:{aa.resnum}"
        wb_to_aa = lambda wblist, state: list(map(format_res, (atom.residue
                                                               for residue in state.mdau.select_atoms("protein").residues
                                                               for atom in residue.atoms
                                                               if is_wb(atom) and bfac(atom) in wblist)))
        
        for state in self.states:
            state.sources_subset, state.targets_subset = wb_to_aa(nodes_subset["sources"], state), wb_to_aa(nodes_subset["targets"], state) 
        
        return nodes_subset    
    
    
    
    
    
    def calculate(self, pkg="all", cores=1): # , ow=False, filterby="incontact"
        pkgs = pkgsl if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
        
        if any(["COM" in pkg for pkg in pkgs]):
            for state in self.states: state._add_comtrajs()
        
        mypool = multiprocess.get_context("fork").Pool(cores)
        utils.pool = mypool
        print(utils.pool)
        
        for state in self.states:
            for pkg in pkgs:
                self._set_pkgclass(state, pkg)
        
        mypool.close()
        mypool.join()
        
        
        
#         pkgs = pkgsl if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
        
#         if filterby == "incontact":
#             if not hasattr(self, "_incontact"):
#                 self.get_incontact(cores)
#             if "getcontacts" in pkgs: pkgs.remove("getcontacts")
#             # with Pool(cores) as pool:
#             #     for state in self.states:
#             #         if not hasattr(state.data, "raw") or not hasattr(state.data.raw, "getcontacts") or ow:
#             #             state._send_calc("getcontacts", ow, filterby, pool)
#             #     pool.close()
#             #     pool.join()
#             # [state._add_pkg("getcontacts") for state in self.states]
#             # pkgs.remove("getcontacts")
        
#         with Pool(cores) as pool:
#             for pkg in pkgs:
#                 for state in self.states:
#                     if not rhasattr(state.data, "raw", pkg) or ow:
#                         #print(f"sending {state.name} {pkg}")
#                         state._send_calc(pkg, ow, filterby, pool)
#             pool.close()
#             pool.join()
        
#         [state._add_pkg(pkg) for pkg in pkgs for state in self.states]
        
#         return
    
    
    
#     def get_incontact(self, cores=None, ow=False):
#         with Pool(cores) as pool:
#             for state in self.states:
#                 if not rhasattr(state.data, "raw", "getcontacts") or ow:
#                     #print(f"sending {state.name} {pkg}")
#                     state._send_calc("getcontacts", ow, "incontact", pool)
#             pool.close()
#             pool.join()
        
#         [state._add_pkg("getcontacts") for state in self.states]
        
#         edges = set((idx for idx in state.data.raw.getcontacts.index for state in self.states))
#         setattr(self, "_incontact", edges)
#         return edges
    
    
    
    
    def analyze(self, pkg="all", metrics="all", filterby="incontact", normalize=True): # ow
        pkgs = pkgsl if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
        metrics = metricsl if metrics=="all" else metrics if isinstance(metrics, list) else [metrics]
        filterbys = filterbyl if filterby=="all" else filterby if isinstance(filterby, list) else [filterby]
        
        mypool = multiprocess.get_context("fork").Pool(cores)
        utils.pool = mypool
        print(utils.pool)
        
        for state in self.states:
            for filterby in filterbys:
                for pkg in pkgs:
                    self._set_anaclass(state, pkg, metrics, filterby, normalize)
        
        mypool.close()
        mypool.join()
        
        
        
        
#         pkgs = pkgsl if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
#         metrics = metricsl if metrics=="all" else metrics if isinstance(metrics, list) else [metrics]
        
#         for state in self.states: state.data._check_raw(pkgs, normalize, filterby)
        
#         ps = []
#         for state in self.states:
#             p = Process(target=state._send_anal, args=(pkgs, metrics, normalize, ow, filterby))
#             #print(f"sent analysis of {state.name} {filterby} data of {pkgs} packages for {metrics} metrics: {p.name}")
#             ps.append(p)
#             p.start()
        
#         [p.join() for p in ps]
        
#         [state.data._add_anal(pkg, metrics, normalize, ow, filterby) for pkg in pkgs for state in self.states]
        
#         return


    
    
    def view(self, pkg, metric, normalize=True, filterby="incontact", num=20):
        # norm = self.state1._get_norm_str(normalize)
        # self.state1.pkg.filterby._norm
        if not hasattr(self, "delta"): self.get_delta()
        if not rhasattr(self.delta, capitalize(pkg), capitalize(filterby), norm(normalize)): self.analyze(pkg, filterby=filterby, normalize=normalize)
        return rgetattr(self.delta, capitalize(pkg), capitalize(filterby), norm(normalize)).view(metric, num)




class Store:
    pass

class State:#(Entity):    
    def __init__(self, pdb='', trajs:list=[], path='', psf=None, parameters=None, GPCRmdID=None):
        if GPCRmdID:
            self._gpcrmdid = GPCRmdID
            self._path = f"{GPCRmdID}"
            os.makedirs(self._path, exist_ok=True)
            if not any([re.search("(pdb$|psf$|xtc$|parameters$)", file) for file in os.listdir(self._path)]):
                self._download_files()
                
            files = os.listdir(self._path)
            self._pdbf = self._path + '/' + next(file for file in files if re.search("pdb$", file))
            self._psff = self._path + '/' + next(file for file in files if re.search("psf$", file))
            self._paramf = self._path + '/' + next(file for file in files if re.search("parameters$", file))
            self._trajs = dict(enumerate(sorted( f"{self._path}/{traj}" for traj in files if re.search("^(?!\.).*\.xtc$", traj) ), 1))
        
        else:
            psf_params = parameters is not None and psf is not None
            
            files_to_check = trajs + [pdb] if not psf_params else trajs + [pdb, psf, parameters]
            files_exist = {file: os.path.isfile(file) for file in files_to_check}
            if any([not file_exist for file_exist in files_exist.values()]):
                raise FileNotFoundError(f"Some of the files could not be found: {files_exist}")
            
            self._path = path
            os.makedirs(self._path, exist_ok=True)
            self._pdbf = pdb
            self._trajs = dict(enumerate(trajs, 1))
            if psf_params:
                self._psff = psf
                self._paramf = parameters
                
        self._datadir = f"{self._path}/data"
        os.makedirs(self._datadir, exist_ok=True)
        self.mdau = self._get_mdau()
        
    
    
    def __sub__(self, other):
        
        
        delta = Store()
        
        for pkg in (key for key in self.__dict__ if key.lower() in [x.lower() for x in pkgsl] and key in other.__dict__):
            setattr(delta, pkg, Store())
            for filterby in (key for key in getattr(self, pkg).__dict__ if key.lower() in filterbyl and key in getattr(other, pkg).__dict__): #if not re.search("(^_|raw)", key)
                setattr(getattr(delta, pkg), filterby, Store())
                for norm in (key for key in rgetattr(self, pkg, filterby).__dict__ if key in ["norm", "no_norm"] and key in rgetattr(self, pkg, filterby).__dict__):
                    dif = rgetattr(self, pkg, filterby, norm) - rgetattr(other, pkg, filterby, norm)
                    setattr(rgetattr(delta, pkg, filterby), norm, Edges(self._pair, dif))
                    
        return delta
        
    
    
    
    
    def _download_files(self):
        from bs4 import BeautifulSoup
        import urllib.request as pwget
        import tarfile, fileinput

        web = "https://submission.gpcrmd.org"
        
        html = requests.get(f"{web}/dynadb/dynamics/id/{self._gpcrmdid}/")
        soup = BeautifulSoup(html.text, features="html.parser").find_all('a')
        links = [link.get('href') for link in soup if re.search("(xtc|pdb|psf|prm)", link.get('href'))]
        
        for link in links:
            fname = f"{self._path}/{link.rsplit('/')[-1]}"
            print(f"downloading {fname}")
            pwget.urlretrieve(f"{web}{link}", fname)
            
            if re.search("prm", fname):
                with tarfile.open(fname) as tar:
                    tar.extractall(self._path)
                os.remove(fname)
                
                for line in fileinput.input(f"{self._path}/parameters", inplace=True):
                    if line.strip().startswith('HBOND'):
                        line = 'HBOND CUTHB 0.5\n'
                    elif line.strip().startswith('END'):
                        line = 'END\n'
                    sys.stdout.write(line)
                
        return
    
    
    

    def _get_mdau(self):
        mdau = mda.Universe(self._pdbf, *self._trajs.values())
        
        if hasattr(self, "_gpcrmdid"):
            prot = mdau.select_atoms("protein")

            prot_numsf = f"{self._datadir}/gpcrdb_gennums.pdb"

            if not os.path.isfile(prot_numsf):
                print(f"retrieving {prot_numsf}")
                with io.StringIO() as protf:
                    with mda.lib.util.NamedStream(protf, "prot.pdb") as f:
                        prot.write(f, file_format="pdb")
                    response = requests.post('https://gpcrdb.org/services/structure/assign_generic_numbers', 
                                             files = {'pdb_file': protf.getvalue()})
                    with open(prot_numsf, "w") as prot_nums:
                        prot_nums.write(response.text)

            nums = mda.Universe(prot_numsf).select_atoms("protein").tempfactors
            prot.tempfactors = nums.round(2)
                    
        return mdau
    
    
    
    
    def _add_comtrajs(self):
        print(f"Making trajectories of residue COM for {self._pdbf}")
        
        compath = f"{self._datadir}/COMtrajs"
        if not os.path.isdir(compath): os.makedirs(compath, exist_ok=True)

        compdb = f"{compath}/ca.pdb"
        if not os.path.isfile(compdb):
            prot = self.mdau.select_atoms("protein and name CA")
            prot.write(compdb)
        setattr(self, "_compdbf", compdb)

        comtrajs = [f"{compath}/{xtc}.xtc" for xtc in self._trajs]
        setattr(self, "_comtrajs", {num: traj for num, traj in enumerate(comtrajs, 1)})

        for xtc, comtraj in self._comtrajs.items():
            if not os.path.isfile(comtraj):
                prot = self.mdau.select_atoms("protein")
                traj = next(traj for traj in self.mdau.trajectory.readers if traj.filename == self._trajs[xtc])
                arr = np.empty((prot.n_residues, traj.n_frames, 3))
                for ts in traj:
                    arr[:, ts.frame] = prot.center_of_mass(compound='residues')

                cau = mda.Universe(compdb, arr, format=mda.coordinates.memory.MemoryReader, order='afc')
                with mda.Writer(comtraj, cau.atoms.n_atoms) as W:
                    for ts in cau.trajectory:
                        W.write(cau.atoms)
        return
    
    
    
    
    def _make_dcds(self):
        import parmed, mdtraj
        print(f"Making dcd trajectories for {self._pdbf}")
        
        dcdpath = f"{self._datadir}/dcds"
        if not os.path.isdir(dcdpath): os.makedirs(dcdpath, exist_ok=True)
        setattr(self, "_protf", lambda ext: f"{dcdpath}/prot.{ext}")
        
        atoms = self.mdau.select_atoms("protein")
        
        dcdpdb = self._protf("pdb")
        if not os.path.isfile(dcdpdb): atoms.write(dcdpdb) # could it be prot_nums? IT HAS THE LIGAND THOUGH
        
        dcdpsf = self._protf("psf")
        if not os.path.isfile(dcdpsf):
            psf = parmed.load_file(self._psff)[atoms.indices]
            psf.title = self._psff
            psf.write_psf(dcdpsf)
        
        
        dcds = [f"{dcdpath}/{xtc}.dcd" for xtc in self._trajs]
        setattr(self, "_dcds", {num: traj for num, traj in enumerate(dcds, 1)})

        for xtc, dcd in self._dcds.items():
            if not os.path.isfile(dcd):
                traj = mdtraj.load(self._trajs[xtc], top=self._pdbf, atom_indices=atoms.indices)
                traj.save_dcd(dcd)
        return
    
    
    
    
    def calculate(self, pkg="all", cores=1): # ow
        pkgs = pkgsl if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
        
        if any(["COM" in pkg for pkg in pkgs]):
            self._add_comtrajs()
            
        if any([re.search("(carma|grinn)", pkg.lower()) for pkg in pkgs]):
            self._make_dcds()
        
        mypool = multiprocess.get_context("fork").Pool(cores)
        utils.pool = mypool
        print(utils.pool)
        
        for pkg in pkgs: self._set_pkgclass(self, pkg)
        
        mypool.close()
        mypool.join()
        
        
    def _set_pkgclass(self, state, pkg):
        pkgclass = eval(f"Pkgs.{capitalize(pkg)}") if isinstance(pkg, str) else pkg
        if not hasattr(state, pkgclass.__name__):
            setattr(state, pkgclass.__name__, pkgclass(state))
    
    
    
    
    def analyze(self, pkg="all", metrics="all", filterby="incontact", normalize=True, cores=1): # ow
        pkgs = pkgsl if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
        metrics = metricsl if metrics=="all" else metrics if isinstance(metrics, list) else [metrics]
        filterbys = filterbyl if filterby=="all" else filterby if isinstance(filterby, list) else [filterby]
        
        mypool = multiprocess.get_context("fork").Pool(cores)
        utils.pool = mypool
        print(utils.pool)
        
        for filterby in filterbys:
            for pkg in pkgs:
                self._set_anaclass(self, pkg, metrics, filterby, normalize)
        
        mypool.close()
        mypool.join()
    
    
    def _set_anaclass(self, state, pkg, metrics, filterby, normalize):
        pkgclass = eval(f"Pkgs.{capitalize(pkg)}") if isinstance(pkg, str) else pkg
        pkgobj = getattr(state, pkgclass.__name__)

        # anaclass = eval(f"Analysis.{filterby[0].upper() + filterby[1:]}") if isinstance(filterby, str) else filterby
        anaclass = eval(f"{capitalize(filterby)}") if isinstance(filterby, str) else filterby
        if not hasattr(pkgobj, anaclass.__name__):
            setattr(pkgobj, anaclass.__name__, anaclass(pkgobj, metrics, normalize))
        
    
    
        
    def view(self, pkg, metric, normalize=True, filterby="incontact", num=20):
        # norm = self._get_norm_str(normalize)
        # if not rhasattr(self.data, filterby, norm, pkg):
        #     self.analyze(pkg, normalize=normalize, filterby=filterby)
            
        return rgetattr(self, capitalize(pkg), capitalize(filterby), norm(normalize)).view(metric, num)






class Edges:
    def __init__(self, parent, df):
        self._parent = parent
        self.df = df
    
    
    def __repr__(self):
        print(self.df.shape)
        return repr(self.df.iloc[:, :1])
    
    
    def __sub__(self, other):
#         if any(col not in other.df.columns for col in data.df.columns):

        cols = [col for col in self.df.columns if col in other.df.columns]
    
        subs = [col for col in cols if "avg" in col]
        sub = pandas.DataFrame.sub(self.df[subs], other.df[subs], fill_value = 0)
    
        adds = [col for col in cols if "std" in col]
        add = pandas.DataFrame.add(self.df[adds], other.df[adds], fill_value = 0)
        
        return pandas.concat([add, sub], axis=1)
    
    
    
    
    def _get_cmap(self):
        if isinstance(self._parent, Pair): 
            # colors = {"inactive": "r", "active": "g",
            #          "Gprotein": "y", "Barr": "b"}
            # color2, color1 = colors[self._parent.state1._pdbf], colors[self._parent.state2._pdbf]
            color2, color1 = ["r", "g"]
        else:
            color1, color2 = "orange", "turquoise"
        
        return matplotlib.colors.LinearSegmentedColormap.from_list('bar', [color1, "w", color2], 2048)
    
    
    
    
    def _get_colors(self, col):
        cmap = self._get_cmap()
        scale = [0, col.min(), col.max()] if isinstance(self._parent, Pair) else [col.mean(), col.min(), col.max()]
        normdata = matplotlib.colors.TwoSlopeNorm(*scale)
        npdata = col.to_numpy()
        edges = np.nonzero(npdata)
        return cmap( normdata( npdata[edges] ).data )




    def _add_edge(self, nv, prot, edge, color, radius):    
        get_coords = lambda resnum: list( prot.select_atoms(f"resnum {resnum} and name CA").center_of_geometry() )
        coords = [get_coords(res.split(':')[-1]) for res in edge]

        return nv.shape.add_cylinder(coords[0], coords[1], list(color),
                                     np.float64(radius), f"{edge[0]}_{edge[1]}")
    
    
    
    def _show_cbar(self):
        cmap = self._get_cmap()
        pl.imshow([[0,1],], cmap = cmap)
        pl.gca().set_visible(False)
        if isinstance(self._parent, Pair):
            cbar = pl.colorbar(orientation = "horizontal", ticks = [0,1])
            cbar.ax.set_xticklabels([self._parent.state2._pdbf.capitalize(), self._parent.state1._pdbf.capitalize()])
        else:
            cbar = pl.colorbar(orientation = "horizontal", ticks = [])

        return
    
    
    
    def view(self, metric, num=20):
        metric = f"{metric}_avg" if not re.search("_avg$", metric) else metric
        # if metric not in df.columns etc
        get_data = lambda num: self.df.sort_values(metric, key = abs, ascending = False)[0:num]
        data = get_data(num)
        
        if isinstance(self._parent, Pair):
            while not any(0 < data[metric]) or not any(0 > data[metric]):
                num += 1
                data = get_data(num)
            print(num)

        self._show_cbar()
        
        mdau = self._parent.state1.mdau if isinstance(self._parent, Pair) else self._parent.mdau
        prot = mda.core.universe.Merge(mdau.select_atoms("protein or segid LIG"))
        nv = nglview.show_mdanalysis(prot, default=False)
        nv.add_cartoon('protein', color='white')

        edges = np.nonzero(data[metric].to_numpy())
        colors = self._get_colors(data[metric])
        
        error = "weight_std" if "weight_avg" in metric else metric.replace("avg", "std")
        radii = np.interp(data[error], (data[error].min(), data[error].max()), (1, 0.1))

        for i in range(len(edges[0])):
            self._add_edge(nv, mdau, data.index[edges[0][i]], colors[i], radii[i])

        return nv

    
    
    
    
class Analysis:
    def __init__(self, pkg, metrics, normalize = True):
        self.pkg = pkg
        self.metrics = metrics
        self.normalize = normalize
        self._name = self.__class__.__name__
        
        self._path = lambda norm: f"{self.pkg.state._datadir}/{self.pkg._name}/{self._name}/{norm}" # lambda norm might not be needed
        self._datapq = lambda norm, metric: f"{self._path(norm)}/{metric}.pq"
        
        self._filtdata = self._get_filt_data()
        
        self.add_metrics(metrics, normalize)
    
    
    
    def _get_filt_data(self):
        return self.pkg.raw#[["weight_avg", "weight_std"]]
    
    
        
    
    def add_metrics(self, metrics, normalize=True):
        norm = utils.norm(normalize)
        os.makedirs(self._path(norm), exist_ok=True)
        
        if not hasattr(self, norm):
            data = self._filtdata[["weight_avg", "weight_std"]]
        else:
            data = getattr(self, norm).df
        
        metrics = metricsl if metrics=="all" else metrics if isinstance(metrics, list) else [metrics]
        pqs = [self._datapq(norm, metric) for metric in metrics]
        no_exist = lambda pqs: [not os.path.isfile(pq) for pq in pqs]
        
        if any([f"{metric}_avg" not in data.columns for metric in metrics]): # or ow
            if any(no_exist(pqs)):
                for metric in (metric for metric in self.metrics if no_exist(pqs)[pqs.index(self._datapq(norm, metric))]):
                    self._analyze(metric, normalize, self._datapq(norm, metric))
                
        
        def wait_analyze(pqs, data):
            while any(no_exist(pqs)):
                time.sleep(5)
            return pqs, data

        def get_data(args):
            pqs, data = args
            print(f"adding analyzed {norm} data of {self._name} for {self.pkg.state._pdbf}")
            
            for pq in pqs:
                metric = pq.rsplit("/", 1)[-1].split(".")[0]
                df = pandas.read_parquet(pq)

                cols = [f"{metric}_{num}" for num in self.pkg.state._trajs]
                df[f"{metric}_avg"] = df[cols].fillna(0).mean(axis=1)
                df[f"{metric}_std"] = df[cols].fillna(0).std(axis=1)
                out = df.drop(cols, axis=1)
                data = data.join(out) if data is not None else out
            return data
        
        add_data = lambda args: setattr(self, norm, Edges(self.pkg.state, get_data(args)))
        utils.get_pool().apply_async(wait_analyze,
                               args=(pqs, data),
                               callback=add_data)
        return
    

    
    
    def _analyze(self, metric, normalize, pq):
        pool = utils.get_pool()
        rawdata = self._filtdata.drop("weight_std", axis=1)
        
        
        def save_pq(df):
            newcolnames = {name: f"{metric}_{name}" for name in df.columns}
            df.rename(columns=newcolnames, inplace=True)
            df.to_parquet(pq)
            return
        
        pool.apply_async(lambda args: rawdata.apply(self._networkx_analysis, args=args),
                         args=((metric, normalize, pq),),
                         callback=save_pq)
        
        
        def _calculate_empty(self, pqf):
            print("sleeping", pqf, os.getpid())
            while not os.path.isfile(pqf):
                time.sleep(5)
            return
        
        for _ in range(len(rawdata.columns)-1): pool.apply_async(_calculate_empty, args=(pq,))
    
    
    
    
    def _networkx_analysis(self, column, metric, normalize, pq):
        weights = column.reset_index().rename(columns={f"{column.name}": "weight"}).dropna()
        network = networkx_from_pandas(weights, "level_0", "level_1", "weight")
        nodes = {}
        
        if callable(metric):
            metricf = metric
        else:        
            if "cfb" in metric:
                metricf = edge_current_flow_betweenness_centrality
            elif "btw" in metric:
                metricf = edge_betweenness_centrality

            if "subset" in metric:
                metricf = eval(f"{metricf.__name__}_subset")
                nodes = {"sources": self._state.sources_subset, "targets": self._state.targets_subset}            
        
        analyzed = metricf(network, normalized=normalize, weight="weight", **nodes)
        edges = {tuple(sorted(k, key = lambda x: int(x.split(":")[-1]))): analyzed[k] for k in analyzed}
        
        return pandas.Series(edges)#, pq
    
    
    
    
class Incontact(Analysis):
    def __init__(self, pkg, metrics, normalize = True):
        super().__init__(pkg, metrics, normalize)
    
    
    def _get_filt_data(self):
        df = super()._get_filt_data()
        
        try:
            indexes = self.pkg.state.Getcontacts.raw.index
        except:
            print("Getcontacts calculation is needed first")
            raise
            
        return df.filter(indexes, axis=0) # maybe pass indexes in class creation
        
        
class Intercontact(Incontact):
    def __init__(self, pkg, metrics, normalize = True):
        super().__init__(pkg, metrics, normalize)

        
    def _get_filt_data(self):

        def get_intercontacts(indexl):
            get_resnum = lambda res: int(res.rsplit(":")[-1])
            return [idx for idx in indexl if abs(get_resnum(idx[0]) - get_resnum(idx[1]) ) >= 4]
        
        df = super()._get_filt_data()

        return df.filter(get_intercontacts(df.index), axis=0)

        
class Whole(Analysis):
    pass