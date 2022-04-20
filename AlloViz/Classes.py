import sys, os, io, re, pandas, time, requests

from multiprocess import Pool
import MDAnalysis as mda
import numpy as np

from importlib import import_module
from lazyasd import LazyObject
matplotlib = LazyObject(lambda: import_module('matplotlib'), globals(), 'matplotlib')
nglview = LazyObject(lambda: import_module('nglview'), globals(), 'nglview')
pl = LazyObject(lambda: import_module('pyplot', package='matplotlib'), globals(), 'matplotlib')
# import matplotlib, nglview#, ipywidgets, matplotlib.cm
# from matplotlib import pyplot as pl

from networkx import from_pandas_edgelist as nx_from_pandas
# import networkx.algorithms.centrality as nx_centrality#import edge_betweenness_centrality, edge_betweenness_centrality_subset # edge_betweenness
from networkx.algorithms.centrality import edge_betweenness_centrality, edge_betweenness_centrality_subset
from networkx.algorithms.centrality import edge_current_flow_betweenness_centrality, edge_current_flow_betweenness_centrality_subset
from networkx.algorithms.centrality import betweenness_centrality, betweenness_centrality_subset
from networkx.algorithms.centrality import current_flow_betweenness_centrality, current_flow_betweenness_centrality_subset

from . import Pkgs
from . import trajutils

from . import utils
rgetattr = utils.rgetattr
rhasattr = utils.rhasattr
capitalize = utils.capitalize


pkgsl = ["MDTASK", "getcontacts", "pyInteraph", "pyInteraphEne", "dynetan", "dynetanCOM", "pytrajCA", "pytrajCB",
         "corrplus", "corrplusLMI", "corrplusCOM", "corrplusCOMLMI", "corrplusPsi", "corrplusPhi", "corrplusOmega", "corrplusDihs",
        "gRINN", "gRINNcorr", "g_corrCAMI", "g_corrCOMMI", "g_corrCALMI", "g_corrCOMLMI",
        "AlloVizPsi", "AlloVizPhi", "AlloVizOmega", "AlloVizDihs",
        "MDEntropyContacts", "MDEntropyDihs", "MDEntropyAlphaAngle"]
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
        
        # nodes_subset = self._add_nodes_subset()
        # self.sources_subset, self.targets_subset = nodes_subset.values()
    
    
    
    def get_delta(self):
        delta = self.state1 - self.state2
        setattr(self, "delta", delta)
        setattr(self.delta, "pair", self)
        return delta
        
     
    
    
#     def _add_nodes_subset(self):
#         bfac = lambda atom: f"{atom.tempfactor:.2f}"
#         is_wb = lambda atom: atom.name == "N" and -8.1 <= atom.tempfactor <= 8.1 and (bfac(atom) != "1.00" and bfac(atom) != "0.00")
        
        
#         targets = [3.50, 6.30, 7.49, 7.50, 7.51, 7.52, 7.53]
#         # In Ballesteros-Weinstein: Ionic lock 3.50 and 6.30; NPxxY 7.49-7.53
    
    
#         def get_sources(states):
#             resnum_to_wb = lambda nums, state: (bfac(atom)
#                                                 for residue in nums
#                                                 for atom in state.mdau.select_atoms(f"resnum {residue}")
#                                                 if is_wb(atom))
            
#             has_lig = [state if ("LIG" in (seg.segid for seg in state.mdau.segments)) else False for state in self.states]
#             if any(has_lig):
#                 state = next(item for item in has_lig if item != False)
#                 aas = state.mdau.select_atoms("(same residue as around 4 segid LIG) and protein").residues
#                 return list(resnum_to_wb(aas.resnums, state))
#             else:
#                 return [3.32] # Conserved position
      
    
#         nodes_subset = {"sources": [val for val in get_sources(self.states) if val not in targets],
#                         "targets": list(map(lambda x: f"{x:.2f}", targets))}
        
        
        
#         format_res = lambda aa: f"A:{aa.resname}:{aa.resnum}"
#         wb_to_aa = lambda wblist, state: list(map(format_res, (atom.residue
#                                                                for residue in state.mdau.select_atoms("protein").residues
#                                                                for atom in residue.atoms
#                                                                if is_wb(atom) and bfac(atom) in wblist)))
        
#         for state in self.states:
#             state.sources_subset, state.targets_subset = wb_to_aa(nodes_subset["sources"], state), wb_to_aa(nodes_subset["targets"], state) 
        
#         return nodes_subset    
    
    
    
    
    
    def calculate(self, pkg="all", cores=1, **kwargs): # , ow=False, filterby="incontact"
        pkgs = pkgsl if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
        
        if any(["COM" in pkg for pkg in pkgs]):
            for state in self.states: state._add_comtrajs()
        
        if cores>1:
            mypool = Pool(cores)
            utils.pool = mypool
        print(utils.pool)
        
        for state in self.states:
            for pkg in pkgs:
                self._set_pkgclass(state, pkg, **kwargs)
        
        if cores>1:
            mypool.close()
            mypool.join()
            utils.pool = utils.dummypool()

    
    
    
    def analyze(self, pkg="all", metrics="all", filterby="incontact", element:list=["edges"], normalize=True): # ow
        pkgs = pkgsl if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
        metrics = metricsl if metrics=="all" else metrics if isinstance(metrics, list) else [metrics]
        filterbys = filterbyl if filterby=="all" else filterby if isinstance(filterby, list) else [filterby]
        
        if cores>1:
            mypool = Pool(cores)
            utils.pool = mypool
        print(utils.pool)
        
        for state in self.states:
            for filterby in filterbys:
                for pkg in pkgs:
                    self._set_anaclass(state, pkg, metrics, filterby, element, normalize)
        
        if cores>1:
            mypool.close()
            mypool.join()
            utils.pool = utils.dummypool()


    
    
    def view(self, pkg, metric, normalize=True, filterby="incontact", num=20):
        # norm = self.state1._get_norm_str(normalize)
        # self.state1.pkg.filterby._norm
        if not hasattr(self, "delta"): self.get_delta()
        if not rhasattr(self.delta, capitalize(pkg), capitalize(filterby)): self.analyze(pkg, filterby=filterby, normalize=normalize)
        return rgetattr(self.delta, capitalize(pkg), capitalize(filterby)).view(metric, num)




class Store:
    pass

class State:    
    def __init__(self, pdb='', trajs:list=[], path='', psf=None, parameters=None, GPCR=False):
        self.GPCR = GPCR
        
        if not isinstance(self.GPCR, bool) and isinstance(self.GPCR, int):
            self._gpcrmdid = self.GPCR
            self._path = f"{self.GPCR}"
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
        self.mdau, self._dihedral_residx = self._get_mdau()
        
        
    _download_files = trajutils._download_files
    _get_mdau = trajutils._get_mdau
    _add_comtrajs = trajutils._add_comtrajs
    _make_dcds = trajutils._make_dcds
    
    
    
#     def __getnewargs_ex__(self):
#         if not isinstance(self.GPCR, bool) and isinstance(self.GPCR, int):
#             return ((), {"GPCR": self.GPCR})#self.state,
        
#         else:
#             if hasattr(self, "_psff") and hasattr(self, "_paramf"):
#                 extra = {"psf": self._psff,
#                          "parameters": self._paramf}
#             else:
#                 extra = {}
            
#             return ((), {"pdb": self._pdbf,
#                         "trajs": self._trajs.values(),
#                         "path": self._path,
#                         "GPCR": self.GPCR}.update(extra))
        
    
#     def __getstate__(self):
#         return self.__dict__

#     def __setstate__(self, statedict):
#         self.__dict__.update(statedict)
        
        
#     def __copy__(self):
#         return self

#     def __deepcopy__(self, memo):
#         return self
        
    
    
    def __sub__(self, other):
        delta = Store()
        
        for pkg in (key for key in self.__dict__ if key.lower() in [x.lower() for x in pkgsl] and key in other.__dict__):
            setattr(delta, pkg, Store())
            for filterby in (key for key in getattr(self, pkg).__dict__ if key.lower() in filterbyl and key in getattr(other, pkg).__dict__): #if not re.search("(^_|raw)", key)
                # setattr(getattr(delta, pkg), filterby, Store())
                # for norm in (key for key in rgetattr(self, pkg, filterby).__dict__ if key in ["norm", "no_norm"] and key in rgetattr(self, pkg, filterby).__dict__):
                    # dif = rgetattr(self, pkg, filterby, norm) - rgetattr(other, pkg, filterby, norm)
                    # setattr(rgetattr(delta, pkg, filterby), norm, Edges(self._pair, dif))
                dif = rgetattr(self, pkg, filterby) - rgetattr(other, pkg, filterby)
                setattr(rgetattr(delta, pkg), filterby, _Edges(self._pair, dif))
                    
        return delta
        
    
    
    
    
    def calculate(self, pkg="all", cores=1, **kwargs): # ow
        pkgs = pkgsl if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
        
        if any(["COM" in pkg for pkg in pkgs]):
            self._add_comtrajs()
            
        if any([re.search("(carma|grinn)", pkg.lower()) for pkg in pkgs]):
            self._make_dcds()
        
        utils.pool = utils.dummypool()
        if cores>1:
            mypool = Pool(cores)
            utils.pool = mypool
        print(utils.pool)
        taskcpus = kwargs.pop("taskcpus") if "taskcpus" in kwargs else cores
        
        for pkg in pkgs: self._set_pkgclass(self, pkg, taskcpus=taskcpus, **kwargs)
        
        if cores>1:
            mypool.close()
            mypool.join()
            utils.pool = utils.dummypool()
            
        
        
    def _set_pkgclass(self, state, pkg, **kwargs):
        pkgclass = eval(f"Pkgs.{capitalize(pkg)}") if isinstance(pkg, str) else pkg
        if not hasattr(state, pkgclass.__name__):
            setattr(state, pkgclass.__name__, pkgclass(state, **kwargs))
    
    
    
    
    def analyze(self, pkg="all", metrics="all", filterby="incontact", element:list=["edges"], normalize=True, cores=1): # ow
        pkgs = pkgsl if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
        metrics = metricsl if metrics=="all" else metrics if isinstance(metrics, list) else [metrics]
        filterbys = filterbyl if filterby=="all" else filterby if isinstance(filterby, list) else [filterby]
        
        utils.pool = utils.dummypool()
        if cores>1:
            mypool = Pool(cores)
            utils.pool = mypool
        print(utils.pool)
        
        for filterby in filterbys:
            for pkg in pkgs:
                self._set_anaclass(self, pkg, metrics, filterby, element, normalize)
        
        if cores>1:
            mypool.close()
            mypool.join()
            utils.pool = utils.dummypool()
    
    
    def _set_anaclass(self, state, pkg, metrics, filterby, element, normalize):
        pkgclass = eval(f"Pkgs.{capitalize(pkg)}") if isinstance(pkg, str) else pkg
        pkgobj = getattr(state, pkgclass.__name__)

        anaclass = eval(f"{capitalize(filterby)}") if isinstance(filterby, str) else filterby
        if not hasattr(pkgobj, anaclass.__name__):
            setattr(pkgobj, anaclass.__name__, anaclass(pkgobj, metrics, element, normalize))
        
    
    
        
    def view(self, pkg, metric, normalize=True, filterby="incontact", num=20):
        # norm = self._get_norm_str(normalize)
        # if not rhasattr(self.data, filterby, norm, pkg):
        #     self.analyze(pkg, normalize=normalize, filterby=filterby)
            
        return rgetattr(self, capitalize(pkg), capitalize(filterby)).view(metric, num)






class Element:
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
    
    
    
    
    def _get_cmap(self): # could be a class attr; even an Analysis or even State/Pair attr
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
    
    
    

    

    
class Edges(Element):
    def __init__(self, *args):
        super().__init__(*args)
        
        

    def _add_edge(self, nv, prot, edge, color, radius):
        get_coords = lambda resnum: list( prot.select_atoms(f"resnum {resnum} and name CA").center_of_geometry() )
        coords = [get_coords(res.split(':')[-1]) for res in edge]

        return nv.shape.add_cylinder(coords[0], coords[1], list(color),
                                     np.float64(radius), f"{edge[0]}_{edge[1]}")
    

    
class Nodes(Element):
     def __init__(self, *args):
        super().__init__(*args)
    
    
    
    
class Analysis: #(_Edges)
    def __init__(self, pkg, metrics, element, normalize = True):
        self.pkg = pkg
        self._parent = self.pkg.state
        self.metrics = metrics
        self.normalize = normalize
        self._name = self.__class__.__name__
        
        self._path = f"{self.pkg.state._datadir}/{self.pkg._name}/{self._name}" # lambda norm might not be needed
        self._datapq = lambda element, metric: f"{self._path}/{element}_{metric}.pq"
        
        self._filtdata = self._get_filt_data()
        # self._graph = nx_from_pandas(df=self._filtdata.reset_index(), 
        #                                    source="level_0", target="level_1", 
        #                                    edge_attr=list(self._filtdata.drop("weight_std", axis=1).columns))
        
        self.add_metrics(metrics, element, normalize)
    
    
    
    def _get_filt_data(self):
        return self.pkg.raw#[["weight_avg", "weight_std"]]
    
    
        
    
    def add_metrics(self, metrics, element, normalize=True):
        # norm = utils.norm(normalize)
        os.makedirs(self._path, exist_ok=True)
        metrics = metricsl if metrics=="all" else metrics if isinstance(metrics, list) else [metrics]
        
        for elem in element:
            if not rhasattr(self, elem, "df"):
                if elem == "edges":
                    data = self._filtdata[["weight_avg", "weight_std"]]
                elif elem == "nodes":
                    data = None
            else:
                data = getattr(self, elem, "df")


            pqs = [self._datapq(elem, metric) for metric in metrics]
            no_exist = lambda pqs: [not os.path.isfile(pq) for pq in pqs]
            is_in_df = any([f"{metric}_avg" not in data.columns for metric in metrics]) if data != None else False
            
            if not is_in_df and any(no_exist(pqs)): # or ow
                for metric in (metric for metric in self.metrics if no_exist(pqs)[pqs.index(self._datapq(elem, metric))]):
                    self._analyze(metric, elem, normalize, self._datapq(elem, metric))


            def wait_analyze(pqs, data):
                while any(no_exist(pqs)):
                    time.sleep(5)
                return pqs, data

            def get_data(args):
                pqs, data = args
                print(f"adding analyzed {elem} {self.pkg} {self._name} data of for {self.pkg.state._pdbf}")

                for pq in pqs:
                    metric = pq.rsplit("/", 1)[-1].split("_", 1)[-1].split(".")[0]
                    df = pandas.read_parquet(pq)

                    cols = [f"{metric}_{num}" for num in self.pkg.state._trajs]
                    df[f"{metric}_avg"] = df[cols].fillna(0).mean(axis=1)
                    df[f"{metric}_std"] = df[cols].fillna(0).std(axis=1)
                    out = df.drop(cols, axis=1)
                    data = pandas.concat([data, out], axis=1)#data.join(out) if data is not None else out
                return data
            
            elemclass = eval(elem.capitalize())
            add_data = lambda args: setattr(self, elem, elemclass(self.pkg.state, get_data(args)))
            utils.get_pool().apply_async(wait_analyze,
                                   args=(pqs, data),
                                   callback=add_data)
        return
    

    
    
    def _analyze(self, metric, elem, normalize, pq):
        pool = utils.get_pool()
        rawdata = self._filtdata.drop("weight_std", axis=1)
        nodes = {}   
        
        if callable(metric):
            metricf = metric
        else:        
            if "cfb" in metric:
                metricf = current_flow_betweenness_centrality
            elif "btw" in metric:
                metricf = betweenness_centrality
            
            if elem == "edges":
                metricf = eval(f"edge_{metricf.__name__}")
            
            if "subset" in metric:
                metricf = eval(f"{metricf.__name__}_subset")
                nodes = {"sources": self._state.sources_subset, "targets": self._state.targets_subset}  
                
                
                
        def save_pq(df):
            newcolnames = {name: f"{metric}_{name}" for name in df.columns}
            df.rename(columns=newcolnames, inplace=True)
            df.to_parquet(pq)
            return
        
        # pool.apply_async(lambda args: rawdata.apply(self._networkx_analysis, args=args),
        #                  args=((metricf, elem, normalize, pq),),
        #                  callback=save_pq)
        pool.apply_async(lambda: rawdata.apply(self._networkx_analysis, args=(metricf, elem, normalize, nodes)), #pq
                         #args=(,),
                         callback=save_pq)
        
        
        def _calculate_empty(pqf):
            print("sleeping", pqf, os.getpid())
            while not os.path.isfile(pqf):
                time.sleep(5)
            return
        
        for _ in range(len(rawdata.columns)-1): pool.apply_async(_calculate_empty, args=(pq,))
    
    
    
    
    def _networkx_analysis(self, column, metricf, elem, normalize, nodes):#pq
        weights = column.reset_index().fillna(0)#.dropna() # this dropna is problematic; nonexistent should be 0? # .rename(columns={f"{column.name}": "weight"})
        network = nx_from_pandas(weights, "level_0", "level_1", column.name)#"weight")         
        
        sort_index = lambda result: {tuple(sorted(k, key = lambda x: int(x.split(":")[-1]))): result[k] for k in result}
        
        try:
            analyzed = metricf(network, normalized=normalize, weight=column.name, **nodes)
        except: # LinAlgError
            print("Singular matrix!", self.pkg._name, self._name, elem, metricf.__name__)
            analyzed = {k: 0 for k in eval(f"network.{elem.lower()}")}#{tuple(sorted(k, key = lambda x: int(x.split(":")[-1]))): 0 for k in network.edges()}
            
        result = sort_index(analyzed) if elem == "edges" else analyzed# if elem == "nodes"
        
        return pandas.Series(result)#, pq
    
    
    
    
class Whole(Analysis):
    pass    
    
    
    
class Incontact(Analysis):
    def __init__(self, *args):
        super().__init__(*args)
    
    
    def _get_filt_data(self):
        df = super()._get_filt_data()
        
        if not rhasattr(self.pkg.state, "Getcontacts", "raw"):
            print("Getcontacts results are needed; sending calculation first...")
            
            pool = utils.get_pool()
            
            self.pkg.state._set_pkgclass(self.pkg.state, "Getcontacts", taskcpus = int(np.ceil(pool._processes/2)))
            
            
            gc = self.pkg.state.Getcontacts
            pqs = [gc._rawpq(xtc) for xtc in gc.state._trajs]
            no_exist = lambda pqs: [not os.path.isfile(pq) for pq in pqs]

            while any(no_exist(pqs)):
                print("Waiting for pq files to be created")
                time.sleep(5)
                
            while not rhasattr(gc, "raw"):
                print("Waiting for raw data to be added to object")
                time.sleep(5)
                
                
        indexes = self.pkg.state.Getcontacts.raw.index
            
        return df.filter(indexes, axis=0) # maybe pass indexes in class creation
        
        
class Intercontact(Incontact):
    def __init__(self, *args):
        super().__init__(*args)

        
    def _get_filt_data(self):
        df = super()._get_filt_data()
        
        
        def get_intercontacts(indexl):
            resnum = lambda res: int(res.rsplit(":")[-1])
            return [idx for idx in indexl if abs(resnum(idx[0]) - resnum(idx[1]) ) >= 4]
        
        return df.filter(get_intercontacts(df.index), axis=0)

        
