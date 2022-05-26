import sys, os, io, re, pandas, time, requests

from multiprocess import Pool
import MDAnalysis as mda
import numpy as np
from Bio import AlignIO

from importlib import import_module
from lazyasd import LazyObject
matplotlib = LazyObject(lambda: import_module('matplotlib'), globals(), 'matplotlib')
nglview = LazyObject(lambda: import_module('nglview'), globals(), 'nglview')
pl = LazyObject(lambda: import_module('matplotlib.pyplot'), globals(), 'pl')
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
filterbyl = ["whole", "incontact", "intercontact"]
#filterbyl = ["whole"]





class Pair:
    def __init__(self, refstate, state2, **kwargs):
        self.state1, self.state2 = refstate, state2
        self.states = [self.state1, self.state2]
        for state in self.states:
            setattr(state, "_pair", self)
            
            # translate_ix = lambda ix: tuple(to_aln_pos[_] for _ in ix) if isinstance(ix, tuple) else to_aln_pos[ix]
            
            
        self._aln = self._make_struct_aln(**kwargs)
        
        # nodes_subset = self._add_nodes_subset()
        # self.sources_subset, self.targets_subset = nodes_subset.values()
    
    
    
    def get_delta(self):
        delta = self.state1 - self.state2
        setattr(self, "delta", delta)
        setattr(self.delta, "pair", self)
        return delta
    
    
    
    def _make_struct_aln(self, aln_method="TMalign_pair"):
        # ****** Pairwise Structural Alignment Methods:
        # --------------------------------------------
        # align_pdbpair        built_in                                                   [pg:        t_coffee is  Installed][built_in]
        # lalign_pdbpair       built_in                                                   [pg:        t_coffee is  Installed][built_in]
        # extern_pdbpair       built_in                                                   [pg:        t_coffee is  Installed][built_in]
        # thread_pair          built_in                                                   [pg:        t_coffee is  Installed][built_in]
        # fugue_pair           http://mizuguchilab.org/fugue/                             [pg:        fugueali is NOT Installed][]
        # pdb_pair             built_in                                                   [pg:        t_coffee is  Installed][built_in]
        # sap_pair             https://mathbio.crick.ac.uk/wiki/Software#SAP              [pg:             sap is  Installed][/gpcr/users/frann/networks/tcoffee/bin/sap]
        # sara_pair            built_in                                                   [pg:        t_coffee is  Installed][built_in]
        # daliweb_pair         built_in                                                   [pg:     dalilite.pl is  Installed][built_in]
        # dali_pair            built_in                                                   [pg:     dalilite.pl is  Installed][built_in]
        # mustang_pair         http://lcb.infotech.monash.edu.au/mustang/                 [pg:         mustang is  Installed][/gpcr/users/frann/networks/tcoffee/bin/mustang]
        # TMalign_pair         http://zhanglab.ccmb.med.umich.edu/TM-align/TMalign.f      [pg:         TMalign is  Installed][/gpcr/users/frann/networks/tcoffee/bin/TMalign]
        
        from subprocess import Popen, PIPE
        from Bio.SeqUtils import seq1
        
        tcoffee = Popen(f"t_coffee -pdb={self.state1._protpdb},{self.state2._protpdb} -method {aln_method} -outfile no -template_file no -output clustalw -align".split(" "),
                        stdout=PIPE, stderr=PIPE, encoding="utf-8")
        
        stderr = tcoffee.stderr.read()
        if "error" not in stderr.lower():
            alignment = AlignIO.read(io.StringIO(tcoffee.stdout.read()), "clustal")
        else:
            print(stderr)
            raise Exception("Structural alignment of the structures couldn't be performed.")
        
        for state, aln in zip(self.states, alignment):
            aas = [f"{res.resname}:{res.resid}" for res in state._prot.residues]
            aln_mapper = dict( (aas.pop(0), pos) for pos, res in enumerate(aln.seq) if res != "-" and res == seq1(aas[0].split(":")[0], custom_map = state._res_dict) )
            aln_mapper.update( dict(map(reversed, aln_mapper.items())) )
            setattr(state, "_aln_mapper", aln_mapper)
        
        return alignment
        
     
    
    
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
    
    
    
    
    
#     def calculate(self, pkg="all", cores=1, **kwargs): # , ow=False, filterby="incontact"
#         pkgs = pkgsl if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
        
#         if any(["COM" in pkg for pkg in pkgs]):
#             for state in self.states: state._add_comtrajs()
        
#         if cores>1:
#             mypool = Pool(cores)
#             utils.pool = mypool
#         print(utils.pool)
        
#         for state in self.states:
#             for pkg in pkgs:
#                 self._set_pkgclass(state, pkg, **kwargs)
        
#         if cores>1:
#             mypool.close()
#             mypool.join()
#             utils.pool = utils.dummypool()

    
    
    
#     def analyze(self, pkg="all", metrics="all", filterby="incontact", element:list=["edges"], normalize=True): # ow
#         pkgs = pkgsl if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
#         metrics = metricsl if metrics=="all" else metrics if isinstance(metrics, list) else [metrics]
#         filterbys = filterbyl if filterby=="all" else filterby if isinstance(filterby, list) else [filterby]
        
#         if cores>1:
#             mypool = Pool(cores)
#             utils.pool = mypool
#         print(utils.pool)
        
#         for state in self.states:
#             for filterby in filterbys:
#                 for pkg in pkgs:
#                     self._set_anaclass(state, pkg, metrics, filterby, element, normalize)
        
#         if cores>1:
#             mypool.close()
#             mypool.join()
#             utils.pool = utils.dummypool()


    
    
    def view(self, pkg, metric, normalize=True, filterby="incontact", num=20):
        # norm = self.state1._get_norm_str(normalize)
        # self.state1.pkg.filterby._norm
        if not hasattr(self, "delta"): self.get_delta()
        if not rhasattr(self.delta, capitalize(pkg), capitalize(filterby)): self.analyze(pkg, filterby=filterby, normalize=normalize)
        return rgetattr(self.delta, capitalize(pkg), capitalize(filterby)).view(metric, num)

    
    
    
    def view(self, pkg, metric, filterby="incontact", element:list=["edges"], num=20, colors=["orange", "turquoise"]):
        # norm = self._get_norm_str(normalize)
        # if not rhasattr(self.data, filterby, norm, pkg):
        #     self.analyze(pkg, normalize=normalize, filterby=filterby)
        
        get_element = lambda element: rgetattr(self, capitalize(pkg), capitalize(filterby), element.lower())
        
        nv = get_element(element[0]).view(metric, num, colors)
        
        if len(element) == 2:
            nv = get_element(element[1]).view(metric, num, colors, nv)
            
        return nv
    
    
    


class Store:
    pass

class State:    
    def __init__(self, pdb='', trajs:list=[], path='', psf=None, parameters=None, GPCR=False, **kwargs):
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
        
        self._res_dict = self._get_res_dict(**kwargs)
        self.mdau, self._prot = self._get_mdau()
        self._dihedral_residx = self._get_dihedral_residx(self._prot)
        self._translate_ix = lambda mapper: lambda ix: tuple(mapper[_] for _ in ix) if isinstance(ix, tuple) else mapper[ix]
        
        
        
    _download_files = trajutils.download_files
    _get_mdau = trajutils.get_mdau
    _add_comtrajs = trajutils.add_comtrajs
    _make_dcds = trajutils.make_dcds
    _get_bonded_cys = trajutils.get_bonded_cys
    _get_res_dict = trajutils.get_res_dict
    _get_dihedral_residx = staticmethod(trajutils.get_dihedral_residx)
    
        
    
    
    def __sub__(self, other):
        delta = Store()
        
        for pkg in (key for key in self.__dict__ if key.lower() in [x.lower() for x in pkgsl] and key in other.__dict__):
            setattr(delta, pkg, Store())
            for filterby in (key for key in getattr(self, pkg).__dict__ if key.lower() in ["whole", "incontact", "intercontact", "raw"] and key in getattr(other, pkg).__dict__): #if not re.search("(^_|raw)", key)
                setattr(getattr(delta, pkg), filterby, Store())
                for elem in (key for key in rgetattr(self, pkg, filterby).__dict__ if key.lower() in ["nodes", "edges"] and key in rgetattr(other, pkg, filterby).__dict__): 
                # setattr(getattr(delta, pkg), filterby, Store())
                # for norm in (key for key in rgetattr(self, pkg, filterby).__dict__ if key in ["norm", "no_norm"] and key in rgetattr(self, pkg, filterby).__dict__):
                    # dif = rgetattr(self, pkg, filterby, norm) - rgetattr(other, pkg, filterby, norm)
                    # setattr(rgetattr(delta, pkg, filterby), norm, Edges(self._pair, dif))
                    dif = rgetattr(self, pkg, filterby, elem) - rgetattr(other, pkg, filterby, elem)
                    elemclass = eval(elem.capitalize())
                    setattr(rgetattr(delta, pkg, filterby), elem, elemclass(self._pair, dif))
                    
        return delta
        
        
        
    
    def calculate(self, pkg="all", cores=1, **kwargs): # ow
        pkgs = pkgsl if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
        
        if any(["COM" in pkg for pkg in pkgs]):
            self._add_comtrajs()
            
        if any([re.search("(carma|grinn)", pkg.lower()) for pkg in pkgs]):
            self._make_dcds()
        
        d = self.__dict__.copy()
        d.update(kwargs)
        
        utils.pool = utils.dummypool()
        if cores>1:
            mypool = Pool(cores)
            utils.pool = mypool
        print(utils.pool)
        # taskcpus = kwargs.pop("taskcpus") if "taskcpus" in kwargs else cores
        
        for pkg in pkgs: #self._set_pkgclass(self, pkg, d) #**kwargs)
            pkgclass = eval(f"Pkgs.{capitalize(pkg)}") if isinstance(pkg, str) else pkg
            if not hasattr(self, pkgclass.__name__):
                setattr(self, pkgclass.__name__, pkgclass(self, d))#**kwargs))
        
        if cores>1:
            mypool.close()
            mypool.join()
            utils.pool = utils.dummypool()
            
        
        
    # def _set_pkgclass(self, state, pkg, d): #**kwargs):
    #     pkgclass = eval(f"Pkgs.{capitalize(pkg)}") if isinstance(pkg, str) else pkg
    #     if not hasattr(state, pkgclass.__name__):
    #         setattr(state, pkgclass.__name__, pkgclass(state, d))#**kwargs))
    
    
    
    
    def analyze(self, pkg="all", filterby="incontact", element="edges", metrics="all", normalize=True, cores=1): # ow
        pkgs = [pkg for pkg in self.__dict__ if pkg.lower() in (pkg.lower() for pkg in pkgl)] if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
        metrics = metricsl if metrics=="all" else metrics if isinstance(metrics, list) else [metrics]
        filterbys = filterbyl if filterby=="all" else filterby if isinstance(filterby, list) else [filterby]
        elements = element if isinstance(element, list) else [element]
        
        utils.pool = utils.dummypool()
        if cores>1:
            mypool = Pool(cores)
            utils.pool = mypool
        print(utils.pool)
        
        for pkg in pkgs:
            pkgclass = eval(f"Pkgs.{capitalize(pkg)}") if isinstance(pkg, str) else pkg
            pkg = rgetattr(self, pkgclass.__name__)
            if not pkg:
                print(f"{pkgclass.__name__} calculation results are needed first")
                continue
                
            for filterby in filterbys:
                anaclass = eval(f"{capitalize(filterby)}") if isinstance(filterby, str) else filterby
                filterby = lambda: rgetattr(pkg, anaclass.__name__)
                if not filterby():
                     setattr(pkg, anaclass.__name__, anaclass(pkg))
                filterby()._add_metrics(elements, metrics, normalize)
                
                
        # for filterby in filterbys:
        #     filterby = rgetattr(self, pkgclass.__name__)
        #     for pkg in pkgs:
        #         # self._set_anaclass(self, pkg, metrics, filterby, element, normalize)
        #         pkgclass = eval(f"Pkgs.{capitalize(pkg)}") if isinstance(pkg, str) else pkg
        #         pkgobj = getattr(self, pkgclass.__name__)
        #         anaclass = eval(f"{capitalize(filterby)}") if isinstance(filterby, str) else filterby
        #         if not hasattr(pkgobj, anaclass.__name__):
        #             setattr(pkgobj, anaclass.__name__, anaclass(pkgobj, metrics, element, normalize))
        
        if cores>1:
            mypool.close()
            mypool.join()
            utils.pool = utils.dummypool()
    
    
#     def _set_anaclass(self, state, pkg, metrics, filterby, element, normalize):
#         pkgclass = eval(f"Pkgs.{capitalize(pkg)}") if isinstance(pkg, str) else pkg
#         pkgobj = getattr(state, pkgclass.__name__)

#         anaclass = eval(f"{capitalize(filterby)}") if isinstance(filterby, str) else filterby
#         if not hasattr(pkgobj, anaclass.__name__):
#             setattr(pkgobj, anaclass.__name__, anaclass(pkgobj, metrics, element, normalize))
        
    
    
        
    def view(self, pkg, metric, filterby="incontact", element:list=["edges"], num=20, colors=["orange", "turquoise"]):
        # norm = self._get_norm_str(normalize)
        # if not rhasattr(self.data, filterby, norm, pkg):
        #     self.analyze(pkg, normalize=normalize, filterby=filterby)
        
        get_element = lambda element: rgetattr(self, capitalize(pkg), capitalize(filterby), element.lower())
        
        nv = get_element(element[0]).view(metric, num, colors)
        
        if len(element) == 2:
            nv = get_element(element[1]).view(metric, num, colors, nv)
            
        return nv






class Element:
    def __init__(self, parent, df):
        self._parent = parent
        self.df = df
    
    
    
    def __repr__(self):
        print(self.df.shape)
        return repr(self.df.iloc[:, :1])
    
    
    def __sub__(self, other):
#         if any(col not in other.df.columns for col in data.df.columns):
        
        selfreindex = self.df.index.map(self._parent._translate_ix(self._parent._aln_mapper))
        selfdf = self.df.abs().reset_index(drop=True).assign(aln_pos=selfreindex.to_numpy()).set_index("aln_pos")
        
        otherreindex = other.df.index.map(other._parent._translate_ix(other._parent._aln_mapper))
        otherdf = other.df.abs().reset_index(drop=True).assign(aln_pos=otherreindex.to_numpy()).set_index("aln_pos")
        

        cols = [col for col in selfdf.columns if col in otherdf.columns]
    
        subs = [col for col in cols if "avg" in col]
        sub = pandas.DataFrame.sub(selfdf[subs], otherdf[subs], axis=0, level="aln_pos").dropna() #fill_value = 0, level="aln_pos"
    
        adds = [col for col in cols if "std" in col]
        add = pandas.DataFrame.add(selfdf[adds], otherdf[adds], axis=0, level="aln_pos").dropna() #fill_value = 0, axis=0, level="aln_pos"
        
        return pandas.concat([add, sub], axis=1)
    
    
    
    
#     def _get_cmap(self): # could be a class attr; even an Analysis or even State/Pair attr
#         if isinstance(self._parent, Pair): 
#             # colors = {"inactive": "r", "active": "g",
#             #          "Gprotein": "y", "Barr": "b"}
#             # color2, color1 = colors[self._parent.state1._pdbf], colors[self._parent.state2._pdbf]
#             color2, color1 = ["r", "g"]
#         else:
#             color1, color2 = "orange", "turquoise"
        
#         return matplotlib.colors.LinearSegmentedColormap.from_list('bar', [color1, "w", color2], 2048)
    
    
    
    
    def _get_colors(self, col, cmap):
        # scale = [0, col.min(), col.max()] if isinstance(self._parent, Pair) else [col.mean(), col.min(), col.max()]
        scale = [0, col.min(), col.max()] if col.min() < 0 else [col.mean(), col.min(), col.max()]
        normdata = matplotlib.colors.TwoSlopeNorm(*scale)
        return cmap( normdata( col.to_numpy() ).data )
    
    
    
    def _show_cbar(self, cmap, minv, maxv):
        pl.imshow([[minv,maxv],], cmap = cmap)
        pl.gca().set_visible(False)
        if isinstance(self._parent, Pair):
            cbar = pl.colorbar(orientation = "horizontal", ticks = [minv,maxv])
            cbar.ax.set_xticklabels([self._parent.state2._pdbf.capitalize(), self._parent.state1._pdbf.capitalize()])
        else:
            cbar = pl.colorbar(orientation = "horizontal", ticks = [minv,maxv])

        return
    
    
    def _get_nv(self, nv):
        mdau_parent = self._parent.state1 if isinstance(self._parent, Pair) else self._parent
        
        if nv is None:
            prot = mda.core.universe.Merge(mdau_parent.mdau.select_atoms("(same segid as protein) or segid LIG"))
            nv = nglview.show_mdanalysis(prot, default=False)
            nv.add_cartoon('protein', color='white')
        
        return nv, mdau_parent
    
    
    
    
    def view(self, metric, num=20, colors=["orange", "turquoise"], nv=None):
        metric = f"{metric}_avg" if not re.search("_avg$", metric) else metric
        # if metric not in df.columns etc
        data = self.df.sort_values(metric, key = abs, ascending = False)
        get_subset = lambda num: data[0:num]
        subset = get_subset(num)
        
        if isinstance(self._parent, Pair):
            while not any(0 < subset[metric]) or not any(0 > subset[metric]):
                num += 1
                subset = get_subset(num)
            print(num)
            
        color1, color2 = colors
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list('bar', [color1, "w", color2], 2048)

        self._show_cbar(cmap, data[metric].min(), data[metric].max())
        colors = self._get_colors(data[metric], cmap)[:num]

        error = "weight_std" if "weight_avg" in metric else metric.replace("avg", "std")
        sizes = np.interp(subset[error], (subset[error].min(), subset[error].max()), (1, 0.1))
        
        nv, mdau_parent = self._get_nv(nv)
        mdau = mdau_parent.mdau
        
        indices = subset.index if isinstance(subset.index[0][0], str) else subset.index.map(mdau_parent._translate_ix(mdau_parent._aln_mapper))
        
        for i in range(len(subset[metric])):
            self._add_element(nv, mdau, indices[i], colors[i], sizes[i])

        return nv
    
    
    

    

    
class Edges(Element):
    def __init__(self, *args):
        super().__init__(*args)
        
        

    def _add_element(self, nv, prot, edge, color, size):
        get_coords = lambda res: list( prot.select_atoms(f"resnum {res.split(':')[-1]} and name CA").positions[0] )
        coords = [get_coords(res) for res in edge]

        return nv.shape.add_cylinder(coords[0], coords[1], list(color),
                                     np.float64(size), f"{edge[0]}_{edge[1]}")
    

    
class Nodes(Element):
    def __init__(self, *args):
        super().__init__(*args)
        
        
       
    def _add_element(self, nv, prot, node, color, size):
        
        get_coords = lambda res: list( prot.select_atoms(f"resnum {res.split(':')[-1]} and name CA").positions[0] )
        coords = get_coords(node)

        return nv.shape.add_sphere(coords, list(color),
                                   np.float64(size)*2.5, f"{node}")
    
    
    
    
class Analysis: #(_Edges)
    def __init__(self, pkg):
        self._pkg = pkg
        # self._parent = self.pkg.state
        # self.metrics = metrics
        # self.normalize = normalize
        self._name = self.__class__.__name__
        
        self._path = f"{self._pkg.state._datadir}/{self._pkg._name}/{self._name}" # lambda norm might not be needed
        os.makedirs(self._path, exist_ok=True)
        self._datapq = lambda element, metric: f"{self._path}/{element}_{metric}.pq"
        
        self._filtdata = self._get_filt_data()
        # self._graph = nx_from_pandas(df=self._filtdata.reset_index(), 
        #                                    source="level_0", target="level_1", 
        #                                    edge_attr=list(self._filtdata.drop("weight_std", axis=1).columns))
        
        # self.add_metrics(metrics, element, normalize)
    
    
    
    def _get_filt_data(self):
        return self._pkg.raw#[["weight_avg", "weight_std"]]
    
    
        
    
    def _add_metrics(self, element, metrics, normalize=True):
        metrics = metricsl if metrics=="all" else metrics if isinstance(metrics, list) else [metrics]
        elements = element if isinstance(element, list) else [element]
        
        for elem in elements:
            if not rhasattr(self, elem, "df"):
                if elem == "edges":
                    data = self._filtdata[["weight_avg", "weight_std"]]
                elif elem == "nodes":
                    data = pandas.DataFrame()
            else:
                data = rgetattr(self, elem, "df")


            pqs = lambda elem: [self._datapq(elem, metric) for metric in metrics]
            no_exist = lambda pqs: [not os.path.isfile(pq) for pq in pqs]
            # is_in_df = any([f"{metric}_avg" not in data.columns for metric in metrics])
            
            # if not is_in_df and 
            if any(no_exist(pqs(elem))): # or ow
                for metric in (metric for metric in metrics if no_exist(pqs(elem))[pqs(elem).index(self._datapq(elem, metric))]):
                    self._analyze(metric, elem, normalize, self._datapq(elem, metric))


            def wait_analyze(elem, data):
                while any(no_exist(pqs(elem))):
                    time.sleep(5)
                return elem, data
                
            
            def add_data(args):
                elem, data = args
                print(f"adding analyzed {elem} {self._pkg} {self._name} data of for {self._pkg.state._pdbf}")

                for pq in pqs(elem):
                    metric = pq.rsplit("/", 1)[-1].split("_", 1)[-1].split(".")[0]
                    df = pandas.read_parquet(pq)

                    cols = [f"{metric}_{num}" for num in self._pkg.state._trajs]
                    df[f"{metric}_avg"] = df[cols].fillna(0).mean(axis=1)
                    df[f"{metric}_std"] = df[cols].fillna(0).std(axis=1)
                    out = df.drop(cols, axis=1)
                    data = pandas.concat([data, out], axis=1)#data.join(out) if data is not None else out

                elemclass = eval(elem.capitalize())
                setattr(self, elem, elemclass(self._pkg.state, data))

                
            utils.get_pool().apply_async(wait_analyze,
                                   args=(elem, data),
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
                nodes = {"sources": self._pkg.state.sources_subset, "targets": self._pkg.state.targets_subset}  
                
                
                
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
        weights = column[column != 0].dropna().abs().rename("weight").reset_index() # btw calculations fail with 0 value weightsa and cfb prob with negative
        #it doesn't make sense either to keep them for others# .rename(columns={f"{column.name}": "weight"})
        # print(column.name, type(column.name), sum(weights["weight"].isna()), weights["weight"].max())
        network = nx_from_pandas(weights, "level_0", "level_1", "weight")
        
        try:
            analyzed = metricf(network, normalized=normalize, weight="weight", **nodes)
            # print(column.name, pandas.Series(analyzed).max())
        except: # LinAlgError
            print("Singular matrix!", self._pkg._name, self._name, elem, metricf.__name__)
            analyzed = {k: 0 for k in eval(f"network.{elem.lower()}")}#{tuple(sorted(k, key = lambda x: int(x.split(":")[-1]))): 0 for k in network.edges()}
        
        sort_index = lambda result: {tuple(sorted(k, key = lambda x: int(x.split(":")[-1]))): result[k] for k in result}
        result = sort_index(analyzed) if elem == "edges" else analyzed# if elem == "nodes"
        
        return pandas.Series(result)#, pq
    
    
    
    
class Whole(Analysis):
    pass    
    
    
    
class Incontact(Analysis):
    def __init__(self, *args):
        super().__init__(*args)
    
    
    def _get_filt_data(self):
        df = super()._get_filt_data()
        
        try:
            indices = self._pkg.state.Getcontacts.raw.index
        except:
            raise Exception("Getcontacts results are needed first")
#         if not rhasattr(self._pkg.state, "Getcontacts", "raw"):
#             print("Getcontacts results are needed; sending calculation first...")
            
#             pool = utils.get_pool()
            
#             self._pkg.state._set_pkgclass(self._pkg.state, "Getcontacts", taskcpus = int(np.ceil(pool._processes/2)))
            
            
#             gc = self._pkg.state.Getcontacts
#             pqs = [gc._rawpq(xtc) for xtc in gc.state._trajs]
#             no_exist = lambda pqs: [not os.path.isfile(pq) for pq in pqs]

#             while any(no_exist(pqs)):
#                 print("Waiting for pq files to be created")
#                 time.sleep(5)
                
#             while not rhasattr(gc, "raw"):
#                 print("Waiting for raw data to be added to object")
#                 time.sleep(5)
                
                
#         indexes = self._pkg.state.Getcontacts.raw.index
            
        return df.filter(indices, axis=0) # maybe pass indexes in class creation
        
        
class Intercontact(Incontact):
    def __init__(self, *args):
        super().__init__(*args)

        
    def _get_filt_data(self):
        df = super()._get_filt_data()
        
        
        def get_intercontacts(indexl):
            resnum = lambda res: int(res.rsplit(":")[-1])
            return [idx for idx in indexl if abs(resnum(idx[0]) - resnum(idx[1]) ) > 5]
        
        return df.filter(get_intercontacts(df.index), axis=0)

        
