from . import Pkgs
import sys, os, io, re, pandas, time
import MDAnalysis as mda
import numpy as np
# import .Pkgs
from .utils import *

pkgsl = ["getcontacts", "pyinteraph", "pyinteraphEne", "dynetan", "dynetanCOM", "pytrajCA", "pytrajCB", "corrplus", "corrplusLMI", "corrplusCOM", "corrplusCOMLMI"]
metricsl = ["cfb", "cfb_subset", "btw", "btw_subset"]





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
            resnum_to_wb = lambda nums, state: (bfac(atom)                                                 for residue in nums                                                 for atom in state.mdau.select_atoms(f"resnum {residue}")                                                 if is_wb(atom))
            
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
        wb_to_aa = lambda wblist, state: list(map(format_res, (atom.residue                                                       for residue in state.mdau.select_atoms("protein").residues                                                       for atom in residue.atoms                                                       if is_wb(atom) and bfac(atom) in wblist)))
        for state in self.states:
            state.sources_subset, state.targets_subset = wb_to_aa(nodes_subset["sources"], state), wb_to_aa(nodes_subset["targets"], state)       
        
        
        return nodes_subset    
    
    
    
    
    
    def calculate(self, pkg="all", cores=None, ow=False, filterby="incontact"):  
        pkgs = pkgsl if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
        
        if filterby == "incontact":
            if not hasattr(self, "_incontact"):
                self.get_incontact(cores)
            if "getcontacts" in pkgs: pkgs.remove("getcontacts")
            # with Pool(cores) as pool:
            #     for state in self.states:
            #         if not hasattr(state.data, "raw") or not hasattr(state.data.raw, "getcontacts") or ow:
            #             state._send_calc("getcontacts", ow, filterby, pool)
            #     pool.close()
            #     pool.join()
            # [state._add_pkg("getcontacts") for state in self.states]
            # pkgs.remove("getcontacts")
        
        with Pool(cores) as pool:
            for pkg in pkgs:
                for state in self.states:
                    if not rhasattr(state.data, "raw", pkg) or ow: #hasattr(state.data, "raw") or not hasattr(state.data.raw, pkg) or ow:
                        #print(f"sending {state.name} {pkg}")
                        state._send_calc(pkg, ow, filterby, pool)
            pool.close()
            pool.join()
        
        [state._add_pkg(pkg) for pkg in pkgs for state in self.states]
        
        return
    
    
    
    def get_incontact(self, cores=None, ow=False):
        with Pool(cores) as pool:
            for state in self.states:
                if not rhasattr(state.data, "raw", "getcontacts") or ow: #hasattr(state.data, "raw") or not hasattr(state.data.raw, "getcontacts") or ow:
                    #print(f"sending {state.name} {pkg}")
                    state._send_calc("getcontacts", ow, "incontact", pool)
            pool.close()
            pool.join()
        
        [state._add_pkg("getcontacts") for state in self.states]
        
        edges = set((idx for idx in state.data.raw.getcontacts.index for state in self.states))
        setattr(self, "_incontact", edges)
        return edges
    
    
    
    
    def analyze(self, pkg="all", metrics="all", normalize=True, ow=False, filterby="incontact"):
        pkgs = pkgsl if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
        metrics = metricsl if metrics=="all" else metrics if isinstance(metrics, list) else [metrics]
        
        for state in self.states: state.data._check_raw(pkgs, normalize, filterby)
        
        ps = []
        for state in self.states:
            p = Process(target=state._send_anal, args=(pkgs, metrics, normalize, ow, filterby))
            #print(f"sent analysis of {state.name} {filterby} data of {pkgs} packages for {metrics} metrics: {p.name}")
            ps.append(p)
            p.start()
        
        [p.join() for p in ps]
        
        [state.data._add_anal(pkg, metrics, normalize, ow, filterby) for pkg in pkgs for state in self.states]
        
        return


    
    
    def view(self, pkg, metric, normalize=True, filterby="incontact", num=20):
        norm = self.state1._get_norm_str(normalize)
        if not hasattr(self, "delta"): self.get_delta()
        if not rhasattr(self.delta, filterby, norm, pkg): self.analyze(pkg) #hasattr(self.delta, filterby) or not hasattr(getattr(self.delta, filterby), norm) or not hasattr(rgetattr(self.delta, filterby, norm), pkg): self.analyze(pkg)
        return rgetattr(self.delta, filterby, norm, pkg).view(metric, num, filterby)











class State:    
    def __init__(self, name, idx):
        self.name = name
        self.idx = idx
        self.data = Data(self)
        
        if not os.path.isdir(self.name): self._download_files()
        
        self._pdbf = f"{self.name}/{self.name}_1.pdb"
        self._psff = f"{self.name}/{self.name}_1.psf"
        self._trajs = dict(enumerate(sorted( traj for traj in os.listdir(self.name) if re.search("^(?!\.).*\.xtc$", traj) ), 1))
        
        self.mdau = self._get_mdau()
        
    
    
    def __sub__(self, other):        
        return self.data - other.data
    
    
    def _rawpn(self, pkg):
        return f"{self.name}/data/raw/{pkg}"
    
    def _rawfn(self, pkg, xtc):
        filen = xtc if isinstance(xtc, int) else xtc.rsplit('.', 1)[0]
        return f"{self._rawpn(pkg)}/{filen}.pq"
    
    
    
    def _get_norm_str(self, normalize):
        return "norm" if normalize else "no_norm"
    
    def _datapn(self, filterby, normalize, pkg):
        norm = self._get_norm_str(normalize)
        return f"{self.name}/data/{filterby}/{norm}/{pkg}"
    
    def _datafn(self, filterby, normalize, pkg, metric):
        return f"{self._datapn(filterby, normalize, pkg)}/{metric}.pq"
        
    
    
    
    
    def _download_files(self):
        from bs4 import BeautifulSoup
        import urllib.request as pwget
        import requests


        web = "https://submission.gpcrmd.org"

        os.makedirs(self.name, exist_ok=True)

        html = requests.get(f"{web}/dynadb/dynamics/id/{self.idx}/")
        soup = BeautifulSoup(html.text, features="html.parser").find_all('a')
        links = (link.get('href') for link in soup if re.search("(xtc|pdb|psf)", link.get('href')))
        
        for link in links:
            ext = f".{link.rsplit('.', 1)[-1]}"
            num = len([file for file in os.listdir(self.name) if ext in file])+1
            fname = f"{self.name}/{self.name}_{num}{ext}"
            print(f"downloading {fname}")
            pwget.urlretrieve(f"{web}{link}", fname)
            
        return
    
    

    def _get_mdau(self):
        trajs = [f"{self.name}/{self._trajs[xtc]}" for xtc in self._trajs]
        mdau = mda.Universe(self._pdbf, *trajs)
        prot = mdau.select_atoms("protein")
        
        prot_numsf = f"{self.name}/{self.name}_nums.pdb"
        
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
        compath = self._rawpn("COMtrajs")
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
                traj = next(traj for traj in self.mdau.trajectory.readers if traj.filename == f"{self.name}/{self._trajs[xtc]}")
                # newu = mda.Universe(self._pdbf, f"{self.name}/{self._trajs[xtc]}")
                # prot = newu.select_atoms("protein")
                arr = np.empty((prot.n_residues, traj.n_frames, 3)) # newu.trajectory.n_frames
                for ts in traj: # newu.trajectory
                    arr[:, ts.frame] = prot.center_of_mass(compound='residues')

                cau = mda.Universe(compdb, arr, format=mda.coordinates.memory.MemoryReader, order='afc')
                # cau.load_new(arr, format=mda.coordinates.memory.MemoryReader, order='afc')
                with mda.Writer(comtraj, cau.atoms.n_atoms) as W:
                    for ts in cau.trajectory:
                        W.write(cau.atoms)
        
        # To add cau
        # cau = mda.Universe(compdb, *comtrajs)
        # setattr(self, "_cau", cau)

        return #compdb, comtraj
    
    
    
    
    
#     def calculate(self, pkg="all", cores=None, ow=False, filterby="incontact"):
#         pkgs = pkgsl if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
        
#         if filterby == "incontact":
#             if hasattr(self, "_pair"):
#                 self._pair.get_incontact(cores)
#             elif not rhasattr(self.data, "raw", "getcontacts"): #hasattr(self.data, "raw") or not hasattr(self.data.raw, "getcontacts") or o
#                 with Pool(cores) as pool:
#                     self._send_calc("getcontacts", filterby, pool)
#                     pool.close()
#                     pool.join()
#                 self._add_pkg("getcontacts")
#             if "getcontacts" in pkgs: pkgs.remove("getcontacts")
            
        
#         with Pool(cores) as pool:
#             for pkg in pkgs:
#                 if not rhasattr(self.data, "raw", pkg) or ow:#hasattr(self.data, "raw") or not hasattr(self.data.raw, pkg) or ow:
#                     #print(f"sending {self.name} {pkg}")
#                     self._send_calc(pkg, ow, filterby, pool)
#             pool.close()
#             pool.join()
        
#         [self._add_pkg(pkg) for pkg in pkgs]
    def calculate(self, pkg="all"):
        pkgs = pkgsl if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
        
        if any(["COM" in pkg for pkg in pkgs]):
            self._add_comtrajs()
        
        
        # pool = get_pool()
        
        for pkg in pkgs:
            pkgclass = eval(f"Pkgs.{pkg[0].upper() + pkg[1:]}") if isinstance(pkg, str) else pkg
            setattr(self.data, pkgclass.__name__, pkgclass(self))
            print(pkg)
    
        
        
        
    
    
    
#     def _send_calc(self, pkg, ow, filterby, pool):
#         def send_empty(self, pqf):
#             while not os.path.isfile(pqf):
#                 time.sleep(5)
#             return
        
#         os.makedirs(self._rawpn(pkg), exist_ok=True)        
#         pkgf = self._getcontactsf if "getcontacts" in pkg else self._corrsf
#         # pkgf = getcontactsf if "getcontacts" in pkg else corrsf
#         if "COM" in pkg and not hasattr(self, "_comtrajs"):#or not hasattr(self, "_cau")
#             self._add_comtrajs()
        
#         for xtc in self._trajs:
#             pqf = self._rawfn(pkg, xtc)
            
#             if not os.path.isfile(pqf) or ow:
#                 print(f"sending calculation of {pkg} for {self.name} {xtc}")
#                 pool.apply_async(pkgf, (xtc, pkg, ow, filterby))
#                 # if "getcontacts" in pkg or "dynetan" in pkg:
#                 #     [pool.apply_async(send_empty, (pqf,))]*2
            
#         return
    
#     def _getcontactsf(self, xtc, pkg, ow, filterby):
#         getcontactsf(self, xtc, pkg, ow)
#         return
    
#     def _corrsf(self, xtc, pkg, ow, filterby):
#         corrsf(self, xtc, pkg, ow, filterby)
#         return
    
    
#     def _add_pkg(self, pkg):
#         if not hasattr(self.data, "raw"): setattr(self.data, "raw", Store())
#         print(f"adding raw data of {pkg} for {self.name}")
        
#         pqs = [self._rawfn(pkg, xtc) for xtc in self._trajs]
#         flist = map(lambda pq: pandas.read_parquet(pq), pqs)
#         df = pandas.concat(flist, axis=1)
#         cols = [f"{num}" for num in self._trajs]
#         df["weight_avg"] = df[cols].fillna(0).mean(axis=1)
#         df["weight_std"] = df[cols].fillna(0).std(axis=1)
        
#         setattr(self.data.raw, pkg, df)
        
#         return
    
    
    
    
    
    
    def analyze(self, pkg="all", metrics="all", normalize=True, ow=False, filterby="incontact"):
        pkgs = pkgsl if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
        metrics = metricsl if metrics=="all" else metrics if isinstance(metrics, list) else [metrics]
        
        self.data._check_raw(pkgs, normalize, filterby)
        self._send_anal(pkgs, metrics, normalize, ow, filterby)
        
        [self.data._add_anal(pkg, metrics, normalize, ow, filterby) for pkg in pkgs]
    
    
    def _send_anal(self, pkgs, metrics, normalize, ow, filterby):
        if not os.path.isdir(self._datapn(filterby, normalize, "")) or         any([not os.path.isfile(self._datafn(filterby, normalize, pkg, metric)) 
             for pkg in pkgs for metric in metrics]) or ow:
            ps = []
            print(f"sending analysis of {self.name} {filterby} data of {pkgs} packages for {self._get_norm_str(normalize)} {metrics} metrics")
            for pkg in pkgs:
                p = Process(target=self.data._analyze, args=(pkg, metrics, normalize, ow, filterby))
                ps.append(p)
                p.start()
            [p.join() for p in ps]
        
        
        return  
    
    
    
    
    def view(self, pkg, metric, normalize=True, filterby="incontact", num=20):
        norm = self._get_norm_str(normalize)
        if not rhasattr(self.data, filterby, norm, pkg):#hasattr(self.data, filterby) or not hasattr(getattr(self.data, filterby), norm) or not hasattr(rgetattr(self.data, filterby, norm), pkg): 
            self.analyze(pkg, normalize=normalize, filterby=filterby)
            
        return rgetattr(self.data, filterby, norm, pkg).view(metric, num, filterby)









class Data:
    def __init__(self, state):
        self._state = state
        
        from networkx import from_pandas_edgelist as networkx_from_pandas
        from networkx.algorithms.centrality import edge_betweenness, edge_betweenness_centrality, edge_betweenness_centrality_subset, edge_current_flow_betweenness_centrality, edge_current_flow_betweenness_centrality_subset
    
    
    
    def __sub__(self, other):
        delta = Store()
        
        for filterby in (key for key in self.__dict__ if (not re.search("(^_|raw)", key) and key in other.__dict__)):
            setattr(delta, filterby, Store())
            for norm in (key for key in getattr(self, filterby).__dict__ if key in getattr(other, filterby).__dict__): #if not re.search("(^_|raw)", key)
                setattr(getattr(delta, filterby), norm, Store())
                get_data = lambda obj: rgetattr(obj, filterby, norm)
                
                subdelta = get_data(delta)

                selfdata = get_data(self)
                otherdata = get_data(other)

                pkgs = [key for key in selfdata.__dict__ if key in otherdata.__dict__ and not re.search("^_", key)]            
                for pkg in pkgs:
                    print(f"calculating delta for {filterby} {norm} {pkg} data")
                    dif = getattr(selfdata, pkg) - getattr(otherdata, pkg)
                    setattr(subdelta, pkg, Edges(self._state._pair, dif))
        
        return delta
        
        
    
    def _check_raw(self, pkgs, normalize, filterby):
        norm = self._state._get_norm_str(normalize)
        
        if not hasattr(self, filterby): 
            setattr(self, filterby, Store())
        if not rhasattr(self, filterby, norm):#hasattr(getattr(self, filterby), norm): 
            setattr(getattr(self, filterby), norm, Store())
            
        if any([not hasattr(self.raw, pkg) for pkg in pkgs]):
            self._state.calculate(pkgs)
            
            
        if filterby == "incontact" or filterby == "intercontact":
            if not rhasattr(self, "raw", "getcontacts"):#hasattr(self, "raw") or not hasattr(self.raw, "getcontacts"): 
                self._state.calculate("getcontacts")
            
        for pkg in (pkg for pkg in pkgs if not rhasattr(self, filterby, norm, pkg)):#hasattr(getattr(self.incontact, norm), pkg)):
            raw = getattr(self.raw, pkg)
            if filterby == "incontact" or filterby == "intercontact":
                raw = raw.filter(self.raw.getcontacts.index, axis=0)
                if filterby == "intercontact":
                    raw = raw.filter(get_intercontacts(raw.index), axis=0)
            setattr(rgetattr(self, filterby, norm), pkg, Edges(self._state, raw[["weight_avg", "weight_std"]]))
        
        # elif filterby == "whole":
        #     for pkg in (pkg for pkg in pkgs if not rhasattr(self, filterby, norm, pkg)):
        #         setattr(rgetattr(self, filterby, norm), pkg, Edges(self._state, getattr(self.raw, pkg)[["weight_avg", "weight_std"]]))
            
            
            

    def _analyze(self, pkg, metrics, normalize, ow, filterby):
        norm = self._state._get_norm_str(normalize)
        
        
        os.makedirs(self._state._datapn(filterby, normalize, pkg), exist_ok=True)
        raw = getattr(self.raw, pkg)
        if filterby == "incontact" or filterby == "intercontact":
            raw = raw.filter(self.raw.getcontacts.index, axis=0)
            if filterby == "intercontact":
                get_resnum = lambda res: int(res.rsplit(":")[-1])
                raw = raw.filter(get_intercontacts(raw.index), axis=0)
        if filterby == "whole":
            metrics = [metric for metric in metrics if "subset" not in metric]
            
        
        data = rgetattr(self, filterby, norm, pkg).df
        
        if any([f"{metric}_avg" not in data.columns for metric in metrics]) or ow:
            #print(f"analyzing {metrics} from {self._state.name} {pkg} data for {filterby} residues")
            raw_data = raw[[col for col in raw.columns if col != "weight_std"]]
            self._send_anal(raw_data, pkg, metrics, normalize, ow, filterby)
            # setattr(self.incontact, pkg, Edges(self._state, data.join(analyzed)))
            
        return
        
        
    def _send_anal(self, rawdata, pkg, metrics, normalize, ow, filterby):
        data = None
        
        if filterby == "whole":
            metrics = [metric for metric in metrics if "subset" not in metric]
        
        for metric in metrics:
            pqf = self._state._datafn(filterby, normalize, pkg, metric)
            
            if not os.path.isfile(pqf) or ow:
                print(f"\tmaking {pqf}")
                newcolnames = {name: f"{metric}_{name}" for name in rawdata.columns}
                df = rawdata.apply(self._networkx_analysis, args=(metric, normalize)).rename(columns=newcolnames)
                df.to_parquet(pqf)
            
        return   
    
    
    
    def _networkx_analysis(self, column, metric, normalize):
        weights = column.reset_index().rename(columns={f"{column.name}": "weight"}).dropna()
        network = networkx_from_pandas(weights, "level_0", "level_1", "weight")
        
        if callable(metric):
            metricf = metric
        else:        
            if "cfb" in metric:
                metricf = edge_current_flow_betweenness_centrality
            elif "btw" in metric:
                metricf = edge_betweenness_centrality

            if "subset" in metric:
                metricf = eval(f"{metricf.__name__}_subset")
        
        if "_subset" in metricf.__name__:
            nodes = {"sources": self._state.sources_subset, "targets": self._state.targets_subset}
        else:
            nodes = {}
        
        analyzed = metricf(network, normalized=normalize, weight="weight", **nodes)
        edges = {tuple(sorted(k, key = lambda x: int(x.split(":")[-1]))): analyzed[k] for k in analyzed}
        return pandas.Series(edges)
    
    
    

    def _add_anal(self, pkg, metrics, normalize, ow, filterby):
        norm = self._state._get_norm_str(normalize)
        print(f"adding analyzed {norm} data of {pkg} for {self._state.name}")
        
        data = None
        for metric in metrics:
            pqf = self._state._datafn(filterby, normalize, pkg, metric)
            df = pandas.read_parquet(pqf)

            cols = [f"{metric}_{num}" for num in self._state._trajs]
            df[f"{metric}_avg"] = df[cols].fillna(0).mean(axis=1)
            df[f"{metric}_std"] = df[cols].fillna(0).std(axis=1)
            out = df.drop(cols, axis=1)
            data = data.join(out) if data is not None else out
        
        exec(f"prev = self.{filterby}.{norm}.{pkg}.df")
        exec(f"self.{filterby}.{norm}.{pkg} = Edges(self._state, prev.join(data))")
        # setattr(self.incontact, pkg, Edges(self._state, prev.join(data)))
        
        return






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
            colors = {"inactive": "r", "active": "g",
                     "Gprotein": "y", "Barr": "b"}
            color2, color1 = colors[self._parent.state1.name], colors[self._parent.state2.name]
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
            cbar.ax.set_xticklabels([self._parent.state2.name.capitalize(), self._parent.state1.name.capitalize()])
        else:
            cbar = pl.colorbar(orientation = "horizontal", ticks = [])

        return
    
    
    
    def view(self, metric, num=20, filterby="incontact"):
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
