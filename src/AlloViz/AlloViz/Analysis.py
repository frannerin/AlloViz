import os
import time

import pandas
import numpy as np
import networkx as nx
from importlib import import_module

from .Visualization import Edges, Nodes
from .utils import rgetattr, rhasattr, capitalize
from . import utils



class Analysis: #(_Edges)
    def __init__(self, pkg):
        self._pkg = pkg
        # self._parent = self.pkg.state
        # self.metrics = metrics
        # self.normalize = normalize
        self._name = self.__class__.__name__
        
        self._path = f"{self._pkg.protein._datadir}/{self._pkg._name}/{self._name}" # lambda norm might not be needed
        os.makedirs(self._path, exist_ok=True)
        self._datapq = lambda element, metric: f"{self._path}/{element}_{metric}.pq"
        
        self._filtdata = self._get_filt_data()
        # self._graph = nx_from_pandas(df=self._filtdata.reset_index(), 
        #                                    source="level_0", target="level_1", 
        #                                    edge_attr=list(self._filtdata.drop("weight_std", axis=1).columns))
        
        # self.add_metrics(metrics, element, normalize)
    
    
    
    def _get_filt_data(self):
        return self._pkg.raw#[["weight_avg", "weight_std"]]
    
    
        
    
    def _add_metrics(self, element, metrics, normalize, nodes_dict, edges_dict):
        # metrics = metricsl if metrics=="all" else metrics if isinstance(metrics, list) else [metrics]
        elements = element if isinstance(element, list) else [element]
        
        for elem in elements:            
            if not rhasattr(self, elem, "df"):
                if elem == "edges":
                    cols = ["weight" in col for col in self._filtdata.columns]
                    data = self._filtdata.loc[:,cols]
                elif elem == "nodes":
                    data = pandas.DataFrame()
            else:
                data = rgetattr(self, elem, "df")
            
            elem_dict = eval(f"{elem}_dict")
            elem_metrics = [metric for metric in metrics if metric in elem_dict]
            pqs = lambda elem: [self._datapq(elem, metric) for metric in elem_metrics]
            no_exist = lambda pqs: [not os.path.isfile(pq) for pq in pqs]
            # is_in_df = any([f"{metric}_avg" not in data.columns for metric in metrics])
            
            # if not is_in_df and 
            if any(no_exist(pqs(elem))): # or ow
                for metric in (metric for metric in elem_metrics if no_exist(pqs(elem))[pqs(elem).index(self._datapq(elem, metric))]):
                    self._analyze(metric, elem_dict[metric], elem, normalize, self._datapq(elem, metric))


            def wait_analyze(elem, data):
                while any(no_exist(pqs(elem))):
                    time.sleep(5)
                return elem, data
                
            
            def add_data(args):
                elem, data = args
                print(f"adding analyzed {elem} {self._pkg} {self._name} data of for {self._pkg.protein._pdbf}")

                for pq in pqs(elem):
                    metric = pq.rsplit("/", 1)[-1].split("_", 1)[-1].split(".")[0]
                    df = pandas.read_parquet(pq)
                    
                    if len(self._pkg.protein._trajs) > 1:
                        cols = [f"{metric}_{num}" for num in self._pkg.protein._trajs]
                        df[f"{metric}"] = df[cols].fillna(0).mean(axis=1)
                        df[f"{metric}_std"] = df[cols].fillna(0).std(axis=1)
                        out = df.drop(cols, axis=1)
                    else:
                        out = df
                        
                    data = pandas.concat([data, out], axis=1)#data.join(out) if data is not None else out

                elemclass = eval(elem.capitalize())
                setattr(self, elem, elemclass(self._pkg.protein, data))

                
            utils.get_pool().apply_async(wait_analyze,
                                   args=(elem, data),
                                   callback=add_data)
        return
    

    
    
    def _analyze(self, metric, metric_import, elem, normalize, pq):
        pool = utils.get_pool()
        cols = ["std" not in col for col in self._filtdata.columns]
        rawdata = self._filtdata.loc[:,cols]
        nodes = {}
        
        module, f = metric_import.rsplit(".", 1)
        metricf = eval(f"import_module('{module}').{f}")
        
        
#         if callable(metric):
#             metricf = metric
#         else:        
#             if "cfb" in metric:
#                 metricf = current_flow_betweenness_centrality
#             elif "btw" in metric:
#                 metricf = betweenness_centrality
            
#             if elem == "edges":
#                 metricf = eval(f"edge_{metricf.__name__}")
            
#             if "subset" in metric:
#                 metricf = eval(f"{metricf.__name__}_subset")
#                 nodes = {"sources": self._pkg.protein.sources_subset, "targets": self._pkg.protein.targets_subset}  
                
                
                
        def save_pq(df):
            if len(df.columns) > 1:
                newcolnames = {name: f"{metric}_{name}" for name in df.columns}
            else:
                newcolnames = {df.columns[0]: metric}
                
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
        # network = nx_from_pandas(weights, "level_0", "level_1", "weight")
        network = nx.from_pandas_edgelist(weights, "level_0", "level_1", "weight")
        largest_component = max(nx.connected_components(network), key=len)
        
        if len(largest_component) < network.number_of_nodes():
            print(f"WARNING! Unconnected network ({network.number_of_nodes()} nodes):", self._pkg._name, self._name, column.name, elem, metricf.__name__, "\n",
            	f"Largest network component will be used for analysis. Sizes (number of nodes) of all components: {[len(comp) for comp in nx.connected_components(network)]}")
            network = network.subgraph(largest_component)
        
        try:
            analyzed = metricf(network, normalized=normalize, weight="weight", **nodes)
            # print(column.name, pandas.Series(analyzed).max())
        # from scipy.linalg import LinAlgError
        except Exception as e: # LinAlgError
            print("ERROR:", self._pkg._name, self._name, elem, metricf.__name__, "\n", e)
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
            indices = self._pkg.protein.GetContacts.raw.index
        except:
            raise Exception("GetContacts results are needed first")
#         if not rhasattr(self._pkg.state, "GetContacts", "raw"):
#             print("GetContacts results are needed; sending calculation first...")
            
#             pool = utils.get_pool()
            
#             self._pkg.state._set_pkgclass(self._pkg.state, "GetContacts", taskcpus = int(np.ceil(pool._processes/2)))
            
            
#             gc = self._pkg.state.GetContacts
#             pqs = [gc._rawpq(xtc) for xtc in gc.state._trajs]
#             no_exist = lambda pqs: [not os.path.isfile(pq) for pq in pqs]

#             while any(no_exist(pqs)):
#                 print("Waiting for pq files to be created")
#                 time.sleep(5)
                
#             while not rhasattr(gc, "raw"):
#                 print("Waiting for raw data to be added to object")
#                 time.sleep(5)
                
                
#         indexes = self._pkg.state.GetContacts.raw.index
            
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