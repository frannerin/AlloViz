from networkx import from_pandas_edgelist as networkx_from_pandas
from networkx.algorithms.centrality import *


class Data:
    def __init__(self, state):
        self._state = state
    
    
    
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
