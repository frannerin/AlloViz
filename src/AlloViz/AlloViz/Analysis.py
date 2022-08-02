import os
import time

import pandas
import numpy as np
import networkx as nx
from importlib import import_module

from .Visualization import Edges, Nodes
from .utils import rgetattr, rhasattr
from . import utils



class Analysis:
    def __init__(self, pkg, **kwargs):
        self._pkg = pkg
        self._name = self.__class__.__name__
        
        # Establish the paths and filenames
        self._path = f"{self._pkg.protein._datadir}/{self._pkg._name}/{self._name}"
        os.makedirs(self._path, exist_ok=True)
        self._datapq = lambda element, metric: f"{self._path}/{element}_{metric}.pq"
        
        # Get the filtered thata according to the filterby/Filtering scheme
        self._filtdata = self._get_filt_data()
        
        # For each column in the filtered data that is not a standard error, create an analyzable NetworkX's graph and save it as an attribute
        for col in [c for c in self._filtdata if "std" not in c]:
            setattr(self, f"_{col}_G", self._get_G(self._filtdata[col]))
            
    
    
    
    def _get_filt_data(self, **kwargs):
        return self._pkg.raw
    
    
    def _get_G(self, column):
        # Transform the column into a 3-column dataframe with the nodes' name and the value as weight. Drop 0s, NAs and use absolute value
        weights = column[column != 0].dropna().abs().rename("weight").reset_index() # btw calculations fail with 0 value weights and cfb prob with negative
        # Create the NetworkX's Graph
        network = nx.from_pandas_edgelist(weights, "level_0", "level_1", "weight")
        
        # Check that the largest connected component has the same size as the total number of nodes, else select the largest component to return as the Graph to be analyzed
        largest_component = max(nx.connected_components(network), key=len)

        if len(largest_component) < network.number_of_nodes():
            print(f"WARNING! Unconnected network ({network.number_of_nodes()} nodes):", self._pkg._name, self._name, column.name, "\n",
                f"Largest network component will be used for analysis. Sizes (number of nodes) of all components: {[len(comp) for comp in nx.connected_components(network)]}")
            network = network.subgraph(largest_component)
        
        return network
    
        
    
    def add_metrics(self, element, metrics, normalize, nodes_dict, edges_dict):
        elements = element if isinstance(element, list) else [element]
        
        for elem in elements:            
            # If the Element's attribute doesn't exist yet, use as initial data the raw edge weights from the filtered data (or an empty DF for nodes)
            if not rhasattr(self, elem, "df"):
                if elem == "edges":
                    cols = ["weight" in col for col in self._filtdata.columns]
                    data = self._filtdata.loc[:,cols]
                elif elem == "nodes":
                    data = pandas.DataFrame()
            # Else, retrieve the Element's attribute DataFrame to add columns to it
            else:
                data = rgetattr(self, elem, "df")
                
            # Retrieve the element's dictionary
            elem_dict = eval(f"{elem}_dict")
            # Create a list of the metrics that (i) have been passed, (ii) are also present in the element's dictionary and (iii) aren't already in the Element's attribute
            elem_metrics = [metric for metric in metrics if (metric in elem_dict and metric not in data)]
            
            # Define the list of .pq files that we expect are going to be saved (or be retrieved) and a function to check which of them already exist
            pqs = lambda elem: [self._datapq(elem, metric) for metric in elem_metrics]
            no_exist = lambda pqs: [not os.path.isfile(pq) for pq in pqs]
            
            # If any of the .pq files don't exist, send the analysis calculations for them
            if any(no_exist(pqs(elem))):
                for metric in (metric for metric in elem_metrics if no_exist(pqs(elem))[pqs(elem).index(self._datapq(elem, metric))]):
                    utils.get_pool().apply_async(self._analyze, args = (metric, elem_dict[metric], elem, normalize, self._datapq(elem, metric)))

                    
            # Function to wait for the analyses to finish in the background; returns the data to be added as attributes when they do
            def wait_analyze(elem, data):
                while any(no_exist(pqs(elem))):
                    time.sleep(5)
                return elem, data
                
            # Function to add the newly calculated (or retrieved) data to the Element's attribute after analysis (wait_analyze) finishes and passes the data
            def add_data(args):
                elem, data = args
                print(f"adding analyzed {elem} {self._pkg} {self._name} data of for {self._pkg.protein._pdbf}")

                for pq in pqs(elem):
                    # Retrieve the metric name from the .pq filename and read it into a DataFrame
                    metric = pq.rsplit("/", 1)[-1].split("_", 1)[-1].split(".")[0]
                    df = pandas.read_parquet(pq)
                    
                    # If there are more than 1 trajectory, calculate the metric average and standard error from the trajectories' analyzed data
                    if len(self._pkg.protein._trajs) > 1:
                        cols = [f"{metric}_{num}" for num in self._pkg.protein._trajs]
                        df[f"{metric}"] = df[cols].fillna(0).mean(axis=1)
                        df[f"{metric}_std"] = df[cols].fillna(0).std(axis=1)
                        # After calculating the mean and std, drop the individual trajectories' columns, retaining the mean and std and also the analysis of the filtered raw weights averages
                        out = df.drop(cols, axis=1)
                    else:
                    # If there is only 1 trajectory, only the filtered raw weights are analyzed, outputing a single analyzed column
                        out = df
                        
                    data = pandas.concat([data, out], axis=1)
                
                # Retrieve the Element's class from the Visualization module and (re-)set it with the data
                elemclass = eval(elem.capitalize())
                setattr(self, elem, elemclass(self._pkg.protein, data))

            # Wait asynchronously for analysis to end and then add the data
            utils.get_pool().apply_async(wait_analyze,
                                         args = (elem, data),
                                         callback = add_data)
        return
    

    
    
    def _analyze(self, metric, metric_import, elem, normalize, pq):
        # Process the NetworkX's module-function that will be used for analysis and import it into the metricf variable
        module, f = metric_import.rsplit(".", 1)
        metricf = eval(f"import_module('{module}').{f}")
        
        # Create an empty DataFrame to store analysis results and the list of columnnames (Graph attributes) to be analyzed
        df = pandas.DataFrame()
        cols = [c for c in self._filtdata if "std" not in c]
        nodes = {} # Temporary fix for future use of source-sink network analyses
        
        for col in cols:
            # Try to apply the NetworkX's analysis function to the selected Graph
            try:
                analyzed = metricf(rgetattr(self, f"_{col}_G"), normalized=normalize, weight="weight", **nodes)
            # If it throws an error, print it along with the performed analysis' information and create fake data with all 0s
            except Exception as e: # from scipy.linalg import LinAlgError
                print("ERROR:", self._pkg._name, self._name, elem, metricf.__name__, "\n", e)
                analyzed = {k: 0 for k in rgetattr(self, f"_{col}_G", elem)}

            # Sort the result's indices if they are a MultiIndex (edges), so that each pair has the residue with the lowest resnum as the first element of the tuple
            sort_index = lambda result: {tuple(sorted(k, key = lambda x: int(x.split(":")[-1]))): result[k] for k in result}
            result = sort_index(analyzed) if elem == "edges" else analyzed
            
            # Add the analyzed data to the df with the appropriate column name (e.g., if there is more than 1 trajectory, there will be "weight", "metric_weight", "metric_1", "metric_"...)
            colname = f"{metric}_{col}" if len(cols) > 1 else metric 
            df = pandas.concat([df, pandas.Series(result, name = colname)], axis=1)
        
        df.to_parquet(pq)
        
        
    
    
    
    
class Whole(Analysis):
    pass    
    
    
    
    
class Incontact(Analysis):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    
    def _get_filt_data(self, **kwargs):
        df = super()._get_filt_data()
        
        # Get GetContacts raw data if its calculation is available, else raise an Exception
        try:
            gc = self._pkg.protein.GetContacts.raw
            # If GetContacts_threshold kwarg is passed, use it to filter the present data without affecting the saved GetContacts' Protein attribute
            if "GetContacts_threshold" in kwargs:
                gc = self._pkg.protein.GetContacts._filter_raw(df, kwargs["GetContacts_threshold"])
        except:
            raise Exception("GetContacts results are needed first")
        
        # Return the data filtered retaining only the indices that are in GetContacts data
        return df.filter(gc.index, axis=0)
        
        
        
        
class Intercontact(Incontact):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        
    def _get_filt_data(self, **kwargs):
        # Get GetContacts-filtered/Incontact data to further filter it with the Intercontact_distance restriction
        df = super()._get_filt_data()
        
        
        def get_intercontacts(indexl, **kwargs):
            # Define Intercontact_distance if it is in kwargs or with the default value
            Intercontact_distance = kwargs["Intercontact_distance"] if "Intercontact_distance" in kwargs else 5
            # Define a function to retrieve the residue number from the residue nomenclature used (chainID:)RES:resnum
            resnum = lambda res: int(res.rsplit(":")[-1])
            # Return the indices whose residue pair satisfies the Intercontact_distance
            # The absolute value of the subtraction between the two residue numbers must be greater than the parameter to indicate that they are at a greater distance in the sequence than it
            return [idx for idx in indexl if abs(resnum(idx[0]) - resnum(idx[1]) ) > Intercontact_distance]
        
        # Return the data filtered according to Incontact + get_intercontacts function
        return df.filter(get_intercontacts(df.index, **kwargs), axis=0)