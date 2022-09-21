"""Module with functions to analyze filtered networks

Main function :func:`~AlloViz.AlloViz.Analysis.analyze` manages the analysis of filtered
networks. It calls :func:`~AlloViz.AlloViz.Analysis.single_analysis` for the analysis of
single element-metric combinations and uses the NetworkX's functions defined in
:data:`~AlloViz.AlloViz.Analysis.nodes_dict` and
:data:`~AlloViz.AlloViz.Analysis.edges_dict`.

"""

import os
import sys
import time

import pandas
import numpy as np
import networkx
from concurrent.futures import ProcessPoolExecutor as Pool
from functools import partial

# from importlib import import_module

# from .Visualization import Edges, Nodes
from . import Elements #import Edges, Nodes
from . import utils
from .utils import rgetattr, rhasattr



nodes_dict = {
    "btw": "networkx.algorithms.centrality.betweenness_centrality",
    "cfb": "networkx.algorithms.centrality.current_flow_betweenness_centrality",
}
"""
Dictionary that maps nodes network metrics custom names (e.g., betweenness centrality,
"btw") with their corresponding NetworkX function (e.g., 
"networkx.algorithms.centrality.betweenness_centrality").
"""

edges_dict = {
    "btw": "networkx.algorithms.centrality.edge_betweenness_centrality",
    "cfb": "networkx.algorithms.centrality.edge_current_flow_betweenness_centrality",
}
"""
Dictionary that maps edges network metrics custom names (e.g., betweenness centrality,
"btw") with their corresponding NetworkX function (e.g., 
"networkx.algorithms.centrality.edge_betweenness_centrality").
"""



# def analyze(filtered, elements="edges", metrics="all", normalize=True, **kwargs):
#     r"""Analyze the filtered network

#     Send the analyses of the passed filtered network for the specified combinations of
#     elements-metrics. Each combination is analyzed independently with 
#     :func:`~AlloViz.AlloViz.Analysis.single_analysis` using NetworkX' functions and
#     results are stored as new instances of classes from the
#     :mod:`AlloViz.AlloViz.Elements` module, which extend the :class:`pandas.DataFrame`
#     class.
    
#     Parameters
#     ----------
#     filtered : :class:`~AlloViz.AlloViz.Filtering.Filtering` object
#         Filtered network object.
#     elements : str or list, {"edges", "nodes"}
#         Network element for which to perform the analysis.
#     metrics : str or list
#         Network metrics to compute, which must be keys in the `nodes_dict` or
#         `edges_dict` dictionaries. Default is "all" and it sends the computation for all
#         the metrics defined in the corresponding dictionary of the selected elements in
#         `element`.
#     normalize : bool
#         Passed to the NetworkX functions that calculate the metrics, to output
#         normalized results or not.

#     Other Parameters
#     ----------------
#     nodes_dict, edges_dict : dict
#         Optional kwarg(s) of the dictionary(ies) that maps network metrics custom names
#         (e.g., betweenness centrality, "btw") with their corresponding NetworkX
#         function (e.g., "networkx.algorithms.centrality.betweenness_centrality").
#         Functions strings must be written as if they were absolute imports, and must
#         return a dictionary of edges or nodes, depending on the element dictionary in
#         which they are. The keys of the dictionaries will be used to name the columns
#         of the analyzed data that the functions produce. Defaults are
#         :data:`~AlloViz.AlloViz.Analysis.nodes_dict` and
#         :data:`~AlloViz.AlloViz.Analysis.edges_dict`.
#     **kwargs
#         `GetContacts_threshold` kwarg can be passed to specify the minimum contact
#         frequency (0-1, default 0) threshold, which will be used to filter out
#         contacts with a frequency (average) lower than it before analysis. Make sure
#         to delete/have deleted all previous analysis attributes and files of any
#         network construction method. `Intercontact_distance` kwarg can be passed to
#         specify the minimum number of sequence positions/distance between residues of
#         a pair to retain in Intercontact filtering, which defaults to 5.
#     """
#     elements = elements if isinstance(elements, list) else [elements]

#     for elem in elements:
#         # If the Element's attribute doesn't exist yet, use as initial data the raw edge weights from the filtered data (or an empty DF for nodes)
#         if not rhasattr(filtered, elem):
#             if elem == "edges":
#                 cols = ["weight" in col for col in filtered._filtdata.columns]
#                 data = filtered._filtdata.loc[:, cols]
#             elif elem == "nodes":
#                 data = pandas.DataFrame()
#         # Else, retrieve the Element's attribute DataFrame to add columns to it
#         else:
#             data = rgetattr(filtered, elem)
            
#         # Retrieve the element's dictionary and if necessary update it with the one passed in kwargs
#         d = eval(f"{elem}_dict").copy()
#         if f"{elem}_dict" in kwargs:
#             d.update(kwargs[f"{elem}_dict"])
#         # Create a list of the desired metrics to calculate
#         metrics = utils.make_list(metrics, if_all = list(d.keys()))
#         # Create a list of the metrics that (i) have been passed, (ii) are also present in the element's dictionary and (iii) aren't already in the Element's attribute
#         elem_metrics = [
#             metric for metric in metrics if (metric in d and metric not in data)
#         ]
        
#         # Define the list of .pq files that we expect are going to be saved (or be retrieved) and a function to check which of them already exist
#         pqs = lambda elem: [filtered._datapq(elem, metric) for metric in elem_metrics]
#         no_exist = lambda pqs: [not os.path.isfile(pq) for pq in pqs]
#         print(utils.get_pool())
#         # If any of the .pq files don't exist, send the analysis calculations for them
#         if any(no_exist(pqs(elem))):
#             for metric in (
#                 metric
#                 for metric in elem_metrics
#                 if no_exist(pqs(elem))[pqs(elem).index(filtered._datapq(elem, metric))]
#             ):
#                 utils.get_pool().apply_async(
#                     single_analysis,
#                     args=(
#                         filtered,
#                         metric,
#                         d[metric],
#                         elem,
#                         normalize,
#                         filtered._datapq(elem, metric),
#                     ),
#                 )

#         # Function to wait for the analyses to finish in the background; returns the data to be added as attributes when they do
#         def wait_analyze(elem, data):
#             while any(no_exist(pqs(elem))):
#                 time.sleep(5)
#             return elem, data

#         # Function to add the newly calculated (or retrieved) data to the Element's attribute after analysis (wait_analyze) finishes and passes the data
#         def add_data(args):
#             elem, data = args
#             print(
#                 f"adding analyzed {elem} {filtered._pkg} {filtered._name} data of for {filtered._pkg.protein._pdbf}"
#             )

#             for pq in pqs(elem):
#                 # Retrieve the metric name from the .pq filename and read it into a DataFrame
#                 metric = pq.rsplit("/", 1)[-1].split("_", 1)[-1].split(".")[0]
#                 df = pandas.read_parquet(pq)

#                 # If there are more than 1 trajectory, calculate the metric average and standard error from the trajectories' analyzed data
#                 if len(filtered._pkg.protein._trajs) > 1:
#                     cols = [f"{metric}_{num}" for num in filtered._pkg.protein._trajs]
#                     df[f"{metric}"] = df[cols].fillna(0).mean(axis=1)
#                     df[f"{metric}_std"] = df[cols].fillna(0).std(axis=1)
#                     # After calculating the mean and std, drop the individual trajectories' columns, retaining the mean and std and also the analysis of the filtered raw weights averages
#                     out = df.drop(cols, axis=1)
#                 else:
#                     # If there is only 1 trajectory, only the filtered raw weights are analyzed, outputing a single analyzed column
#                     out = df

#                 data = pandas.concat([data, out], axis=1)

#             # Retrieve the Element's class from the Visualization module and (re-)set it with the data
#             elemclass = eval(elem.capitalize())
#             setattr(filtered, elem, elemclass(data))
#             getattr(filtered, elem)._parent = filtered._pkg.protein

#         # Wait asynchronously for analysis to end and then add the data
#         utils.get_pool().apply_async(
#             wait_analyze, args=(elem, data), callback=add_data
#         )

#     #return lambda: getattr(filtered, elements[0]) if len(elements) == 1 else None

# def single_analysis(filtered, metric, metricf, elem, normalize, pq):
#     r"""Analyze raw data with a single element-metric

#     Analyze stored, filtered raw data with the passed combination of
#     element-metric/NetworkX' analysis function and save the results.

#     Parameters
#     ----------
#     filtered : :class:`~AlloViz.AlloViz.Filtering.Filtering` object
#         Filtered network object.
#     metric : str
#         Network metric to compute, which must be a key in the `nodes_dict` or
#         `edges_dict` dictionaries.
#     metricf : str
#         NetworkX function to analyze data. It must be written as if it were an
#         absolute import
#         (e.g., "networkx.algorithms.centrality.betweenness_centrality"). 
#     element : str or list, {"edges", "nodes"}
#         Network element for which the analysis is performed.
#     normalize : bool
#         Passed to the NetworkX functions that calculate the metrics, to output
#         normalized results or not.
#     pq : str
#         Name of the parquet (.pq) file in which to save the analysis results.
#     """
#     # Process the NetworkX's module-function that will be used for analysis and import it into the metricf variable
#     # module, f = metric_import.rsplit(".", 1)
#     metricf = eval(metricf)  # eval(f"import_module('{module}').{f}")

#     # Create an empty DataFrame to store analysis results and the list of columnnames (Graph attributes) to be analyzed
#     df = pandas.DataFrame()
#     cols = [c for c in filtered._filtdata if "std" not in c]
#     nodes = {}  # Temporary fix for future use of source-sink network analyses

#     for col in cols:
#         # Try to apply the NetworkX's analysis function to the selected Graph
#         try:
#             analyzed = metricf(
#                 rgetattr(filtered, f"_{col}_G"),
#                 normalized=normalize,
#                 weight="weight",
#                 **nodes,
#             )
#         # If it throws an error, print it along with the performed analysis' information and create fake data with all 0s
#         except Exception as e:  # from scipy.linalg import LinAlgError
#             print(
#                 "ERROR:",
#                 filtered._pkg._name,
#                 filtered._name,
#                 elem,
#                 metricf.__name__,
#                 "\n",
#                 e,
#             )
#             analyzed = {k: 0 for k in rgetattr(filtered, f"_{col}_G", elem)}

#         # Sort the result's indices if they are a MultiIndex (edges), so that each pair has the residue with the lowest resnum as the first element of the tuple
#         sort_index = lambda result: {
#             tuple(sorted(k, key=lambda x: int(x.split(":")[-1]))): result[k]
#             for k in result
#         }
#         result = sort_index(analyzed) if elem == "edges" else analyzed

#         # Add the analyzed data to the df with the appropriate column name (e.g., if there is more than 1 trajectory, there will be "weight", "metric_weight", "metric_1", "metric_"...)
#         colname = f"{metric}_{col}" if len(cols) > 1 else metric
#         df = pandas.concat([df, pandas.Series(result, name=colname)], axis=1)

#     df.to_parquet(pq)

def analyze_graph(args):
    r"""Analyze a graph/column from raw filtered data with an element-metric

    Analyze a stored, filtered raw data column with the passed combination of
    element-metric/NetworkX' analysis function and return the results.

    Parameters
    ----------
    graph : :external:ref:`Graph <graph>` object
        Single column/graph object to analyze.
    metricf : str
        NetworkX function to analyze data. It must be written as if it were an
        absolute import
        (e.g., "networkx.algorithms.centrality.betweenness_centrality").
    normalize : bool
        Passed to the NetworkX functions that calculate the metrics, to output
        normalized results or not.
    colname : str
        Name of the analyzed column that it will have in the final DataFrame for saving.
    """
    graph, metricf, normalize, colname = args
    nodes = {} # Temporary fix for future use of source-sink network analyses

    # Try to apply the NetworkX's analysis function to the selected Graph
    try:
        analyzed = metricf(
            graph,
            normalized=normalize,
            weight="weight",
            **nodes,
        )
    # If it throws an error, return False to print it along with the performed analysis' information and create fake data with all 0s
    except Exception as e:
        return e

    # Sort the result's indices if they are a MultiIndex (edges), so that each pair has the residue with the lowest resnum as the first element of the tuple
    sort_index = lambda result: {
        tuple(sorted(k, key=lambda x: int(x.split(":")[-1]))): result[k]
        for k in result
    }
    result = sort_index(analyzed) if len(list(analyzed.keys())[0]) == 2 else analyzed

    # Return a series
    return pandas.Series(result, name=colname)




def single_analysis(graphs, metricf, metric, elem, normalize, pq):
    r"""Analyze raw data with a single element-metric

    Analyze stored, filtered raw data with the passed combination of
    element-metric/NetworkX' analysis function and save the results.

    Parameters
    ----------
    graphs : dict of external:ref:`Graph <graph>` objects
        Graphs to analyze.
    metricf : str
        NetworkX function to analyze data. It must be written as if it were an
        absolute import
        (e.g., "networkx.algorithms.centrality.betweenness_centrality"). 
    metric : str
        Network metric to compute, which must be a key in the `nodes_dict` or
        `edges_dict` dictionaries.
    elem : str or list, {"edges", "nodes"}
        Network element for which the analysis is performed.
    normalize : bool
        Passed to the NetworkX functions that calculate the metrics, to output
        normalized results or not.
    pq : str
        Name of the parquet (.pq) file in which to save the analysis results.
    """
    # Process the NetworkX's module-function that will be used for analysis and import it into the metricf variable
    metricf = eval(metricf)
    
    # Function to gget the columns' new names (e.g., if there is more than 1 trajectory, there will be "weight", "metric_weight", "metric_1", "metric_"...)
    get_colname = lambda metric, col: f"{metric}_{col}" if len(graphs) > 1 else metric
    # Analyze all columns in parallel, returning a Series for each (or False if it couldn't be analyzed)
    with Pool(len(graphs)) as p:
        args = [(graph, metricf, normalize, get_colname(metric, col)) for col, graph in graphs.items()]
        results = list(p.map(analyze_graph, args))
        p.shutdown()
    
    # Check if any of the analyses failed and print the information
    fails = [(args[i][-1], results.pop(i)) for i in range(len(results)) if not isinstance(results[i], pandas.Series)]
    if len(fails) > 0:
        for f in fails:
            print(
                "ERROR:",
                pq,
                f[0],
                elem,
                metricf.__name__,
                "\n",
                f[-1],
            )
    
    pandas.concat(results, axis=1).to_parquet(pq)
    
    

# Function to wait for the analyses to finish in the background; returns the data to be added as attributes when they do
def wait_analyze(pqs):
    while any([not os.path.isfile(pq) for pq in pqs]):
        time.sleep(5)
    return pqs

# Function to add the newly calculated (or retrieved) data to the Element's attribute after analysis (wait_analyze) finishes and passes the data
def add_data(pqs, elem, data, filtered):
    print(
        f"adding analyzed {elem} {filtered._pkg} {filtered._name} data of for {filtered._pkg.protein._pdbf}"
    )
    sys.stdout.flush()

    for pq in pqs:
        # Retrieve the metric name from the .pq filename and read it into a DataFrame
        metric = pq.rsplit("/", 1)[-1].split("_", 1)[-1].split(".")[0]
        df = pandas.read_parquet(pq)

        # If there are more than 1 trajectory, calculate the metric average and standard error from the trajectories' analyzed data
        if len(filtered._pkg.protein._trajs) > 1:
            cols = [f"{metric}_{num}" for num in filtered._pkg.protein._trajs]
            df[f"{metric}"] = df[cols].fillna(0).mean(axis=1)
            df[f"{metric}_std"] = df[cols].fillna(0).std(axis=1)
            # After calculating the mean and std, drop the individual trajectories' columns, retaining the mean and std and also the analysis of the filtered raw weights averages
            out = df.drop(cols, axis=1)
        else:
            # If there is only 1 trajectory, only the filtered raw weights are analyzed, outputing a single analyzed column
            out = df

        data = pandas.concat([data, out], axis=1)

    # Retrieve the Element's class from the Visualization module and (re-)set it with the data
    elemclass = eval(f"Elements.{elem.capitalize()}")
    setattr(filtered, elem, elemclass(data))
    getattr(filtered, elem)._parent = filtered._pkg.protein
    
    


def analyze(filtered, elements, metrics, normalize, **kwargs):
    r"""Analyze the filtered network

    Send the analyses of the passed filtered network for the specified combinations of
    elements-metrics. Each combination is analyzed independently with 
    :func:`~AlloViz.AlloViz.Analysis.single_analysis` using NetworkX' functions and
    results are stored as new instances of classes from the
    :mod:`AlloViz.AlloViz.Elements` module, which extend the :class:`pandas.DataFrame`
    class.
    
    Parameters
    ----------
    filtered : :class:`~AlloViz.AlloViz.Filtering.Filtering` object
        Filtered network object.
    elements : str or list, {"edges", "nodes"}
        Network element for which to perform the analysis.
    metrics : str or list
        Network metrics to compute, which must be keys in the `nodes_dict` or
        `edges_dict` dictionaries.
    normalize : bool
        Passed to the NetworkX functions that calculate the metrics, to output
        normalized results or not.

    Other Parameters
    ----------------
    nodes_dict, edges_dict : dict
        Optional kwarg(s) of the dictionary(ies) that maps network metrics custom names
        (e.g., betweenness centrality, "btw") with their corresponding NetworkX
        function (e.g., "networkx.algorithms.centrality.betweenness_centrality").
        Functions strings must be written as if they were absolute imports, and must
        return a dictionary of edges or nodes, depending on the element dictionary in
        which they are. The keys of the dictionaries will be used to name the columns
        of the analyzed data that the functions produce. Defaults are
        :data:`~AlloViz.AlloViz.Analysis.nodes_dict` and
        :data:`~AlloViz.AlloViz.Analysis.edges_dict`.
    """
    elements = elements if isinstance(elements, list) else [elements]

    for elem in elements:
        # If the Element's attribute doesn't exist yet, use as initial data the raw edge weights from the filtered data (or an empty DF for nodes)
        if not rhasattr(filtered, elem):
            if elem == "edges":
                cols = ["weight" in col for col in filtered._filtdata.columns]
                data = filtered._filtdata.loc[:, cols]
            elif elem == "nodes":
                data = pandas.DataFrame()
        # Else, retrieve the Element's attribute DataFrame to add columns to it
        else:
            data = rgetattr(filtered, elem)
            
        # Retrieve the element's dictionary and if necessary update it with the one passed in kwargs
        d = eval(f"{elem}_dict").copy()
        if f"{elem}_dict" in kwargs:
            d.update(kwargs[f"{elem}_dict"])
        # Create a list of the desired metrics to calculate
        metrics = utils.make_list(metrics, if_all = list(d.keys()))
        # Create a list of the metrics that (i) have been passed, (ii) are also present in the element's dictionary and (iii) aren't already in the Element's attribute
        elem_metrics = [
            metric for metric in metrics if (metric in d and metric not in data)
        ]
        
        # Define the list of .pq files that we expect are going to be saved (or be retrieved) and a function to check which of them already exist
        pqs = lambda elem: [filtered._datapq(elem, metric) for metric in elem_metrics]
        no_exist = lambda pqs: [not os.path.isfile(pq) for pq in pqs]
        # If any of the .pq files don't exist, send the analysis calculations for them
        if any(no_exist(pqs(elem))):
            for metric in (
                metric
                for metric in elem_metrics
                if no_exist(pqs(elem))[pqs(elem).index(filtered._datapq(elem, metric))]
            ):
                utils.get_pool().apply_async(
                    single_analysis,
                    args=(
                        filtered.graphs,
                        d[metric],
                        metric,
                        elem,
                        normalize,
                        filtered._datapq(elem, metric),
                    ),
                )
                
        # Wait asynchronously for analysis to end and then add the data
        utils.get_pool().apply_async(
            wait_analyze, args=(pqs(elem),), callback=partial(add_data, elem=elem, data=data, filtered=filtered)
        )

    #return lambda: getattr(filtered, elements[0]) if len(elements) == 1 else None


    
    
    

    
    
    
    
