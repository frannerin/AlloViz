import AlloViz

pdbfile = "117/11159_dyn_117.pdb"
dcdfile = "117/11156_trj_117_r.dcd"
cache_path="deleteme"

if __name__ == "__main__":
    prot = AlloViz.Protein(pdb=pdbfile, trajs=dcdfile, path=cache_path)

    method="pytraj_CA"
    prot.calculate(method)
    calc_result = getattr(prot, method) 

    flist=["All"]
    fargs={}
    # flist = ['No_Sequence_Neighbors', 'Spatially_distant']
    # fargs = {'Sequence_Neighbor_distance': 5, 'Interresidue_distance': 10.0}
    prot.filter(method, filterings=[flist], **fargs)   
    flist_as_string = "_".join(flist) # :( (((
    filter_result = getattr(calc_result, flist_as_string)

    el="edges"
    # el = "nodes"

    # met="btw"
    met = "raw"

    prot.analyze(method, elements=el, metrics=met) 

    if met == "raw":
        analysis_result = filter_result._filtdata
    else:
        _ = getattr(filter_result, el)
        analysis_result = getattr(_, met)


    # result in prot.pytraj_CA.All.edges.btw
