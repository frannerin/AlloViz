import AlloViz

pdbfile = "117/11159_dyn_117.pdb"
dcdfile = "117/11156_trj_117.dcd"
cache_path="deleteme"

if __name__ == "__main__":
    prot = AlloViz.Protein(pdb=pdbfile, trajs=dcdfile, path=cache_path)

    method="pytraj_CA"
    prot.calculate(method)

    flist=["All"]
    fargs={}
    prot.filter(method, filterings=[flist], **fargs)   

    el="edges"
    met="btw"
    prot.analyze(method, elements=el, metrics=met) 

    # result in prot.pytraj_CA.All.edges.btw
    