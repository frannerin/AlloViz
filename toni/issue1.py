import AlloViz

pdbfile = "117/11159_dyn_117.pdb"
dcdfile = "117/11156_trj_117.dcd"
cache_path="deleteme"

prot = AlloViz.Protein(pdb=pdbfile, trajs=dcdfile, path=cache_path)

method="pytraj_CA"
prot.calculate(method)

flist=[]
fargs={}
prot.filter(method, filterings=[flist], **fargs)  # "all"?

el="edges"
met="btw"
prot.analyze(method, elements=el, metrics=met) # "all"?
