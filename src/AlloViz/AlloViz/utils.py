pkgsl = ["MDTASK", "Getcontacts", "PyInteraph", "PyInteraphEne", "Dynetan", "DynetanCOM", "PytrajCA", "PytrajCB",
         "Corrplus", "CorrplusLMI", "CorrplusCOM", "CorrplusCOMLMI", "CorrplusPsi", "CorrplusPhi", "CorrplusOmega", "CorrplusDihs",
        "gRINN", "gRINNcorr", "G_corrCAMI", "G_corrCOMMI", "g_corrCALMI", "g_corrCOMLMI",
        "AlloVizPsi", "AlloVizPhi", "AlloVizOmega", "AlloVizDihs",
        "MDEntropyContacts", "MDEntropyDihs", "MDEntropyAlphaAngle"]

metricsl = ["cfb", "btw"]

filterbyl = ["whole", "incontact", "intercontact"]




def rgetattr(obj, *attrs):
    attl = list(attrs)
    while len(attl) >= 1:
        attr = attl.pop(0)
        try:
            obj = getattr(obj, attr)
        except AttributeError:
            return False
    return obj

def rhasattr(obj, *attrs):
    result = rgetattr(obj, *attrs)
    return True if not isinstance(result, bool) else False



class dummypool:
    def apply_async(self, function, args=[], callback=None):
        if callback is not None:
            return callback(function(*args))
        else:
            return function(*args)


pool = dummypool()
def get_pool():
    global pool
    # print(pool)
    return pool





def capitalize(string):
    return string[0].upper() + string[1:]

# def norm(normalize):
#     return "norm" if normalize else "no_norm"






from importlib import import_module
from lazyasd import LazyObject

def lazy_import(key, val):
    extra_arg = "package='AlloViz'" if 'Packages' in val else ''
    return LazyObject(lambda: import_module(val, extra_arg(val)), globals(), key)






# def get_intercontacts(indexl):
#     get_resnum = lambda res: int(res.rsplit(":")[-1])
#     return [idx for idx in indexl if abs(get_resnum(idx[0]) - get_resnum(idx[1]) ) >= 4]






# pdict = {}

# def update_pdict(name, p):
#     global pdict
#     pdict[name] = p
#     print("adding to pdict", pdict)
#     p.start()
    
# def get_p(name):
#     global pdict
#     print("getting from pdict", pdict)
#     p = pdict[name]
#     print(p)
#     p.get()
#     print(p)