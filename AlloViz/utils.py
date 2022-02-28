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
    def apply_async(self, function, args, callback=None):
        if callback is not None:
            return callback(function(*args))
        else:
            return function(*args)


pool = dummypool()
def get_pool():
    global pool
    # print(pool)
    return pool








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