import os, time

import pandas
import numpy as np

from contextlib import redirect_stdout, redirect_stderr
from importlib import import_module
from lazyasd import LazyObject

from ..AlloViz.utils import get_pool, rgetattr, rhasattr




def lazy_import(key, val):
    extra_arg = {"package": 'AlloViz'} if 'Packages' in val else {}
    return LazyObject(lambda: import_module(val, **extra_arg), globals(), key)







class Base:
    def __new__(cls, protein, d):#**kwargs):
        #print("new", os.getpid(), protein, dir(protein), dir())
        new = super().__new__(cls)
        new._name = new.__class__.__name__
        new.protein = protein
        new._d = d
        
        #new._selection = d["_protein_sel"]
        new._pdbf = d["_pdbf"]
        new._trajs = d["_trajs"]
        
        new._path = f"{d['_datadir']}/{new._name}/raw"
        os.makedirs(new._path, exist_ok=True)
        new._rawpq = lambda xtc: f"{new._path}/{xtc if isinstance(xtc, int) else xtc.rsplit('.', 1)[0]}.pq"
        
        return new
    
    
    def __getnewargs__(self):
        return self.protein, self._d
        
        
    def __init__(self, *args):
        pqs = [self._rawpq(xtc) for xtc in self._trajs]
        no_exist = lambda pqs: [not os.path.isfile(pq) for pq in pqs]
        
        if any(no_exist(pqs)):
            for xtc in (xtc for xtc in self._trajs if no_exist(pqs)[xtc-1]):
                self._calculate(xtc)
                # time.sleep(5)
                
        
        def wait_calculate(pqs):
            while any(no_exist(pqs)):
                time.sleep(5)
            return pqs
                
        def get_raw(pqs):
            print(f"adding raw data of {self._name} for {self._pdbf}: ", pqs)
            flist = map(lambda pq: pandas.read_parquet(pq), pqs)
            df = pandas.concat(flist, axis=1)
            
            if len(self._trajs) > 1:
                cols = [f"{num}" for num in self._trajs]
                df["weight"] = df[cols].fillna(0).mean(axis=1)
                df["weight_std"] = df[cols].fillna(0).std(axis=1)
            else:
                df.rename(columns={"1": "weight"}, inplace=True)
                
            return df
        
        add_raw = lambda pqs: setattr(self, "raw", get_raw(pqs))
        get_pool().apply_async(wait_calculate,
                               args=(pqs,),
                               callback=add_raw)
        
        
    
    def _calculate(self, xtc, *args):
        def send_and_log(xtc, *args):
            #print(f"sending {xtc}", os.getpid())
            with open(f"{self._path}/{self._name}.log", "a+") as f:
                with redirect_stdout(f), redirect_stderr(f):
                    return self._computation(xtc, *args)
                
        get_pool().apply_async(send_and_log,
                         args=(xtc, *args),
                         callback=self._save_pq)
        
    # def _computation
    
    
    
    def _save_pq(self, args):
        corr, xtc, *resl = args
        
        if len(resl) != 0:
            resl = resl[0]
            if corr.shape != (len(resl), len(resl)):
                corr = corr[np.ix_(resl, resl)]
        elif len(resl) == 0:
            resl = slice(0, corr.shape[0])
            
        resnames = [f"{aa.resname}:{aa.resid}" for aa in self._d["protein"].residues[resl]]
        
        df = pandas.DataFrame(corr, columns=resnames, index=resnames)
        df = df.where( np.triu(np.ones(df.shape), k=1).astype(bool) )
        df = pandas.DataFrame({f"{xtc}": df.stack()})
        # if not len(df[f"{xtc}"].unique()) == 1:
        df.to_parquet(self._rawpq(xtc))

    
    
    


class Multicore(Base):
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        
        if "taskcpus" not in new._d:
            new.taskcpus = int(np.ceil(os.cpu_count()/2))
        else:
            new.taskcpus = new._d["taskcpus"]
        
        new._empties = new.taskcpus-1
            
        return new
    
    
    
    def _calculate(self, xtc, *args):
        super()._calculate(xtc, *args)
        
        
        def calculate_empty(pq):
            # print("sleeping", pq, os.getpid())
            while not os.path.isfile(pq):
                time.sleep(5)
            return
        
        for _ in range(self._empties): get_pool().apply_async(calculate_empty, args=(self._rawpq(xtc),))
    
    
    



        
class Combined_Dihs(Base):
    
    def _calculate(self, xtc):
        pkg = self._name.replace("Dihs", "")
        Dihl = ["Phi", "Psi", "Omega"]
        get_rawpq = lambda Dih: rgetattr(self, "protein", f"{pkg}{Dih}", "_rawpq")(xtc)
        no_exist = lambda Dihl: [not rhasattr(self, "protein", f"{pkg}{Dih}") for Dih in Dihl]
        
        if any(no_exist(Dihl)):
            for Dih in (Dih for Dih in Dihl if no_exist(Dihl)[Dihl.index(Dih)]):
                #pkgclass = eval(f"{pkg}{Dih}") if isinstance(pkg, str) else pkg
                pkgclass = eval(f"self._{Dih}")
                setattr(self.protein, pkgclass.__name__, pkgclass(self.protein, self._d))

        
        def wait_calculate(Dihl):
            not_finished = lambda Dihl: [not os.path.isfile(get_rawpq(Dih)) for Dih in Dihl]
            while any(not_finished(Dihl)):
                time.sleep(5)
            return Dihl
                
        def save_pq(Dihl):
            dfs = [pandas.read_parquet(get_rawpq(Dih))[f"{xtc}"].abs() for Dih in Dihl]
            
            final = None
            for df in dfs:
                if final is None:
                    final = df
                else:
                    final = final + df
            df = final / len(Dihl) # average of the absolute number

            # df = (final - final.min()) / (final.max() - final.min()) # This is done column-wise # This would be needed for absolute number sum; we are doing averaging
            pandas.DataFrame(df).to_parquet(self._rawpq(xtc))
            return
            
            
        get_pool().apply_async(wait_calculate,
                         args=(Dihl,),
                         callback=save_pq)
        
        




class Use_COM(Base):
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        
        new._pdbf = new._d["_compdbf"]
        new._trajs = new._d["_comtrajs"]
        
        return new
    
    

    
    
# class Use_dcd(Base):
#     def __new__(cls, protein, d):
#         new = super().__new__(cls, protein, d)
        
        
        
#         return new