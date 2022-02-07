class Pkg: #abc.ABC
    def __init__(self, state):#, ¿pkg?): #metrics="all", filterby="whole", normalize=True, cores=None, ow=False
        args = locals()
        del args["self"]
        if "args" in args: del args["args"]
        self.__dict__.update(args)
        # self.pool = get_pool()
        
        self._name = self.__class__.__name__
        self._path = lambda filterby: f"{self.state.name}/data/{self._name}/{filterby}"
        self._rawpq = lambda xtc: f"{self._path('raw')}/{xtc if isinstance(xtc, int) else xtc.rsplit('.', 1)[0]}.pq"
        
        self.raw = self._initialize()
        
        
#     def _get_norm_str(self, normalize):
#         return "norm" if normalize else "no_norm"
    
#     def _datapn(self, filterby, normalize, pkg):
#         norm = self._get_norm_str(normalize)
#         return f"{self.name}/data/{filterby}/{norm}/{pkg}"
    
#     def _datafn(self, filterby, normalize, pkg, metric):
#         return f"{self._datapn(filterby, normalize, pkg)}/{metric}.pq"
        
        
        
        
    
    def _initialize(self):        
        pqs = [self._rawpq(xtc) for xtc in self.state._trajs]
        no_exist = [not os.path.isfile(pq) for pq in pqs]
        
        if any(no_exist):
            os.makedirs(self._path("raw"), exist_ok=True)
            for xtc in (xtc for xtc in self.state._trajs if no_exist[xtc-1]):
                self._calculate(xtc)
        
        
        
        def wait_calculate(pqs):
            while any([not os.path.isfile(pq) for pq in pqs]):
                time.sleep(5)
                
        def get_raw(*args):
            print(f"adding raw data of {self._name} for {self.state.name}")
            flist = map(lambda pq: pandas.read_parquet(pq), pqs)
            df = pandas.concat(flist, axis=1)
            cols = [f"{num}" for num in self.state._trajs]
            df["weight_avg"] = df[cols].fillna(0).mean(axis=1)
            df["weight_std"] = df[cols].fillna(0).std(axis=1)
            return df
            
        add_raw = lambda _: setattr(self, "raw", get_raw())
        get_pool().apply_async(wait_calculate,
                               args=(pqs,),
                               callback=add_raw)
        return
        

        
    
    def _calculate(self, xtc):
        pool = get_pool()
        pdb = self.state._pdbf
        traj = f"{self.state.name}/{self.state._trajs[xtc]}"
        pq = self._rawpq(xtc)
        print(f"making {pq}")
        return pool, pdb, traj, pq






class Multicorepkg(Pkg):
    def __init__(self, state, taskcpus = 3):
        self.taskcpus = taskcpus
        super().__init__(state)
    
    def _calculate_empty(self, pqf):
        print("sleeping", pqf, os.getpid())
        while not os.path.isfile(pqf):
            time.sleep(5)
        return






class Correlationpkg(Pkg):
    def __init__(self, state):
        self.selection = "protein"
        super().__init__(state)
    
    
    def _save_pq(self, args):
        corr, xtc, pq = args
        
        resl = [f"A:{aa.resname}:{aa.resid}" for aa in self.state.mdau.select_atoms(self.selection).residues]

        df = pandas.DataFrame(corr, columns=resl, index=resl)
        df = df.where( np.triu(np.ones(df.shape), k=1).astype(bool) )
        df = pandas.DataFrame({f"{xtc}": df.stack()})
        df.to_parquet(pq)






class COMpkg(Pkg):
    def __init__(self, state):
        if not hasattr(state, "_comtrajs"):
            self.state._add_comtrajs()
            
        super().__init__(state)
    
    
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        pdb = self.state._compdbf
        traj = self.state._comtrajs[xtc]
        
        return pool, pdb, traj, pq
