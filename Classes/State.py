class State:    
    def __init__(self, name, idx):
        self.name = name
        self.idx = idx
        self.data = Data(self)
        
        if not os.path.isdir(self.name): self._download_files()
        
        self._pdbf = f"{self.name}/{self.name}_1.pdb"
        self._psff = f"{self.name}/{self.name}_1.psf"
        self._trajs = dict(enumerate(sorted( traj for traj in os.listdir(self.name) if re.search("^(?!\.).*\.xtc$", traj) ), 1))
        
        self.mdau = self._get_mdau()
        
    
    
    def __sub__(self, other):        
        return self.data - other.data
    
    
    def _rawpn(self, pkg):
        return f"{self.name}/data/raw/{pkg}"
    
    def _rawfn(self, pkg, xtc):
        filen = xtc if isinstance(xtc, int) else xtc.rsplit('.', 1)[0]
        return f"{self._rawpn(pkg)}/{filen}.pq"
    
    
    
    def _get_norm_str(self, normalize):
        return "norm" if normalize else "no_norm"
    
    def _datapn(self, filterby, normalize, pkg):
        norm = self._get_norm_str(normalize)
        return f"{self.name}/data/{filterby}/{norm}/{pkg}"
    
    def _datafn(self, filterby, normalize, pkg, metric):
        return f"{self._datapn(filterby, normalize, pkg)}/{metric}.pq"
        
    
    
    
    
    def _download_files(self):
        from bs4 import BeautifulSoup
        import urllib.request as pwget
        import requests


        web = "https://submission.gpcrmd.org"

        os.makedirs(self.name, exist_ok=True)

        html = requests.get(f"{web}/dynadb/dynamics/id/{self.idx}/")
        soup = BeautifulSoup(html.text, features="html.parser").find_all('a')
        links = (link.get('href') for link in soup if re.search("(xtc|pdb|psf)", link.get('href')))
        
        for link in links:
            ext = f".{link.rsplit('.', 1)[-1]}"
            num = len([file for file in os.listdir(self.name) if ext in file])+1
            fname = f"{self.name}/{self.name}_{num}{ext}"
            print(f"downloading {fname}")
            pwget.urlretrieve(f"{web}{link}", fname)
            
        return
    
    

    def _get_mdau(self):
        trajs = [f"{self.name}/{self._trajs[xtc]}" for xtc in self._trajs]
        mdau = mda.Universe(self._pdbf, *trajs)
        prot = mdau.select_atoms("protein")
        
        prot_numsf = f"{self.name}/{self.name}_nums.pdb"
        
        if not os.path.isfile(prot_numsf):
            print(f"retrieving {prot_numsf}")
            with io.StringIO() as protf:
                with mda.lib.util.NamedStream(protf, "prot.pdb") as f:
                    prot.write(f, file_format="pdb")
                response = requests.post('https://gpcrdb.org/services/structure/assign_generic_numbers', 
                                         files = {'pdb_file': protf.getvalue()})
                with open(prot_numsf, "w") as prot_nums:
                    prot_nums.write(response.text)
        
        nums = mda.Universe(prot_numsf).select_atoms("protein").tempfactors
        prot.tempfactors = nums.round(2)
                    
        return mdau
    
    
    def _add_comtrajs(self):
        compath = self._rawpn("COMtrajs")
        if not os.path.isdir(compath): os.makedirs(compath, exist_ok=True)

        compdb = f"{compath}/ca.pdb"
        if not os.path.isfile(compdb):
            prot = self.mdau.select_atoms("protein and name CA")
            prot.write(compdb)
        setattr(self, "_compdbf", compdb)

        comtrajs = [f"{compath}/{xtc}.xtc" for xtc in self._trajs]
        setattr(self, "_comtrajs", {num: traj for num, traj in enumerate(comtrajs, 1)})

        for xtc, comtraj in self._comtrajs.items():
            if not os.path.isfile(comtraj):
                prot = self.mdau.select_atoms("protein")
                traj = next(traj for traj in self.mdau.trajectory.readers if traj.filename == f"{self.name}/{self._trajs[xtc]}")
                # newu = mda.Universe(self._pdbf, f"{self.name}/{self._trajs[xtc]}")
                # prot = newu.select_atoms("protein")
                arr = np.empty((prot.n_residues, traj.n_frames, 3)) # newu.trajectory.n_frames
                for ts in traj: # newu.trajectory
                    arr[:, ts.frame] = prot.center_of_mass(compound='residues')

                cau = mda.Universe(compdb, arr, format=mda.coordinates.memory.MemoryReader, order='afc')
                # cau.load_new(arr, format=mda.coordinates.memory.MemoryReader, order='afc')
                with mda.Writer(comtraj, cau.atoms.n_atoms) as W:
                    for ts in cau.trajectory:
                        W.write(cau.atoms)
        
        # To add cau
        # cau = mda.Universe(compdb, *comtrajs)
        # setattr(self, "_cau", cau)

        return #compdb, comtraj
    
    
    
    
    
#     def calculate(self, pkg="all", cores=None, ow=False, filterby="incontact"):
#         pkgs = pkgsl if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
        
#         if filterby == "incontact":
#             if hasattr(self, "_pair"):
#                 self._pair.get_incontact(cores)
#             elif not rhasattr(self.data, "raw", "getcontacts"): #hasattr(self.data, "raw") or not hasattr(self.data.raw, "getcontacts") or o
#                 with Pool(cores) as pool:
#                     self._send_calc("getcontacts", filterby, pool)
#                     pool.close()
#                     pool.join()
#                 self._add_pkg("getcontacts")
#             if "getcontacts" in pkgs: pkgs.remove("getcontacts")
            
        
#         with Pool(cores) as pool:
#             for pkg in pkgs:
#                 if not rhasattr(self.data, "raw", pkg) or ow:#hasattr(self.data, "raw") or not hasattr(self.data.raw, pkg) or ow:
#                     #print(f"sending {self.name} {pkg}")
#                     self._send_calc(pkg, ow, filterby, pool)
#             pool.close()
#             pool.join()
        
#         [self._add_pkg(pkg) for pkg in pkgs]
    def calculate(self, pkg="all"):
        pkgs = pkgsl if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
        
        if any(["COM" in pkg for pkg in pkgs]):
            self._add_comtrajs()
        
        
        # pool = get_pool()
        
        for pkg in pkgs:
            pkgclass = eval(pkg[0].upper() + pkg[1:]) if isinstance(pkg, str) else pkg
            setattr(self.data, pkgclass.__name__, pkgclass(self))
            print(pkg)
    
        
        
        
    
    
    
#     def _send_calc(self, pkg, ow, filterby, pool):
#         def send_empty(self, pqf):
#             while not os.path.isfile(pqf):
#                 time.sleep(5)
#             return
        
#         os.makedirs(self._rawpn(pkg), exist_ok=True)        
#         pkgf = self._getcontactsf if "getcontacts" in pkg else self._corrsf
#         # pkgf = getcontactsf if "getcontacts" in pkg else corrsf
#         if "COM" in pkg and not hasattr(self, "_comtrajs"):#or not hasattr(self, "_cau")
#             self._add_comtrajs()
        
#         for xtc in self._trajs:
#             pqf = self._rawfn(pkg, xtc)
            
#             if not os.path.isfile(pqf) or ow:
#                 print(f"sending calculation of {pkg} for {self.name} {xtc}")
#                 pool.apply_async(pkgf, (xtc, pkg, ow, filterby))
#                 # if "getcontacts" in pkg or "dynetan" in pkg:
#                 #     [pool.apply_async(send_empty, (pqf,))]*2
            
#         return
    
#     def _getcontactsf(self, xtc, pkg, ow, filterby):
#         getcontactsf(self, xtc, pkg, ow)
#         return
    
#     def _corrsf(self, xtc, pkg, ow, filterby):
#         corrsf(self, xtc, pkg, ow, filterby)
#         return
    
    
#     def _add_pkg(self, pkg):
#         if not hasattr(self.data, "raw"): setattr(self.data, "raw", Store())
#         print(f"adding raw data of {pkg} for {self.name}")
        
#         pqs = [self._rawfn(pkg, xtc) for xtc in self._trajs]
#         flist = map(lambda pq: pandas.read_parquet(pq), pqs)
#         df = pandas.concat(flist, axis=1)
#         cols = [f"{num}" for num in self._trajs]
#         df["weight_avg"] = df[cols].fillna(0).mean(axis=1)
#         df["weight_std"] = df[cols].fillna(0).std(axis=1)
        
#         setattr(self.data.raw, pkg, df)
        
#         return
    
    
    
    
    
    
    def analyze(self, pkg="all", metrics="all", normalize=True, ow=False, filterby="incontact"):
        pkgs = pkgsl if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
        metrics = metricsl if metrics=="all" else metrics if isinstance(metrics, list) else [metrics]
        
        self.data._check_raw(pkgs, normalize, filterby)
        self._send_anal(pkgs, metrics, normalize, ow, filterby)
        
        [self.data._add_anal(pkg, metrics, normalize, ow, filterby) for pkg in pkgs]
    
    
    def _send_anal(self, pkgs, metrics, normalize, ow, filterby):
        if not os.path.isdir(self._datapn(filterby, normalize, "")) or         any([not os.path.isfile(self._datafn(filterby, normalize, pkg, metric)) 
             for pkg in pkgs for metric in metrics]) or ow:
            ps = []
            print(f"sending analysis of {self.name} {filterby} data of {pkgs} packages for {self._get_norm_str(normalize)} {metrics} metrics")
            for pkg in pkgs:
                p = Process(target=self.data._analyze, args=(pkg, metrics, normalize, ow, filterby))
                ps.append(p)
                p.start()
            [p.join() for p in ps]
        
        
        return  
    
    
    
    
    def view(self, pkg, metric, normalize=True, filterby="incontact", num=20):
        norm = self._get_norm_str(normalize)
        if not rhasattr(self.data, filterby, norm, pkg):#hasattr(self.data, filterby) or not hasattr(getattr(self.data, filterby), norm) or not hasattr(rgetattr(self.data, filterby, norm), pkg): 
            self.analyze(pkg, normalize=normalize, filterby=filterby)
            
        return rgetattr(self.data, filterby, norm, pkg).view(metric, num, filterby)
