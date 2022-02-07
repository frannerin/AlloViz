class Pair:
    def __init__(self, refstate, state2):
        self.state1, self.state2 = refstate, state2
        self.states = [self.state1, self.state2]
        for state in self.states:
            state._pair = self
        
        nodes_subset = self._add_nodes_subset()
        self.sources_subset, self.targets_subset = nodes_subset.values()
    
    
    
    def get_delta(self):
        delta = self.state1 - self.state2
        setattr(self, "delta", delta)
        setattr(self.delta, "pair", self)
        return delta
        
     
    
    
    def _add_nodes_subset(self):
        bfac = lambda atom: f"{atom.tempfactor:.2f}"
        is_wb = lambda atom: atom.name == "N" and -8.1 <= atom.tempfactor <= 8.1 and (bfac(atom) != "1.00" and bfac(atom) != "0.00")
        
        
        targets = [3.50, 6.30, 7.49, 7.50, 7.51, 7.52, 7.53]
        # In Ballesteros-Weinstein: Ionic lock 3.50 and 6.30; NPxxY 7.49-7.53
    
    
        def get_sources(states):
            resnum_to_wb = lambda nums, state: (bfac(atom)                                                 for residue in nums                                                 for atom in state.mdau.select_atoms(f"resnum {residue}")                                                 if is_wb(atom))
            
            has_lig = [state if ("LIG" in (seg.segid for seg in state.mdau.segments)) else False for state in self.states]
            if any(has_lig):
                state = next(item for item in has_lig if item != False)
                aas = state.mdau.select_atoms("(same residue as around 4 segid LIG) and protein").residues
                return list(resnum_to_wb(aas.resnums, state))
            else:
                return [3.32] # Conserved position
      
    
        nodes_subset = {"sources": [val for val in get_sources(self.states) if val not in targets],
                        "targets": list(map(lambda x: f"{x:.2f}", targets))}
        
        
        
        format_res = lambda aa: f"A:{aa.resname}:{aa.resnum}"
        wb_to_aa = lambda wblist, state: list(map(format_res, (atom.residue                                                       for residue in state.mdau.select_atoms("protein").residues                                                       for atom in residue.atoms                                                       if is_wb(atom) and bfac(atom) in wblist)))
        for state in self.states:
            state.sources_subset, state.targets_subset = wb_to_aa(nodes_subset["sources"], state), wb_to_aa(nodes_subset["targets"], state)       
        
        
        return nodes_subset    
    
    
    
    
    
    def calculate(self, pkg="all", cores=None, ow=False, filterby="incontact"):  
        pkgs = pkgsl if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
        
        if filterby == "incontact":
            if not hasattr(self, "_incontact"):
                self.get_incontact(cores)
            if "getcontacts" in pkgs: pkgs.remove("getcontacts")
            # with Pool(cores) as pool:
            #     for state in self.states:
            #         if not hasattr(state.data, "raw") or not hasattr(state.data.raw, "getcontacts") or ow:
            #             state._send_calc("getcontacts", ow, filterby, pool)
            #     pool.close()
            #     pool.join()
            # [state._add_pkg("getcontacts") for state in self.states]
            # pkgs.remove("getcontacts")
        
        with Pool(cores) as pool:
            for pkg in pkgs:
                for state in self.states:
                    if not rhasattr(state.data, "raw", pkg) or ow: #hasattr(state.data, "raw") or not hasattr(state.data.raw, pkg) or ow:
                        #print(f"sending {state.name} {pkg}")
                        state._send_calc(pkg, ow, filterby, pool)
            pool.close()
            pool.join()
        
        [state._add_pkg(pkg) for pkg in pkgs for state in self.states]
        
        return
    
    
    
    def get_incontact(self, cores=None, ow=False):
        with Pool(cores) as pool:
            for state in self.states:
                if not rhasattr(state.data, "raw", "getcontacts") or ow: #hasattr(state.data, "raw") or not hasattr(state.data.raw, "getcontacts") or ow:
                    #print(f"sending {state.name} {pkg}")
                    state._send_calc("getcontacts", ow, "incontact", pool)
            pool.close()
            pool.join()
        
        [state._add_pkg("getcontacts") for state in self.states]
        
        edges = set((idx for idx in state.data.raw.getcontacts.index for state in self.states))
        setattr(self, "_incontact", edges)
        return edges
    
    
    
    
    def analyze(self, pkg="all", metrics="all", normalize=True, ow=False, filterby="incontact"):
        pkgs = pkgsl if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
        metrics = metricsl if metrics=="all" else metrics if isinstance(metrics, list) else [metrics]
        
        for state in self.states: state.data._check_raw(pkgs, normalize, filterby)
        
        ps = []
        for state in self.states:
            p = Process(target=state._send_anal, args=(pkgs, metrics, normalize, ow, filterby))
            #print(f"sent analysis of {state.name} {filterby} data of {pkgs} packages for {metrics} metrics: {p.name}")
            ps.append(p)
            p.start()
        
        [p.join() for p in ps]
        
        [state.data._add_anal(pkg, metrics, normalize, ow, filterby) for pkg in pkgs for state in self.states]
        
        return


    
    
    def view(self, pkg, metric, normalize=True, filterby="incontact", num=20):
        norm = self.state1._get_norm_str(normalize)
        if not hasattr(self, "delta"): self.get_delta()
        if not rhasattr(self.delta, filterby, norm, pkg): self.analyze(pkg) #hasattr(self.delta, filterby) or not hasattr(getattr(self.delta, filterby), norm) or not hasattr(rgetattr(self.delta, filterby, norm), pkg): self.analyze(pkg)
        return rgetattr(self.delta, filterby, norm, pkg).view(metric, num, filterby)
