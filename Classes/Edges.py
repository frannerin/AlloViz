class Edges:
    def __init__(self, parent, df):
        self._parent = parent
        self.df = df
    
    
    def __repr__(self):
        print(self.df.shape)
        return repr(self.df.iloc[:, :1])
    
    
    def __sub__(self, other):
#         if any(col not in other.df.columns for col in data.df.columns):

        cols = [col for col in self.df.columns if col in other.df.columns]
    
        subs = [col for col in cols if "avg" in col]
        sub = pandas.DataFrame.sub(self.df[subs], other.df[subs], fill_value = 0)
    
        adds = [col for col in cols if "std" in col]
        add = pandas.DataFrame.add(self.df[adds], other.df[adds], fill_value = 0)
        
        return pandas.concat([add, sub], axis=1)
    
    
    
    
    def _get_cmap(self):
        if isinstance(self._parent, Pair): 
            colors = {"inactive": "r", "active": "g",
                     "Gprotein": "y", "Barr": "b"}
            color2, color1 = colors[self._parent.state1.name], colors[self._parent.state2.name]
        else:
            color1, color2 = "orange", "turquoise"
        
        return matplotlib.colors.LinearSegmentedColormap.from_list('bar', [color1, "w", color2], 2048)
    
    
    
    
    def _get_colors(self, col):
        cmap = self._get_cmap()
        scale = [0, col.min(), col.max()] if isinstance(self._parent, Pair) else [col.mean(), col.min(), col.max()]
        normdata = matplotlib.colors.TwoSlopeNorm(*scale)
        npdata = col.to_numpy()
        edges = np.nonzero(npdata)
        return cmap( normdata( npdata[edges] ).data )




    def _add_edge(self, nv, prot, edge, color, radius):    
        get_coords = lambda resnum: list( prot.select_atoms(f"resnum {resnum} and name CA").center_of_geometry() )
        coords = [get_coords(res.split(':')[-1]) for res in edge]

        return nv.shape.add_cylinder(coords[0], coords[1], list(color),
                                     np.float64(radius), f"{edge[0]}_{edge[1]}")
    
    
    
    def _show_cbar(self):
        cmap = self._get_cmap()
        pl.imshow([[0,1],], cmap = cmap)
        pl.gca().set_visible(False)
        if isinstance(self._parent, Pair):
            cbar = pl.colorbar(orientation = "horizontal", ticks = [0,1])
            cbar.ax.set_xticklabels([self._parent.state2.name.capitalize(), self._parent.state1.name.capitalize()])
        else:
            cbar = pl.colorbar(orientation = "horizontal", ticks = [])

        return
    
    
    
    def view(self, metric, num=20, filterby="incontact"):
        metric = f"{metric}_avg" if not re.search("_avg$", metric) else metric
        # if metric not in df.columns etc
        get_data = lambda num: self.df.sort_values(metric, key = abs, ascending = False)[0:num]
        data = get_data(num)
        
        if isinstance(self._parent, Pair):
            while not any(0 < data[metric]) or not any(0 > data[metric]):
                num += 1
                data = get_data(num)
            print(num)

        self._show_cbar()
        
        mdau = self._parent.state1.mdau if isinstance(self._parent, Pair) else self._parent.mdau
        prot = mda.core.universe.Merge(mdau.select_atoms("protein or segid LIG"))
        nv = nglview.show_mdanalysis(prot, default=False)
        nv.add_cartoon('protein', color='white')

        edges = np.nonzero(data[metric].to_numpy())
        colors = self._get_colors(data[metric])
        
        error = "weight_std" if "weight_avg" in metric else metric.replace("avg", "std")
        radii = np.interp(data[error], (data[error].min(), data[error].max()), (1, 0.1))

        for i in range(len(edges[0])):
            self._add_edge(nv, mdau, data.index[edges[0][i]], colors[i], radii[i])

        return nv
