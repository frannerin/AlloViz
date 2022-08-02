import pandas
import MDAnalysis as mda
import numpy as np

from importlib import import_module
from lazyasd import LazyObject
matplotlib = LazyObject(lambda: import_module('matplotlib'), globals(), 'matplotlib')
pl = LazyObject(lambda: import_module('matplotlib.pyplot'), globals(), 'pl')
nglview = LazyObject(lambda: import_module('nglview'), globals(), 'nglview')


from . import Classes




class Element:
    def __init__(self, parent, df):
        self._parent = parent
        self.df = df
    
    
    
    def __repr__(self):
        print(self.df.shape)
        return repr(self.df.iloc[:, :1])
    
    
    def __sub__(self, other):
#         if any(col not in other.df.columns for col in data.df.columns):
        
        selfreindex = self.df.index.map(self._parent._translate_ix(self._parent._aln_mapper))
        selfdf = self.df.abs().reset_index(drop=True).assign(aln_pos=selfreindex.to_numpy()).set_index("aln_pos")
        
        otherreindex = other.df.index.map(other._parent._translate_ix(other._parent._aln_mapper))
        otherdf = other.df.abs().reset_index(drop=True).assign(aln_pos=otherreindex.to_numpy()).set_index("aln_pos")
        

        cols = [col for col in selfdf.columns if col in otherdf.columns]
    
        subs = [col for col in cols if "std" not in col]
        sub = pandas.DataFrame.sub(selfdf[subs], otherdf[subs], axis=0, level="aln_pos").dropna() #fill_value = 0, level="aln_pos"
    
        adds = [col for col in cols if "std" in col]
        add = pandas.DataFrame.add(selfdf[adds], otherdf[adds], axis=0, level="aln_pos").dropna() #fill_value = 0, axis=0, level="aln_pos"
        
        return pandas.concat([add, sub], axis=1).dropna()
    
    
    
    
    
    def _get_colors(self, col, cmap):
        # scale = [0, col.min(), col.max()] if isinstance(self._parent, Pair) else [col.mean(), col.min(), col.max()]
        scale = [0, col.min(), col.max()] if col.min() < 0 else [col.mean(), col.min(), col.max()]
        normdata = matplotlib.colors.TwoSlopeNorm(*scale)
        return cmap( normdata( col.to_numpy() ).data )
    
    
    
    def _show_cbar(self, cmap, minv, maxv):
        pl.imshow([[minv,maxv],], cmap = cmap)
        pl.gca().set_visible(False)
        if isinstance(self._parent, Classes.Delta):
            cbar = pl.colorbar(orientation = "horizontal", ticks = [minv,maxv])
            cbar.ax.set_xticklabels([self._parent.state2.name, self._parent.state1.name])
            cbar.ax.set_title("Delta-network")
        else:
            cbar = pl.colorbar(orientation = "horizontal", ticks = [minv,maxv])
            cbar.ax.set_title(self._parent.name)

        return
    
    
    def _get_nv(self, nv):
        mdau_parent = self._parent.state1 if isinstance(self._parent, Classes.Delta) else self._parent
        
        if nv is None:
            # prot = mda.core.universe.Merge(mdau_parent.protein).select_atoms(f"({mdau_parent._protein_sel}) or segid LIG"))
            nv = nglview.show_mdanalysis(mdau_parent.protein, default=False)
            nv.add_cartoon('protein', color='white')
        
        return nv, mdau_parent
    
    
    
    
    def view(self, metric, num=20, colors=["orange", "turquoise"], nv=None):
        # metric = f"{metric}_avg" if not re.search("_avg$", metric) else metric
        # if metric not in df.columns etc
        data = self.df.sort_values(metric, key = abs, ascending = False)
        subset = data[0:num]
            
        color1, color2 = colors
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list('bar', [color1, "w", color2], 2048)

        self._show_cbar(cmap, data[metric].min(), data[metric].max())
        colors = self._get_colors(data[metric], cmap)[:num]
        
        error = "weight_std" if "weight" in metric else f"{metric}_std"
        if error in data.columns:
            sizes = np.interp(subset[error], (subset[error].min(), subset[error].max()), (1, 0.1))
        else:
            sizes = np.ones(len(subset))
        
        nv, mdau_parent = self._get_nv(nv)
        #mdau = mdau_parent.mdau
        
        indices = subset.index.map(mdau_parent._translate_ix(mdau_parent._aln_mapper)) if subset.index.name == "aln_pos" else subset.index
        
        for i in range(len(subset[metric])):
            self._add_element(nv, mdau_parent.protein, indices[i], colors[i], sizes[i])

        return nv
    
    
    

    

    
class Edges(Element):
    def __init__(self, *args):
        super().__init__(*args)
        
        

    def _add_element(self, nv, prot, edge, color, size):
        get_coords = lambda res: list( prot.select_atoms(f"resnum {res.split(':')[-1]} and name CA").positions[0] )
        coords = [get_coords(res) for res in edge]

        return nv.shape.add_cylinder(coords[0], coords[1], list(color),
                                     np.float64(size), f"{edge[0]}_{edge[1]}")
    

    
class Nodes(Element):
    def __init__(self, *args):
        super().__init__(*args)
        
        
       
    def _add_element(self, nv, prot, node, color, size):
        
        get_coords = lambda res: list( prot.select_atoms(f"resnum {res.split(':')[-1]} and name CA").positions[0] )
        coords = get_coords(node)

        return nv.shape.add_sphere(coords, list(color),
                                   np.float64(size)*2.5, f"{node}")