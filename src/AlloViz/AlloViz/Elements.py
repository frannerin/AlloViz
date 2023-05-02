"""Module with classes to store and represent analyzed network elements

Parent base class :class:`~AlloViz.AlloViz.Elements.Element` extends
:class:`pandas.DataFrame` to store data and defines additional methods and private
methods to represent the networks on the protein structure. Children classes differ on
the type of network element (nodes or edges) displayed.

"""

import pandas
import MDAnalysis as mda
import numpy as np

from importlib import import_module
from lazyasd import LazyObject

matplotlib = LazyObject(lambda: import_module("matplotlib"), globals(), "matplotlib")
pl = LazyObject(lambda: import_module("matplotlib.pyplot"), globals(), "pl")
nglview = LazyObject(lambda: import_module("nglview"), globals(), "nglview")


from . import Classes


class Element(pandas.DataFrame):
    """Base class for network storage and representation
    
    This class extends :class:`pandas.DataFrame` with a
    :meth:`~AlloViz.AlloViz.Elements.Element.view` method (and related private methods)
    to visualize the networks on the protein structures, besides storing the analyzed
    data. It uses the private methods
    :meth:`~AlloViz.AlloViz.Elements.Element._get_colors` and
    :meth:`~AlloViz.AlloViz.Elements.Element._show_cbar` to establish the color
    scale and represent it as a colorbar, respectively, and
    :meth:`~AlloViz.AlloViz.Elements.Element._get_nv` to retrieve information about
    the protein structure to show and/or the visualization widget to add the network
    elements to. The different elements, cylinders for
    :class:`~AlloViz.AlloViz.Elements.Edges` and spheres for
    :class:`~AlloViz.AlloViz.Elements.Nodes`, are represented using the corresponding
    child class' `_add_element` private method.
    
    Once a class instance is created, the private attribute `_parent` must be established
    with the :class:`~AlloViz.Protein` or :class:`~AlloViz.Delta` object the information
    belongs to, as it is needed for representation.
    """
    
    
    # # temporary properties
    # _internal_names = pd.DataFrame._internal_names + ["internal_cache"]
    # _internal_names_set = set(_internal_names)

    # normal properties; retained after manipulating df; should new methods be here?
    # _metadata = 
    
    _metadata = ["_parent"]
    
    @property
    def _constructor(self):
        return self.__class__._internal_ctor
    
    @classmethod
    def _internal_ctor(cls, *args, **kwargs):
        kwargs['parent'] = None
        return cls(*args, **kwargs)
    
    
    @property
    def _constructor_sliced(self):
        return pandas.Series
    
    
    
    def __init__(self, data, parent, index=None, columns=None, dtype=None, copy=True):
        super().__init__(data=data,
                          index=index,
                          columns=columns,
                          dtype=dtype,
                          copy=copy)
        self._parent = parent
    
    
    
#     @property
#     def _constructor(self):
#         return self.__class__

    
    
#     def __init__(self, *args, parent):
#         super().__init__(*args)
#         self._parent = parent
    
    
    
    
    def __sub__(self, other):
        # Translate the current dataframe index (residues or residue pairs) to positions of the structural alignment made for Delta calculation
        selfreindex = self.index.map(
            self._parent._translate_ix(self._parent._aln_mapper)
        )
        # Transform data to absolute value and change the index to the alignment positions
        selfdf = (
            self.abs()
            .reset_index(drop=True)
            .assign(aln_pos=selfreindex.to_numpy())
            .set_index("aln_pos")
        )

        # Repeat with the corresponding Element instance of the "other"
        otherreindex = other.index.map(
            other._parent._translate_ix(other._parent._aln_mapper)
        )
        otherdf = (
            other.abs()
            .reset_index(drop=True)
            .assign(aln_pos=otherreindex.to_numpy())
            .set_index("aln_pos")
        )

        # Select the columns that the two have in common for delta-network calculation and that also averages and not from an individual trajectory (e.g., not _1, _2...)
        cols = [col for col in selfdf.columns if col in otherdf.columns and not any([str(i) in col for i in self._parent._trajs])]
        
        # The columns to be subtracted are the ones that are averages and also not standard errors; cells with equivalent "aln_pos" (index) are subtracted and NA are dropped in the end
        subs = [col for col in cols if "std" not in col]
        sub = pandas.DataFrame.sub(
            selfdf[subs], otherdf[subs], axis=0, level="aln_pos"
        ).dropna()  # fill_value = 0, level="aln_pos"

        # Standard errors are added up
        adds = [col for col in cols if "std" in col]
        add = pandas.DataFrame.add(
            selfdf[adds], otherdf[adds], axis=0, level="aln_pos"
        ).dropna()  # fill_value = 0, axis=0, level="aln_pos"

        return pandas.concat([add, sub], axis=1).dropna()

    def _get_colors(self, col, cmap):
        r"""Retrieve the colormapped data

        Scale the data passed as a :class:`pandas.Series` using
        :class:`matplotlib.colors.TwoSlopeNorm`, centered in 0 if there are negative
        values or in the mean, and apply the passed colormap.

        Parameters
        ----------
        col : :class:`pandas.Series`
            Whole data of the element-metric of which a subset is going to be represented.
        cmap : :class:`matplotlib.colors.LinearSegmentedColormap`
            Used colormap function.
        """
        # If the data has negative (it always will have positive) values, set the middle in 0, else use the mean
        scale = (
            [0, col.min(), col.max()]
            if col.min() < 0
            else [col.mean(), col.min(), col.max()]
        )
        # Interpolate the data (if we force the middle to be in 0, the two sides may not be equivalent and TwoSlopeNorm is used)
        normdata = matplotlib.colors.TwoSlopeNorm(*scale)
        # Return the colormapped data
        return cmap(normdata(col.to_numpy()).data)

    def _show_cbar(self, cmap, minv, maxv):
        r"""Plot the colorbar

        Show the colorbar with the colors used for the network, adapting the title and
        the text of the ticks depending on the information shown: the network of a single
        :class:`~AlloViz.Protein` structure or a :class:`~AlloViz.Delta`-network.

        Parameters
        ----------
        cmap : :class:`matplotlib.colors.LinearSegmentedColormap`
            Used colormap function.
        minv, maxv : float
            Minimum and maximum value of the whole data of which a subset is represented.
        """
        # matplotlib stuff to only show the colorbar with the selected cmap and min and max values
        pl.imshow(
            [
                [minv, maxv],
            ],
            cmap=cmap,
        )
        pl.gca().set_visible(False)
        # If a delta-network is represented, use the Proteins' passed names to show in the colorbar ticks and use "Delta-network" as title
        if isinstance(self._parent, Classes.Delta):
            cbar = pl.colorbar(orientation="horizontal", ticks=[minv, maxv])
            cbar.ax.set_xticklabels(
                [self._parent.state2.name, self._parent.refstate.name]
            )
            cbar.ax.set_title("Delta-network")
        # Else, if a single Protein's network is represented, show the values in the ticks and the passed name parameter as title
        else:
            cbar = pl.colorbar(orientation="horizontal", ticks=[minv, maxv])
            cbar.ax.set_title(self._parent.name)

        return

    def _get_nv(self, nv):
        r"""Establish the visualization object to be used for representation

        It selects the reference state of the Delta object as parent if a delta-network
        is visualized, else the parent parameter/attribute itself. If a nv is already
        passed, network elements should be added to it; else one is created with the
        selected Protein's protein attribute.

        Parameters
        ----------
        nv : None or :class:`nglview.NGLWidget`
            Existing representation to which new network elements should be added.
        """
        # Select the reference state of the Delta object as parent if a delta-network is visualized, else the parent parameter/attribute itself
        mdau_parent = (
            self._parent.refstate
            if isinstance(self._parent, Classes.Delta)
            else self._parent
        )

        # If a nv is passed, network elements should be added to it; else create one with the selected Protein's protein attribute
        if nv is None:
            nv = nglview.show_mdanalysis(mdau_parent.protein, default=False)
            nv.add_cartoon("protein", color="white")
            nv.center("protein")

        return nv, mdau_parent

    def view(self, metric, num=20, colors=["orange", "turquoise"], nv=None):
        r"""Represent the selected metric in the structure

        Retrieves the analyzed data corresponding to the present Element (depending on
        the class) and from it the corresponding metric column from the DataFrame. It is
        used to obtain the elements' colors, sizes (inv. proportional to errors, if
        available) and names (resnames); and the parent instance attribute is used to
        retrieve the structure for representation using :class:`nglview.NGLWidget`.

        Data is sorted according to the selected metric in absolute value and descending
        order, a :class:`~matplotlib.colors.LinearSegmentedColormap` is made with the
        passed colors. The colormap is represented in a colorbar through
        :meth:`~AlloViz.AlloViz.Elements.Element._show_cbar` and is used to
        establish the elements' colors through
        :meth:`~AlloViz.AlloViz.Elements.Element._get_colors`.

        Errors, if available (i.e., if more than one trajectory has been used and
        averages were calculated), are used to establish the elements' sizes to be
        inversely proportional to them (interpolated between 1 and 0.1), thus directly
        proportional to the "confidence" in the calculated value in a way.

        Elements are shown on a representation of the selected structure or added to the
        passed :class:`~nglview.NGLWidget` if applicable.

        Parameters
        ----------
        metric : str
            Metric/Name of the column in the object's df attribute to represent.
        num : int, default: 20
            Number of (each of the) network elements to show on the structure.
        colors : list, default: ["orange", "turquoise"]
            List of two colors to assign to the minimum and maximum values of the network
            to be represented, respectively. Middle value is assigned "white" and it will
            be the mean of the network values or 0 if the network has both negative and
            positive values.
        nv : :class:`nglview.NGLWidget`, optional
            A structure representation into which the shapes representing the chosen
            network elements will be added.

        See Also
        --------
        AlloViz.Protein.view
        """
        # Sort the data by the selected metric for visualization, in absolute value and descending order
        data = self.sort_values(metric, key=abs, ascending=False)
        # Select the subset to be represented according to 'num'
        subset = data[0:num]

        # Create the linear colormap to use with the passed colors
        color1, color2 = colors
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            "bar", [color1, "w", color2], 2048
        )

        # Show the colorbar
        self._show_cbar(cmap, data[metric].min(), data[metric].max())
        # Retrieve the colors of the data and store only those of the subset ([:num])
        colors = self._get_colors(data[metric], cmap)[:num]

        # Establish the name of the corresponding standard error column
        ## weight_std if the selected metric is the one resulting from analyzing the averaged "weight" column (f"{metric}_weight")
        ## f"{metric}_std" if the selected metric is the one resulting from averaging the independently-analyzed trajectories' raw data (if there are more than 1)
        error = "weight_std" if "weight" in metric else f"{metric}_std"
        # If there is more than 1 trajectory and there are error columns because averages could be calculated, inversely interpolate the errors to the (1, 0.1) range to use as sizes
        if error in data.columns:
            sizes = np.interp(
                subset[error], (subset[error].min(), subset[error].max()), (1, 0.1)
            )
        # Else, use all size "1"
        else:
            sizes = np.ones(len(subset))

        # Retrieve the NGLWidget to use (the one passed to this function or one made from the mdau_parent) and the parent whose structure is used for representation
        nv, mdau_parent = self._get_nv(nv)

        # Retrieve the names of the elements (nodes or edges) to be shown in the structure
        # If a delta-network is visualized, translate from alignment positions to the residue names of the reference structure
        indices = (
            subset.index.map(mdau_parent._translate_ix(mdau_parent._aln_mapper))
            if subset.index.name == "aln_pos"
            else subset.index
        )

        # Add iteratively the elements to the NGLWidget with the child classes' own method to add cylinders or spheres
        for i in range(len(subset[metric])):
            self._add_element(nv, mdau_parent.protein, indices[i], colors[i], sizes[i])

        return nv
    
    
class Edges(Element):
    """Class for storage and viz of Edges

    See this class' :meth:`~AlloViz.AlloViz.Elements.Edges._add_element`
    """

    def _add_element(self, nv, prot, edge, color, size):
        r"""Add a network element to the representation

        Adds a cylinder representing a network edge to the passed representation,
        retrieving the coordinates according to the participating residues included
        in the element's name from the passed structure and using the desired color
        and size.

        Parameters
        ----------
        nv : :class:`nglview.NGLWidget`
            A structure representation into which the shapes representing the chosen
            network elements will be added.
        prot : :class:`~MDAnalysis.core.universe.Universe`
            Universe of the chosen structure from which to retrieve atom coordinates.
        edge : str
            Name of the edge: contains the name of the participating residues in the
            selected structure.
        color : list
            RGB values of the element's color.
        size : float
            Size of the element.

        See Also
        --------
        AlloViz.AlloViz.Elements.Nodes._add_element
        """
        # Get the tridimensional positions of the CA atom of the passed residue number; nglview needs a list of three floats
        get_coords = lambda res: list(
            prot.select_atoms(f"resnum {res.split(':')[-1]} and name CA").positions[0]
        )
        # Edges/cylinders need two sets of coordinates to define their ends
        coords = [get_coords(res) for res in edge]

        # The result of the adding function must be returned for everything to work
        return nv.shape.add_cylinder(
            coords[0], coords[1], list(color), np.float64(size), f"{edge[0]}_{edge[1]}"
        )


class Nodes(Element):
    """Class for storage and viz of Nodes

    See this class' :meth:`~AlloViz.AlloViz.Elements.Nodes._add_element`
    """

    def _add_element(self, nv, prot, node, color, size):
        r"""Add a network element to the representation

        Adds a sphere representing a network node to the passed representation,
        retrieving the coordinates according to the passed residue/node's name
        from the passed structure and using the desired color and size.

        Parameters
        ----------
        nv : :class:`nglview.NGLWidget`
            A structure representation into which the shapes representing the chosen
            network elements will be added.
        prot : :class:`~MDAnalysis.core.universe.Universe`
            Universe of the chosen structure from which to retrieve atom coordinates.
        node : str
            Name of the node, i.e. of the residues as it is in theselected structure.
        color : list
            RGB values of the element's color.
        size : float
            Size of the element.

        See Also
        --------
        AlloViz.AlloViz.Elements.Edges._add_element
        """
        # Get the tridimensional positions of the CA atom of the passed residue number; nglview needs a list of three floats
        get_coords = lambda res: list(
            prot.select_atoms(f"resnum {res.split(':')[-1]} and name CA").positions[0]
        )
        coords = get_coords(node)

        # The result of the adding function must be returned for everything to work
        return nv.shape.add_sphere(
            coords, list(color), np.float64(size) * 2.5, f"{node}"
        )