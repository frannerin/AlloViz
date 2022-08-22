"""AlloViz package

The main classes are imported into the namespace: :class:`~AlloViz.Protein` and
:class:`~AlloViz.Delta`.

:mod:`AlloViz.AlloViz` contains the rest of modules used by :class:`~AlloViz.Protein`
and :class:`~AlloViz.Delta` to process information and send calculations and analyses
and produce representations.

:mod:`AlloViz.Wrappers` contains the classes that wrap the code of the packages/network
construction methods available in AlloViz, to send their calculations through the AlloViz
code without needing to execute them independently, and process their results to store
them in AlloViz objects.

"""

from .AlloViz.Classes import Protein, Delta

__all__ = ["Protein", "Delta"]