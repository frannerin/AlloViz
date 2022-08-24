"""Wrappers of AlloViz's network construction methods/packages

Each module contains the class(es) that wrap the code of a package used by AlloViz as a
network construction method. They are used to send their calculations through the AlloViz
code without needing to execute them independently, and process their results to store
them in AlloViz objects. They are child classes of :mod:`~AlloViz.Wrappers.Base` module's
:class:`~AlloViz.Wrappers.Base.Base` class and other classes therein (by multiple 
inheritance).

"""

# from .AlloViz import AlloViz_Psi, AlloViz_Phi, AlloViz_Omega, AlloViz_Dihs
# from .correlationplus import correlationplus_CA_Pear, correlationplus_CA_LMI
# from .correlationplus import correlationplus_COM_Pear, correlationplus_COM_LMI
# from .correlationplus import correlationplus_Phi, correlationplus_Psi, correlationplus_Omega, correlationplus_Dihs
# from .dynetan import dynetan, dynetan_COM
# from .g_correlation import g_correlation_CA_MI, g_correlation_COM_MI, g_correlation_CA_LMI, g_correlation_COM_LMI
# from .GetContacts import GetContacts
# from .gRINN import gRINN
# # from .GSAtools import GSAtools
# from .MDEntropy import MDEntropy_Contacts, MDEntropy_Dihs, MDEntropy_AlphaAngle
# from .MDTASK import MDTASK
# from .PyInteraph2 import PyInteraph2_Contacts, PyInteraph2_Energy
# from .pytraj import pytraj_CA, pytraj_CB