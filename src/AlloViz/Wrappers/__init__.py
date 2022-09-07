"""Wrappers of AlloViz's network construction methods/packages

Each module contains the class(es) that wrap the code of a package used by AlloViz as a
network construction method. They are used to send their calculations through the AlloViz
code without needing to execute them independently, and process their results to store
them in AlloViz objects. They are child classes of :mod:`~AlloViz.Wrappers.Base` module's
:class:`~AlloViz.Wrappers.Base.Base` class and other classes therein (by multiple 
inheritance).

"""

from .AlloViz_w import AlloViz_Psi, AlloViz_Phi, AlloViz_Omega
from .AlloViz_w import AlloViz_Backbone_Dihs_Avg, AlloViz_Backbone_Dihs_Max

from .AlloViz_w import AlloViz_Chi1, AlloViz_Chi2, AlloViz_Chi3, AlloViz_Chi4, AlloViz_Chi5
from .AlloViz_w import AlloViz_Sidechain_Dihs_Avg, AlloViz_Sidechain_Dihs_Max

from .AlloViz_w import AlloViz_Dihs_Avg, AlloViz_Dihs_Max



from .correlationplus_w import correlationplus_CA_Pear, correlationplus_CA_LMI
from .correlationplus_w import correlationplus_COM_Pear, correlationplus_COM_LMI

from .correlationplus_w import correlationplus_Phi, correlationplus_Psi, correlationplus_Omega
from .correlationplus_w import correlationplus_Backbone_Dihs_Avg, correlationplus_Backbone_Dihs_Max



from .PyInteraph2_w import PyInteraph2_COM_Contacts, PyInteraph2_COM_Contacts_Corrected
from .PyInteraph2_w import PyInteraph2_Atomic_Contacts_Strength, PyInteraph2_Atomic_Contacts_Occurrence
from .PyInteraph2_w import PyInteraph2_Energy



from .dynetan_w import dynetan, dynetan_COM

# from .g_correlation_w import g_correlation_CA_MI, g_correlation_COM_MI, g_correlation_CA_LMI, g_correlation_COM_LMI

from .GetContacts_w import GetContacts

from .gRINN_w import gRINN

# from .GSAtools_w import GSAtools

from .MDEntropy_w import MDEntropy_Contacts, MDEntropy_Dihs, MDEntropy_AlphaAngle

from .MDTASK_w import MDTASK

from .pytraj_w import pytraj_CA, pytraj_CB