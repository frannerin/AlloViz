"""Wrappers of AlloViz's network construction methods/packages

Each module contains the class(es) that wrap the code of a package used by AlloViz as a
network construction method. They are used to send their calculations through the AlloViz
code without needing to execute them independently, and process their results to store
them in AlloViz objects. They are child classes of :mod:`~AlloViz.Wrappers.Base` module's
:class:`~AlloViz.Wrappers.Base.Base` class and other classes therein (by multiple 
inheritance).

"""

from .CARDS_w import CARDS_MI_Phi, CARDS_Disorder_Phi, CARDS_Disorder_mediated_Phi, CARDS_Holistic_Phi
from .CARDS_w import CARDS_MI_Psi, CARDS_Disorder_Psi, CARDS_Disorder_mediated_Psi, CARDS_Holistic_Psi
from .CARDS_w import CARDS_MI_Backbone_Dihs, CARDS_Disorder_Backbone_Dihs, CARDS_Disorder_mediated_Backbone_Dihs, CARDS_Holistic_Backbone_Dihs
# from .CARDS_w import CARDS_MI_Backbone_Dihs_Max, CARDS_Disorder_Backbone_Dihs_Max, CARDS_Disorder_mediated_Backbone_Dihs_Max, CARDS_Holistic_Backbone_Dihs_Max

from .CARDS_w import CARDS_MI_Chi1, CARDS_Disorder_Chi1, CARDS_Disorder_mediated_Chi1, CARDS_Holistic_Chi1
from .CARDS_w import CARDS_MI_Chi2, CARDS_Disorder_Chi2, CARDS_Disorder_mediated_Chi2, CARDS_Holistic_Chi2
from .CARDS_w import CARDS_MI_Chi3, CARDS_Disorder_Chi3, CARDS_Disorder_mediated_Chi3, CARDS_Holistic_Chi3
from .CARDS_w import CARDS_MI_Chi4, CARDS_Disorder_Chi4, CARDS_Disorder_mediated_Chi4, CARDS_Holistic_Chi4
from .CARDS_w import CARDS_MI_Sidechain_Dihs, CARDS_Disorder_Sidechain_Dihs, CARDS_Disorder_mediated_Sidechain_Dihs, CARDS_Holistic_Sidechain_Dihs
# from .CARDS_w import CARDS_MI_Sidechain_Dihs_Max, CARDS_Disorder_Sidechain_Dihs_Max, CARDS_Disorder_mediated_Sidechain_Dihs_Max, CARDS_Holistic_Sidechain_Dihs_Max

from .CARDS_w import CARDS_MI_Dihs, CARDS_Disorder_Dihs, CARDS_Disorder_mediated_Dihs, CARDS_Holistic_Dihs
# from .CARDS_w import CARDS_MI_Dihs_Max, CARDS_Disorder_Dihs_Max, CARDS_Disorder_mediated_Dihs_Max, CARDS_Holistic_Dihs_Max


from .AlloViz_w import AlloViz_Psi, AlloViz_Phi, AlloViz_Omega
from .AlloViz_w import AlloViz_Backbone_Dihs#, AlloViz_Backbone_Dihs_Max

from .AlloViz_w import AlloViz_Chi1, AlloViz_Chi2, AlloViz_Chi3, AlloViz_Chi4, AlloViz_Chi5
from .AlloViz_w import AlloViz_Sidechain_Dihs#, AlloViz_Sidechain_Dihs_Max

from .AlloViz_w import AlloViz_Dihs#, AlloViz_Dihs_Max



from .correlationplus_w import correlationplus_CA_Pear, correlationplus_CA_LMI
from .correlationplus_w import correlationplus_COM_Pear, correlationplus_COM_LMI

from .correlationplus_w import correlationplus_Phi, correlationplus_Psi, correlationplus_Omega
from .correlationplus_w import correlationplus_Backbone_Dihs#_Avg, correlationplus_Backbone_Dihs_Max



from .PyInteraph2_w import PyInteraph2_COM_Contacts#, PyInteraph2_COM_Contacts_Corrected # would be better to used the corrected as default
from .PyInteraph2_w import PyInteraph2_Atomic_Contacts_Strength, PyInteraph2_Atomic_Contacts_Occurrence
from .PyInteraph2_w import PyInteraph2_Energy



from .dynetan_w import dynetan#, dynetan_COM

from .g_correlation_w import g_correlation_CA_MI, g_correlation_COM_MI, g_correlation_CA_LMI, g_correlation_COM_LMI

from .GetContacts import GetContacts

# from .gRINN_w import gRINN

# from .GSAtools import GSAtools

from .MDEntropy_w import MDEntropy_Phi, MDEntropy_Psi, MDEntropy_Backbone_Dihs
from .MDEntropy_w import MDEntropy_Contacts, MDEntropy_AlphaAngle

from .MDTASK import MDTASK

from .pytraj_w import pytraj_CA, pytraj_CB