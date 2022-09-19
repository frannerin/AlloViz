"""Module containing variables and elements with AlloViz information

The main member is the :data:`~AlloViz.AlloViz.info.wrappers` dictionary, which contains
all the available network construction methods in AlloViz as keys and the kind of
information they use as their respective values.

"""

# Atom/angle names
alpha = "Carbon \u03B1"
beta = "Carbon \u03B2"
COM = "Residue COM"
whole = "Whole residue"
Dihs = "All backbone dihedrals (Phi, psi and omega)"

# Package/Network construction method-level common information
dynetani = ("Atoms' movement correlation", "dynetan", "Mutual Information (MI)")
pytraji = ("Atoms' movement correlation", "pytraj", "Pearson's")
correlationplusi = ("Atoms' movement correlation", "correlationplus")
correlationplusdihsi = ("Dihedrals' movement correlation", "correlationplus", "Pearson's")
AlloVizi = ("Dihedrals' movement correlation", "AlloViz", "MI")
MDEntropyi = ("MDEntropy", "MI")
g_correlationi = ("Atoms' movement correlation", "g_correlation")
pyinteraphi = ("Contacts", "PyInteraph2", "-")



#: Dictionary with all the available network construction methods and their information
#:
#: All network construction methods available in AlloViz are the keys of the dictionary
#: and each of the values is a 4-member tuple. In order, the tuple contains the kind of
#: information it uses for network construction, the main package name (not the same as
#: the AlloViz accession name of the wrapper), the correlaton metric they use if
#: applicable, and the atom/angle whose information it uses.
wrappers = {
    "MDTASK": ("Atoms' movement correlation", "MD-TASK", "Pearson's", alpha),

    "pytraj_CA": pytraji + (alpha,),
    "pytraj_CB": pytraji + (beta,),

    "dynetan": dynetani + (whole,),
    "dynetan_COM": dynetani + (COM,),
    
    "g_correlation_CA_MI": g_correlationi + ("MI", alpha), 
    "g_correlation_COM_MI": g_correlationi + ("MI", COM), 
    "g_correlation_CA_LMI": g_correlationi + ("LMI", alpha), 
    "g_correlation_COM_LMI": g_correlationi + ("LMI", COM),
    
    "correlationplus_CA_Pear": correlationplusi + ("Pearson's", alpha),
    "correlationplus_COM_Pear": correlationplusi + ("Pearson's", COM),
    "correlationplus_CA_LMI": correlationplusi + ("Linear MI (LMI)", alpha),
    "correlationplus_COM_LMI": correlationplusi + ("Linear MI (LMI)", COM),

    "correlationplus_Phi": correlationplusdihsi + ("Phi",),
    "correlationplus_Psi": correlationplusdihsi + ("Psi",),
    "correlationplus_Omega": correlationplusdihsi + ("Omega",),

    "correlationplus_Backbone_Dihs_Avg": correlationplusdihsi + (Dihs + " (average)",),
    "correlationplus_Backbone_Dihs_Max": correlationplusdihsi + (Dihs + " (max. value)",),
    
    "AlloViz_Phi": AlloVizi + ("Phi",),
    "AlloViz_Psi": AlloVizi + ("Psi",),
    "AlloViz_Omega": AlloVizi + ("Omega",),
    "AlloViz_Backbone_Dihs_Avg": AlloVizi + (Dihs + " (average)",),
    "AlloViz_Backbone_Dihs_Max": AlloVizi + (Dihs + " (max. value)",),

    "AlloViz_Chi1": AlloVizi + ("Chi1",),
    "AlloViz_Chi2": AlloVizi + ("Chi2",),
    "AlloViz_Chi3": AlloVizi + ("Chi3",),
    "AlloViz_Chi4": AlloVizi + ("Chi4",),
    "AlloViz_Chi5": AlloVizi + ("Chi5",),
    "AlloViz_Sidechain_Dihs_Avg": AlloVizi + ("All side-chain dihedrals (average)",),
    "AlloViz_Sidechain_Dihs_Max": AlloVizi + ("All side-chain dihedrals (max. value)",),

    "AlloViz_Dihs_Avg": AlloVizi + ("All dihedrals (average)",),
    "AlloViz_Dihs_Max": AlloVizi + ("All dihedrals (max. value)",),
    
    "MDEntropy_Dihs": ("Dihedrals' movement correlation",) + MDEntropyi + (Dihs,),
    "MDEntropy_AlphaAngle": ("Dihedrals' movement correlation",) + MDEntropyi + ("Alpha angle",),
    "MDEntropy_Contacts": ("Contacts",) + MDEntropyi + ("Contact frequency",),
    
    "GetContacts": ("Contacts", "GetContacts", "-", "Contact frequency"),
    
    "PyInteraph2_Atomic_Contacts_Occurrence": pyinteraphi + ("Contact frequency",),
    "PyInteraph2_Atomic_Contacts_Strength": pyinteraphi + ("Contact strength",),
    "PyInteraph2_COM_Contacts": pyinteraphi + ("Residue COM contacts",),
    "PyInteraph2_COM_Contacts_Corrected": ("Contacts", "PyInteraph2 (with Rg correction)", "-", "Residue COM contacts"),

    "PyInteraph2_Energy": ("Interaction energy", "PyInteraph2", "-", whole),
    
    "gRINN": ("Interaction energy", "gRINN", "-", whole),
    "gRINN_corr": ("Interaction energy", "gRINN", "Pearson's", whole),
}


# The inverse dictionary is used to construct and display a table with all the available network construction methods in the package documentation
# In that case, the aim is to show the general, more wide-spread columns first and the AlloViz accession name (keys in 'wrappers') last
inverse = {}

for key, val in wrappers.items():
    # g_correlation-related methods need local compilation of g_correlation and are only available locally; dynetan_COM is prohibitively slow
    if "g_correlation" not in key and "dynetan_COM" not in key:
        inverse[val] = key


# The df variable of the processed inverse dictionary is used when this module is imported by the documentation's conf.py script to print the updated table in the README.rst file
import pandas
df = pandas.DataFrame.from_dict(inverse, orient="index", columns=["Name in AlloViz"])
df.index = pandas.MultiIndex.from_tuples(list(df.index), names=["Residue information extracted from trajectories",
                                                            "Package",
                                                            "Correlation measurement",
                                                            "Atom/angle tracked"])