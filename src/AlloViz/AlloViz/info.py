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
Dihs = "All backbone dihedrals (Phi and psi)"

# Package/Network construction method-level common information
dynetani = ("Atoms' displacements", "dynetan", "MI")
pytraji = ("Atoms' displacements", "pytraj", "Pearson's")
correlationplusi = ("Atoms' displacements", "correlationplus")
correlationplusdihsi = ("Dihedral angles", "correlationplus", "Pearson's")
AlloVizi = ("Dihedral angles", "AlloViz", "MI")
MDEntropyi = ("MDEntropy", "MI")
g_correlationi = ("Atoms' displacements", "g_correlation")
pyinteraphi = ("Contacts", "PyInteraph2", "None")
CARDSi = ("Dihedral angles", "CARDS")



#: Dictionary with all the available network construction methods and their information
#:
#: All network construction methods available in AlloViz are the keys of the dictionary
#: and each of the values is a 4-member tuple. In order, the tuple contains the kind of
#: information it uses for network construction, the main package name (not the same as
#: the AlloViz accession name of the wrapper), the correlaton metric they use if
#: applicable, and the atom/angle whose information it uses.
wrappers = {
    "MDTASK": ("Atoms' displacements", "MD-TASK", "Pearson's", alpha),

    "pytraj_CA": pytraji + (alpha,),
    "pytraj_CB": pytraji + (beta,),

    "dynetan": dynetani + (whole,),
    # "dynetan_COM": dynetani + (COM,),
    
    "g_correlation_CA_MI": g_correlationi + ("MI", alpha), 
    "g_correlation_COM_MI": g_correlationi + ("MI", COM), 
    "g_correlation_CA_LMI": g_correlationi + ("LMI", alpha), 
    "g_correlation_COM_LMI": g_correlationi + ("LMI", COM),
    
    "correlationplus_CA_Pear": correlationplusi + ("Pearson's", alpha),
    # "correlationplus_COM_Pear": correlationplusi + ("Pearson's", COM),
    "correlationplus_CA_LMI": correlationplusi + ("LMI", alpha),
    # "correlationplus_COM_LMI": correlationplusi + ("LMI", COM),

    "correlationplus_Phi": correlationplusdihsi + ("Phi",),
    "correlationplus_Psi": correlationplusdihsi + ("Psi",),
    #"correlationplus_Omega": correlationplusdihsi + ("Omega",),

    "correlationplus_Backbone_Dihs": correlationplusdihsi + (Dihs,),
    # "correlationplus_Backbone_Dihs_Max": correlationplusdihsi + (Dihs + " (max. value)",),
    
    "AlloViz_Phi": AlloVizi + ("Phi",),
    "AlloViz_Psi": AlloVizi + ("Psi",),
    #"AlloViz_Omega": AlloVizi + ("Omega",),
    "AlloViz_Backbone_Dihs": AlloVizi + (Dihs,),
    # "AlloViz_Backbone_Dihs_Max": AlloVizi + (Dihs + " (max. value)",),

    "AlloViz_Chi1": AlloVizi + ("Chi1",),
    "AlloViz_Chi2": AlloVizi + ("Chi2",),
    "AlloViz_Chi3": AlloVizi + ("Chi3",),
    "AlloViz_Chi4": AlloVizi + ("Chi4",),
    #"AlloViz_Chi5": AlloVizi + ("Chi5",),
    "AlloViz_Sidechain_Dihs": AlloVizi + ("All side-chain dihedrals",),
    # "AlloViz_Sidechain_Dihs_Max": AlloVizi + ("All side-chain dihedrals (max. value)",),

    "AlloViz_Dihs": AlloVizi + ("All dihedrals",),
    # "AlloViz_Dihs_Max": AlloVizi + ("All dihedrals (max. value)",),
    
    'CARDS_MI_Phi': CARDSi + ('MI', 'Phi'), 
    'CARDS_Disorder_Phi': CARDSi + ('Pure-disorder MI', 'Phi'), 
    'CARDS_Disorder_mediated_Phi': CARDSi + ('Disorder-mediated MI', 'Phi'), 
    'CARDS_Holistic_Phi': CARDSi + ('Holistic MI', 'Phi'), 
    'CARDS_MI_Psi': CARDSi + ('MI', 'Psi'), 
    'CARDS_Disorder_Psi': CARDSi + ('Pure-disorder MI', 'Psi'), 
    'CARDS_Disorder_mediated_Psi': CARDSi + ('Disorder-mediated MI', 'Psi'), 
    'CARDS_Holistic_Psi': CARDSi + ('Holistic MI', 'Psi'),
    
    'CARDS_MI_Backbone_Dihs': CARDSi + ('MI', Dihs), 
    'CARDS_Disorder_Backbone_Dihs': CARDSi + ('Pure-disorder MI', Dihs),
    'CARDS_Disorder_mediated_Backbone_Dihs': CARDSi + ('Disorder-mediated MI', Dihs), 
    'CARDS_Holistic_Backbone_Dihs': CARDSi + ('Holistic MI', Dihs),
    # 'CARDS_MI_Backbone_Dihs_Max': CARDSi + ('MI', 'All backbone dihedrals (Phi and psi) (max. value)'), 
    # 'CARDS_Disorder_Backbone_Dihs_Max': CARDSi + ('Pure-disorder MI', 'All backbone dihedrals (Phi and psi) (max. value)'), 
    # 'CARDS_Disorder_mediated_Backbone_Dihs_Max': CARDSi + ('Disorder-mediated MI', 'All backbone dihedrals (Phi and psi) (max. value)'), 
    # 'CARDS_Holistic_Backbone_Dihs_Max': CARDSi + ('Holistic MI', 'All backbone dihedrals (Phi and psi) (max. value)'), 
    
    'CARDS_MI_Chi1': CARDSi + ('MI', 'Chi1'), 
    'CARDS_Disorder_Chi1': CARDSi + ('Pure-disorder MI', 'Chi1'), 
    'CARDS_Disorder_mediated_Chi1': CARDSi + ('Disorder-mediated MI', 'Chi1'), 
    'CARDS_Holistic_Chi1': CARDSi + ('Holistic MI', 'Chi1'), 
    'CARDS_MI_Chi2': CARDSi + ('MI', 'Chi2'), 
    'CARDS_Disorder_Chi2': CARDSi + ('Pure-disorder MI', 'Chi2'), 
    'CARDS_Disorder_mediated_Chi2': CARDSi + ('Disorder-mediated MI', 'Chi2'), 
    'CARDS_Holistic_Chi2': CARDSi + ('Holistic MI', 'Chi2'), 
    'CARDS_MI_Chi3': CARDSi + ('MI', 'Chi3'), 
    'CARDS_Disorder_Chi3': CARDSi + ('Pure-disorder MI', 'Chi3'), 
    'CARDS_Disorder_mediated_Chi3': CARDSi + ('Disorder-mediated MI', 'Chi3'), 
    'CARDS_Holistic_Chi3': CARDSi + ('Holistic MI', 'Chi3'), 
    'CARDS_MI_Chi4': CARDSi + ('MI', 'Chi4'), 
    'CARDS_Disorder_Chi4': CARDSi + ('Pure-disorder MI', 'Chi4'), 
    'CARDS_Disorder_mediated_Chi4': CARDSi + ('Disorder-mediated MI', 'Chi4'), 
    'CARDS_Holistic_Chi4': CARDSi + ('Holistic MI', 'Chi4'), 
    
    'CARDS_MI_Sidechain_Dihs': CARDSi + ('MI', 'All side-chain dihedrals'), 
    'CARDS_Disorder_Sidechain_Dihs': CARDSi + ('Pure-disorder MI', 'All side-chain dihedrals'), 
    'CARDS_Disorder_mediated_Sidechain_Dihs': CARDSi + ('Disorder-mediated MI', 'All side-chain dihedrals'), 
    'CARDS_Holistic_Sidechain_Dihs': CARDSi + ('Holistic MI', 'All side-chain dihedrals'),
    # 'CARDS_MI_Sidechain_Dihs_Max': CARDSi + ('MI', 'All side-chain dihedrals (max. value)'), 
    # 'CARDS_Disorder_Sidechain_Dihs_Max': CARDSi + ('Pure-disorder MI', 'All side-chain dihedrals (max. value)'), 
    # 'CARDS_Disorder_mediated_Sidechain_Dihs_Max': CARDSi + ('Disorder-mediated MI', 'All side-chain dihedrals (max. value)'), 
    # 'CARDS_Holistic_Sidechain_Dihs_Max': CARDSi + ('Holistic MI', 'All side-chain dihedrals (max. value)'), 
    
    'CARDS_MI_Dihs': CARDSi + ('MI', 'All dihedrals'), 
    'CARDS_Disorder_Dihs': CARDSi + ('Pure-disorder MI', 'All dihedrals'), 
    'CARDS_Disorder_mediated_Dihs': CARDSi + ('Disorder-mediated MI', 'All dihedrals'), 
    'CARDS_Holistic_Dihs': CARDSi + ('Holistic MI', 'All dihedrals'), 
    # 'CARDS_MI_Dihs_Max': CARDSi + ('MI', 'All dihedrals (max. value)'), 
    # 'CARDS_Disorder_Dihs_Max': CARDSi + ('Pure-disorder MI', 'All dihedrals (max. value)'), 
    # 'CARDS_Disorder_mediated_Dihs_Max': CARDSi + ('Disorder-mediated MI', 'All dihedrals (max. value)'), 
    # 'CARDS_Holistic_Dihs_Max': CARDSi + ('Holistic MI', 'All dihedrals (max. value)'),
    
    "MDEntropy_Phi": ("Dihedral angles",) + MDEntropyi + ("Phi",),
    "MDEntropy_Psi": ("Dihedral angles",) + MDEntropyi + ("Psi",),
    "MDEntropy_Backbone_Dihs": ("Dihedral angles",) + MDEntropyi + (Dihs,),
    "MDEntropy_AlphaAngle": ("Dihedral angles",) + MDEntropyi + ("Alpha angle",),
#     "MDEntropy_Contacts": ("Contacts",) + MDEntropyi + ("Contact frequency",),
    
    "GetContacts": ("Contacts", "GetContacts", "None", "Contact frequency"),
    
    "PyInteraph2_Atomic_Contacts_Occurrence": pyinteraphi + ("Contact frequency",),
    "PyInteraph2_Atomic_Contacts_Strength": pyinteraphi + ("Contact strength",),
    "PyInteraph2_COM_Contacts": pyinteraphi + ("Residue COM contacts",),
    # "PyInteraph2_COM_Contacts_Corrected": ("Contacts", "PyInteraph2 (with Rg correction)", "None", "Residue COM contacts"),
    # it would be better to use the Rg corrected version

    "PyInteraph2_Energy": ("Interaction energy", "PyInteraph2", "None", whole),
    
    # "gRINN": ("Interaction energy", "gRINN", "None", whole),
    # "gRINN_corr": ("Interaction energy", "gRINN", "Pearson's", whole),
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




#: Dictionary with MDAnalysis' selection functions-ready dihedral-participating atomnames
#:
#: Each entry from the dictionary corresponds to a dihedral and is a dictionary itself
#: ready to be passed to MDAnalysis' selection functions (phi_selection, etc) as kwargs.
#: Chi2 and onward can be selected with the chi1_selection function passing custom atom
#: names.
#:
#: Info retrieved from `<https://docs.mdanalysis.org/stable/documentation_pages/core/topologyattrs.html#MDAnalysis.core.topologyattrs.Atomnames.chi1_selection>`_
#: and `<http://www.mlb.co.jp/linux/science/garlic/doc/commands/dihedrals.html>`_.
dihedrals_atoms = {
    "phi": { # C’-N-CA-C
        "c_name": 'C',
        "n_name": 'N',
        "ca_name": 'CA'
    },
    "psi": { # N-CA-C-N’
        "c_name": 'C',
        "n_name": 'N',
        "ca_name": 'CA'
    },
    "omega": { # CA-C-N’-CA’
        "c_name": 'C',
        "n_name": 'N',
        "ca_name": 'CA'
    },
    "chi1": {
        "n_name": 'N',
        "ca_name": 'CA',
        "cb_name": 'CB',
        "cg_name": 'CG CG1 OG OG1 SG'
    },
    "chi2": {
        "n_name": "CA",
        "ca_name": "CB",
        "cb_name": "CG CG1",
        "cg_name": "CD OD1 ND1 CD1 SD"
    },
    "chi3": {
        "n_name": "CB",
        "ca_name": "CG",
        "cb_name": "CD SD",
        "cg_name": "NE OE1 CE"
    },
    "chi4": {
        "n_name": "CG",
        "ca_name": "CD",
        "cb_name": "CE NE",
        "cg_name": "CZ NZ"
    },
    "chi5": {
        "n_name": "CD",
        "ca_name": "NE",
        "cb_name": "CZ",
        "cg_name": "NH1"
    }
}

