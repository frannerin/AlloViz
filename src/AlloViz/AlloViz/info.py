alpha = "Carbon \u03B1"
beta = "Carbon \u03B2"
COM = "Residue COM"
whole = "Whole residue"
Dihs = "Backbone dihedrals\n(Phi, psi and omega)"

dynetani = ("Movement correlation", "dynetan", "Mutual Information (MI)")
pytraji = ("Movement correlation", "pytraj", "Pearson's")
correlationplusi = ("Movement correlation", "correlationplus")
correlationplusdihsi = ("Dihedral correlation", "correlationplus", "Pearson's")
AlloVizi = ("Dihedral correlation", "AlloViz", "MI")
MDEntropyi = ("MDEntropy", "MI")
g_correlationi = ("Movement correlation", "g_correlation")




wrappers = {
    "dynetan": dynetani + (whole,),
    "dynetan_COM": dynetani + (COM,),
    
    "pytraj_CA": pytraji + (alpha,),
    "pytraj_CB": pytraji + (beta,),
    
    "MDTASK": ("Movement correlation", "MD-TASK", "Pearson's", alpha),
    
    "correlationplus_CA_Pear": correlationplusi + ("Pearson's", alpha),
    "correlationplus_COM_Pear": correlationplusi + ("Pearson's", COM),
    "correlationplus_CA_LMI": correlationplusi + ("Linear MI (LMI)", alpha),
    "correlationplus_COM_LMI": correlationplusi + ("Linear MI (LMI)", COM),
    "correlationplus_Phi": correlationplusdihsi + ("Phi",),
    "correlationplus_Psi": correlationplusdihsi + ("Psi",),
    "correlationplus_Omega": correlationplusdihsi + ("Omega",),
    "correlationplus_Dihs": correlationplusdihsi + (Dihs,),
    
    "AlloViz_Phi": AlloVizi + ("Phi",),
    "AlloViz_Psi": AlloVizi + ("Psi",),
    "AlloViz_Omega": AlloVizi + ("Omega",),
    "AlloViz_Dihs": AlloVizi + (Dihs,),
    
    "MDEntropy_Dihs": ("Dihedral correlation",) + MDEntropyi + (Dihs,),
    "MDEntropy_AlphaAngle": ("Dihedral correlation",) + MDEntropyi + ("Alpha angle",),
    "MDEntropy_Contacts": ("Contact frequency",) + MDEntropyi + (whole,),
    
    "GetContacts": ("Contact frequency", "GetContacts", "-", whole),
    
    "PyInteraph2_Contacts": ("Contact frequency", "PyInteraph2", "-", whole),
    "PyInteraph2_Energy": ("Interaction energy", "PyInteraph2", "-", whole),
    
#     "gRINN": ("Interaction energy", "gRINN", NA, whole),
#     "gRINN_corr": ("Interaction energy", "gRINN", "Pearson's", whole),
    
#     "g_correlation_CA_MI": g_correlationi + ("MI", alpha), 
#     "g_correlation_COM_MI": g_correlationi + ("MI", COM), 
#     "g_correlation_CA_LMI": g_correlationi + ("LMI", alpha), 
#     "g_correlation_COM_LMI": g_correlationi + ("LMI", COM), 
}


inverse = {}

for key, val in wrappers.items():
    inverse[val] = key



import pandas
df = pandas.DataFrame.from_dict(inverse, orient="index", columns=["Name"])