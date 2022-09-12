import os, json

running_dir = os.path.dirname(os.path.realpath(__file__))

if "hydra" in os.environ["HOSTNAME"]:
    try:
        from rsync import rsync
        rsync(f"rsync -z --inplace ori:/protwis/sites/files/Precomputed/get_contacts_files/compl_info.json {running_dir}/compl_info.json")
        print("compl_info.json WAS updated")
    except:
        print("compl_info.json dictionary couldn't be updated")

with open(f"{running_dir}/compl_info.json", "r") as f:
    d = json.load(f)

    
    
    
    
    
    
    
ac_options = {
    'cf' : {
        'at' : ['gc', 'py2'],
        'rc' : ['py2', 'py2c'],
        'st' : ['py2'],
    },
    
    'mc' : {
        'wr' : {'mi' : ['dt']},
        'rc' : {
            'mi' : ['g_c'],
            'lmi': ['g_c', 'cp'],
            'pr' : ['cp']
        },
        'ac' : {
            'mi' : ['g_c'],
            'lmi': ['g_c', 'cp'],
            'pr' : ['cp', 'pt', 'mdt']
        },
        'bc' : {'pr' : ['pt']}
    },
    
    'dh' : {
        'bb' : {
            'ph' : {
                'pr' : ['cp'],
                'mi' : ['ad']    },
            'ps' : {
                'pr' : ['cp'],
                'mi' : ['ad']    },
            'om' : {
                'pr' : ['cp'],
                'mi' : ['ad']    },
            'bba' : {
                'pr' : ['cp'],
                'mi' : ['ad', 'mde']    },
            'bbm' : {
                'pr' : ['cp'],
                'mi' : ['ad']    },
            'ap' : {'mi' : ['mde']} },
        'sc' : {
            'ch1' : {'mi' : ['ad']},
            'ch2' : {'mi' : ['ad']},
            'ch3' : {'mi' : ['ad']},
            'ch4' : {'mi' : ['ad']},
            'ch5' : {'mi' : ['ad']},
            'sca' : {'mi' : ['ad']},
            'scm' : {'mi' : ['ad']}     },
        'all' : {
            'alla' : {'mi' : ['ad']},
            'allm' : {'mi' : ['ad']}    },
    },
    
    'en' : ['py2']
}


ac_options_codes = {
    #'cf' : 'Contact frequency', ####
    'gc' : 'GetContacts', ####
    'py2' : 'PyInteraph2',
    'mc' : "Atoms' movement correlation", ####
    'wr' : 'Whole residue',
    'mi' : 'Mutual Information (MI)',
    'dt' : 'dynetan',
    'rc' : 'Residue COM',
    'lmi' : 'Linear Mutual Information (LMI)',
    'pr' : 'Pearson',
    'g_c' : 'g_correlation',
    'cp' : 'correlationplus',
    'pt' : 'pytraj',
    'mdt': 'MD-TASK',
    'ac' : 'Alpha-carbon',
    'bc' : 'Beta-carbon',
    'dh' : "Dihedrals' movement correlation", ####
    'ph' : 'Phi',
    'ps' : 'Psi',
    'om' : 'Omega',
    'all' : 'All dihedrals',
    'ad' : 'AlloViz',
    # 'gr' : 'gRINN',
    # 'grc' : 'gRINN correlation',
    'en' : 'Interaction energy',
    'mde' : 'MDEntropy',
    'ap' : 'Alpha angle',
    
    'bb' : 'Backbone dihedrals',
    'bba' : 'All backbone dihedrals (average)',
    'bbm' : 'All backbone dihedrals (max. value)',
    'sc' : 'Side-chain dihedrals',
    'ch1' : 'Chi1',
    'ch2' : 'Chi2',
    'ch3' : 'Chi3',
    'ch4' : 'Chi4',
    'ch5' : 'Chi5',
    'sca' : 'All side-chain dihedrals (average)',
    'scm' : 'All side-chain dihedrals (max. value)',
    'alla' : 'All dihedrals (average)',
    'allm' : 'All dihedrals (max. value)',
    
    'cf' :  'Contacts',
    'cof' : 'Contact frequency',
    'comc' : 'Residue COM contacts',
    'cs' : 'Contact strength',
    'py2c' : 'PyInteraph2 (with Rg correction)',
    
    # 'wh' : 'Whole', ######
    # 'in' : 'Incontact', ######
    # 'it' : 'Intercontact', ######
    'wh' : 'All',
    'in' : 'GetContacts_edges',
    'it' : 'GPCR_Interhelix',
    'sd' : 'Spatially_distant',
}
    
    
abrvs = {val: key for key, val in ac_options_codes.items()}