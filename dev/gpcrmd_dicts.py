ac_options = {
            'cf' : ['gc', 'py2'],
            'mc' : {
                'wr' : { 'mi' : ['dt']},
                'rc' : {
                    'mi' : ['g_c'],
                    'lmi': ['g_c', 'cp'],
                    'pr' : ['cp']
                },
                'ac' : {
                    'mi' : ['g_c'],
                    'lmi': ['g_c', 'cp'],
                    'pr' : ['cp','pt', 'mdt']
                },
                'bc' : {
                    'pr' : ['pt']
                }
            },
            'dh' : {
                'ph' : {
                    'pr' : ['cp'],
                    'mi' : ['ad']    },
                'ps' : {
                    'pr' : ['cp'],
                    'mi' : ['ad']    },
                'om' : {
                    'pr' : ['cp'],
                    'mi' : ['ad']    },
                'ap' : {'mi' : ['mde']},
                'all' : {
                    'pr' : ['cp'],
                    'mi' : ['ad', 'mde']    }
            },
            'en' : ['py2']
        }
        ac_options_codes = {
            'cf' : 'Contact frequency',
            'gc' : 'getcontacts',
            'py2' : 'PyInteraph2',
            'mc' : 'Movement correlation',
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
            'dh' : 'Dihedral correlation',
            'ph' : 'Phi',
            'ps' : 'Psi',
            'om' : 'Omega',
            'all' : 'All dihedrals',
            'ad' : 'AlloViz',
            'gr' : 'gRINN',
            'grc' : 'gRINN correlation',
            'en' : 'Interaction energy',
            'mde' : 'MDEntropy',
            'ap' : 'Alpha'
        }