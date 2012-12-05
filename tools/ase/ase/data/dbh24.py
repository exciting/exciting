"""
The following contains a database of 24 gas-phase reaction barrier heights for small molecules.
It is the DBH24 (diverse barrier heights) set of the Truhlar group,
with 12 forward and 12 backward barriers.
All geometries are from
Zheng, Zhao and Truhler, "J. Chem. Theo. Comput.", 3:569-582, 2007
while energies are from
Zheng, Zhao and Truhler, "J. Chem. Theo. Comput.", 5:808-821, 2009
"""

from ase.atoms import Atoms

dbh24 = ['dbh24_H', 'dbh24_N2O','dbh24_OH','dbh24_N2','dbh24_tst_H_N2O__OH_N2',
	 'dbh24_HCl','dbh24_tst_H_ClH__HCl_H',
	 'dbh24_CH3','dbh24_FCl','dbh24_CH3F','dbh24_Cl','dbh24_tst_CH3_FCl__CH3F_Cl',
         'dbh24_Cl-ion_CH3Cl','dbh24_tst_Cl-ion_CH3Cl',
         'dbh24_F-ion_CH3Cl','dbh24_Cl-ion_CH3F','dbh24_tst-Cl-ion_CH3F__F_ion_CH3Cl',
         'dbh24_OH-ion','dbh24_CH3OH','dbh24_F-ion','dbh24_tst-OH-ion_CH3F__F_ion_CH3OH',
         'dbh24_HN2','dbh24_tst_H_N2__HN2',
         'dbh24_C2H4','dbh24_CH3CH2','dbh24_tst_H_C2H4__CH3CH2',
         'dbh24_HCN','dbh24_HNC','dbh24_tst_HCN__HNC',
         'dbh24_CH4','dbh24_H2O','dbh24_tst_OH_CH4__CH3_H2O',
         'dbh24_H2','dbh24_O','dbh24_tst_H_OH__O_H2',
         'dbh24_H2S','dbh24_HS','dbh24_tst_H_H2S__H2_HS']

dbh24_reaction_list = {
'dbh24_r1': {
    'description': 'HAT 1',
    'number': 1,
    'initial': ['dbh24_H', 'dbh24_N2O'],
    'final': ['dbh24_OH','dbh24_N2'],
    'tst': 'dbh24_tst_H_N2O__OH_N2'},
'dbh24_r2': {
    'description': 'HAT 2',
    'number': 2,
    'initial': ['dbh24_H', 'dbh24_HCl'],
    'final': ['dbh24_HCl', 'dbh24_H'],
    'tst': 'dbh24_tst_H_ClH__HCl_H'},
'dbh24_r3': {
    'description': 'HAT 3',
    'number': 3,
    'initial': ['dbh24_CH3', 'dbh24_FCl'],
    'final': ['dbh24_CH3F', 'dbh24_Cl'],
    'tst': 'dbh24_tst_CH3_FCl__CH3F_Cl'},
'dbh24_r4': {
    'description': 'NS 1',
    'number': 4,
    'initial': ['dbh24_Cl-ion_CH3Cl'],
    'final': ['dbh24_Cl-ion_CH3Cl'],
    'tst': 'dbh24_tst_Cl-ion_CH3Cl'},
'dbh24_r5': {
    'description': 'NS 2',
    'number': 5,
    'initial': ['dbh24_F-ion_CH3Cl'],
    'final': ['dbh24_Cl-ion_CH3F'],
    'tst': 'dbh24_tst-Cl-ion_CH3F__F_ion_CH3Cl'},
'dbh24_r6': {
    'description': 'NS 3',
    'number': 6,
    'initial': ['dbh24_OH-ion', 'dbh24_CH3F'],
    'final': ['dbh24_CH3OH', 'dbh24_F-ion'],
    'tst': 'dbh24_tst-OH-ion_CH3F__F_ion_CH3OH'},
'dbh24_r7': {
    'description': 'UA 1',
    'number': 7,
    'initial': ['dbh24_H', 'dbh24_N2'],
    'final': ['dbh24_HN2'],
    'tst': 'dbh24_tst_H_N2__HN2'},
'dbh24_r8': {
    'description': 'UA 2',
    'number': 8,
    'initial': ['dbh24_H', 'dbh24_C2H4'],
    'final': ['dbh24_CH3CH2'],
    'tst': 'dbh24_tst_H_C2H4__CH3CH2'},
'dbh24_r9': {
    'description': 'UA 3',
    'number': 9,
    'initial': ['dbh24_HCN'],
    'final': ['dbh24_HNC'],
    'tst': 'dbh24_tst_HCN__HNC'},
'dbh24_r10': {
    'description': 'HT 1',
    'number': 10,
    'initial': ['dbh24_OH', 'dbh24_CH4'],
    'final': ['dbh24_CH3', 'dbh24_H2O'],
    'tst': 'dbh24_tst_OH_CH4__CH3_H2O'},
'dbh24_r11': {
    'description': 'HT 2',
    'number': 11,
    'initial': ['dbh24_H', 'dbh24_OH'],
    'final': ['dbh24_O', 'dbh24_H2'],
    'tst': 'dbh24_tst_H_OH__O_H2'},
'dbh24_r12': {
    'description': 'HT 3',
    'number': 12,
    'initial': ['dbh24_H', 'dbh24_H2S'],
    'final': ['dbh24_H2', 'dbh24_HS'],
    'tst': 'dbh24_tst_H_H2S__H2_HS'}
}

data = {
# reaction 1 = HAT 1
'dbh24_H': {
    'name': 'dbh24_H',
    'symbols': 'H',
    'magmoms': [1.],
    'charge': 0.,
    'positions': [[0. , 0. , 0.]]},
'dbh24_N2O': {
    'name': "dbh24_N2O",
    'symbols': 'NNO',
    'magmoms': None,
    'charge': 0.,
    'positions': [[ 0.      ,  0.      , -1.195674],
                  [ 0.      ,  0.      , -0.075111],
                  [ 0.      ,  0.      ,  1.111937]]},
'dbh24_OH': {
    'name': "dbh24_OH",
    'symbols': 'OH',
    'magmoms': [ 1., 0.],
    'charge': 0.,
    'positions': [[ 0.      ,  0.      ,  0.106894],
                  [ 0.      ,  0.      , -0.855149]]},
'dbh24_N2': {
    'name': "dbh24_N2",
    'symbols': 'NN',
    'magmoms': None,
    'charge': 0.,
    'positions': [[ 0.     ,  0.     ,  0.548555],
                  [ 0.     ,  0.     , -0.548555]]},
'dbh24_tst_H_N2O__OH_N2': {
    'name': "dbh24_tst_H_N2O__OH_N2",
    'Vf': 17.13, # kcal/mol
    'Vb': 82.47, # kcal/mol
    'symbols': 'HONN',
    'magmoms': [1., 0., 0., 0.],
    'charge': 0.,
    'positions': [[ -0.303286,	-1.930712, 0.],
		  [ -0.861006,  -0.621526, 0.],
		  [  0.000000, 	 0.257027, 0.],
                  [  1.027333,   0.729104, 0.]]},
# reaction 2 = HAT 2
'dbh24_HCl': {
    'name': "dbh24_HCl",
    'symbols': 'HCl',
    'magmoms': None,
    'charge': 0.,
    'positions': [[ 0.     ,  0.     , -1.203645],
                  [ 0.     ,  0.     , 0.070803]]},
'dbh24_tst_H_ClH__HCl_H': {
    'name': "dbh24_tst_H_ClH__HCl_H",
    'Vf': 18.00, # kcal/mol
    'Vb': 18.00, # kcal/mol
    'symbols': 'HClH',
    'magmoms': [1., 0., 0.],
    'charge': 0.,
    'positions': [[ 0., 0.,  1.485800],
		  [ 0., 0.,  0.      ],
		  [ 0., 0., -1.485800]]},
# reaction 3 = HAT 3
'dbh24_CH3': {
    'name': "dbh24_CH3",
    'symbols': 'CHHH',
    'magmoms': [1.,0.,0.,0.],
    'charge': 0.,
    'positions': [[  0.,	0., 	  0.],
		  [  1.077317,	0., 	  0.],
		  [ -0.538659,  0.932984, 0.],
		  [ -0.538659, -0.932984, 0.]]},
'dbh24_FCl': {
    'name': "dbh24_FCl",
    'symbols': 'FCl',
    'magmoms': None,
    'charge': 0.,
    'positions': [[  0.,  0., -1.065985],
		  [  0.,  0.,   0.564345]]},
'dbh24_CH3F': {
    'name': "dbh24_CH3F",
    'symbols': 'CFHHH',
    'magmoms': None,
    'charge': 0.,
    'positions': [[ -0.632074,  0.000001,   0.000000],
		  [  0.749117,  0.000002,  -0.000002],
		  [ -0.983182, -0.338489,   0.972625],
		  [ -0.983222,  1.011553,  -0.193172],
		  [ -0.983203, -0.673084,  -0.779437]]},
'dbh24_Cl': {
    'name': 'dbh24_Cl',
    'symbols': 'Cl',
    'magmoms': [1.],
    'charge': 0.,
    'positions': [[0. , 0. , 0.]]},
'dbh24_tst_CH3_FCl__CH3F_Cl': {
    'name': "dbh24_tst_CH3_FCl__CH3F_Cl",
    'Vf': 6.75, # kcal/mol
    'Vb': 60.00, # kcal/mol
    'symbols': 'ClFCHHH',
    'magmoms': [0.,0.,1.,0.,0.,0.],
    'charge': 0.,
    'positions': [[  1.454749, 	-0.001237,	-0.000040],
		  [ -0.323587,	 0.004631,	 0.000124],
		  [ -2.387418,	-0.002147,	-0.000073],
		  [ -2.495086,  -0.855361,	-0.649404],
		  [ -2.497313,  -0.138673,	 1.063139],
		  [ -2.501537,	 0.986269, 	-0.413734]]},
# reaction 4 = NS 1
'dbh24_Cl-ion_CH3Cl': {
    'name': "dbh24_Cl-ion_CH3Cl",
    'symbols': 'ClCHHHCl',
    'magmoms': None,
    'charge': -1.,
    'positions': [[ 0.000000,  0.000000, -2.384735],
                  [ 0.000000,  0.000000, -0.566331],
                  [ 0.000000,  1.025066, -0.224379],
                  [-0.887734, -0.512533, -0.224379],
                  [ 0.887734, -0.512533, -0.224379],
                  [ 0.000000,  0.000000,  2.624213]]},
'dbh24_tst_Cl-ion_CH3Cl': {
    'name': "dbh24_tst_Cl-ion_CH3Cl",
    'Vf': 13.41, # kcal/mol
    'Vb': 13.41, # kcal/mol
    'symbols': 'ClCHHHCl',
    'magmoms': None,
    'charge': -1.,
    'positions': [[ 0.000025,  0.019526,  2.322499],
                  [ 0.000513,  0.000486, -0.000089],
                  [ 0.761278, -0.750733,  0.006377],
                  [-1.030451, -0.282724,  0.002147],
                  [ 0.270728,  1.034927, -0.008697],
                  [-0.000297, -0.019784, -2.322458]]},
# reaction 5 = NS 2
'dbh24_F-ion_CH3Cl': {
    'name': "dbh24_F-ion_CH3Cl",
    'symbols': 'ClCHHHF',
    'magmoms': None,
    'charge': -1.,
    'positions': [[ 0.000000,  0.000000,  1.623138],
                  [ 0.000000,  0.000000, -0.227358],
                  [ 0.000000,  1.026321, -0.555141],
                  [ 0.888820, -0.513160, -0.555141],
                  [-0.888820, -0.513160, -0.555141],
                  [ 0.000000,  0.000000, -2.729308]]},
'dbh24_Cl-ion_CH3F': {
    'name': "dbh24_Cl-ion_CH3F",
    'symbols': 'FCHHHCl',
    'magmoms': None,
    'charge': -1.,
    'positions': [[ 0.000000,  0.000000, -2.648539],
                  [ 0.000000,  0.000000, -1.240170],
                  [ 0.000000,  1.024719, -0.886406],
                  [-0.887432, -0.512359, -0.886406],
                  [ 0.887432, -0.512359, -0.886406],
                  [ 0.000000,  0.000000,  1.996299]]},
'dbh24_tst-Cl-ion_CH3F__F_ion_CH3Cl': {
    'name': "dbh24_tst-Cl-ion_CH3F__F_ion_CH3Cl",
    'Vf': 3.44, # kcal/mol
    'Vb': 29.42, # kcal/mol
    'symbols': 'FCHHHCl',
    'magmoms': None,
    'charge': -1.,
    'positions': [[ 0.000000,  0.000000, -2.537929],
                  [ 0.000000,  0.000000, -0.488372],
                  [ 1.062087,  0.000000, -0.614972],
                  [-0.531044,  0.919794, -0.614972],
                  [-0.531044, -0.919794, -0.614972],
                  [ 0.000000,  0.000000,  1.624501]]},
# reaction 6 = NS 3
'dbh24_OH-ion': {
    'name': "dbh24_OH-ion",
    'symbols': 'OH',
    'magmoms': None,
    'charge': -1.,
    'positions': [[ 0.000000,  0.000000,  0.106894],
                  [ 0.000000,  0.000000, -0.855149]]},
'dbh24_CH3OH': {
    'name': "dbh24_CH3OH",
    'symbols': 'COHHHH',
    'magmoms': None,
    'charge': 0.,
    'positions': [[ -0.046423,  0.663069,  0.000000],
                  [ -0.046423, -0.755063,  0.000000],
                  [ -1.086956,  0.975938,  0.000000],
                  [  0.860592, -1.057039,  0.000000],
                  [  0.438145,  1.071594,  0.889539],
                  [  0.438145,  1.071594, -0.889539]]},
'dbh24_F-ion': {
    'name': "dbh24_F-ion",
    'symbols': 'F',
    'magmoms': None,
    'charge': -1.,
    'positions': [[ 0.0, 0.0, 0.0]]},
'dbh24_tst-OH-ion_CH3F__F_ion_CH3OH': {
    'name': "dbh24_tst-OH-ion_CH3F__F_ion_CH3OH",
    'Vf': -2.44, # kcal/mol
    'Vb': 17.66, # kcal/mol
    'symbols': 'FCHHHOH',
    'magmoms': None,
    'charge': -1.,
    'positions': [[ 1.850614, -0.013179, -0.000128],
                  [ 0.090857,  0.010586,  0.000269],
                  [ 0.040907,  1.079548, -0.011749],
                  [ 0.037163, -0.528013, -0.922944],
                  [ 0.037486, -0.507463,  0.935132],
                  [-1.892801,  0.103266, -0.000118],
                  [-2.173821, -0.815112,  0.000039]]},
# reaction 7 = UA 1
'dbh24_HN2': {
    'name': "dbh24_HN2",
    'symbols': 'NNH',
    'magmoms': [1., 0., 0.],
    'charge': 0.,
    'positions': [[ -0.062442,  0.659491,  0.000000],
                  [ -0.062442, -0.518709,  0.000000],
                  [  0.874194, -0.985478,  0.000000]]},
'dbh24_tst_H_N2__HN2': {
    'name': "dbh24_tst_H_N2__HN2",
    'Vf': 14.36, # kcal/mol
    'Vb': 10.61, # kcal/mol
    'symbols': 'NNH',
    'magmoms': [1., 0., 0.],
    'charge': 0.,
    'positions': [[ 0.084563, -0.642934,  0.000000],
                  [ 0.084563,  0.479877,  0.000000],
                  [-1.183883,  1.141399,  0.000000]]},
# reaction 8 = UA 2
'dbh24_C2H4': {
    'name': "dbh24_C2H4",
    'symbols': 'CCHHHH',
    'magmoms': None,
    'charge': 0.,
    'positions': [[ 0.000000,  0.000000,  0.665593],
                  [ 0.000000,  0.000000, -0.665593],
                  [ 0.000000,  0.921495,  1.231668],
                  [ 0.000000, -0.921495,  1.231668],
                  [ 0.000000,  0.921495, -1.231668],
                  [ 0.000000, -0.921495, -1.231668]]},
'dbh24_CH3CH2': {
    'name': "dbh24_CH3CH2",
    'symbols': 'CCHHHHH',
    'magmoms': [1.,0.,0.,0.,0.,0.,0.],
    'charge': 0.,
    'positions': [[ -0.258719, -0.816829, 0.000000],
                  [ -0.250987,  0.674191, 0.000000],
                  [  0.758830, -1.225939, 0.000000],
                  [ -0.758830, -1.213866, 0.883419],
                  [ -0.758830, -1.213866,-0.883419],
                  [ -0.170021,  1.225939,-0.924320],
                  [ -0.170021,  1.225939, 0.924320]]},
'dbh24_tst_H_C2H4__CH3CH2': {
    'name': "dbh24_tst_H_C2H4__CH3CH2",
    'Vf': 1.72, # kcal/mol
    'Vb': 41.75, # kcal/mol
    'symbols': 'CCHHHHH',
    'magmoms': [1.,0.,0.,0.,0.,0.,0.],
    'charge': 0.,
    'positions': [[ -0.567877,  0.000051, -0.218958],
                  [  0.751139, -0.000036,  0.041932],
                  [ -1.493884, -0.000488,  1.531765],
                  [ -1.101691,  0.920651, -0.408626],
                  [ -1.102022, -0.920234, -0.409110],
                  [  1.299128, -0.922344,  0.173763],
                  [  1.298899,  0.922325,  0.174363]]},
# reaction 9 = UA 3
'dbh24_HCN': {
    'name': "dbh24_HCN",
    'symbols': 'CNH',
    'magmoms': None,
    'charge': 0.,
    'positions': [[ 0.000000,  0.000000,  -0.500365],
                  [ 0.000000,  0.000000,   0.652640],
                  [ 0.000000,  0.000000,  -1.566291]]},
'dbh24_HNC': {
    'name': "dbh24_HNC",
    'symbols': 'CNH',
    'magmoms': None,
    'charge': 0.,
    'positions': [[ 0.000000,  0.000000, -0.737248],
                  [ 0.000000,  0.000000,  0.432089],
                  [ 0.000000,  0.000000,  1.426960]]},
'dbh24_tst_HCN__HNC': {
    'name': "dbh24_tst_HCN__HNC",
    'Vf': 48.07, # kcal/mol
    'Vb': 32.82, # kcal/mol
    'symbols': 'CNH',
    'magmoms': None,
    'charge': 0.,
    'positions': [[ 0.080319,  0.620258,  0.000000],
                  [ 0.080319, -0.568095,  0.000000],
                  [-1.044148,  0.255121,  0.000000]]},
# reaction 10 = HT 1
'dbh24_CH4': {
    'name': "dbh24_CH4",
    'symbols': 'CHHHH',
    'magmoms': None,
    'charge': 0.,
    'positions': [[ 0.000000,  0.000000,  0.000000],
                  [ 0.627837,  0.627837,  0.627837],
                  [-0.627837, -0.627837,  0.627837],
                  [ 0.627837, -0.627837, -0.627837],
                  [-0.627837,  0.627837, -0.627837]]},
'dbh24_H2O': {
    'name': "dbh24_H2O",
    'symbols': 'OHH',
    'magmoms': None,
    'charge': 0.,
    'positions': [[ 0.000000,  0.000000,  0.117145],
                  [ 0.000000,  0.756709, -0.468582],
                  [ 0.000000, -0.756709, -0.468582]]},
'dbh24_tst_OH_CH4__CH3_H2O': {
    'name': "dbh24_tst_OH_CH4__CH3_H2O",
    'Vf': 6.7, # kcal/mol
    'Vb': 19.6, # kcal/mol
    'symbols': 'COHHHHH',
    'magmoms': [0.,1.,0.,0.,0.,0.,0.],
    'charge': 0.,
    'positions': [[ -1.211487,  0.007968,  0.000407],
                  [  1.293965, -0.108694,  0.000133],
                  [  0.009476, -0.118020,  0.002799],
                  [ -1.525529, -0.233250,  1.010070],
                  [ -1.430665,  1.033233, -0.278082],
                  [ -1.552710, -0.710114, -0.737702],
                  [  1.416636,  0.849894, -0.000591]]},
# reaction 11 = HT 2
'dbh24_O': {
    'name': 'dbh24_O',
    'symbols': 'O',
    'magmoms': [2.],
    'charge': 0.,
    'positions': [[0. , 0. , 0.]]},
'dbh24_H2': {
    'name': "dbh24_H2",
    'symbols': 'HH',
    'magmoms': None,
    'charge': 0.,
    'positions': [[ 0.000000,  0.000000,  0.370938],
                  [ 0.000000,  0.000000, -0.370938]]},
'dbh24_tst_H_OH__O_H2': {
    'name': "dbh24_tst_H_OH__O_H2",
    'Vf': 10.7, # kcal/mol
    'Vb': 13.1, # kcal/mol
    'symbols': 'HOH',
    'magmoms': [1.,0.,1.],
    'charge': 0.,
    'positions': [[ 0.000000,  0.000000, -0.860287],
                  [ 0.000000,  0.000000,  0.329024],
                  [ 0.000000,  0.000000, -1.771905]]},
# reaction 12 = HT 3
'dbh24_H2S': {
    'name': "dbh24_H2S",
    'symbols': 'SHH',
    'magmoms': None,
    'charge': 0.,
    'positions': [[ 0.000000,  0.000000,  0.102519],
                  [ 0.000000,  0.966249, -0.820154],
                  [ 0.000000, -0.966249, -0.820154]]},
'dbh24_HS': {
    'name': "dbh24_HS",
    'symbols': 'SH',
    'magmoms': [0.,1.],
    'charge': 0.,
    'positions': [[ 0.000000,  0.000000,  0.078835],
                  [ 0.000000,  0.000000, -1.261367]]},
'dbh24_tst_H_H2S__H2_HS': {
    'name': "dbh24_tst_H_H2S__H2_HS",
    'Vf': 3.6, # kcal/mol
    'Vb': 17.3, # kcal/mol
    'symbols': 'HSHH',
    'magmoms': [0.,1.,0.,0.],
    'charge': 0.,
    'positions': [[ 1.262097, -0.220097,  0.000000],
                  [ 0.000000,  0.223153,  0.000000],
                  [-0.500576, -1.115445,  0.000000],
                  [-0.761521, -2.234913,  0.000000]]},
}

def create_dbh24_system(name, **kwargs):
    """Creates a DBH24 system.
    """
    if name not in data:
        raise NotImplementedError('System %s not in database.' % name)
    d = data[name]
    if 'magmoms' not in kwargs:
        kwargs['magmoms'] = d['magmoms']
    return Atoms(d['symbols'], d['positions'], **kwargs)

def get_dbh24_magmoms(name):
    """Returns the magnetic moments of DBH24 systems.
    """
    if name not in data:
        raise KeyError('System %s not in database.' % name)
    else:
        return data[name]['magmoms']

def get_dbh24_charge(name):
    """ Returns the total charge of DBH24 systems.
    """
    assert name in dbh24
    d = data[name]
    charge = d['charge']
    return charge

def get_dbh24_Vf(name):
    """ Returns forward DBH24 TST barrier in kcal/mol
    """
    assert name in dbh24
    d = data[name]
    Vf = d['Vf']
    return Vf

def get_dbh24_Vb(name):
    """ Returns backward DBH24 TST barrier in kcal/mol
    """
    assert name in dbh24
    d = data[name]
    Vb = d['Vb']
    return Vb

def get_dbh24_initial_states(name):
    """ Returns initial DBH24 states
    """
    assert name in dbh24_reaction_list
    d = dbh24_reaction_list[name]
    initial = d['initial']
    return initial

def get_dbh24_final_states(name):
    """ Returns final DBH24 states
    """
    assert name in dbh24_reaction_list
    d = dbh24_reaction_list[name]
    final = d['final']
    return final

def get_dbh24_tst(name):
    """ Returns DBH24 TST names
    """
    assert name in dbh24_reaction_list
    d = dbh24_reaction_list[name]
    tst = d['tst']
    return tst
