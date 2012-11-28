"Test the ase.utils.geometry module"

import numpy as np

from ase.lattice.spacegroup import crystal
from ase.utils.geometry import get_layers, cut, stack

np.set_printoptions(suppress=True)

al = crystal('Al', [(0,0,0)], spacegroup=225, cellpar=4.05)


# Cut out slab of 5 Al(001) layers
al001 = cut(al, nlayers=5)
correct_pos = np.array([[ 0. ,  0. ,  0. ],
                        [ 0. ,  0.5,  0.2],
                        [ 0.5,  0. ,  0.2],
                        [ 0.5,  0.5,  0. ],
                        [ 0. ,  0. ,  0.4],
                        [ 0. ,  0.5,  0.6],
                        [ 0.5,  0. ,  0.6],
                        [ 0.5,  0.5,  0.4],
                        [ 0. ,  0. ,  0.8],
                        [ 0.5,  0.5,  0.8]])
assert np.allclose(correct_pos, al001.get_scaled_positions())

# Check layers along 001
tags, levels = get_layers(al001, (0, 0, 1))
assert np.allclose(tags, [0, 1, 1, 0, 2, 3, 3, 2, 4, 4])
assert np.allclose(levels, [ 0., 2.025, 4.05, 6.075, 8.1])

# Check layers along 101
tags, levels = get_layers(al001, (1, 0, 1))
assert np.allclose(tags, [0, 1, 5, 3, 2, 4, 8, 7, 6, 9])
assert np.allclose(levels, [0.000, 0.752, 1.504, 1.880, 2.256,
                            2.632, 3.008, 3.384, 4.136, 4.888], atol=0.001)

# Check layers along 111
tags, levels = get_layers(al001, (1, 1, 1))
assert np.allclose(tags, [0, 2, 2, 4, 1, 5, 5, 6, 3, 7])
assert np.allclose(levels, [0.000, 1.102, 1.929, 2.205, 
                            2.756, 3.031, 3.858, 4.960], atol=0.001)


# Cut out slab of three Al(111) layers
al111 = cut(al, (1,-1,0), (0,1,-1), nlayers=3)
correct_pos = np.array([[ 0.5       ,  0.        ,  0.        ],
                        [ 0.        ,  0.5       ,  0.        ],
                        [ 0.5       ,  0.5       ,  0.        ],
                        [ 0.        ,  0.        ,  0.        ],
                        [ 1/6.      ,  1/3.      ,  1/3.      ],
                        [ 1/6.      ,  5/6.      ,  1/3.      ],
                        [ 2/3.      ,  5/6.      ,  1/3.      ],
                        [ 2/3.      ,  1/3.      ,  1/3.      ],
                        [ 1/3.      ,  1/6.      ,  2/3.      ],
                        [ 5/6.      ,  1/6.      ,  2/3.      ],
                        [ 5/6.      ,  2/3.      ,  2/3.      ],
                        [ 1/3.      ,  2/3.      ,  2/3.      ]])
assert np.allclose(correct_pos, al111.get_scaled_positions())

# Cut out cell including all corner and edge atoms (non-periodic structure)
al = cut(al, extend=1.1)
correct_pos = np.array([[ 0.   ,  0.   ,  0.   ],
                        [ 0.   ,  2.025,  2.025],
                        [ 2.025,  0.   ,  2.025],
                        [ 2.025,  2.025,  0.   ],
                        [ 0.   ,  0.   ,  4.05 ],
                        [ 2.025,  2.025,  4.05 ],
                        [ 0.   ,  4.05 ,  0.   ],
                        [ 2.025,  4.05 ,  2.025],
                        [ 0.   ,  4.05 ,  4.05 ],
                        [ 4.05 ,  0.   ,  0.   ],
                        [ 4.05 ,  2.025,  2.025],
                        [ 4.05 ,  0.   ,  4.05 ],
                        [ 4.05 ,  4.05 ,  0.   ],
                        [ 4.05 ,  4.05 ,  4.05 ]])
assert np.allclose(correct_pos, al.positions)

# Create an Ag(111)/Si(111) interface
ag = crystal(['Ag'], basis=[(0,0,0)], spacegroup=225, cellpar=4.09)
si = crystal(['Si'], basis=[(0,0,0)], spacegroup=227, cellpar=5.43)

ag111 = cut(ag, a=(4, -4, 0), b=(4, 4, -8), nlayers=5)
si111 = cut(si, a=(3, -3, 0), b=(3, 3, -6), nlayers=5)
#
#interface = stack(ag111, si111)
#assert len(interface) == 1000
#assert np.allclose(interface.positions[::100],
#                   [[  4.08125   ,  -2.040625  ,   -2.040625  ],
#                    [  8.1625    ,   6.121875  ,  -14.284375  ],
#                    [ 10.211875  ,   0.00875   ,    2.049375  ],
#                    [ 24.49041667,  -4.07833333,  -16.32208333],
#                    [ 18.37145833,  14.29020833,  -24.48166667],
#                    [ 24.49916667,  12.25541667,  -20.39458333],
#                    [ 18.36854167,  16.32791667,  -30.60645833],
#                    [ 19.0575    ,   0.01166667,    5.45333333],
#                    [ 23.13388889,   6.80888889,    1.36722222],
#                    [ 35.3825    ,   5.45333333,  -16.31333333]])
#
