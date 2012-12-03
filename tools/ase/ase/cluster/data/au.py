"""Element data - 79 Au Gold"""

name = 'Gold'
symbol = 'Au'
symmetry = 'fcc'

energy_slope = -0.090
energy_intersection = -1.537

energy_types = {1: -2.618, #bulk
                2: -2.237, #100 surface
                3: -2.343, #110 surface (top) or 111-111 edge
                4: -2.369, #110 surface (bottom)
                5: -2.352, #111 surface
                6: -2.028, #100-110 edge
                7: -2.215, #100-111 edge
                }

data = {'name': name,
        'symbol': symbol,
        'symmetry': symmetry,
        'energy_slope': energy_slope,
        'energy_intersection': energy_intersection,
        'energy_types': energy_types,
        }

