# creates: spacegroup-al.png spacegroup-fe.png spacegroup-rutile.png spacegroup-cosb3.png spacegroup-mg.png spacegroup-skutterudite.png spacegroup-diamond.png spacegroup-nacl.png

import ase.io

execfile('spacegroup-al.py')
execfile('spacegroup-mg.py')
execfile('spacegroup-fe.py')
execfile('spacegroup-diamond.py')
execfile('spacegroup-nacl.py')
execfile('spacegroup-rutile.py')
execfile('spacegroup-skutterudite.py')

for name in ['al', 'mg', 'fe', 'diamond', 'nacl', 'rutile', 'skutterudite']:
    atoms = globals()[name]
    ase.io.write('spacegroup-%s.pov'%name, 
                 atoms, 
                 transparent=False, 
                 display=False, 
                 run_povray=True,
                 #canvas_width=128,
                 show_unit_cell=2,
                 rotation='10x,-10y', 
                 #celllinewidth=0.02,
                 celllinewidth=0.05,
                 )

execfile('spacegroup-cosb3.py')
    
