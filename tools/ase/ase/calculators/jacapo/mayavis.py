'''
mayavi interface to plot atoms, unit cells, and volumetric data
'''
import numpy as np
from enthought.mayavi import mlab
mlab.figure(1, bgcolor=(1,1,1), size=(350, 350))
mlab.clf()

def plot_cylinder(start,end,tube_radius=0.1,color=(0,0,0)):
    
    mlab.plot3d([start[0],end[0]],[start[1],end[1]],[start[2],end[2]],
                tube_radius=tube_radius,color=color)
    
    
def plot_atoms(atoms):

    for atom in atoms:
        pos = atom.position
        mlab.points3d([pos[0]],[pos[1]],[pos[2]],
                      scale_factor=4,
                      resolution=20,
                      color=(1,0,0),        #this should get species specifuc
                      scale_mode='none')

    (u0,u1,u2) = atoms.get_cell()
    origin = np.array([0.0,0.0,0.0])
    
    plot_cylinder(origin,u0)
    plot_cylinder(origin,u1)
    plot_cylinder(origin,u2)

    plot_cylinder(u0,u0+u1)
    plot_cylinder(u0,u0+u2)

    plot_cylinder(u1,u1+u0)
    plot_cylinder(u1,u1+u2)

    plot_cylinder(u2,u2+u0)
    plot_cylinder(u2,u2+u1)

    plot_cylinder(u0+u1,u0+u1+u2)
    plot_cylinder(u1+u2,u0+u1+u2)
    plot_cylinder(u0+u2,u0+u1+u2)
    mlab.show()


    



if __name__ == '__main__':
    from ase.lattice.cubic import *

    from ase.lattice.bravais import cross
    import numpy as np

    a = np.array([0.5,0,0])
    c = np.array([0,1,0],dtype=np.float)
    b1 = c - a

    a = np.array([0,1,0],np.float)
    c = np.array([0,0.5,0.5])
    b2 = c - a

    a3 = np.array([2,1,1],np.float)

    a1 = cross(b1,a3)
    a2 = cross(b2,a3)
    v211 = FaceCenteredCubic(directions=[a1,a2,a3],
                             miller=(None,None,[2,1,1]),
                             symbol='Pd',    
                             size=(1,1,2),  
                             debug=0)

    uc = v211.get_cell()
    uc[2][2] += 10.0
    v211.set_cell(uc)

    plot_atoms(v211.repeat((2,2,1)))
    
