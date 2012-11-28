from ase import Atoms
from ase.optimize import QuasiNewton
from ase.neb import NEB
from ase.optimize.mdmin import MDMin
try:
    from asap3 import EMT
except ImportError:
    pass
else:

    a = 3.6
    b = a / 2
    initial = Atoms('Cu4',
                    positions=[(0, 0, 0),
                               (0, b, b),
                               (b, 0, b),
                               (b, b, 0)],
                    cell=(a, a, a),
                    pbc=True)
    initial *= (4, 4, 4)
    del initial[0]
    images = [initial] + [initial.copy() for i in range(6)]
    images[-1].positions[0] = (0, 0, 0)
    for image in images:
        image.set_calculator(EMT())
        #image.set_calculator(ASAP())

    for image in [images[0], images[-1]]:
        QuasiNewton(image).run(fmax=0.01)
    neb = NEB(images)
    neb.interpolate()

    for a in images:
        print a.positions[0], a.get_potential_energy()

    dyn = MDMin(neb, dt=0.1, trajectory='mep1.traj')
    #dyn = QuasiNewton(neb)
    print dyn.run(fmax=0.01, steps=25)
    for a in images:
        print a.positions[0], a.get_potential_energy()
