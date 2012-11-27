#/usr/bin/env python
#PBS -l nodes=4:ppn=8
#PBS -l walltime=13:00:00
from ase.optimize import QuasiNewton
from ase.constraints import FixAtoms
from ase.calculators.emt import EMT
from ase.neb import NEB
from ase.lattice.surface import fcc100, add_adsorbate
from ase.optimize.test import run_test
from gpaw import GPAW, Mixer, PoissonSolver
import gpaw.mpi as mpi


name = 'neb'

def get_atoms():
    # 2x2-Al(001) surface with 3 layers and an
    # Au atom adsorbed in a hollow site:
    slab = fcc100('Al', size=(2, 2, 3))
    add_adsorbate(slab, 'Au', 1.7, 'hollow')
    slab.center(axis=2, vacuum=4.0)

    # Fix second and third layers:
    mask = [atom.tag > 1 for atom in slab]
    slab.set_constraint(FixAtoms(mask=mask))

    # Use EMT potential:
    slab.set_calculator(EMT())

    # Initial state:
    qn = QuasiNewton(slab, logfile=None)
    qn.run(fmax=0.05)
    initial = slab.copy()

    # Final state:
    slab[-1].x += slab.get_cell()[0, 0] / 2
    qn = QuasiNewton(slab, logfile=None)
    qn.run(fmax=0.05)
    final = slab.copy()

    # Setup a NEB calculation
    constraint = FixAtoms(mask=[atom.tag > 1 for atom in initial])

    images = [initial]
    for i in range(3):
        image = initial.copy()
        image.set_constraint(constraint)
        images.append(image)

    images.append(final)

    neb = NEB(images, parallel=mpi.parallel)
    neb.interpolate()

    def set_calculator(calc):
        i = 0
        for image in neb.images[1:-1]:
            if not mpi.parallel or mpi.rank // (mpi.size // 3) == i:
                image.set_calculator(calc)
            i += 1
    neb.set_calculator = set_calculator

    return neb

def get_calculator_emt():
    return EMT()

def get_calculator_gpaw():
    if mpi.parallel:
        assert mpi.size % 3 == 0
        s = mpi.size // 3
        r0 = mpi.rank // s * s
        comm = range(r0, r0 + s)
    else:
        comm = mpi.world
    calc = GPAW(h=0.25,
                kpts=(2, 2, 1),
                communicator=comm,
                txt='neb-%d.txt' % r0)
    return calc

run_test(get_atoms, get_calculator_emt, name + '-emt')
run_test(get_atoms, get_calculator_gpaw, name + '-gpaw')
