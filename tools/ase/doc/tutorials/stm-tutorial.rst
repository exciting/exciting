.. _stm-tutorial:

==============================
Tutorial: STM images - Al(100)
==============================

The STM is a revolutionary experimental surface probe that has
provided direct local insight into the surface electronic
structure. Sometimes the interpretation of STM topographs are not
straightforward and therefore theoretically modeled STM images may
resolve conflicting possibilities and point to an underlying atomistic
model. The CAMPOS code includes python modules for generating
Tersoff-Hamann STM topographs. The STM code is illustrated here for a
Al(100) in a 2x2 unit cell.

Let's make the Al(100) fcc surface by using the :mod:`lattice` module::

  from ase.lattice.surface import *
  atoms = fcc100('Al', size=(2,2,2))

Now a calculator must be defined, in this tutorial we will make a STM
image from the GPAW calculator.

For the GPAW code the calculator for the Al(100) surface can be
defined like this::

  from gpaw import GPAW
  calc = GPAW(gpts=(28,28,20),nbands=28,
  	kpts=(4,4,1),txt='Al100.out')
  atoms.set_calculator(calc)
  energy = atoms.get_potential_energy() 
  calc.write('Al100.gpw', 'all')


Linescans
=========

In this section we will make simulated STM linescans and contour plot
using matplotlib. First initialize the :class:`STM` object and get the
averaged current along the z-direction::

  from ase import STM
  stm = STM(atoms, symmetries=[0, 1, 2])
  z = 2.5
  c = stm.get_averaged_current(z)

From the current we make a scan to get a 2D array of constant current
heights::

  h = stm.scan(c)

Finally we make a contour plot::

  import pylab as p
  p.contourf(h, 40)
  p.hot()
  p.colorbar()
  p.show()	
  

