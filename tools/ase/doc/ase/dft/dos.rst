.. module:: dft.dos
   :synopsis: Density of states

=================
Density of states
=================

Example::

  calc = ...
  dos = DOS(calc, width=0.2)
  d = dos.get_dos()
  e = dos.get_energies()

You can plot the result like this::

  import pylab as plt
  plt.plot(e, d)
  plt.xlabel('energy [eV]')
  plt.ylabel('DOS')
  plt.show()

.. image:: dos.png

Calculations involving moments of a DOS distribution may be
facilitated by the use of :func:`~ase.dft.get_distribution_moment`
method, as in the following example::

  from ase.dft import get_distribution_moment
  volume = get_distribution_moment(e,d)
  center, width = get_distribution_moment(e,d,(1,2))


More details
------------

.. autoclass:: ase.dft.dos.DOS
   :members: get_energies, get_dos

.. automodule:: ase.dft
   :members: get_distribution_moment
