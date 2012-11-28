.. module:: transport
   :synopsis: Electron transport

==================
Electron transport
==================

.. default-role:: math

The :mod:`transport` module of ASE assumes the generic setup of the system in
question sketched below:

. . . |setup| . . .

.. |setup| image:: transport_setup.png
   :align: middle

There is a central region (blue atoms plus the molecule) connected to
two semi-infinite leads constructed by infinitely repeated *principal
layers* (red atoms). The entire structure may be periodic in the
transverse direction, which can be effectively sampled using
**k**-points (yellowish atoms).

The system is described by a Hamiltonian matrix which must be
represented in terms of a localized basis set such that each element
of the Hamiltonian can be ascribed to either the left, central, or
right region, *or* the coupling between these.

The Hamiltonian can thus be decomposed as:

.. math::

    H = \begin{pmatrix}
      \ddots      & V_L         &             &             &     \\
      V_L^\dagger & H_L         & V_L         &             &     \\
                  & V_L^\dagger & H_C         & V_R         &     \\
                  &             & V_R^\dagger & H_R         & V_R \\
                  &             &             & V_R^\dagger & \ddots
    \end{pmatrix}

where `H_{L/R}` describes the left/right principal layer, and `H_C`
the central region. `V_{L/R}` is the coupling between principal
layers, *and* from the principal layers into the central region.  The
central region must contain at least one principal layer on each side,
and more if the potential has not converged to its bulk value at this
size. The central region is assumed to be big enough that there is no
direct coupling between the two leads. The principal layer must be so
big that there is only coupling between nearest neighbor layers.

Having defined `H_{L/R}`, `V_{L/R}`, and `H_C`, the elastic
transmission function can be determined using the Nonequilibrium
Green Function (NEGF) method.  This is achieved by the class:
:class:`~ase.transport.calculators.TransportCalculator` (in
ase.transport.calculators) which makes no requirement on the origin of
these five matrices.

.. class:: ase.transport.calculators.TransportCalculator(energies, h, h1, h2, s=None, s1=None, s2=None, align_bf=False)

  Determine transport properties of device sandwiched between
  semi-infinite leads using nonequillibrium Green function methods.

  energies is the energy grid on which the transport properties should
  be determined.

  h1 (h2) is a matrix representation of the Hamiltonian of two
  principal layers of the left (right) lead, and the coupling between
  such layers.

  h is a matrix representation of the Hamiltonian of the scattering
  region. This must include at least on lead principal layer on each
  side. The coupling in (out) of the scattering region is assumed to
  be identical to the coupling between left (right) principal layers.

  s, s1, and s2 are the overlap matrices corresponding to h, h1, and
  h2. Default is the identity operator.

  If align_bf is True, the onsite elements of the Hamiltonians will be
  shifted to a common fermi level.


This module is stand-alone in the sense that it makes no requirement
on the origin of these five matrices. They can be model Hamiltonians
or derived from different kinds of electronic structure codes.

For an example of how to use the :mod:`transport` module, see the GPAW
exercise on `electron transport`_

.. _electron transport: http://wiki.fysik.dtu.dk/gpaw/exercises/transport/transport.html

.. default-role::
