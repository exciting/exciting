.. module::  emt
   :synopsis: Effective Medium Theory

==========================
Pure Python EMT calculator
==========================

The EMT potential is included in the ASE package in order to have a
simple calculator that can be used for quick demonstrations and
tests.

.. warning::

   If you want to do a real application using EMT, you should used the
   *much* more efficient implementation in the ASAP_ calculator.

.. class:: EMT()

Right now, the only supported elements are: H, C, N, O, Al, Ni, Cu,
Pd, Ag, Pt and Au.  The EMT parameters for the metals are quite
realistic for many purposes, whereas the H, C, N and O parameters are
just for fun!


.. _ASAP: http://wiki.fysik.dtu.dk/asap
