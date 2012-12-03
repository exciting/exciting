.. module:: structure

.. seealso:: 

   * The :mod:`lattice` module
   * The :mod:`~lattice.spacegroup` module
   * The :mod:`~lattice.surface` module


.. _bulk-crystal-section:

Common bulk crystals
====================

.. autofunction:: ase.structure.bulk

examples:

>>> from ase.structure import bulk
>>> a1 = bulk('Cu', 'fcc', a=3.6)
>>> a2 = bulk('Cu', 'fcc', a=3.6, orthorhombic=True)
>>> a3 = bulk('Cu', 'fcc', a=3.6, cubic=True)
>>> a1.cell
array([[ 0. ,  1.8,  1.8],
       [ 1.8,  0. ,  1.8],
       [ 1.8,  1.8,  0. ]])
>>> a2.cell
array([[ 2.54558441,  0.        ,  0.        ],
       [ 0.        ,  2.54558441,  0.        ],
       [ 0.        ,  0.        ,  3.6       ]])
>>> a3.cell
array([[ 3.6,  0. ,  0. ],
       [ 0. ,  3.6,  0. ],
       [ 0. ,  0. ,  3.6]])

|a1| |a2| |a3|

.. |a1| image:: a1.png
.. |a2| image:: a2.png
.. |a3| image:: a3.png


.. _nanotubes-section:

Nanotubes
=========

.. autofunction:: ase.structure.nanotube

examples:

>>> from ase.structure import nanotube
>>> cnt1 = nanotube(6, 0, length=4)
>>> cnt2 = nanotube(3, 3, length=6, bond=1.4, symbol='Si')

|cnt1| |cnt2|

.. |cnt1| image:: cnt1.png
.. |cnt2| image:: cnt2.png

.. _nanoribbons-section:

Graphene nanoribbons
====================

.. autofunction:: ase.structure.graphene_nanoribbon

examples:

>>> from ase.structure import graphene_nanoribbon
>>> gnr1 = graphene_nanoribbon(3, 4, type='armchair')
>>> gnr2 = graphene_nanoribbon(2, 6, type='zigzag', saturated=True,
>>>                             C_H=1.1, C_C=1.4, vacc=6.0, 
>>>                            magnetic=True,initial_mag=1.12)                     

|gnr1| |gnr2|

.. |gnr1| image:: gnr1.png
.. |gnr2| image:: gnr2.png

Special points in the Brillouin zone
====================================

You can find the psecial points in the Brillouin zone:

>>> from ase.structure import bulk
>>> from ase.dft.kpoints import ibz_points, get_bandpath
>>> si = bulk('Si', 'diamond', a=5.459)
>>> points = ibz_points['fcc']
>>> G = points['Gamma']
>>> X = points['X']
>>> W = points['W']
>>> K = points['K']
>>> L = points['L']
>>> kpts, x, X = get_bandpath([W, L, G, X, W, K], si.cell)
>>> print len(kpts), len(x), len(X)
50 50 6
