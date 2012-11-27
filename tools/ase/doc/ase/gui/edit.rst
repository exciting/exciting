.. module:: edit

====
Edit
====

Add atoms
---------

Allows to add single atoms or a group of atoms to an existing atoms
object. If the description is an atom or a known molecule from the g1,
g2, or g3 set (e.g. CH4), then the structure from the data molecule is
used. In addition, a molecule can also be imported from file via the
``load molecule`` button. 

The specified position can either be absolute, or determined
automatically via 

::
   
   auto+<dist>

where auto is the centre of mass of the currently selected atoms, and
``<dist>`` is the distance toward the viewing plane. 

The molecule-to-be is rotated into the current viewing plane before
addition into the system. Two options exist for chosing the origin
within the new atoms, it can be either the centre of mass or the
origin of the loaded geometry. 

Modify
------

Menu to allow modification of the atomic symbol, an attached tag, or
its magnetic moment.

Copy/Paste
----------

Allows to copy parts of an existing system to a clipboard, and pastes
via the same infrastructure as the ``Add atoms`` functionality. Note
that the on-screen orientation of the pasted atoms will be the same as
at the time of copying. Note also that, by default, the origin of the
pasted system is taken to be the atom furthest away from the viewing
point. 

