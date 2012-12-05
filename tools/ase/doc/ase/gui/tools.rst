.. module:: tools

=====
Tools
=====

Graphs
------
Allows to graph different quantities for a given trajectory. A 'save' button also gives the opportunity to save the data to file.

This example plots the maximal force for each image i and could help in investigating the convergence properties for relaxations:

::

  i, e-min(E), fmax

These are the symbols that can be used:

================ ==================================================
 Symbol	         Interpretation
================ ==================================================
e                total energy
epot        	 potential energy
ekin        	 kinetic energy
fmax        	 maximum force
fave        	 average force
d(n1,n2)    	 distance between two atoms
R[n,0-2]    	 position of atom number n
i           	 current image number
E[i]        	 energy of image number i
F[n,0-2]    	 force on atom number n
M[n]        	 magnetic moment of atom number n
A[0-2,0-2]  	 unit-cell basis vectors 
s           	 path length
a(n1,n2,n3) 	 tangle between atoms n1, n2 and n3, centered on n2
dih(n1,n2,n3,n4) dihedral angle between n1, n2, n3, and n4
T	         temperature (requires velocity)
================ ==================================================

Movie
-----

Allows to play the current trajectory as a movie using a number of
different settings. Default duration is 5 s. 

Expert mode
-----------
Python interface to all ag functions, with numerous extra commands
defined that help to modify and visualize a system. The commands come
in two flavors, the first is interpreted on an atom-by-atom basis
(e.g. operates on position, color, etc) and the second is based on the
entire frame. The flavor of a given line is determined from the first
command. Note that the frame-based commands can be used in atom-based
operations, but not vice versa. See below for some examples. 

Regular python syntax applies to the commands and ``numpy`` has been
imported as ``np``.

Two buttons allow to reduce the operation to a given frame:

======================== =============================================
Only selected atoms (sa) Restricts operation only to the selected
     	      	    	 atoms. The text command ``sa`` activates or
			 deactivates this button. 
Only current frame (cf)  Restricts operation only to the current
     	     	   	 frame, useful for long trajectories. The text
			 command ``cf`` activates or deactivates this
			 button.  
======================== =============================================

List of atom-based commands with a few examples:

================ ===================================================== 
Command	         Interpretation
================ ===================================================== 
x,y,z            Cartesian coordinates. Example:
		 ``x += A[0][0]``
r,g,b		 Color components, invoking the expert mode changes
		 the color mode to manual and allows to address all
		 colors individually. Example:
		 ``r = (z-min(R[:,2]))/(max(R[:,2])-min(R[:,2]))``
rad		 atomic display radius
s		 Boolean to control the selection of an atom. Example:
		 ``s = Z == 6 and x > 5`` or ``s = d == False``
f		 force on an atom
Z		 atomic number
m		 magnetic moment
d		 dynamic, e.g. d = False fixes an atom
================ =====================================================

List of frame-based and global commands and global objects with
examples: 

================ =====================================================
Command	         Interpretation
================ =====================================================
e		 total energy
fmax 		 maximal force
A		 unit cell
E		 total energy array of all frames. Example:
		 ``e-min(E)`` 
F		 all forces in one frame
M		 all magnetic moments
R		 all atomic positions
S		 boolean array of the entire selection
D		 boolean array of dynamic atoms (False = atom is fixed)
del S		 deletes current selection
sa,cf		 toggles the selected-atoms-only or the
		 current-frame-only buttons 
frame		 provides and edits the frame number in a trajectory
center		 centers system in its unit cell
cov		 array of original covalent radii
self		 expert mode window
gui		 ag GUI object, this controls the entire ag session
img		 ag images object, all physical data is stored here
================ =====================================================

To save data between commands, one has to assign variables to parent
objects in the gui, e.g. via ``self.temp_var =
R-img.P[0,:]``. DISCLAIMER: Doing so might risk the functionality
ofthe entire ag session if you accidentally overwrite basic
functionality of the gui or the image objects stored within.  

Finally, recurring selections of commands collected as scripts can be
executed as 

::
	
	exec <filename>

If the file in question is saved in the directory ``~/.ase/`` then
just the filename will also do. 

Constraints
-----------
Allows to set (or remove) constraints based on the currently selected atoms. 

Render scene
------------
Graphical interface to the ASE povray interface, ideally it requires
that povray is installed on your computer to function, but it also can
be used just to export the complete set of povray files. 

The texture of each atom is adjustable: The default texture is applied
to all atoms, but then additional textures can be defined based on
selections (``Create new texture from current selection``). These can
be obtained either from selecting atoms by hand or by defining a
selection with a boolean expression, for example ``Z==6 and x>5 and
y<0`` will select all carbons with coordinates x>5 and y<0. The
available commands are listed in the ``Help on textures`` window. 

A movie-making mode (``render all N frames``) is also available. After
rendering, the frames can be stitched together using the ``convert``
unix program e.g. 

::

	localhost:doc hanke$ convert -delay 4.17 temp.*.png temp.gif

For this particular application it might be a good idea to use a white
background instead of the default transparent option. 

Move atoms
----------
Allows selected atoms to be moved using the arrow keys. The direction
is always parallel to the plane of the screen. Two possible movements
are available: Just pressing the arrow keys will move by 0.1
Angstrom, ``shift`` + arrow keys will move by 0.01 Angstrom. 

Rotate atoms
------------
Allows sets of atoms to be rotated using the arrow keys. Different
rotation modes are available depending on the number of selected
atoms. Again, two modes are available. Just the arrow keys will rotate
by 2.5 degrees, and ``shift`` + arrow keys will rotate by 0.5 deg.  

============================== =======================================
number of atoms labeled	       rotation mode
============================== =======================================
0 atoms, 1, 3, 5 or more atoms uses the centre of mass of the atoms to
  			       be rotated as the rotation centre. 
2 atoms			       Defines the vector connecting the two
  			       atoms as rotation axis. 
4 atoms, selected sequentially Defines the vector connecting the two
  	 	  	       atoms as rotation axis. This mode has
			       the advantage that the dihedral angle
			       is measured at the same time, thus
			       allowing one to monitor the degree of
			       rotation. 
============================== =======================================

Orient atoms
------------
stub


NEB
---
stub

Bulk Modulus
------------
stub
