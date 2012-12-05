.. _download_and_install:

=================================
Installation of ASE: Requirements
=================================

The following packages are required for basic ASE functionality:

1) Python_ 2.4 - 2.7.
2) NumPy_.

.. _Python: http://www.python.org
.. _NumPy: http://www.scipy.org/NumPy

It is highly recommended (but not required) to install also these two:

3) matplotlib_.
4) pygtk_.

Matplotlib is needed for :mod:`writing png and eps files <io>`, and
both packages are needed for ASE's simple GUI (called **ag**, see :mod:`gui`).
Some of these packages may already be installed on your system. OSX
users, please see below for additional details to install the
prerequisite libraries. 

.. _matplotlib: http://matplotlib.sourceforge.net
.. _pygtk: http://www.pygtk.org

========
Download
========

.. highlight:: bash

.. _latest_stable_release:

Latest stable release
=====================

The latest stable release can be obtained from ``svn`` or as a ``tarball``.

.. note::

   The recommended installation path is :envvar:`$HOME`.

When using svn please set the following variable:

- bash::

   export ASE_TAGS=https://svn.fysik.dtu.dk/projects/ase/tags/

- csh/tcsh::

   setenv ASE_TAGS https://svn.fysik.dtu.dk/projects/ase/tags/

======= =========== ============================================ =============================
Release Date        Retrieve as svn checkout                     Retrieve as tarball
======= =========== ============================================ =============================
 3.6.0_ Feb 24 2012 ``svn co -r 2515 $ASE_TAGS/3.6.0 ase-3.6.0`` python-ase-3.6.0.2515.tar.gz_
 3.5.1_ May 24 2011 ``svn co -r 2175 $ASE_TAGS/3.5.1 ase-3.5.1`` python-ase-3.5.1.2175.tar.gz_
 3.4.1_ Aug 11 2010 ``svn co -r 1765 $ASE_TAGS/3.4.1 ase-3.4.1`` python-ase-3.4.1.1765.tar.gz_
 3.4.0_ Apr 23 2010 ``svn co -r 1574 $ASE_TAGS/3.4.0 ase-3.4.0`` python-ase-3.4.0.1574.tar.gz_
 3.3.1_ Jan 20 2010 ``svn co -r 1390 $ASE_TAGS/3.3.1 ase-3.3.1`` python-ase-3.3.1.1390.tar.gz_
 3.2.0_ Sep 4 2009  ``svn co -r 1121 $ASE_TAGS/3.2.0 ase-3.2.0`` python-ase-3.2.0.1121.tar.gz_
 3.1.0_ Mar 27 2009 ``svn co -r 846 $ASE_TAGS/3.1.0 ase-3.1.0``  python-ase-3.1.0.846.tar.gz_
 3.0.0_ Nov 13 2008 ``svn co -r 657 $ASE_TAGS/3.0.0 ase-3.0.0``  python-ase-3.0.0.657.tar.gz_
======= =========== ============================================ =============================

.. _3.6.0:
    https://trac.fysik.dtu.dk/projects/ase/browser/tags/3.6.0

.. _python-ase-3.6.0.2515.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.6.0.2515.tar.gz

.. _3.5.1:
    https://trac.fysik.dtu.dk/projects/ase/browser/tags/3.5.1

.. _python-ase-3.5.1.2175.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.5.1.2175.tar.gz

.. _3.4.1:
    https://trac.fysik.dtu.dk/projects/ase/browser/tags/3.4.1

.. _python-ase-3.4.1.1765.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.4.1.1765.tar.gz

.. _3.4.0:
    https://trac.fysik.dtu.dk/projects/ase/browser/tags/3.4.0

.. _python-ase-3.4.0.1574.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.4.0.1574.tar.gz

.. _3.3.1:
    https://trac.fysik.dtu.dk/projects/ase/browser/tags/3.3.1

.. _python-ase-3.3.1.1390.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.3.1.1390.tar.gz

.. _3.2.0:
    https://trac.fysik.dtu.dk/projects/ase/browser/tags/3.2.0

.. _python-ase-3.2.0.1121.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.2.0.1121.tar.gz

.. _3.1.0:
    https://trac.fysik.dtu.dk/projects/ase/browser/tags/3.1.0

.. _python-ase-3.1.0.846.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.1.0.846.tar.gz

.. _3.0.0:
    https://trac.fysik.dtu.dk/projects/ase/browser/tags/3.0.0

.. _python-ase-3.0.0.657.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.0.0.657.tar.gz



.. _latest_development_release:

Latest development release
==========================

The latest revision can be obtained like this::

  $ svn checkout https://svn.fysik.dtu.dk/projects/ase/trunk ase

or from the daily snapshot: `<python-ase-snapshot.tar.gz>`_.

.. note::

   The recommended checkout path is :envvar:`$HOME`.

.. _installation:

============
Installation
============

After downloading create the link to the requested version, e.g.:

- if retrieved from ``svn``::

   $ cd $HOME
   $ ln -s ase-3.5.1 ase
    
- if retrieved as ``tarball``::

   $ cd $HOME
   $ tar zxf python-ase-3.5.1.2175.tar.gz
   $ ln -s python-ase-3.5.1.2175 ase

It is sufficient to
put the directory :file:`$HOME/ase` in your :envvar:`PYTHONPATH`
environment variable, and the directory :file:`$HOME/ase/tools` in
your :envvar:`PATH` environment variable.  Do this permanently in
your :file:`~/.bashrc` file::

  export PYTHONPATH=$HOME/ase:$PYTHONPATH
  export PATH=$HOME/ase/tools:$PATH

or your :file:`~/.cshrc` file::

  setenv PYTHONPATH ${HOME}/ase:${PYTHONPATH}
  setenv PATH ${HOME}/ase/tools:${PATH}

Instead of :envvar:`HOME`, you may use any other directory.

.. index:: test

If you have root-permissions, you can install ASE system-wide::

  $ cd ase
  $ sudo python setup.py install

.. _running_tests:

Installation of prerequisites for OSX
=====================================

For Apple users, the MacPorts_ Project provides a straight-forward
route to obtain all necessary requirements. Unfortunately, MacPorts
does not install the `gtk` bindings to matplotlib_ by default, which
are required to open the GUI. To get all the ASE prerequisites for
python 2.7 in one single command anyway, install MacPorts and then
run::

  $ port install py27-matplotlib +gtk2

Use the `sudo` command if you have root access and if you require 
a system-wide install. Once finished, please follow the instructions
above to install ASE. 

.. _MacPorts: http://www.macports.org/

Run the tests
=============

Make sure that everything works by running the :mod:`test
suite <test>`.  This will create many files, so run the tests in a new
directory (preferably using bash)::

  $ bash
  $ mkdir /tmp/testase.$$; cd /tmp/testase.*
  $ python ~/ase/tools/testase.py 2>&1 | tee testase.log

.. note:: 

   In the development version of ASE, and in future stable versions,
   the test script is just named :file:`testase`.

.. note::

   The last test :trac:`ase/test/COCu111.py` requires closing
   the graphics windows to terminate the whole test-suite.

If any of the tests fail,
then please send us :file:`testase.log` (see :ref:`bugs`).

.. note::

   If matplotlib_ or pygtk_ is not installed, one of the tests will
   fail - avoid this with::

     $ testase.py --no-display

