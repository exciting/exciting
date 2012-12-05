.. _versioncontrol:

===========================
Using version control (SVN)
===========================

The version control system used in ASE development is subversion. A thorough
subversion manual can be found at http://svnbook.red-bean.com/, here
is a brief overview of the most basic features needed when developing ASE.

* perform svn checkout of ase.

  Place it, for example, in ``ase-svn`` directory::
 
     cd
     svn checkout https://svn.fysik.dtu.dk/projects/ase/trunk ase-svn

  This retrieves the code tree from the subversion repository.
  Prepend :envvar:`PYTHONPATH` and :envvar:`PATH` environment variables
  as described at :ref:`installation`.

* Updating the working copy of the code (in the directory ``ase-svn``)::

    svn update

  After each (important) update, remove ``ase/svnrevision.py*`` files,
  and run::

    python setup.py sdist

  to keep the ``ase/svnrevision.py`` file up-to-date.

* Checking the status of the working copy (in the directory ``ase-svn``)::

    svn stat

  The status about the files which are not in version control can be
  surpassed with the ``-q`` flag, and the status with respect to latest
  additions in server can be checked with the ``-u`` flag.

* Committing the changes to the repository

  Before sending the changes in the working copy to the repository, working
  copy should be updated. After that, the changes can be send with::

    svn commit -m "Message to describe the committed changes"

  If the ``-m`` option is omitted, an editor is opened for writing the
  log message.

* Adding files or directories to version control::

    svn add filename

  If ``filename`` is directory, also all the files within the
  directory are added. Note that ``svn commit`` is required before the
  new files are actually included in version control.

* ASE documentation resides under ``doc`` directory,
  and is also under subversion control.
  See more at :ref:`writing_documentation_ase`.
