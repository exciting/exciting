.. _writing_documentation_ase:

=====================
Writing documentation
=====================

We use the Sphinx_ tool to generate the documentation (both HTML
and PDF_).  The documentation is stored in SVN as text files in the
:trac:`doc` directory using the reStructuredText_ markup language.

.. _reStructuredText: http://docutils.sf.net/rst.html
.. _Sphinx: http://sphinx.pocoo.org
.. _PDF: ../ase-manual.pdf


Installing Docutils and Sphinx
==============================

The reStructuredText_ parser that Sphinx needs, is part of the Docutils_
project.  So, we need to install docutils and sphinx (version>= 0.5).

.. _Docutils: http://docutils.sf.net


Other requirements
==================

When building the documentation, a number of png-files are generated.
For that to work, you need the following installed:

* matplotlib
* povray
* dvipng
* pdflatex
* bibtex
* AUCTex
* convert (ImageMagick)

.. _using_sphinx:

Using Sphinx
============

.. highlight:: bash


First, you should take a look at the documentation for Sphinx_ and
reStructuredText_.

If you don't already have your own copy of the ASE package, then get
the :ref:`latest_development_release` and install it.

Then :command:`cd` to the :file:`doc` directory and build the html-pages::

  $ cd ~/ase/doc
  $ sphinx-build . _build

.. Note::

   Make sure that you build the Sphinx documentation using the corresponding ASE version
   by setting the environment variables :envvar:`$PYTHONPATH` and :envvar:`$PATH`.

Make your changes to the ``.rst`` files, run the
:command:`sphinx-build` command again, check the results and if things
looks ok, commit::

  $ emacs index.rst
  $ sphinx-build . _build
  $ firefox _build/index.html
  $ svn ci -m "..." index.rst

To build a pdf-file, you do this::

  $ sphinx-build -b latex . _build
  $ cd _build
  $ make ase-manual.pdf



Extensions to Sphinx
====================

.. highlight:: rest

We have a couple of extensions to Sphinx:

**:mol:**

   Use ``:mol:`CH_3OH``` to get :mol:`CH_3OH`.

**:svn:**

   A role for creating a link to a file in SVN.  If you write
   ``:svn:`ase/atoms.py```, you
   will get: :svn:`ase/atoms.py`.

**:trac:**

   A role for creating a link to a file in Trac.  If you write
   ``:trac:`ase/atoms.py```, you
   will get: :trac:`ase/atoms.py`.

**:epydoc:**

   A role for creating a link to the API-documentation generated with
   epydoc_.  If you
   write ``:epydoc:`ase.atoms.Atoms```, you will get:
   :epydoc:`ase.atoms.Atoms`.

**:math:**

   This role is for inline LaTeX-style math.  Example:
   ``:math:`\sin(x_n^2)``` gives you :math:`\sin(x_n^2)`.

**.. math::**

   Write displayed LaTeX-style math.  Example::

     .. math:: \frac{1}{1+x^2}

   gives you:

   .. math:: \frac{1}{1+x^2}


If you add the line ``.. default-role:: math``, then you can leave out
the ``:math:`` part like here: ```\sin(x_n^2)```.

The implementation of these roles is here: :svn:`doc/ext.py`.  Our
custom, obsolete, implementation of the math role and directive is
here: :svn:`doc/mathpng.py`.  With sphinx >= 0.5 please use
``sphinx.ext.pngmath``.


.. _epydoc:  http://epydoc.sf.net



reStructedText in emacs
=======================

.. highlight:: common-lisp

For people using emacs, the `reStructuredText extension`_ is highly
recommended. The intallation procedure is described in the top of the
file, but for most people, it is enough to place it in your emacs
load-path (typically ``.emacs.d/``) and add the lines::

  (add-to-list 'load-path "~/.emacs.d")
  (require 'rst)

somewhere in your ``.emacs`` file.

To make the mode auto load for relevant file extension, you can write
something like::

  (setq auto-mode-alist
        (append '(("\\.rst$" . rst-mode)
                  ("\\.rest$" . rst-mode)) auto-mode-alist))

In your ``.emacs`` file.

.. _reStructuredText extension: http://docutils.sourceforge.net/tools/editors/emacs/rst.el

Updating Sphinx
===============

Starting a new project with sphinx requires an initial configuration.
This is achieved by running :command:`sphinx-quickstart`.
When updating from a very old sphinx you may consider
generating new configuration files and updating the old files accordingly.

**Note** that the current project is configured up-to-date,
so if you are "simply" writing the documentation
you **must** skip the :command:`sphinx-quickstart` step
and focus on :ref:`using_sphinx`.

Here is how do you setup the GPAW project with sphinx:

 - :command:`cd` to the :file:`doc` directory,

 - run :command:`sphinx-quickstart`
   and answer the questions (example given for GPAW)::

    > Root path for the documentation [.]:

    > Separate source and build directories (y/N) [n]:

    > Name prefix for templates and static dir [.]: _

    > Project name: GPAW
    > Author name(s): 2008, CAMd et al.
  
    > Project version: 0.5
    > Project release [0.5]:

    > Source file suffix [.rst]:

    > Name of your master document (without suffix) [index]: contents

    > autodoc: automatically insert docstrings from modules (y/N) [n]: y
    > doctest: automatically test code snippets in doctest blocks (y/N) [n]:
    > intersphinx: link between Sphinx documentation of different projects (y/N) [n]: y

   This will create :file:`doc/conf.py` and :file:`doc/contents.rst`.
   Both these files need to be edited further
   (:file:`doc/conf.py` may for example include
   options for ``sphinx.ext.pngmath``)

