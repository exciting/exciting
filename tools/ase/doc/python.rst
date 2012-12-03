.. _python_info:

---------------
What is Python?
---------------

This section will give a very brief introduction to the Python
language.

.. seealso::

   * The Python_ home page.
   * Python Recipes_.
   * Try a `Python quick reference card`_ or a `different reference card`_.


.. _Recipes: http://code.activestate.com/recipes/langs/python
.. _Python quick reference card: http://www.limsi.fr/Individu/pointal/python/pqrc
.. _different reference card: http://rgruet.free.fr/
.. _Python: http://www.python.org


Executing Python code
---------------------

You can execute Python code interactively by starting the interpreter
like this::

  $ python
  Python 2.5.1 (r251:54863, Mar  7 2008, 04:10:12) 
  [GCC 4.1.3 20070929 (prerelease) (Ubuntu 4.1.2-16ubuntu2)] on linux2
  Type "help", "copyright", "credits" or "license" for more information.
  >>> print 'hello'
  hello

You can also put the ``print 'hello'`` line in a file (``hello.py``)
and execute it as a Python script::

  $ python hello.py
  hello

Or like this::

  $ python -i hello.py
  hello
  >>> print 'hi!'
  hi!

Finally, you can put ``#!/usr/bin/env python`` in the first line of
the ``hello.py`` file, make it executable (``chmod +x hello.py``) and
execute it like any other executable.

.. tip::

   For interactive Python sessions, it is very convenient to have a
   personal ``.pythonrc`` file::

     import rlcompleter
     import readline
     readline.parse_and_bind("tab: complete")
     from ase import *

   and point the :envvar:`PYTHONSTARTUP` environment variable at it (see
   here_ for details).


   .. _here: http://www.python.org/doc/current/lib/module-rlcompleter.html


.. tip::

   For an even better interactive experience, use ipython_.

   .. _ipython: http://ipython.scipy.org



Types
-----

Python has the following predefined types:

===========  =====================  ==========================
type         description            example
===========  =====================  ==========================
``bool``     boolean                ``False``
``int``       integer                ``117``
``float``    floating point number  ``1.78``
``complex``  complex number         ``0.5 + 2.0j``
``str``      string                 ``'abc'``
``tuple``    tuple                  ``(1, 'hmm', 2.0)``
``list``     list                   ``[1, 'hmm', 2.0]``
``dict``     dictionary             ``{'a': 7.0, 23: True}``
``file``     file                   ``open('stuff.dat', 'w')``
===========  =====================  ==========================

A ``dict`` object is mapping from keys to values:

>>> d = {'s': 0, 'p': 1}
>>> d['d'] = 2
>>> d
{'p': 1, 's': 0, 'd': 2}
>>> d['p']
1

A ``list`` object is an ordered collection of arbitrary objects:

>>> l = [1, ('gg', 7), 'hmm']
>>> l[1]
('gg', 7)
>>> l.append(1.2)
>>> l[-2]
'hmm'

A ``tuple`` behaves like a ``list`` - except that it can't be modified
inplace.  Objects of types ``list`` and ``dict`` are *mutable* - all
the other types listed in the table are *immutable*, which means that
once an object has been created, it can not change.

.. note::

   List and dictionary objects *can* change.  Variables in
   Python are references to objects.  This is demonstrated here:

   >>> a = ['q', 'w']
   >>> b = a
   >>> a.append('e')
   >>> a
   ['q', 'w', 'e']
   >>> b
   ['q', 'w', 'e']


.. note::

   Another very important type is the :term:`ndarray` type described
   here: :ref:`numpy`.



Loops
-----

A loop in Python can be done like this:

>>> things = ['a', 7]
>>> for x in things:
...     print x
...
a
7

The ``things`` object could be any sequence.  Strings, tuples, lists,
dictionaries, ndarrays and files are sequences.


Functions and classes
---------------------

A function is defined like this:

>>> def f(x, m=2, n=1):
...     y =  x + n
...     return y**m

Here ``f`` is a function, ``x`` is an argument, ``m`` and ``n`` are keywords with default values ``2`` and ``1`` and ``y`` is a variable.

A :term:`class` is defined like this:

>>> class A:
...     def __init__(self, b):
...         self.c = b
...     def m(self):
...         return f(self.c, n=0)

The ``__init__()`` function is called a :term:`constructor`.  You can think
of a class as a template for creating user defined objects.

In the class ``A`` ``__init__`` is a constructor, ``c`` is an attribute and ``m`` is a method.

>>> a = A(7)
>>> a.c
7
>>> a.m()
49
>>> g = a.m
>>> g()
49

Here we make an :term:`instance`/object ``a`` of type ``A`` and ``g`` is a :term:`method` bound to ``a``.


Importing modules
-----------------

If you put the definitions of the function ``f`` and the class ``C``
in a file ``stuff.py``, then you can use that code from another piece
of code::

  from stuff import f, C
  print f(1, 2)
  print C(1).m(2)

or::

  import stuff
  print stuff.f(1, 2)
  print stuff.C(1).m(2)

or::

  import stuff as st
  print st.f(1, 2)
  print st.C(1).m(2)


Python will look for ``stuff.py`` in these directories:

1) current working directory
2) directories listed in your :envvar:`PYTHONPATH`
3) Python's own system directory (typically :file:`/usr/lib/python2.5`)

and import the first one found.
