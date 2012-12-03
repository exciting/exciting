.. module:: test

================
Testing the code
================

All additions and modifications to ASE should be tested.

.. index:: testase

Test scripts should be put in the :trac:`ase/test` directory.  The
scripts in this directory may be run by using the script
:svn:`tools/testase` or by using the function:

.. function:: test.test(verbosity=1, dir=None)
    
    Runs the test scripts in :trac:`ase/test`.


.. important::

  When you fix a bug, add a test to the test suite checking that it is
  truly fixed.  Bugs sometimes come back, do not give it a second
  chance!


How to fail successfully
========================

The test suite provided by :func:`test.test` automatically runs all test
scripts in the :trac:`ase/test` directory and summarizes the results.

.. note::

  Test scripts are run from within Python using the :func:`execfile` function.
  Among other things, this provides the test scripts with an specialized global
  :term:`namespace`, which means they may fail or behave differently if you try
  to run them directly e.g. using :command:`python testscript.py`.

If a test script causes an exception to be thrown, or otherwise terminates
in an unexpected way, it will show up in this summary. This is the most
effective way of raising awareness about emerging conflicts and bugs during
the development cycle of the :ref:`latest revision <Latest_development_release>`.


Remember, great tests should serve a dual purpose:

**Working interface**
    To ensure that the :term:`class`'es and :term:`method`'s in ASE are
    functional and provide the expected interface. Empirically speaking, code
    which is not covered by a test script tends to stop working over time.

**Replicable results**
    Even if a calculation makes it to the end without crashing, you can never
    be too sure that the numerical results are consistent. Don't just assume 
    they are, :func:`assert` it!

.. function:: assert(expression)
    
    Raises an `AssertionError` if the `expression` does not evaluate to `True`.

Example::

  from ase import molecule
  atoms = molecule('C60')
  atoms.center(vacuum=4.0)
  result = atoms.get_positions().mean(axis=0)
  expected = 0.5*atoms.get_cell().diagonal()
  tolerance = 1e-4
  assert (abs(result - expected) < tolerance).all()


Using functions to repeat calculations with different parameters::

  def test(parameter):
      # setup atoms here...
      atoms.set_something(parameter)
      # calculations here...
      assert everything_is_going_to_be_alright

  if __name__ in ['__main__', '__builtin__']:
      test(0.1)
      test(0.3)
      test(0.7)
          
.. important::

  Unlike normally, the module `__name__` will be set to ``'__builtin__'``
  when a test script is run by the test suite.

