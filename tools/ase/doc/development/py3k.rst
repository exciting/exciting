Python 3 strategy
=================

* One codebase for both 2 and 3.

* Use "print(...)" if possible::

    print 'bla bla'   # no
    print('bla bla')  # yes
    print 'bla bla:', x       # no
    print('bla bla: %s' % x)  # yes

* Don't do this: ``print >> f, ...``.  Use ``f.write(... + '\n')`` or
  ``ase.utils.prnt(..., file=f)``.

* More help here: http://packages.python.org/six/
