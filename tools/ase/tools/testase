#!/usr/bin/env python

from optparse import OptionParser

description = """Run ASE test suite.  ***WARNING***: This will leave a
large number of files in current working directory, so be sure to do
it in a new directory!"""

p = OptionParser(usage='%prog [OPTION]', description=description)
p.add_option('--no-display', action='store_true',
             help='do not open graphical windows')

opts, args = p.parse_args()
if len(args) != 0:
    raise p.error('Unexpected arguments: %s' % args)

from ase.test import test
test(2, display=not opts.no_display)
