# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

version_base = '3.6.1'

try:
    from ase.svnversion import svnversion
except ImportError:
    version = version_base
else:
    version = version_base + '.' + svnversion
