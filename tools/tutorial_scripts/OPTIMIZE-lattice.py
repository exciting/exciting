#!/usr/bin/env python
#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#
#%!%!% -------------------------------- OPTIMIZE-lattice.py -------------------------------- %!%!%#
#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#
#
# AUTHOR:
#    Rostam Golesorkhtabar
#    r.golesorkhtabar@gmail.com
# 
# DATE:
#    Wed Jan 01 00:00:00 2014
#
# SYNTAX:
#    python OPTIMIZE-lattice.py
#           OPTIMIZE-lattice.py
# 
# EXPLANATION:
# 
#__________________________________________________________________________________________________

import os
import sys
import glob

#%!%!%--- Reading the optimization INFO file ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%

INFOlist = sorted(glob.glob('INFO_*'))
if (len(INFOlist) == 0):
    os.system('OPTIMIZE-setup.py')
elif(len(INFOlist) == 1):
    os.system('OPTIMIZE-analyze.py')
else:
    sys.exit('\n.... Oops ERROR: There are more than one INFO file in this directory !?!?!?\n')