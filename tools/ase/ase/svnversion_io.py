# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

from os import path
import sys

ON_POSIX = 'posix' in sys.builtin_module_names

try:
    from subprocess import Popen
except ImportError:
    from os import popen3
else:
    def popen3(cmd):
        from subprocess import PIPE
        p = Popen(cmd, shell=True, close_fds=ON_POSIX,
                  stdin=PIPE, stdout=PIPE, stderr=PIPE)
        return p.stdin, p.stdout, p.stderr

def write_svnversion(svnversion, dir):
    svnversionfile = path.join(dir, 'svnversion.py')
    f = open(svnversionfile,'w')
    f.write('svnversion = "%s"\n' % svnversion)
    f.close()
    print 'svnversion = ' +svnversion+' written to '+svnversionfile
    # assert svn:ignore property if the installation is under svn control
    # because svnversion.py has to be ignored by svn!
    cmd = popen3('svn propset svn:ignore svnversion.py '+dir)[1]
    output = cmd.read()
    cmd.close()

def get_svnversion_from_svn(dir):
    # try to get the last svn version number from svnversion
    cmd = popen3('svnversion -n '+dir)[1] # assert we are in the project dir
    output = cmd.read().strip()
    cmd.close()
    if output.startswith('exported'):
        # we build from exported source (e.g. rpmbuild)
        output = None
    return output

svnversion = get_svnversion_from_svn(dir='ase')
if svnversion:
    write_svnversion(svnversion, dir='ase')
