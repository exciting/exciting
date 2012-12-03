import sys

msg = "\nThe matplotlib python module is missing or not installed properly.\n"
msg += "Is the PYTHONPATH environment variable set correctly?\n"
msg += "Please verify your installation by running on the command line:\n"
msg += "python -c 'import matplotlib'\n"
msg += "\n"
msg += "This module is optional and required in order to use "
msg += "ASE's simple GUI (ag).\n"
msg += "If you don't wish to use ag ignore this error, otherwise\n"
msg += "please install the package using "
msg += "your distribution package manager, i.e.:\n"
msg += "\n"
msg += "  Debian/Ubuntu: sudo apt-get python-matplotlib\n"
msg += "\n"
msg += "  OpenSUSE: yast -i python-matplotlib\n"
msg += "\n"
msg += "  Red Hat/Fedora: yum install python-matplotlib\n"
msg += "\n"
msg += "or perform manual installation, preferably as non-root user,\n"
msg += "following http://matplotlib.sourceforge.net/users/installing.html."

if locals().get('display'):
    try:
        import matplotlib
    except ImportError:
        print >> sys.stderr, msg
        raise
