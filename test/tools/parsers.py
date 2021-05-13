"""
Wrapper module to expose exciting parsers
"""
import sys
import os
import warnings


def install_excitingtools():
    """
    Install excitingtools to provide the exciting parsers
    """
    import subprocess

    if 'excitingtools' in sys.modules:
        return
    else:
        print("Running pip install for exciting_tools. "
              "See <EXCITINGROOT>/test/tools/parsers.py for more info")
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-e", "../tools/exciting_tools/"])
        # Restart script such that package is found following installation
        os.execv(sys.executable, [sys.executable] + sys.argv)

try:
    __import__('excitingtools')
except ModuleNotFoundError:
    # Give instructions to user rather than automatically installing:
    message = """excitingtools must be installed for the test suite to run.
To install excitingtools, from the test directory, type:

  pip3 install -e ../tools/exciting_tools

excitingtools can be uninstalled by typing:

  pip3 uninstall excitingtools  
    """
    warnings.warn(message)
finally:
    from excitingtools.parser.ErroneousFileError import ErroneousFileError
    from excitingtools.parser.initParser import parseInit, getInitFile
    from excitingtools.parser.parserChooser import parser_chooser
