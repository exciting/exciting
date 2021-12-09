"""
Check the installed python version. 
Raises exception when the version is lower than python 3.5
"""
import sys

if sys.version_info <= (3, 5):
    raise Exception('Running the test suite requires python 3.5. See test/README for more detailed information.')