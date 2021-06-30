"""
Clean all test directories.
"""
import os

from tools.runner.clean import clean_tests
from tools.constants import settings

tests = next(os.walk(settings.test_farm))[1]
clean_tests(settings.test_farm, 
            tests,
            settings.run_dir,
            settings.ref_dir,
            settings.ignored_output)
