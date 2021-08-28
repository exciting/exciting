"""
Clean all test directories.
Compatible with python 2 and 3
"""
import os
import shutil

run_dirs = [d[0] for d in os.walk('test_farm') if os.path.basename(d[0]) == 'run']

for run_dir in run_dirs:
    shutil.rmtree(run_dir)

print('Test directories cleaned.')
