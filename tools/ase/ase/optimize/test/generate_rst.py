import os
import re
dirlist = os.listdir('.')
name = '.*\.csv'
filterre = re.compile(name)
dirlist = filter(filterre.search, dirlist)
namelist = [d.strip('.csv') for d in dirlist]

f = open('testoptimize.rst', 'w')
f.write(
""".. _optimizer_tests:

===============
Optimizer tests
===============
This page shows benchmarks of optimizations done with our different optimizers.
Note that the iteration number (steps) is not the same as the number of force
evaluations. This is because some of the optimizers uses internal line searches
or similar.
"""
)

for name in namelist:
    lines = open(name + '.csv', 'r').read().split('\n')
    firstline = lines.pop(0)
    f.write(
        '\n' + 
        name + '\n' + \
        '=' * len(name) + '\n'
        'Calculator used: %s\n' % firstline.split(',')[-1] + \
        '\n' + \
        '=============== ===== ================= ========== ===============\n' + \
        'Optimizer       Steps Force evaluations Energy     Note           \n' + \
        '=============== ===== ================= ========== ===============\n'
    )
    for line in lines:
        if len(line):
            print line.split(',')
            f.write(
                '%-15s %5s %17s %10s %s\n' % tuple(line.split(','))
            )
    f.write(
        '=============== ===== ================= ========== ===============\n'
    )
f.close()
