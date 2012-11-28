#!/usr/bin/python
"""Simple shallow test of the CASTEP interface"""

import os
import shutil
import tempfile
import traceback

from ase.test import NotAvailable

# check if CASTEP_COMMAND is set a environment variable

if not os.environ.has_key('CASTEP_COMMAND'):
    print("WARNING: Environment variable CASTEP_COMMAND is not set")
    print("Will set CASTEP_COMMAND  = castep for the sake of this test")
    print("Please change it if this does not run castep in your environment")
    os.environ['CASTEP_COMMAND'] = 'castep'


if not (os.system('which %s' % os.environ['CASTEP_COMMAND']) == 0):
    raise NotAvailable("""Could not find CASTEP. If you have it
                          installed make sure, you set the CASTEP_COMMAND
                          environment variable correctly""")


# check if we can import everything
ase_castep_dir = "ase"

try:
    castep_calc = __import__(ase_castep_dir + ".calculators.castep", globals(), locals(), ["Castep", "CastepParam", "create_castep_keywords"])
    Castep = castep_calc.Castep
    CastepParam = castep_calc.CastepParam
    create_castep_keywords = castep_calc.create_castep_keywords

except Exception, e:
    traceback.print_exc()
    print(e)
    assert False, 'Castep calculator module could not be loaded'

try:
    __import__(ase_castep_dir + ".io.castep")
except Exception, e:
    assert False, 'Castep io module could not be loaded'


tmp_dir = tempfile.mkdtemp()
cwd = os.getcwd()

from ase.calculators.castep import Castep

try:
    c = Castep(directory=tmp_dir, label='test_label')
except Exception, e:
    traceback.print_exc()
    print(e)
    assert False, 'Could not instantiate castep calculator'


try:
    c.xc_functional = 'PBE'
except Exception, e:
    traceback.print_exc()
    print(e)
    assert False, 'Setting xc_functional  failed'

import ase.lattice.cubic
lattice = ase.lattice.cubic.BodyCenteredCubic('Li' )

print('For the sake of evaluating this test, warnings')
print('about auto-generating pseudo-potentials are')
print('normal behavior and can be safely ignored')

try:
    lattice.set_calculator(c)
except Exception, e:
    traceback.print_exc()
    print(e)
    assert False, 'Setting the calculator %s failed' % c



try:
    create_castep_keywords(
        castep_command=os.environ['CASTEP_COMMAND'],
        path=tmp_dir,
        fetch_only=20)
except Exception, e:
    traceback.print_exc()
    print(e)
    assert  False, "Cannot create castep_keywords, this usually means a  bug"\
    + " in the interface or the castep binary cannot be called"


param_fn = os.path.join(tmp_dir, 'myParam.param')
param = open(param_fn,'w')
param.write('XC_FUNCTIONAL : PBE #comment\n')
param.write('XC_FUNCTIONAL : PBE #comment\n')
param.write('#comment\n')
param.write('CUT_OFF_ENERGY : 450.\n')
param.close()
try:
    c.merge_param(param_fn)
except Exception, e:
    traceback.print_exc()
    print(e)
    assert False,"Error in merge_param_filename, go figure"


# check if the CastepOpt, CastepCell comparison mechanism works

p1 = CastepParam()
p2 = CastepParam()
assert p1._options == p2._options, "Print two newly created CastepParams are not the same"

p1._options['xc_functional'].value = 'PBE'
p1.xc_functional = 'PBE'

assert not p1._options == p2._options, "Changed one CastepParam, but the still look the same"

assert c.calculation_required(lattice), 'Calculator does not fetch that a calculation is required'

if not c.dryrun_ok():
    print(c._error)
    assert False, "Dryrun_ok does not work, where it should"
else:
    print("Dryrun is ok")

c.prepare_input_files(lattice)

os.chdir(cwd)
shutil.rmtree(tmp_dir)


print("Test finished without errors")