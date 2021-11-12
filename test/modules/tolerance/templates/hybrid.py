"""
Module containing default tolerances for hybrid XC output.
This is essentially an extension of the ground state output for INFO.OUT
and a copy of all other tested ground state outputs.

Tolerances in default_tolerances cover files:
  * INFO.OUT
  * eigval.xml
  * atoms.xml
  * geometry.xml
  * evalcore.xml
"""
from tolerance.tol_classes import DefaultTolerances, Tol, TolWithMessage
from templates.groundstate import info_out as gs_info_out, eigval, evalcore, atoms, geometry

default = DefaultTolerances(integer=Tol(0),
                            float=Tol(1.e-8),
                            str=Tol(''),
                            )


hybrid_message = """Note, not all keys are present in each hybrid calculation.
We suggest the developer reviews which of these keys are required after generating a tolerance file: 
* 'Correlation type'. Only required for xctype="EXX"
* 'Exchange type'
* 'Mixing coefficient for exact exchange'
* 'Screening parameter (omega)'. Only currently required if Exchange-correlation type="408" (HSE)
"""

# Extend ground state tolerances
info_out = gs_info_out
info_out.update({'Correlation type': default.integer,
                 'Exchange type': default.str,
                 'Mixing coefficient for exact exchange': default.float,
                 'Screening parameter (omega)': default.float
                 }
                )

# Single dictionary for all tested hybrid XC outputs, for dumping to JSON
hybrid_tolerances = {'files_under_test': ['INFO.OUT', 'eigval.xml', 'evalcore.xml', 'atoms.xml', 'geometry.xml'],
                     'INFO.OUT': info_out,
                     'eigval.xml': eigval,
                     'evalcore.xml': evalcore,
                     'atoms.xml': atoms,
                     'geometry.xml': geometry
                     }
