"""
Module containing default tolerances for TDDFT output

Each entry has a unique key and value = {'tolerance':value, 'unit':Unit.unit}
where Unit is an enum class.  If a quantity does not have an associated unit,
unit = null.

Keys were generated by using the tolerance.py script. See that for command line
details.

Tolerances in default_tolerances cover files (** is one of 11, 22, 33):
  * EPSILON_NAR_FXCMB1_OC??_QMT001.OUT
  * EPSILON_NAR_NLF_FXCMB1_OC??_QMT001.OUT
  * LOSS_NAR_FXCMB1_OC??_QMT001.OUT
  * LOSS_NAR_NLF_FXCMB1_OC??_QMT001.OUT

TODO(Alex/Bene/Keith) Issue #110. Increase files under test in TDDFT
Files that require test coverage:
  * SIGMA_NAR_??.OUT
  * SUMRULES_NAR_FXCMB1_??.OUT
  * DIELTENS0_??.OUT
  * SCREEN_Q??.OUT
  * FXC_BSE_HEAD_NAR_QMT001.OUT
  * GQPOINTS_SCR_Q00001.OUT
  * EIGVAL_SCR.OUT
"""

from excitingtools.units import Unit

from tol_classes import DefaultTolerances, Tol


default = DefaultTolerances(integer=Tol(0),
                            float=Tol(1.e-8),
                            str=Tol(''),
                            energy=Tol(1.e-8, Unit.ev),
                            frequency=Tol(1.e-8, Unit.ev),
                            oscillator_strength=Tol(1.e-8, Unit.null),
                            loss_function=Tol(1.e-8, Unit.inv_ev), 
                            structure_factor=Tol(1.e-1, Unit.null)
                            )

tddft_epsilon_tols = {
    'frequency': default.frequency,
    'real_oscillator_strength': default.oscillator_strength,
    'imag_oscillator_strength': default.oscillator_strength,
    'real_oscillator_strength_kkt': default.oscillator_strength
}

tddft_loss_tols = {
    'frequency': default.frequency,
    'loss_function': default.loss_function,
    'structure_factor': default.structure_factor
}

tddft_tolerances = {'files_under_test': ['EPSILON_??.OUT',
                                         'LOSS_??.OUT'],
                    'EPSILON_??.OUT': tddft_epsilon_tols,
                    'LOSS_??.OUT': tddft_loss_tols}
