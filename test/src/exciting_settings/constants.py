"""Environment and directory variables
"""
from collections import namedtuple
import os
from abc import ABC, abstractmethod

from ..runner.profile import BuildType, build_type_enum_to_str

Defaults = namedtuple('Defaults', ['max_time',        # Time after which a test is killed (in seconds)
                                   'test_farm',       # Test farm directory
                                   'input_file',      # Input file for exciting
                                   'run_dir',         # Run directory for tests
                                   'ref_dir',         # Directory for test reference data
                                   'exe_dir',         # Location of exciting executable
                                   'ignored_output',  # Output files to not reference
                                   'binary_names',    # List of exciting executables
                                   'binary_mpismp',   # Exciting executable with smp and mpi parallelisation
                                   'binary_purempi',  # Exciting executable with mpi parallelisation
                                   'binary_smp',      # Exciting executable with smp parallelisation
                                   'binary_serial',   # Serial exciting executable
                                   'default_np',      # Dict of default MPI processes per executable
                                   'default_threads', # Dict of default threads per executable
                                   'exe_ref',         # Executable for running the reference calculations
                                   'SAFE_MODE',       # Run the code with extra checks/turn warnings to errors
                                   'tests_file',      # Yaml configuration file, defining all tests
                                   'defaults_file'    # Yaml file with test configuration defaults
                                   ])


_binary_names    = [binary for binary in build_type_enum_to_str.values()]
_default_np      = {BuildType.serial: 1, BuildType.puresmp: 1, BuildType.purempi: 2, BuildType.mpiandsmp: 2}
_default_threads = {BuildType.serial: 1, BuildType.puresmp: 2, BuildType.purempi: 1, BuildType.mpiandsmp: 2}


# Define an immutable instance of the default settings
# Access like settings.max_time
settings = Defaults(max_time       = 1800,
                    test_farm      = 'test_farm',
                    input_file     = 'input.xml',
                    run_dir        = 'run',
                    ref_dir        = 'ref',
                    exe_dir        =  os.path.normpath(os.path.join(os.getcwd(), '../bin')),
                    ignored_output = ['STATE.OUT', 'OCC', 'EVEC', 'EVALSV', 'EVALFV', 'APWCMT', 'SYM',
                                      'PMAT', 'FERMISURF', 'RMSDVEFF', 'LOCMT', 'EXCLI', 'SCCLI'],
                    binary_names    = _binary_names,
                    binary_mpismp   = BuildType.mpiandsmp,
                    binary_purempi  = BuildType.purempi,
                    binary_smp      = BuildType.puresmp,
                    binary_serial   = BuildType.serial,
                    default_np      = _default_np,
                    default_threads = _default_threads,
                    exe_ref         = BuildType.serial,
                    SAFE_MODE       = bool(os.getenv('SAFE_MODE')) if os.getenv('SAFE_MODE') is not None else False,
                    tests_file      = 'tests_config.yml',
                    defaults_file   = 'defaults_config.yml'
                    )

species_files = ['Ni.xml', 'La.xml', 'K.xml', 'Xe.xml', 'Ag.xml', 'Bk.xml', 'Co.xml', 'Md.xml', 'Lu.xml', 'Ar.xml',
                 'Bi.xml', 'Cm.xml', 'H.xml', 'Yb.xml', 'Zn.xml', 'Te.xml', 'I.xml', 'Cl.xml', 'As.xml', 'Mg.xml',
                 'No.xml', 'Ta.xml', 'N.xml', 'Ac.xml', 'Y.xml', 'At.xml', 'Tb.xml', 'Tc.xml', 'Au.xml', 'O.xml',
                 'Lr.xml', 'In.xml', 'Ge.xml', 'Re.xml', 'Pm.xml', 'Gd.xml', 'Kr.xml', 'Po.xml', 'Sc.xml', 'Rf.xml',
                 'Sb.xml', 'Rb.xml', 'Ru.xml', 'Dy.xml', 'Ho.xml', 'Ra.xml', 'Se.xml', 'Sr.xml', 'Fr.xml', 'Ga.xml',
                 'Fe.xml', 'Es.xml', 'Si.xml', 'Pr.xml', 'Pd.xml', 'Er.xml', 'Rn.xml', 'Ir.xml', 'He.xml', 'Eu.xml',
                 'Pt.xml', 'Pu.xml', 'Sn.xml', 'Pb.xml', 'Hf.xml', 'Fm.xml', 'Rh.xml', 'Sm.xml', 'Pa.xml', 'Hg.xml',
                 'Os.xml', 'B.xml', 'U.xml', 'Zr.xml', 'Cf.xml', 'C.xml', 'Na.xml', 'Li.xml', 'Mo.xml', 'Cs.xml',
                 'Al.xml', 'V.xml', 'Cd.xml', 'Tm.xml', 'Tl.xml', 'Ba.xml', 'Ce.xml', 'W.xml', 'Am.xml', 'Cr.xml',
                 'Nb.xml', 'Mn.xml', 'S.xml', 'Ca.xml', 'Be.xml', 'Br.xml', 'Th.xml', 'Ti.xml', 'Np.xml', 'Ne.xml',
                 'P.xml', 'Cu.xml', 'F.xml', 'Nd.xml']

# Keys to be removed from test and reference dictionaries (i.e. that we do not want to test)
keys_to_remove = {'INFO.OUT': [['scl', 'Wall time (seconds)']]}


def main_output(method: str) -> str:
    """Given a method, return the main output file.

    :param method: Method
    :return Output file name
    """
    if method.lower() == 'gw':
        return 'GW_INFO.OUT'
    else:
        return "INFO.OUT"


class RunProperties(ABC):
    """Base class for properties used by a run."""
    executable_cmd: str
    max_time: int
    ref_dir: str
    run_dir: str

    @abstractmethod
    def calculation_completed(self, terminated_cleanly: bool) -> bool:
        """ Has a calculation successfully completed.

        :param terminated_cleanly: Subprocess returned cleanly.
        """
        ...


class ExcitingRunProperties(RunProperties):
    """Container grouping properties associated with a program execution

    Note, this could be further specified for GS and excited state calculations,
    whether the measure of whether a calculation successfully completed will
    differ.
    """
    def __init__(self, executable_cmd: str, max_time: int, ref_dir, run_dir, outputs=None, main_output=None):
        self.executable_cmd = executable_cmd
        self.max_time = max_time
        self.ref_dir = ref_dir
        self.run_dir = run_dir
        self.outputs = outputs
        self.main_output = main_output

    def calculation_completed(self, terminated_cleanly: bool) -> bool:
        """Has an exciting calculation successfully completed.

        Defined as the calculation cleanly exiting, and the
        main output file being written. This should be sufficient for GS.
        Should not need to check INFO.OUT or warnings file to confirm max SCF reached:
        - This is responsibility of the caller, and will be caught be the reference comparison

        TODO(Bene). Issue 133. Extend ExcitingRunProperties for Excited States
        To confirm an XS calculation finished successfully:
          * One needs to consider what to check in this case.
          * May need a new class or method

        :param terminated_cleanly: Subprocess returned cleanly.
        :return run_success: Whether the run succeeded.
        """
        main_output = os.path.join(self.run_dir, self.main_output)
        run_success = terminated_cleanly and os.path.isfile(main_output)
        return run_success
