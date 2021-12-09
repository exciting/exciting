"""
Environment and directory variables
"""

from collections import namedtuple
import os

from ..runner.profile import Build_type, build_type_enum_to_str

Defaults = namedtuple('Defaults', ['max_time',        # Time after which a test is killed (in seconds)
                                   'test_farm',       # Test farm directory
                                   'input_file',      # Input file for exciting
                                   'main_output',     # Main output from exciting
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
                                   'exe_ref'          # Executable for running the reference calculations
                                   ])


_binary_names    = [binary for binary in build_type_enum_to_str.values()]
_exe_ref         = Build_type.puresmp
_default_np      = {Build_type.serial: 1, Build_type.puresmp: 1, Build_type.purempi: 2, Build_type.mpiandsmp: 2}
_default_threads = {Build_type.serial: 1, Build_type.puresmp: 2, Build_type.purempi: 1, Build_type.mpiandsmp: 2}


# Define an immutable instance of the default settings
# Access like settings.max_time
settings = Defaults(max_time       = 1800,
                    test_farm      = 'test_farm',  
                    input_file     = 'input.xml',
                    main_output    = 'INFO.OUT',    
                    run_dir        = 'run',       
                    ref_dir        = 'ref',        
                    exe_dir        =  os.path.normpath(os.path.join(os.getcwd(), '../bin')),
                    ignored_output = ['STATE.OUT', 'OCC', 'EVEC', 'EVALSV', 'EVALFV', 'APWCMT', 'SYM',
                                      'PMAT', 'FERMISURF', 'RMSDVEFF', 'LOCMT', 'EXCLI', 'SCCLI'], 
                    binary_names    = _binary_names,
                    binary_mpismp   = Build_type.mpiandsmp,
                    binary_purempi  = Build_type.purempi,
                    binary_smp      = Build_type.puresmp,
                    binary_serial   = Build_type.serial,
                    default_np      = _default_np,
                    default_threads = _default_threads,
                    exe_ref         = _exe_ref
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
