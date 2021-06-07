"""
 Environment and directory variables
"""

from collections import namedtuple

_binary_names    = ['exciting_serial', 'exciting_smp', 'exciting_purempi', 'exciting_mpismp']
_exe_ref         = 'exciting_smp'
_default_np      = {'exciting_serial':1, 'exciting_smp':1, 'exciting_purempi':2, 'exciting_mpismp':2}
_default_threads = {'exciting_serial':1, 'exciting_smp':2, 'exciting_purempi':1, 'exciting_mpismp':2}
_action_choices  = ['run', 'ref', 'clean']

Defaults = namedtuple('Defaults', ['max_time',       # Time after which a test is killed (in seconds)
                                   'test_farm',      # Test farm directory
                                   'species',        # Species file directory
                                   'input_file',     # Input file for exciting
                                   'main_output',    # Main output from exciting
                                   'run_dir',        # Run directory for tests
                                   'ref_dir',        # Directory for test reference data
                                   'exe_dir',        # Location of exciting executable
                                   'init_default',   # Template for init xml
                                   'ignored_output', # Output files to not reference
                                   'not_clean',      # Files protected from make clean
                                   'binary_names',   # List of exciting executables
                                   'exe_ref',        # Executable for running the reference calculations
                                   'default_np',     # Dict of default MPI processes per executable
                                   'default_threads',# Dict of default threads per executable 
                                   'action_choices'  # Action choices for test script 
                                   ])

# Define an immutable instance of the default settings
# Access like settings.max_time
settings = Defaults(max_time      = 1800,
                    test_farm     = 'test_farm',  
                    species       = '../species',
                    input_file    = 'input.xml',
                    main_output   = 'INFO.OUT',    
                    run_dir       = 'run',       
                    ref_dir       = 'ref',        
                    exe_dir       = '../../../../bin/',        
                    init_default  = 'xml/init_templates/init_default.xml' ,   
                    ignored_output = ['STATE.OUT', 'OCC', 'EVEC', 'EVALSV', 'EVALFV', 'APWCMT', 'SYM', 
                                      'PMAT', 'FERMISURF', 'RMSDVEFF', 'LOCMT', 'EXCLI', 'SCCLI'], 
                    not_clean     = ['input.xml'],
                    binary_names  = _binary_names,
                    exe_ref       = _exe_ref,
                    default_np    = _default_np,
                    default_threads = _default_threads,
                    action_choices  = _action_choices
                    )
