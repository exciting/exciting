"""
 Environment and directory variables
"""

from collections import namedtuple

_binary_names    = ['excitingser', 'excitingmpi', 'excitingmpismp']
_default_np      = {'excitingser':1, 'excitingmpi':2, 'excitingmpismp':2}
_default_threads = {'excitingser':1, 'excitingmpi':1, 'excitingmpismp':2}
_action_choices  = ['run', 'ref', 'clean', 'report']

Defaults = namedtuple('Defaults', ['max_time',       # Time after which a test is killed (in seconds)
                                   'test_farm',      # Test farm directory
                                   'species',        # Species file directory
                                   'main_output',    # Main output from exciting
                                   'run_dir',        # Run directory for tests
                                   'ref_dir',        # Directory for test reference data
                                   'exe_dir',        # Location of exciting executable
                                   'init_default',   # Template for init xml
                                   'ignored_output', # Output files to not reference
                                   'not_clean',      # Files protected from make clean
                                   'binary_names',   # List of exciting executables
                                   'default_np',     # Dict of default MPI processes per executable
                                   'default_threads',# Dict of default threads per executable 
                                   'action_choices'  # Action choices for test script 
                                   ])

# Define an immutable instance of the default settings
# Access like settings.max_time
settings = Defaults(max_time      = 1800,
                    test_farm     = 'test_farm',  
                    species       = '../species',        
                    main_output   = 'INFO.OUT',    
                    run_dir       = 'run',       
                    ref_dir       = 'ref',        
                    exe_dir       = '../../../../bin/',        
                    init_default  = 'xml/init_templates/init_default.xml' ,   
                    ignored_output = ['input.xml', 'STATE.OUT', 'OCC', 'EVEC', 'EVALSV', 'EVALFV'], 
                    not_clean     = ['input.xml'],
                    binary_names  = _binary_names,
                    default_np    = _default_np,
                    default_threads = _default_threads,
                    action_choices  = _action_choices
                    )
