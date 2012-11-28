'Some ag tests.'

import os
import sys

# Make sure ag can run in terminal mode without $DISPLAY and gtk:
sys.argv = ['ag', '--terminal']
display = os.environ.pop('DISPLAY', None)
error = False
try:
    from ase.gui.ag import main
    main()
    assert 'gtk' not in sys.modules
except:
    error = True

if display is not None:
    os.environ['DISPLAY'] = display

if error:
    raise
