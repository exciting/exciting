import os
import sys
import tempfile
import textwrap
import traceback

from ase.tasks.task import Task
from ase.tasks.bulk import BulkTask
from ase.tasks.molecule import MoleculeTask
from ase.tasks.calcfactory import calcnames


usage = """\
Usage: ase [calculator] [task] [options] system(s)

%s
task:       'molecule', 'bulk' or the name of Python script that instantiates
            a Task object.  Default value is 'molecule'.
systems:    chemical formulas or filenames of files containing the atomic
            structure.

Try "ase molecule --help" or "ase bulk --help".
"""


def run(args=sys.argv[1:], calcname='emt', task=None):

    if isinstance(args, str):
        args = args.split(' ')

    argsoriginal = args[:]

    if len(args) > 0 and args[0] in calcnames:
        calcname = args.pop(0)

    if task is None:
        taskname = 'molecule'
        if (len(args) > 0 and
            (args[0] in ['molecule', 'bulk'] or args[0].endswith('.py'))):
            taskname = args.pop(0)

        if taskname.endswith('.py'):
            locals = {}
            execfile(taskname, locals, locals)
            tasks = [task for task in locals.values() if isinstance(task, Task)]
            assert len(tasks) == 1
            task = tasks[0]
        elif taskname == 'bulk':
            task = BulkTask()
        else:
            task = MoleculeTask()

    if len(args) == 0 and task.collection is None:
        sys.stderr.write(
            usage % textwrap.fill(', '.join(calcnames[:-1]) +
                                  ' or ' + calcnames[-1] +
                                  '.  Default value is emt.',
                                  initial_indent='calculator: ',
                                  subsequent_indent=' ' * 12))
        return
    
    task.set_calculator_factory(calcname)

    args = task.parse_args(args)

    if task.interactive_python_session:
        if '-i' in argsoriginal:
            argsoriginal.remove('-i')
        if '--interactive-python-session' in argsoriginal:
            argsoriginal.remove('--interactive-python-session')
        file = tempfile.NamedTemporaryFile()
        file.write('import os\n')
        file.write('if "PYTHONSTARTUP" in os.environ:\n')
        file.write('    execfile(os.environ["PYTHONSTARTUP"])\n')
        file.write('from ase.tasks.main import run\n')
        file.write('atoms, task = run(%r, %r)\n' % (argsoriginal, calcname))
        file.flush()
        os.system('python -i %s' % file.name)
        return

    atoms = task.run(args)

    return atoms, task


def main():
    try:
        run()
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception:
        traceback.print_exc()
        sys.stderr.write("""
An exception occurred!  Please report the issue to
ase-developer@listserv.fysik.dtu.dk - thanks!  Please also report this
if it was a user error, so that a better error message can be provided
next time.""")
        raise
