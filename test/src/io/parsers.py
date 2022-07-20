"""
Wrapper module to expose exciting parsers

For other codes, please replace with your own parser module/s.
"""
import sys
import os
import warnings
from typing import Union

from excitingtools.parser_utils.grep_parser import grep

from ..runner.profile import compiler_enum_map, Compiler
from ..tester.failure import Failure, Failure_code


def install_excitingtools():
    """
    Install excitingtools to provide the exciting parsers
    """
    import subprocess

    if 'excitingtools' in sys.modules:
        return
    else:
        print("Running pip install for exciting_tools. "
              "See <EXCITINGROOT>/test/tools/parsers.py for more info")
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-e", "../tools/exciting_tools/"])
        # Restart script such that package is found following installation
        os.execv(sys.executable, [sys.executable] + sys.argv)

try:
    __import__('excitingtools')
except ModuleNotFoundError:
    # Give instructions to user rather than automatically installing:
    message = """excitingtools must be installed for the test suite to run.
To install excitingtools, from the test directory, type:

  pip3 install -e ../tools/exciting_tools

excitingtools can be uninstalled by typing:

  pip3 uninstall excitingtools  
    """
    warnings.warn(message)
finally:
    from excitingtools.parser_utils.erroneous_file_error import ErroneousFileError
    from excitingtools.exciting_dict_parsers.parser_factory import parser_chooser


def read_output_file(file_name: str) -> Union[dict, Failure]:
    """
    Read an exciting output file

    :param str file_name: File name prepended by directory path
    :return Union[dict, Failure] ref_data: Reference data or Failure object
    """
    # Assume file_name of the form test_farm/method/test_case/ref/file.out
    sub_dir = file_name.split('/')[-2]
    if not (sub_dir in ['ref', 'run']):
        raise ValueError(f"Subdirectory in which output file exists is not 'ref' or 'run' \a"
                         "which is expected by the test suite: {file_name}")

    try:
        data = parser_chooser(file_name)
        return data
    except OSError:
        failure_code = {'ref': Failure_code.REFERENCE, 'run': Failure_code.RUN}
        return Failure(test_name=file_name, failure_code=failure_code[sub_dir])
    except ErroneousFileError:
        return Failure(test_name=file_name, failure_code=Failure_code.ERRORFILE)


def get_exciting_root() -> Union[str, None]:
    """
    Get the absolute path to exciting's root directory.

    This function only works if executed from with the app-test directory,
    such that the directory structure is:
    some/path/exciting_root/test/...

    :return str root: exciting's root directory.
    """
    exciting_test_dir = 'test'
    sub_dirs = os.getcwd().split('/')

    root = '/'
    for sub_dir in sub_dirs:
        if sub_dir == exciting_test_dir:
            return root
        root = os.path.join(root, sub_dir)

    return None


def get_compiler_type() -> Union[Compiler, None]:
    """
    Get the compiler used to build exciting from the make.inc file

    Uses relative directories => Assumes test suite always ran from
    the test directory

    :return compiler: enum compiler type
    """
    exciting_root = get_exciting_root()
    make_inc = os.path.join(exciting_root, 'build/make.inc')
    result = grep('F90', make_inc).splitlines()

    for line in result:
        file_comment = line[0] == '#'

        if not file_comment:
            for compiler in compiler_enum_map.keys():
                if compiler in line.lower():
                    return compiler_enum_map[compiler]

    return None
