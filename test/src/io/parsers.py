"""
Wrapper module to expose exciting parsers

For other codes, please replace with your own parser module/s.
"""
import os
import warnings
from typing import Union
import subprocess

from excitingtools.parser_utils.grep_parser import grep

from ..runner.profile import compiler_version_identifier_map, Compiler
from ..tester.failure import Failure, Failure_code

try:
    __import__('excitingtools')
except ModuleNotFoundError:
    # Give instructions to user rather than automatically installing:
    message = """excitingtools must be installed for the test suite to run.
To install excitingtools, from the test directory, type:

  python3 -m pip install ../tools/exciting_tools

excitingtools can be uninstalled by typing:

  pip uninstall excitingtools  
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
    # TODO(Alex) Move `exciting_test_dir` to exciting_settings
    # Move this function to excitingtools
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
    Get the compiler command used to build exciting from the make.inc file 
    and determines the version via cmd line.

    Uses relative directories => Assumes test suite always ran from
    the test directory

    :return compiler: enum compiler type
    """
    exciting_root = get_exciting_root()
    make_inc = os.path.join(exciting_root, 'build/make.inc')
    result = grep('F90', make_inc, options = {'w': ' '}).splitlines()

    for line in result:
        file_comment = line[0] == '#'
        f77_line = 'F77' in line 

        if not file_comment and not f77_line:
            compiler_version_cmd_str = line.split()[2] + " --version"
            compiler_version_output = subprocess.check_output(compiler_version_cmd_str, shell=True).decode()

            for compiler_id in compiler_version_identifier_map.keys():
                if compiler_id in compiler_version_output:
                    return compiler_version_identifier_map[compiler_id]

    return None