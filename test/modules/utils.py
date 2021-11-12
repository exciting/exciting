import enum
import subprocess
import os
import re
from typing import List
import sys


class ExcitingCalculation(enum.Enum):
    """
    Broad categories of calculations performed by exciting
    """
    groundstate = enum.auto()
    gw = enum.auto()
    tddft = enum.auto()
    bse = enum.auto()
    hybrid = enum.auto()
    phonon = enum.auto()
    band_structure = enum.auto()
    dos = enum.auto()
    plot = enum.auto()
    wannier = enum.auto()
    transport = enum.auto()
    optical_properties = enum.auto()
    electric_properties = enum.auto()
    core_properties = enum.auto()
    spin_properties = enum.auto()


def get_calculation_types(input_calcs: List[str]) -> List[ExcitingCalculation]:
    """
    Given a list of strings, determine which calculation types are
    a match.

    :param List[str] input_calcs: Input strings for calculation names
    :return  List[Calculation] List of Calculation enums
    """
    all_calculations_str = "\n".join(calc for calc in ExcitingCalculation._member_names_)

    matched_calculations = []
    for calc in input_calcs:
        matched_calculations += re.findall("^.*" + calc + ".*$", all_calculations_str, re.MULTILINE)

    if len(matched_calculations) != len(input_calcs):
        unmatched_calcs = set(input_calcs) - set(matched_calculations)
        print("Some calculation inputs did not match any valid method choices: ", unmatched_calcs)
        print("Here is a complete list of valid method strings (substring matches are also valid):")
        all_calculations_pretty_str = "\n".join(" * " + calc for calc in ExcitingCalculation._member_names_)
        print(all_calculations_pretty_str)
        sys.exit()

    names_to_enums = {calc.name: calc for calc in ExcitingCalculation}
    return [names_to_enums[name] for name in matched_calculations]


class Compiler(enum.Enum):
    """
    Only include compilers that exciting supports 
    """
    intel = enum.auto()
    gcc = enum.auto()
    all = enum.auto()


compiler_enum_map = {'ifort': Compiler.intel, 'gfortran': Compiler.gcc}


class Build_type(enum.Enum):
    """
    """
    debug_serial = enum.auto()
    debug_mpiandsmp = enum.auto()
    serial = enum.auto()
    mpiandsmp = enum.auto()
    puresmp = enum.auto()
    purempi = enum.auto()
    all = enum.auto()


build_type_str_to_enum = {'exciting_serial': Build_type.serial,
                          'exciting_mpismp': Build_type.mpiandsmp,
                          'exciting_smp': Build_type.puresmp,
                          'exciting_purempi': Build_type.purempi
                          }

build_type_enum_to_str = {value: key for key, value in build_type_str_to_enum.items()}


class CompilerBuild:
    def __init__(self, compiler: Compiler, build: Build_type):
        self.compiler = compiler
        self.build = build


def grep(string: str, fname: str) -> str:
    """
    Wrapper for grep 

    :param str string: string to grep for in file fname
    :param str fname:  file name to search
    :return str output: Grep'ed string
    """

    opts = ''
    grep_str = "grep " + opts + " '" + string + "' " + fname

    try:
        output = subprocess.check_output(grep_str, shell=True).decode("utf-8")
    except subprocess.CalledProcessError as grepexc:
        print("subprocess error:", grepexc.returncode, "grep found:", grepexc.output)
        output = grepexc.output

    return output


def get_exciting_root() -> str:
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


def get_compiler_type():
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
