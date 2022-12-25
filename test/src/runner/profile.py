"""
Build profile and exciting methods
"""
import enum
import re
from typing import List, Union
import sys


class Compiler(enum.Enum):
    """
    Only include compilers that exciting supports
    """
    intel = enum.auto()
    gcc = enum.auto()
    all = enum.auto()


compiler_enum_map = {'ifort': Compiler.intel, 'gfortran': Compiler.gcc}


class BuildType(enum.Enum):
    """
    """
    debug_serial = enum.auto()
    debug_mpiandsmp = enum.auto()
    serial = enum.auto()
    mpiandsmp = enum.auto()
    smp = enum.auto()
    purempi = enum.auto()
    all = enum.auto()


build_type_str_to_enum = {'exciting_serial': BuildType.serial,
                          'exciting_mpismp': BuildType.mpiandsmp,
                          'exciting_smp': BuildType.smp,
                          'exciting_purempi': BuildType.purempi
                          }

build_type_enum_to_str = {value: key for key, value in build_type_str_to_enum.items()}


class CompilerBuild:
    def __init__(self, compiler: Union[Compiler, str], build: Union[BuildType, str]):
        self.compiler = Compiler[compiler] if isinstance(compiler, str) else compiler
        self.build = BuildType[build] if isinstance(build, str) else build


# TODO(Alex) Move this to exciting_settings
class ExcitingCalculation(enum.Enum):
    """
    Broad categories of calculations performed by exciting

    To generalise the test framework, this should inherit from an abstract class
    """
    groundstate = enum.auto()
    gw = enum.auto()
    tddft = enum.auto()
    rt_tddft = enum.auto()
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
