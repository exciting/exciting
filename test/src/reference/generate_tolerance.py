"""
Module to generate tolerance files for an application test.

A complete set of reference data that's test_suite-independent
requires all tolerances to be explicitly defined independently.
Tolerances are written in JSON as it's the simplest structured
file format to work with.

To generate default tolerances files:

```
python3 generate_tolerance.py -for groundstate gw -write-to path/2/test_farm/ref
```

for a specific file. To generate ALL tolerance files in the directory
of execution, simply omit the options:

```
python3 tolerance.py
```

"""
import argparse
import json
import copy
import os
import sys

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(os.path.dirname(parentdir))


from src.tolerance.tol_classes import tol_file_name
from src.tolerance.templates.groundstate import ground_state_tolerances
from src.tolerance.templates.gw import gw_tolerances
from src.tolerance.templates.hybrid import hybrid_tolerances, hybrid_message
from src.tolerance.templates.bse import bse_tolerances
from src.tolerance.templates.tddft import tddft_tolerances
from src.tolerance.templates.rt_tddft import rt_tddft_tolerances
from src.tolerance.templates.properties import optical_properties_tolerances, core_properties_tolerances, \
    electric_field_properties_tolerances, spin_properties_tolerances
from src.tolerance.templates.transport import transport_tolerances
from src.tolerance.templates.plotting import plotting_tolerances
from src.tolerance.templates.bandstructure_dos import bandstructure_tolerances, dos_tolerances
from src.tolerance.templates.wannier import wannier_tolerances
from src.runner.profile import ExcitingCalculation, get_calculation_types


def parse_input_args() -> dict:
    """
    Parse command line inputs

    :return dict: Dictionary of parsed inputs
    """
    parser = argparse.ArgumentParser(description="Module to generate tolerance files for an application test.")
    all_calculation_types = list(tol_file_name.keys())

    parser.add_argument('-for',
                        help='Calculation type to generate tolerance file for',
                        dest='calculation_type',
                        type=str,
                        nargs='+',
                        default=all_calculation_types
                        )

    parser.add_argument('-write-to',
                        help='Directory to write tolerance file to',
                        dest='destination',
                        default='./',
                        type=str
                        )

    return parser.parse_args().__dict__


def serialise_tolerance_values(tol_dict: dict) -> dict:
    """
    JSON can't serialise objects, therefore iterate through the
    tolerance dictionary and convert value objects to dictionaries.

    NOTE: This routine expects dictionary values of type Tol(value, Unit)
    where Tol is a class and Unit is an enum class.

    Structure of tol_dict should be guaranteed:

    {'key': default.value} which is equivalent to {'key': Tol(value, Unit)}

    and return:

    {'key': {'tol': numerical_tolerance, 'unit': 'unit_string'}
    """
    output_dict = {}
    for file_name, file_tols in tol_dict.items():
        serialised_tol_dict = {}
        for key, value in file_tols.items():
            serialised_tol_dict[key] = value.to_dict()

        output_dict[file_name] = copy.deepcopy(serialised_tol_dict)
        serialised_tol_dict.clear()

    return output_dict


def write_tolerance_with_json(tolerances: dict, file_name: str):
    """
    Wrapper for specifically writing tolerance files with JSON.

    :param dict, tolerances: Tolerances for writing to file
    :param str, file_name: File
    """
    with open(file_name, "w") as fid:
        json.dump(serialise_tolerance_values(tolerances), fid, indent=2)


def generate_tolerance_file(calculation: ExcitingCalculation, file_path: str):
    """
    Generate default tolerance files for exciting outputs.

    Write a tolerance file associated with a calculation of type
    'calculation' to file_path/tol_file_name[calculation]

    For example for ground_state, it writes:
      test_farm/groundstate/test/ref/tolerance_ground_state.json

    :param ExcitingCalculation calculation: exciting calculation type.
    :param str file_path: Directory to write tolerances to.
    """
    tolerances = {}

    full_file_name = os.path.join(file_path, tol_file_name[calculation.name])

    if calculation == ExcitingCalculation.groundstate:
        tolerances = ground_state_tolerances

    if calculation == ExcitingCalculation.gw:
        tolerances = gw_tolerances

    if calculation == ExcitingCalculation.hybrid:
        print(hybrid_message)
        tolerances = hybrid_tolerances

    if calculation == ExcitingCalculation.tddft:
        tolerances = tddft_tolerances

    if calculation == ExcitingCalculation.rt_tddft:
        tolerances = rt_tddft_tolerances

    if calculation == ExcitingCalculation.bse:
        tolerances = bse_tolerances

    # Properties
    if calculation == ExcitingCalculation.band_structure:
        tolerances = bandstructure_tolerances

    if calculation == ExcitingCalculation.dos:
        tolerances = dos_tolerances

    if calculation == ExcitingCalculation.plot:
        tolerances = plotting_tolerances

    if calculation == ExcitingCalculation.wannier:
        tolerances = wannier_tolerances

    if calculation == ExcitingCalculation.transport:
        tolerances = transport_tolerances

    if calculation == ExcitingCalculation.optical_properties:
        tolerances = optical_properties_tolerances

    if calculation == ExcitingCalculation.electric_properties:
        tolerances = electric_field_properties_tolerances

    if calculation == ExcitingCalculation.core_properties:
        tolerances = core_properties_tolerances

    if calculation == ExcitingCalculation.spin_properties:
        tolerances = spin_properties_tolerances

    write_tolerance_with_json(tolerances, full_file_name)


if __name__ == "__main__":
    inputs = parse_input_args()
    calculation_types = get_calculation_types(inputs['calculation_type'])

    for calculation in calculation_types:
        print('Generating ', calculation)
        generate_tolerance_file(calculation, inputs['destination'])
        print("\n\n")
