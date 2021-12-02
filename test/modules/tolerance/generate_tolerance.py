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

from tol_classes import tol_file_name
from templates.groundstate import ground_state_tolerances
from templates.gw import gw_tolerances
from templates.hybrid import hybrid_tolerances, hybrid_message
from templates.bse import bse_tolerances
from templates.tddft import tddft_tolerances
from utils import ExcitingCalculation, get_calculation_types

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
    output_dict = {"files_under_test": tol_dict.pop("files_under_test")}

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

    """
    full_file_name = os.path.join(file_path, tol_file_name[calculation.name])

    # TODO(A/H/B/E) Issue 94. Fill in templates in subsequent merge requests
    if calculation == ExcitingCalculation.groundstate:
        write_tolerance_with_json(ground_state_tolerances, full_file_name)

    if calculation == ExcitingCalculation.gw:
        write_tolerance_with_json(gw_tolerances, full_file_name)

    if calculation == ExcitingCalculation.hybrid:
        print(hybrid_message)
        write_tolerance_with_json(hybrid_tolerances, full_file_name)

    if calculation == ExcitingCalculation.tddft:
        write_tolerance_with_json(tddft_tolerances, full_file_name)

    if calculation == ExcitingCalculation.bse:
        write_tolerance_with_json(bse_tolerances, full_file_name)

    if calculation == ExcitingCalculation.phonon:
        # TODO(Ignacio) Add phonon tolerances
        sys.exit('TODO(Ignacio). Add phonon tolerances')
        # write_tolerance_with_json(phonon_tolerances, full_file_name)

    # Properties
    if calculation in [ExcitingCalculation.band_structure, ExcitingCalculation.dos]:
        sys.exit('tolerance template needs to be defined:' + calculation.name)
        # write_tolerance_with_json(bandstructure_dos_tolerances, full_file_name)

    if calculation == ExcitingCalculation.plot:
        sys.exit('tolerance template needs to be defined:' + calculation.name)
        # write_tolerance_with_json(plotting_tolerances, full_file_name)

    if calculation == ExcitingCalculation.wannier:
        sys.exit('tolerance template needs to be defined:' + calculation.name)
        # write_tolerance_with_json(wannier_tolerances, full_file_name)

    if calculation == ExcitingCalculation.transport:
        sys.exit('tolerance template needs to be defined:' + calculation.name)
        # write_tolerance_with_json(transport_tolerances, full_file_name)

    if calculation == ExcitingCalculation.optical_properties:
        sys.exit('tolerance template needs to be defined:' + calculation.name)
        # write_tolerance_with_json(optical_properties_tolerances, full_file_name)

    if calculation == ExcitingCalculation.electric_properties:
        sys.exit('tolerance template needs to be defined:' + calculation.name)
        # write_tolerance_with_json(electric_field_properties_tolerances, full_file_name)

    if calculation == ExcitingCalculation.core_properties:
        sys.exit('tolerance template needs to be defined:' + calculation.name)
        # write_tolerance_with_json(core_properties_tolerances, full_file_name)

    if calculation == ExcitingCalculation.spin_properties:
        sys.exit('tolerance template needs to be defined:' + calculation.name)
        # write_tolerance_with_json(spin_properties_tolerances, full_file_name)


if __name__ == "__main__":
    inputs = parse_input_args()
    calculation_types = get_calculation_types(inputs['calculation_type'])

    for calculation in calculation_types:
        print('Generating ', calculation)
        generate_tolerance_file(calculation, inputs['destination'])
        print("\n\n")

