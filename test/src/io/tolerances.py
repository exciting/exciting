"""
Parse and preprocess JSON tolerance files
"""
from typing import List, Tuple
import re
import os
import json
from copy import deepcopy

from ..utilities.wildcard_processor import wildcard_processor, wildcard


def update_wildcard_files_under_test(files_under_test: List[str], file_names: List[str]) -> List[str]:
    """
    Create copies of 'files_under_test' entries for every file name with a common name prefix.

    TODO(Alex) Ultimately we can remove this and have the config file define files_under_test.

    For files for which there can be multiple extensions (FILE_0.OUT, FILE_1.OUT, etc), it is expected that
    a single tolerance will be defined with the key 'FILE_<WILDCARD>' in the tolerances dictionary and therefore
    only a single entry 'FILE_<WILDCARD>' in the 'files_under_test' list

    This routine will copy this element for every reference file that matches the prefix.

    For example:
      On input, files_under_test = ['FILE_<WILDCARD>', 'OTHER_XY.OUT']
      Reference files = ['FILE_ABC_0.OUT', 'FILE_DEF_1.OUT', 'OTHER_XY.OUT']
      On output, files_under_test_updated = ['FILE_ABC_0.OUT', 'FILE_DEF_1.OUT', 'OTHER_XY.OUT']

    :param List[str] files_under_test: File names under regression testing, which may or may not contain wildcards in their
     names.
    :param List[str] file_names: List of file names containing (but not exclusively containing) all reference output
    files to test i.e. ['INFO.OUT', 'GW_INFO.OUT', 'random_file.out'].
    :return List[str] files_under_test_updated: File names under regression testing, with wildcards replaced according
    to reference files listed in file_names.
    """
    assert files_under_test, "No 'files_under_test' are specified"
    files_under_test_updated = []
    for file_under_test in files_under_test:
        regex = re.compile(wildcard_processor(file_under_test))

        for file_name in file_names:
            if regex.match(file_name):
                files_under_test_updated.append(file_name)

    return files_under_test_updated


def update_wildcard_tolerances(tolerances: dict, files_under_test: List[str]) -> dict:
    """
    Create copies of tolerance entries for every file with a common name prefix.

    For files for which there can be multiple extensions (FILE_0.OUT, FILE_1.OUT, etc), it is expected that
    a single tolerance will be defined with the key 'FILE_<WILDCARD>' in the tolerances dictionary.

    This routine will copy this key:value for every reference file that matches the prefix.

    For example:
      On input, tolerances = {'FILE_<WILDCARD>': tol_dict}
      Reference file_names = [FILE_0.OUT, FILE_1.OUT]
      On output, tolerances = {'FILE_0.OUT': tol_dict, 'FILE_1.OUT': tol_dict,}

    If a key is already present, the initial tolerances will be retained.

    :param dict tolerances: Tolerances for all files under regression testing.
    :param List[str] files_under_test: List of all file names under testing (should not contain wildcards)
    :return dict tolerances_updated: Tolerances with for all files under regression testing.
    """
    for file in files_under_test:
        if any([True for w in wildcard if w in file]):
            raise ValueError('Wild card string present in file name: ' + file)

    tolerances_updated = {}
    unmatched_file_names = []

    for file_name, file_tolerances in tolerances.items():
        regex = re.compile(wildcard_processor(file_name))
        matched = False

        for reference_file in files_under_test:
            if regex.match(reference_file):
                tolerances_updated[reference_file] = file_tolerances
                matched = True

        if not matched:
            unmatched_file_names.append(file_name)

    if unmatched_file_names:
        raise KeyError(f"One or more keys in the tolerances JSON does not correspond to specified files_under_test "
                       f"or cannot be matched with the wildcard symbols defined in wildcards.py: "
                       f"{unmatched_file_names}")

    return tolerances_updated


def list_tolerance_files_in_directory(directory: str) -> List[str]:
    """
    List the tolerances from any files in 'directory' of the form `*tol*.json`

    :param str directory: Directory containing files
    :return List[str] tolerance_files: List of tolerance files
    """
    files_in_directory = next(os.walk(directory))[2]
    r = re.compile("tolerance*.*json")
    tolerance_files = list(filter(r.match, files_in_directory))
    return tolerance_files


def load_tolerances(directory: str, tolerance_files: List[str]) -> dict:
    """
    Load JSON tolerance files

    :param str directory: Directory containing `*tol*.json` file
    :param List[str] tolerance_files:
    :return dict tolerances: Dictionary of tolerances
    """
    if not tolerance_files:
        raise FileNotFoundError('No tolerance files listed in `tolerance_files` argument')

    tolerances = {}
    for file in tolerance_files:
        with open(os.path.join(directory, file)) as fid:
            tolerances.update(json.load(fid))

    return tolerances


def strip_tolerance_units(json_tolerance: dict) -> dict:
    """
    Strip the units from the json_tolerances:

    json_tolerance[file_name][key] =  {'tol': 1e-08, 'unit': 'Bohr', 'comment': 'text'}
    to
    just_tolerance[file_name][key] =  1e-08
    """

    just_tolerances = {}
    for file_name, tolerances in json_tolerance.items():
        tmp = {}
        for key, entry in tolerances.items():
            tmp[key] = entry['tol']
        just_tolerances[file_name] = deepcopy(tmp)
        del tmp
    return just_tolerances


def get_json_tolerances(full_ref_dir: str, files_under_test: List[str]) -> Tuple[dict, List[str]]:
    """
    Load and preprocess JSON tolerance files.

    TODO(A/B/H) Issue 100. We currently discard the tolerance units.
       Update ErrorFinder class to use tol AND units in data comparison,
       Could achieve this by passing the comparison function as an argument
       See function documentation for a description of the issue

    :param str full_ref_dir: Full path to reference directory.
    :param List[str] files_under_test: Files under test.
    :return Tuple[dict, List[str]]: tolerances_without_units, reference_outputs:
       Tolerances dictionary of the form {'method': tols}, and output files to compare.
    """
    tolerance_files: List[str] = list_tolerance_files_in_directory(full_ref_dir)
    tolerances: dict = load_tolerances(full_ref_dir, tolerance_files)

    files_in_directory = next(os.walk(full_ref_dir))[2]
    reference_outputs = update_wildcard_files_under_test(files_under_test, files_in_directory)

    tolerances = update_wildcard_tolerances(tolerances, reference_outputs)
    tolerances_without_units = strip_tolerance_units(tolerances)

    return tolerances_without_units, reference_outputs
