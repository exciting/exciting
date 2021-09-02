"""
Functions for operating on the test suite's JSON-formatted tolerance files
"""
import os
from typing import List, Callable, NoReturn, Union
import json
import pathlib


def extend_keys(tolerances: dict, new_entries: dict) -> dict:
    """
    Add new entries to existing tolerances (mutates input)

    :param dict tolerances: Tolerances associated with a given output file
    :param dict new_entries: New entries to add to the tolerances

    :return: dict tolerances: Input tolerances, with new entries added
    """
    tolerances.update(new_entries)
    return tolerances


def remove_keys(tolerances: dict, entries_to_remove: List[str]) -> dict:
    """
    Remove existing entries from existing tolerances (mutates input)

    :param dict tolerances: Tolerances associated with a given output file
    :param List[str] entries_to_remove: Entries to remove from the tolerances
    Each element of the list is a key.

    :return: dict tolerances: Input tolerances, with entries_to_remove removed
    """
    for entry_key in entries_to_remove:
        del tolerances[entry_key]
    return tolerances


def modify_tolerance_file_keys(files: List[str],
                               modified_entries: dict,
                               key_modifier: Callable[[dict, Union[list, dict]], dict],
                               treat_original: Callable[[str], NoReturn],
                               ) -> list:
    """
    Modify key:value pairs in a list of tolerance `files`, according to `key_modifier`

    If an exciting file gets a new output value, all existing tests/tolerance files should (probably) have that key:tol
    pair added to them. This is a drawback of pure regression testing, which compares an output to a reference output.
    If any of these things occurs:

     * addition of new information
     * removal of old information
     * output labels/keys change

    either the files containing the tolerances must be updated or the test framework must be robust enough to skip
    over missing or depreciated tolerance keys (which may not even be desirable behaviour).

    This routine expects each tolerance file to have two levels of nesting, of the form:

    ground_state_tol.json =
       {'info_out':
          {'Lattice vectors (cartesian)': 1.e-8,
           'Reciprocal lattice vectors (cartesian)': 1.e-8
          }
        'some_other_gs_output_file':
          {'some other data': 1.e-7}
       }

    :param List[str] files: List of files to parse and modify, defined with relative paths.

    :param dict modified_entries: keys to modify in existing tolerance files, with the expected form:

        new_entries[file_name_tag] = {'some_new_quantity': 1.e-8} when adding a key:tolerance
        or
        new_entries[file_name_tag] = ['some_depreciated_quantity'] when removing a key:tolerance

    :param Callable[[dict, list], dict] key_modifier: Function for key modification

    :param Callable[[str], NoReturn] treat_original: Function defining how existing tolerance files should be treated

    :return List[str] missed: List of files that were specified in 'files' but were not modified because they could not
    be found/do not exist
    """

    missed = []
    for file in files:

        try:
            full_file_path = pathlib.Path(file).resolve()
        except FileNotFoundError:
            missed.append(file)
            continue

        with open(full_file_path) as fid:
            tolerance_groups = json.load(fid)

        assert tolerance_groups.keys() == modified_entries.keys(), \
            'Top-level keys of input tolerance file inconsistent with top-level keys of modified_entries'

        treat_original(full_file_path.name)

        for file_name_tag, tolerances in tolerance_groups.items():
            tolerance_groups[file_name_tag] = key_modifier(tolerances, modified_entries[file_name_tag])

        with open(full_file_path, "w") as fid:
            json.dump(tolerance_groups, fid, indent=2)

    return missed


def extend_tolerance_file_keys(files: List[str],
                               new_entries: dict,
                               rename_original=lambda name: os.rename(name, 'orig_' + name)
                               ) -> list:
    """
    Extend tolerance files adding key:value pairs to tolerance files, specified in `files`

    :param List[str] files: List of files to parse and extend, defined with relative paths.

    :param dict new_entries: key:tol pairs to add to the existing tolerance files, with the expected form:

        new_entries = {'file_name_tag':
                        {'some_new_quantity': 1.e-8}
                      }
        when adding a key:tolerance

    :param function rename_original: Behaviour for dealing with the original tolerance file in the file system.
     Default behaviour creates a copy of the original tolerances file named `'orig_' + name`.
     To overwrite the original tolerance file, pass `rename_original = lambda name: name`, which
     does nothing.

    :return List[str] missed: List of files that were specified in 'files' but were not modified because they could not
    be found/do not exist.
    """
    assert all(isinstance(value, dict) for value in new_entries.values()), \
        "All values in new_entries should be dictionaries"

    return modify_tolerance_file_keys(files, new_entries, extend_keys, rename_original)


def remove_tolerance_file_keys(files: List[str],
                               entries_to_remove: dict,
                               rename_original=lambda name: os.rename(name, 'orig_' + name)
                               ) -> list:
    """
    Remove depreciated keys from the tolerance files, specified in `files`

    :param List[str] files: List of files to parse and remove keys from, defined with relative paths.

    :param dict entries_to_remove: keys to remove from the existing tolerance files, with the expected form:

        entries_to_remove = {'file_name_tag':
                               ['some_depreciated_quantity']
                            }

    :param function rename_original: Behaviour for dealing with the original tolerance file in the file system.
     Default behaviour creates a copy of the original tolerances file named `'orig_' + name`.
     To overwrite the original tolerance file, pass `rename_original = lambda name: name`, which
     does nothing.

    :return List[str] missed: List of files that were specified in 'files' but were not modified because they could not
    be found/do not exist.
    """
    assert all(isinstance(value, list) for value in entries_to_remove.values()), \
        "All values in entries_to_remove should be lists"

    return modify_tolerance_file_keys(files, entries_to_remove, remove_keys, rename_original)
