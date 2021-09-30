import os
import shutil
from xml.etree.ElementTree import ParseError
from typing import List, Optional

from .parsers import parseInit


def copy_exciting_input(source: str, destination: str, species_files: List[str],
                        input_file: Optional[str] = 'input.xml'):
    """
    Copy input files for an exciting calculation from source to destination.

    The input files are the input.xml and the species files.

    :param str source: Path to the exciting calculation where to copy the input files from
    :param str destination: Path to the directory where the input files shall copied to.
    :param List[str] species_files: Exciting species files
    :param Optional[str] input_file: Exciting input file
    """
    if not os.path.isdir(destination):
        raise NotADirectoryError('Target directory does not exist:', destination)

    species_required = has_species_files(source, species_files)

    files_to_copy = [input_file] + species_required

    for file in files_to_copy:
        shutil.copyfile(os.path.join(source, file), os.path.join(destination, file))


def has_species_files(directory: str, species_files: List[str]) -> List[str]:
    """
    Raise exception if no species file is in directory.

    :param str directory: Path to directory to check
    :param List[str] species_files: List of possible species files
    :return List[str]: Species files present in directory
    """
    files_in_directory = next(os.walk(directory))[2]

    species_present = list(set(files_in_directory) & set(species_files))

    if len(species_present) == 0:
        raise FileNotFoundError
    else:
        return species_present


def create_run_dir(path_to_test_case: str, run_dir: str):
    """
    Create run directory for a test case with the name test_name in test_farm.
    If a run directory already exists, it will be deleted.

    :param str path_to_test_case:  path to the test case
    :param str run_dir:    name of the run directory
    """

    # check if the test case exists
    try:
        os.path.isdir(path_to_test_case)
    except Exception:
        raise NotADirectoryError('%s does not exist.', path_to_test_case)

    if os.path.isdir(os.path.join(path_to_test_case, run_dir)):
        try:
            shutil.rmtree(os.path.join(path_to_test_case, run_dir))
        except Exception:
            raise OSError('%s could not be removed.' % os.path.join(path_to_test_case, run_dir))

    os.mkdir(os.path.join(path_to_test_case, run_dir))


def get_test_from_init(path_to_test_case: str, init_file: str) -> dict:
    """
    Read files for testing, and corresponding tolerances, from init.xml.

    :param path_to_test_case: path to the test case
    :param init_file:  Name of the init file. By default init.xml.
    """
    try:
        init = parseInit(os.path.join(path_to_test_case, 'init.xml'))
    except FileNotFoundError:
        init = parseInit(init_file)
    except Exception:
        raise ParseError('Could not parse init file.')

    return init


def flatten_directory(path: str):
    """
    Flatten the file structure by factor 1 for the directory at path.
    :param path: path to directory which gets flattened
    """
    try:
        os.path.isdir(path)
    except Exception:
        raise NotADirectoryError('%s does not exist.' % path)

    dirs = next(os.walk(path))[1]

    for dir in dirs:
        files = next(os.walk(os.path.join(path, dir)))[2]
        for file in files:
            try:
                shutil.move(os.path.join(path, dir, file), path)
            except shutil.Error:
                os.remove(os.path.join(path, file))
                shutil.move(os.path.join(path, dir, file), path)

        shutil.rmtree(os.path.join(path, dir))


def remove_ignored_files(path: str, ignored_output: list):
    """
    Remove files that are ignored for the tests in a directory.
    :param path:              Path to the directpry 
    :param ignored_output:    Files that are ignored by the tests.
    """
    try:
        files = next(os.walk(path))[2]
    except Exception:
        return NotADirectoryError('%s is not a directory.' % path)

    for file in files_to_remove(files, ignored_output):
        os.remove(os.path.join(path, file))


def files_to_remove(files: List[str], remove_list: List[str]) -> List[str]:
    """
    Return the intersection between a list of files and a list of files that shall be removed.

    :param List[str] files: List of files
    :param List[str] remove_list: List of files to be removed
    :return List[str] to_remove: List of files to be removed
    """
    to_remove = []
    for remove in remove_list:
        for file in files:
            if remove in file:
                to_remove.append(file)

    return to_remove


def create_init(test_location: str, description: str, init_name: str):
    """
    Copy init_default.xml from xml/init_templates to the directory
    of the new test case.

    :param str test_location: Test name
    :param str description: Test description
    :param str init_name: Name of init template file
    """
    init_template = os.path.join("xml/init_templates", init_name)
    init_file = os.path.join(test_location, "init.xml")

    shutil.copy(init_template, init_file)

    with open(init_file, 'r') as file:
        lines = file.read()
    lines = lines.replace('test_name', os.path.basename(test_location))
    lines = lines.replace('test_description', description)

    with open(init_file, 'w') as file:
        file.write(lines)
