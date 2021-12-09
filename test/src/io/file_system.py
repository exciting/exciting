"""
Reading, writing and interacting with the file system.
"""
import os
import shutil
from typing import List


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
    except NotADirectoryError:
        raise NotADirectoryError('%s does not exist.', path_to_test_case)

    if os.path.isdir(os.path.join(path_to_test_case, run_dir)):
        try:
            shutil.rmtree(os.path.join(path_to_test_case, run_dir))
        except OSError:
            raise OSError('%s could not be removed.' % os.path.join(path_to_test_case, run_dir))

    os.mkdir(os.path.join(path_to_test_case, run_dir))


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


def remove_ignored_files(path: str, ignored_output: list):
    """
    Remove files that are ignored for the tests in a directory.

    :param str path:              Path to the directory
    :param list ignored_output:    Files that are ignored by the tests.
    """
    try:
        files = next(os.walk(path))[2]
    except NotADirectoryError:
        raise NotADirectoryError('%s is not a directory.' % path)

    for file in files_to_remove(files, ignored_output):
        os.remove(os.path.join(path, file))


def copy_calculation_inputs(source: str, destination: str, input_files: List[str]):
    """
    Copy input files for a calculation from source to destination.
    Performs checking to ensure all inputs exist.

    :param str source: Path to the exciting calculation where to copy the input files from
    :param str destination: Path to the directory where the input files shall copied to.
    :param List[str] input_files: Input files required to perform a calculation
    """
    inputs_present = input_files_in_directory(source, input_files)
    missing_files = set(input_files) - set(inputs_present)

    if missing_files:
        raise FileNotFoundError(f'Required input file(s) are missing from the source directory: {missing_files}')

    if not os.path.isdir(destination):
        raise NotADirectoryError('Target directory does not exist:', destination)

    for file in input_files:
        shutil.copyfile(os.path.join(source, file), os.path.join(destination, file))


def input_files_in_directory(directory: str, files: List[str]) -> List[str]:
    """
    Return files which are present in files argument and directory.

    :param str directory: Path to directory to check
    :param List[str] files: List of possible files
    :return List[str]:  Files present in directory
    """
    files_in_directory = next(os.walk(directory))[2]
    if not files_in_directory:
        raise FileNotFoundError(f'Directory "{directory}" contains no files')

    species_files = set(files_in_directory) & set(files)
    return list(species_files)
