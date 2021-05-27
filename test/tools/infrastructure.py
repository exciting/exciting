import os 
import shutil
from xml.etree.ElementTree import ParseError
from typing import List, Optional

from .parsers import parseInit

def copy_exciting_input(source:str, destination:str, species_files:List[str], input_file:Optional[str]='input.xml'):
    """
    Copy input files for an exciting calculation from source to destination. 
    The input files are the input.xml and the species files. 
    :param source:          path to the exciting calculation where to copy the input files from
    :param destination:     path to the directory where the input files shall copied to.
    :param species_files:   exciting species files
    :param input_file:      exciting input file
    """
    
    # check if input.xml exists in source
    try:
        os.path.isfile(os.path.join(source, input_file))
    except Exception:
        raise FileNotFoundError('No %s in %s.'%(input_file, source))
    
    try:
        os.path.isdir(destination)
    except Exception:
        raise NotADirectoryError('%s does not exist.'%source)

    try:
        has_species_files(source, species_files)
    except Exception:
        raise FileNotFoundError('%s has no species files.'%source)

    files = [input_file] + species_files
    files_in_reference = next(os.walk(source))[2]
    for file in files:
        if file in files_in_reference:
            shutil.copyfile(os.path.join(source, file), os.path.join(destination, file))

def has_species_files(directory:str, species_files:list):
    """
    Raise exception if no species file is in directory.
    :param directory:       path to directory to check
    :param species_fiels:   list of possible species files
    """

    files_in_directory = next(os.walk(directory))[2]

    num_of_species_in_directory = len(list(set(files_in_directory) & set(species_files)))

    if num_of_species_in_directory == 0:
        raise FileNotFoundError


def create_run_dir(path_to_test_case:str, run_dir:str):
    """
    Create run directory for a test case with the name test_name in test_farm. If a run directory already exists, it will be deleted.
    :param test_farm:  path to the test case
    :param test_name:  name of the test
    :param run_dir:    name of the run directory
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
            raise OSError('%s could not be removed.'%os.path.join(path_to_test_case, run_dir))
    
    os.mkdir(os.path.join(path_to_test_case, run_dir))

def get_test_from_init(path_to_test_case:str, init_file:str)->dict:
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

def flatten_directory(path:str):
    """
    Flatten the file structure by factor 1 for the directory at path.
    :param path: path to directory which gets flattened
    """
    try:
        os.path.isdir(path)
    except Exception:
        raise NotADirectoryError('%s does not exist.'%path)

    dirs = next(os.walk(path))[1] 

    for dir in dirs:
        files = next(os.walk(os.path.join(path, dir)))[2]
        for file in files:
            shutil.move(os.path.join(path, dir, file), path)
        
        shutil.rmtree(os.path.join(path, dir))

def remove_ignored_files(path:str, ignored_output:list):
    """
    Remove files that are ignored for the tests in a directory.
    :param path:              Path to the directpry 
    :param ignored_output:    Files that are ignored by the tests.
    """
    try:
        files = next(os.walk(path))[2]
    except Exception:
        return NotADirectoryError('%s is not a directory.'%path)

    files_for_deleting = set(files) & set(ignored_output)
    for file in files_for_deleting:
        os.remove(os.path.join(path, file))

def create_test_directory(test_farm:str, name:str, ref_dir:str):
    """
    Create a directory for a new test case.
    :param testFarm:    location of the test farm
    :param name:        name of the new test case
    :param runDir:      name of the run directory of the test case
    :param refDir:      name of the ref directory of the test case
    """
    try:
        not os.path.isdir(os.path.join(test_farm, name))
    except Exception:
        raise FileExistsError

    path = os.path.join(test_farm, name)

    os.makedirs(path)
    os.makedirs(os.path.join(path, ref_dir))

def create_init(testFarm:str, name:str, description:str, init:str):
    """
    Copy init_default.xml from xml/init_templates to the directory of the new test case.
    :param testFarm:     location of the test farm   
    :param name:         name of the new test case
    :param description:  description in the init file
    :param init:         init file template
    """
    shutil.copy(os.path.join("xml/init_templates", init), os.path.join(testFarm, name, "init.xml"))
    with open(os.path.join(testFarm, name, "init.xml"), 'r') as file :
        lines = file.read()
    lines = lines.replace('test_name', name)
    lines = lines.replace('test_description', description)
    file.close
    with open(os.path.join(testFarm, name, "init.xml"), 'w') as file :
        file.write(lines)
    file.close