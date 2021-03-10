'''
This script creates a new test case for the test suite.
'''
import sys
import os
import argparse as ap
from collections import namedtuple
import warnings

sys.path.insert(1, 'tools')
from procedures import newTest_dir, copyInputFiles, runSingleReference, create_init
from constants import settings

def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return '%s: %s \n' % (category.__name__, message)
warnings.formatwarning = warning_on_one_line

def optionParser(exedir:str):
    p = ap.ArgumentParser(description="Usage: python3 newtestcase.py -n <name> -i <reference calc> -e <executable> -np <NP>")
    p.add_argument ('-n',
                   metavar = '--name',
                   help = "Name of the new test case.",
                   type = str,
                   default=None)
    p.add_argument ('-d', 
                   metavar='--description', 
                   help="Description of the new test case.",
                   type=str, 
                   default=None)
    p.add_argument ('-r', 
                    metavar='--reference-input', 
                    help="Path to an existing caclculation, that will be the reference.",
                    type=str, 
                    default=None)
    p.add_argument ('-i', 
                    metavar='--init_file', 
                    help="Init file for the new test case. The available init files can be found at '/xml/init_templates/'. Default is init_groundstate.xml.",\
                    type=str, 
                    default='init_groundstate.xml')

    args = p.parse_args()

    input_options = {'name': args.n, 
                     'description': args.d, 
                     'reference_location': args.r, 
                     'init_file': args.i
                     }
    
    return input_options

def interactiveUserInterface(args:dict):
    """
    Interactive user interface for creating a new test case. 
    Interaction is only triggered if the elements for 
    input_options "name", "reference_input", "description"
    are not specified or not compatible.

    :param args: parsed command line arguments.
    """
    input_options = args
    if args['name']==None:
        name = input("Please enter a name for the new test case: \n")
    else:
        name = args['name']
    while os.path.exists(os.path.join(settings.test_farm, name)):
        name_old = name
        name = input("A test case with the same name exists already. " + \
                     "Please choose a different name, replace the existing test case (enter 'replace', this will remove the existing test case!!)" + \
                     " or quit (enter 'exit').\n")
        if name=='replace':
            os.system('rm -r %s/%s'%(settings.test_farm,name_old))
            name = name_old
            break
        elif name.lower() =='exit':
            sys.exit()
    input_options['name'] = name

    if args['description']==None:
        input_options['description'] = input("Please enter a description for the new test case:\n")
    else:
        input_options['description'] = args['description']

    reenter = True
    if args['reference_location']==None:
        reference_location = input("If you want to set up a test case with an existing calculation as reference, " + \
                                   "please enter the path to the corresponding directory. \n" +\
                                   "Make sure that a valid %s and all necessary species files are contained.\n"%settings.input_file + \
                                   "If you want to set up an empty test case, please enter 'no'.\n")
        if reference_location.lower() in ["no","n"]:
            reference_location = None
            reenter = False
    else:
        reference_location = args['reference_location']
    if reenter:
        while not os.path.exists(os.path.join(reference_location, settings.input_file)):
            reference_location = input("The reference path you entered does not exists. Please reenter, create a blank test case (type 'no') " + \
                                    "or exit (type 'exit').\n")
            if reference_location.lower() in ["no","n"]:
                reference_location = None
                break
            elif reference_location.lower() == "exit":
                sys.exit()
    input_options['reference_location'] = reference_location

    return input_options
    

def main(settings:namedtuple, input_options:dict):
    """
    Create a new test case from an input file and run the reference.

    :stettings:        default settings
    :input_options:    definitions of user input or parsed command line arguments
    """
    speciesFiles = next(os.walk(settings.species))[2]
    notClean = settings.not_clean + speciesFiles

    newTest_dir(settings.test_farm, 
                input_options['name'], 
                settings.run_dir, 
                settings.ref_dir)
    create_init(settings.test_farm, 
                input_options['name'], 
                input_options['description'], 
                input_options['init_file'])
    if input_options['reference_location']==None:
        print("Create blank test case %s. Please put the input files for your test case in the directories "%input_options["name"] + \
              "%s/ and %s/ in the test directory %s/%s. \n"%(settings.run_dir, settings.ref_dir, settings.test_farm, input_options['name']))
        warnings.warn("Before you can run the test case, you must generate the reference data.")
        print("Run the reference: \n\n" + \
              "    python3 runtest.py -a ref -t %s\n"%input_options['name'])
    else:
        print('Create new test case %s.'%(input_options['name']))        
        print("Take reference input from %s."%input_options['reference_location'])
        copyInputFiles(input_options['reference_location'], 
                       settings.test_farm, 
                       input_options['name'], 
                       settings.run_dir, 
                       settings.ref_dir, 
                       settings.input_file, 
                       speciesFiles)
            
        runSingleReference(settings.test_farm, 
                           settings.main_output, 
                           input_options['name'], 
                           settings.ref_dir,
                           os.path.join(settings.exe_dir, settings.exe_ref),
                           settings.ignored_output, 
                           settings.max_time)
    
    print("Run the test case: \n\n" + \
          "    python3 runtest.py -t %s\n"%input_options['name'])
if __name__ == "__main__":
    args = optionParser(settings.exe_dir)
    input_options = interactiveUserInterface(args)
    main(settings, input_options)  
