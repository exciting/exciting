'''
This script creates a new test case for the test suite.
'''
import sys
import os
import argparse as ap

sys.path.insert(1, 'tools')
from procedures import newTest_dir, copyInputFiles, runSingleReference, create_init

def optionParser():
    p = ap.ArgumentParser(description="Usage: python3 newtestcase.py -n <name> -i <reference calc> -e <executable> -np <NP>")
    p.add_argument ('-n', metavar='--name', help="Name of the new test case.",\
                type=str, default=None)
    p.add_argument ('-d', metavar='--description', help="Description of the new test case.",\
                type=str, default=None)
    p.add_argument ('-r', metavar='--reference', help="Path to the reference caclculation.",\
                type=str, default=None)
    p.add_argument ('-i', metavar='--init_file', help="Init file for the new test case. The available init files can be found at '/xml/init_templates'. Default is init_groundstate.xml.",\
                type=str, default=None)
    p.add_argument('-e', metavar='--executable', help="Executable for test run or rerunning references.", \
                type=str, default=None, choices=['excitingser', 'excitingmpi'])
    p.add_argument('-np', metavar='--NP', help='Number of cores for MPI run. Can only be used in combination with excitingmpi as executable.', \
                type=int, default=None)
    
    args = p.parse_args()
    name = args.n
    description = args.d
    inputLoc = args.r
    init = args.i
    executable = args.e
    np = args.np

    if init == None:
        init = 'init_groundstate.xml'
    
    if np != None:
        if executable=='excitingmpi':
            executable = 'mpirun -np %i %s'%(np,executable)
        else:
            raise ValueError('np can not be set in combination with excitingser as executable')
    return name, description, inputLoc, init, executable

    
        



def main():
    ## some global definitions ####################################
    testFarm = 'test_farm' 
    species = '../species'
    mainOut  = 'INFO.OUT'                                                                                      #directory of the testfarm
    runDir = 'run'                                                                                                
    refDir = 'ref'
    inputDef = 'input.xml'
    notRef = ['input.xml', 'STATE.OUT', 'OCCSV.OUT', 'EVECSV.OUT', 'EVECFV.OUT', 'EVALSV.OUT', 'EVALFV.OUT']     #Files that will not be saved as reference
    notClean = ['input.xml']                                                                                     #Files that will not be removed while cleaning
    ###############################################################
    speciesFiles = next(os.walk(species))[2]
    notClean = notClean+speciesFiles

    name, description, inputLoc, init, executable = optionParser()
    if description==None:
        raise ValueError("Please specify a description for the new test case.")
        return
    if name==None:
        raise ValueError("Please specify a name for the new test case.")
        return
    path = os.path.join(testFarm, name)
    if (os.path.exists(path)):
        raise OSError('This directory is already existing. Please choose an other name.')
        return
    
    if inputLoc==None:
        print('A new test case will be created without input files.')
        newTest_dir(testFarm, name, runDir, refDir)
        create_init(testFarm, name, description, init)
    else:
        if not os.path.exists(os.path.join(inputLoc, inputDef)):
            raise OSError('%s does not exist.'%os.path.join(inputLoc, inputDef))
            return
        newTest_dir(testFarm, name, runDir, refDir)
        create_init(testFarm, name, description, init)
        copyInputFiles(inputLoc, testFarm, name, runDir, refDir, inputDef, speciesFiles)
        
        while(True):
            answer = input("Do you want to run the reference calculation? (This is only possible if speciespath is set correct) (y/n)")
            if answer=='y':
                break
            elif answer=='n':
                sys.exit()
        runSingleReference(testFarm, mainOut, name, refDir, executable, notRef)

if __name__ == "__main__":
    main()  
