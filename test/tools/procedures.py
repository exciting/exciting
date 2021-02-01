'''
Procedures for testing
'''

import sys
import os
import shutil
import xml.etree.ElementTree as ET 
from subprocess import PIPE, CalledProcessError, check_call, Popen, TimeoutExpired
import glob
import time

sys.path.insert(1, 'tools/parser')
from ErrornousFileError import ErrornousFileError
from initParser import parseInit, getInitFile
from parserChooser import parserChooser

sys.path.insert(1, 'tools/tester')
from test import Test, fromInit
from report import Report, indent, test_suite_summary, skipped_test_summary, timing_summary
from failure import *

#sys.path.insert(1, 'xml/schema')        If we want to validat init.xml, the stuff is there.
#from validate import validate


def collectSingleReport(testFarm, testDir, root):
    '''
    Collects content from report.xml for a single test case.
    Input:
        testFarm        string              location of the test farm
        testDir         string              test case for that report.xml will be collected
        root            ET root element     root element for reports.xml in report/
    '''
    try:
        rootT = ET.parse(os.path.join(testFarm, testDir, 'report.xml')).getroot()
        root.insert(1, rootT)
        return
    except FileNotFoundError:
        return

def collectReports(testFarm, testList):
    '''
    Collects content from report.xml for test cases in test list (see collectSingleReport)
    Input:
        testFarm        string              location of the test farm
        testList        list of string      test cases for that report.xml will be collected
    '''
    root = ET.Element("reports")
    for testDir in testList:
        collectSingleReport(testFarm, testDir, root)
    indent(root)
    tree = ET.ElementTree(root)
    tree.write('report/reports.xml')

    os.system('xsltproc report/reports2html.xsl report/reports.xml > report/report.html')

def exciting_run(executable, mainOut, maxTime):
    '''
    Executes a exciting run, checks if it was successfull.
    Input:
        executable  string    program that will be executed (excitingser, excitingmpi, mpirun -np NP excitingmpi)
        mainOut     string    main output file (INFO.OUT)
        maxTime     int       maximum time for an exciting run in seconds
    Output:
        success     bool                true if run was successfull, false else
        errMess     list of strings     terminal output of exciting
        run_time    float     Run time of job 
    '''
    t_start = time.time()
    exciting_run = Popen(executable.split(), stdout = PIPE)
    
    try:
        errMess = exciting_run.communicate(timeout=maxTime)[0]
        runSucc = os.path.isfile(mainOut)
    except TimeoutExpired:
        exciting_run.kill()
        errMess = 'Time expired.'
        runSucc = False
    
    exciting_run.wait()
    t_end = time.time()
    
    return runSucc, errMess, t_end-t_start


def runSingleTest(testFarm:str, mainOut:str, testDir:str, runDir:str,
                  refDir:str, init_default:str, executable:str, maxTime:str, 
                  timing:dict, handle_errors:bool):
    '''
    Runs a singel test.
    Input:
        testFarm        string      location of the test farm
        mainOut         string      main output file of the exciting calculation
        testDir         string      test case that will be ran
        runDir          string      name of the run directory of the test case
        refDir          string      name of the ref directory of the test case
        init_default    string      location of the default init.xml
        execuatable     string      executable for the exciting run
        timing          dict        test run times in seconds
        handle_errors   bool        Whether or not failures and passes are allowed to propagate

    Output:
        report          object      report instance of the test     
    '''
    print('Run test %s:'%testDir)
    os.chdir(os.path.join(testFarm,testDir))
    try:
        init = parseInit(getInitFile(testDir))
    except FileNotFoundError:
        init = parseInit(os.path.join('../..', init_default))
    except OSError:
        os.chdir('../../')
        return
        
    name = init['name']
    description = init['description']
    tests = init['tests']

    # Initialize report object:
    report = Report(name, description)

    # The actual test run is performed. If the run fails, then an error message is added
    # and all other tests are skipped.
    os.chdir(runDir)
    runSuc, errMess, timing[name] = exciting_run(executable, mainOut, maxTime)
    
    if runSuc==False:
        report.runFailed(testDir, errMess)
        os.chdir('../')
    else:
        print('Run succeeded')
        os.chdir('../')

        # test all the files specified in init.xml:
            
        for testInit in tests:
            testFile = testInit['file']
            runPath = os.path.join(runDir, testFile)
            # Handles the case, if files are saved in the run directory or in a sub directory of the run directory
            if len(os.path.split(testFile)[0])>0:
                testDirSub = os.path.split(testFile)[0]
                testFile = os.path.split(testFile)[1]
                refDirSub = os.path.join(refDir, '%s_REF'%testDirSub)
                refPath = os.path.join(refDirSub, '%s.ref'%testFile)
            else:
                refPath = os.path.join(refDir, '%s.ref'%testFile)

            test = fromInit(testInit)
            # Check if reference file exists
            if not os.path.isfile(refPath):
                test.append(Failure(Failure_code.REFERENCE, err_msg=testInit['file']))           
                report.collectTest(test)
                continue
            # Checks if the file exists and is not broken
            try:
                runData = parserChooser(runPath)
            except OSError:
                test.append(Failure(Failure_code.FILENOTEXIST, err_msg=testInit['file']))
                report.collectTest(test)
                continue
            except ErrornousFileError:
                test.append(Failure(Failure_code.ERRORFILE, err_msg=testInit['file']))
                report.collectTest(test)
                continue

            if '_REF' in refPath:
                os.rename(refPath, os.path.join(refDirSub, testFile))
                refPath = os.path.join(refDirSub, testFile)
                refData = parserChooser(refPath)
                os.rename(refPath, os.path.join(refDirSub, '%s.ref'%testFile))
            else:
                os.rename(refPath, os.path.join(refDir, testFile))
                refPath = os.path.join(refDir, testFile)
                refData = parserChooser(refPath)
                os.rename(refPath, os.path.join(refDir, '%s.ref'%testFile))

            if 'info.xml' in testFile:
                test.evaluate_info(runData.data, refData.data)
            elif 'INFO.OUT' in testFile and 'WANNIER' not in testFile:
                test.evaluate_INFO(runData.data, refData.data)
            elif 'eigval.xml' in testFile:
                test.evaluate_eigval(runData.data, refData.data)
            else:                                                                               
                test.evaluate(runData.data, refData.data)

            report.collectTest(test)

    print('Time (s): %.1f' % timing[name])
    report.writeToTerminal()
    report.writeToXML(testDir)
    report.assert_errors(handle_errors)

    os.chdir('../../')
    return report 

def runTests(testFarm:str, mainOut:str, testList:list, runDir:str, refDir:str,
             init_default:str, executable:str, np:int, omp:int, maxTime:int,
             skipped_tests:list, handle_errors:bool):
    '''
    Runs tests in testList (see runSingleTest).
    Input:
        testFarm        string              location of the test farm
        mainOut         string              main output file of the exciting calculation
        testList        list of string      test cases that will be ran
        runDir          string              name of the run directory of the test case
        refDir          string              name of the ref directory of the test case
        init_default    string              location of the default init.xml
        executable      string              executable command for the exciting run
        np              int                 number of MPI processes
        omp             int                 number of OMP threads
        maxTime         int                 max time before a job is killed
        skipped_tests   list                list of tests to skip
        handle_errors   bool                Whether or not failures and passes are allowed to propagate
    '''
    if 'excitingser' in executable:
        print('Run tests with excitingser.')
    elif 'excitingmpismp' in executable:
        print('Run tests with excitingmpismp with %i open MP threads and %i MPI processes.'%(omp, np))
    elif 'excitingmpi' in executable:
        print('Run tests with excitingmpi %i MPI processes.'%np)

    testList = remove_tests_to_skip(testList, skipped_tests)
    timing = {}
    test_suite_report = []
    
    for testDir in testList:
        test_suite_report.append(
            runSingleTest(testFarm, 
                          mainOut, 
                          testDir, 
                          runDir, 
                          refDir, 
                          init_default, 
                          executable, 
                          maxTime, 
                          timing,
                          handle_errors)
                                )

    all_asserts_succeeded = test_suite_summary(test_suite_report)
    skipped_test_summary(skipped_tests)
    timing_summary(timing, verbose=True)
    assert all_asserts_succeeded, "Some test suite assertions failed or were passed over"

    return 
        

def isNotForbiddenFile(f, forbiddenFiles):
    '''
    Returns True if a file is not in forbidden files.
    Input:
        f               string          file name
        forbiddenFiles  list of string  list of forbidden files
    Output:
        out             bool
    '''
    out = True
    for fF in forbiddenFiles:
        out = out and (fF not in f)
    return out

def runSingleReference(testFarm, mainOut, testDir, refDir, executable, forbiddenFiles, maxTime):
    '''
    Reference run for single test case.
    Input:
        testFarm        string              location of the test farm
        mainOut         string              main output file of the exciting calculation
        testDir         string              test case for that the reference will be calculated
        refDir          string              name of the ref directory of the test case
        execuatable     string              executable for the exciting run
        forbiddenFiles  list of strings     files that will not be saved as reference
    '''

    os.chdir(os.path.join(testFarm,testDir,refDir))

    dirs, files = next(os.walk('.'))[1:]
    for f in files:
        if '.ref' in f:
            os.remove(f)
    for d in dirs:
        if '_REF' in d:
            shutil.rmtree(d)
    runSuc, errMess = exciting_run(executable, mainOut, maxTime)
    dirs, files = next(os.walk('.'))[1:]
    if runSuc:
        for f in files:
            if f not in forbiddenFiles and 'OCCSV' not in f and 'EVEC' not in f and 'EVALFV' and f not in 'EVALSV':
                os.rename(f, '%s.ref'%f)
        for d in dirs:
            if isNotForbiddenFile(f, forbiddenFiles):
                os.chdir(d)
                files = next(os.walk('.'))[-1]
                for f in files:
                    if isNotForbiddenFile(f, forbiddenFiles):
                        os.rename(f, '%s.ref'%f)
                os.chdir('..')
                os.rename(d, '%s_REF'%d)
                
        print('  %s: Run succeeded'%testDir)
    else:
        print('  %s: Run failed'%testDir)
        for err in errMess:
            print(err)
    os.chdir('../../..')

def runReferences(testFarm, mainOut, testList, refDir, executable, forbiddenFiles, np, omp, maxTime):
    '''
    Reference run for all tests (see runSingleReference).
    Input:
        testFarm        string              location of the test farm
        mainOut         string              main output file of the exciting calculation
        testList        string              list of test cases for that the reference will be calculated
        runDir          string              name of the run directory of the test case
        refDir          string              name of the ref directory of the test case
        execuatable     string              executable for the exciting run
        forbiddenFiles  list of strings     files that will not be saved as reference
    '''
    if 'excitingser' in executable:
        print('Run references with excitingser.')
    elif 'excitingmpi' in executable:
        print('Run references with excitingmpi with %i open MP threads and %i MPI processes.'%(omp, np))
    
    if testList == next(os.walk(testFarm))[1]:
        print('Rerun references for all tests in %s.'%testFarm)
    else:
        print('Rerun references for tests {} in {}.'.format(testList, testFarm))
    for testDir in testList:
        runSingleReference(testFarm, mainOut, testDir, refDir, executable, forbiddenFiles, maxTime)

def cleanSingleTest(testFarm, testDir, runDir, refDir, forbiddenFiles):
    '''
    Removes all files from exciting calculation in a single test case except files stored as reference and files defined in forbiddenFiles.
    Input:
        testFarm        string              location of the test farm
        testDir         string              test case that will be cleaned
        runDir          string              name of the run directory of the test case
        refDir          string              name of the ref directory of the test case
        forbiddenFiles  list of strings     files that will not be removed
    '''
    os.chdir(os.path.join(testFarm, testDir, runDir))
    dirs, files = next(os.walk('.'))[1:]
    for f in files:
        if f not in forbiddenFiles:
            os.remove(f)
    for d in dirs:
        if f not in forbiddenFiles:
            shutil.rmtree(d)
    os.chdir(os.path.join('..', refDir))
    dirs, files = next(os.walk('.'))[1:]
    for f in files:
        if f not in forbiddenFiles and '.ref' not in f:
            os.remove(f)
    for d in dirs:
        if f not in forbiddenFiles and '_REF' not in d:
            shutil.rmtree(d)
    os.chdir('../../..')

def cleanTests(testFarm, testList, runDir, refDir, forbiddenFiles):
    '''
    Cleans the tests in testList (see cleanSingleTest).
    Input:
        testFarm        string              location of the test farm
        testList        list of strings     test cases that will be cleaned
        runDir          string              name of the run directory of the test case
        refDir          string              name of the ref directory of the test case
        forbiddenFiles  list of strings     files that will not be removed
    '''
    print('Clean test directories.')
    
    for testDir in testList:
        cleanSingleTest(testFarm, testDir, runDir, refDir, forbiddenFiles)

def newTest_dir(testFarm, name, runDir, refDir):
    '''
    Creates new test case at location path.
    Input:
        testFarm    string      location of the test farm
        name        string      name of the new test case
        runDir      string      name of the run directory of the test case
        refDir      string      name of the ref directory of the test case
    '''
    path = os.path.join(testFarm, name)
    os.makedirs(path)
    os.makedirs(os.path.join(path, runDir))
    os.makedirs(os.path.join(path, refDir))

def copyInputFiles(src, testFarm, name, runDir, refDir, inputInd, species):
    '''
    Copies input.xml and species files from a reference exciting calculation to a test case, located at path.
    Input:
        src         string          source of the reference calculation
        testFarm    string      location of the test farm
        name        string      name of the new test case
        refDir      string          name of the ref directory of the test case
        inputInd    string          indicator of the input file. In case of input.xml, "input" will work.
        species     list od strings list of possible species files
    '''
    files_src = next(os.walk(src))[2]
    files_copy = []
    for f in files_src:
        if inputInd in f:
            files_copy.append(f)
        elif f in species:
            files_copy.append(f)
        else:
            pass
    path = os.path.join(testFarm, name)
    for f in files_copy:
        shutil.copy(os.path.join(src,f), os.path.join(path,runDir,f))
        shutil.copy(os.path.join(src,f), os.path.join(path,refDir,f))

def create_init(testFarm, name, description, init):
    '''
    Copies init_default.xml from xml/init_templates to the directory of the new test case.
    Input:
        testFarm    string      location of the test farm   
        name        string      name of the new test case
    '''
    shutil.copy(os.path.join("xml/init_templates", init), os.path.join(testFarm, name, "init.xml"))
    with open(os.path.join(testFarm, name, "init.xml"), 'r') as file :
        lines = file.read()
    lines = lines.replace('test_name', name)
    lines = lines.replace('test_description', description)
    file.close
    with open(os.path.join(testFarm, name, "init.xml"), 'w') as file :
        file.write(lines)
    file.close


def remove_tests_to_skip(all_tests:list, skipped_tests:list) -> list:
    """
    Remove tests given in 'skipped_tests' from the test suite,
    for a specific executable choice. 

    This is useful if a particular test crashes or hangs, and needs to be
    debugged BUT shoudn't cause the test suite to report a failure.

    :param all_tests:     list of all test names
    :param skipped_tests: list of tests to skip. Each entry is a dict

    :return tests_to_run: list of tests to run (with skipped_tests removed) 
    """

    tests_to_skip = [test['name'] for test in skipped_tests]    
    # Using sets is probably faster but they won't preserve ordering 
    tests_to_run =  []
    for test in all_tests:
        if test not in tests_to_skip:
            tests_to_run.append(test)
    
    return tests_to_run