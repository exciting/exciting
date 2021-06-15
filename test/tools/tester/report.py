import xml.etree.ElementTree as ET

from .test import Test
from .failure import *

def indent(elem, level=0):
    i = "\n" + level*"  "
    j = "\n" + (level-1)*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for subelem in elem:
            indent(subelem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = j
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = j
    return elem        

class Report(dict):
    """
    Collects the fail lists of a test run and writes a report.
    In:
        name            string      Gives the thing a name. Important writing a report.
        description     string      Description of the tests that are performed.
    """
    def __init__(self, name, description):
        self.name = name
        self.description = description
        self.numberOfFails = None
        self.failedTests = 0
        self.passedTests = 0
        self.succeededTests = 0
        self.numberOfTests = 0
        super(Report, self).__init__()

    def updateNumberOfFails(self, newFailNumbers):
        """
        Updates the number of failures.
        """
        if self.numberOfFails==None:
            self.numberOfFails = newFailNumbers
        else:
            for key in self.numberOfFails.keys():
                self.numberOfFails[key] += newFailNumbers[key]

    def collectTest(self, Test):
        """
        Adds a new Test to the report. Automatically updates the attributes.
        """
        self[Test.name]=Test
        self.updateNumberOfFails(Test.countFailures())
        # TODO(Alex) This should be assertions, for a given test case 
        self.numberOfTests = self.numberOfTests+1
        if Test.succeeded:
            self.succeededTests = self.succeededTests+1
        elif Test.passed:
            self.passedTests = self.passedTests+1
        else:
            self.failedTests = self.failedTests+1

    def runFailed(self, name, err_msg):
        test = Test(name=name)
        test.append(Failure(Failure_code.RUN, err_msg=err_msg))
        self.collectTest(test)

    def writeToTerminal(self):
        """
        Write a report to stdout
        """
        print('Report:', end = ' ')
        print_color('SUCCESS %i/%i'%(self.succeededTests, self.numberOfTests), 'green', end = '')
        print(',', end=' ') 
        print_color('PASS %i/%i'%(self.passedTests, self.numberOfTests), 'yellow', end = '')
        print(',', end=' ')
        print_color('FAIL %i/%i' %(self.failedTests, self.numberOfTests), 'red', end = '')
        print('.')
        
        for key in self.keys():
            self[key].printFailList()
    
    def writeToXML(self, directory, FileName='report.xml'):
        root = ET.Element('report')
        
        name = ET.SubElement(root, 'name')
        name.text = self.name

        description = ET.SubElement(root, 'description')
        description.text = self.description
        
        for key in self.keys():
            self[key].XMLoutput(root, directory)
        
        indent(root)
        tree = ET.ElementTree(root)
        tree.write(FileName)
        
    def assert_errors(self, handle_errors:bool):
        """
        Throws if the any assertions for a test got passed over or failed
        (where passed over means some reference data was missing)

        :param handle_errors  If true, allow any errors to propagate 
                              to the end of the test suite 
        """
        if not handle_errors:
            assert self.failedTests == 0, "Test suite assertion(s) failed"
            assert self.passedTests == 0, "Test suite reference data requires regenerating"
        return 


def skipped_test_summary(skipped_tests:list):
    """
    Summarises skipped tests 

    :param skipped_tests: list of skipped tests. 
           Each element is a dictionary of the skipped test properties 
    :return None  
    """

    if len(skipped_tests) == 0:
        return
    else: 
        print('Summary of SKIPPED tests:')
        for test in skipped_tests:
            print(' ', test['name'], '. ', test['comment'])


def test_suite_summary(test_suite_report:list):
    """
    Summarises the test suite assertions.
    Writes summary to stdout
    Skipped test cases are not included in the report

    :param test_suite_report: list of test reports, each of class Report
    :return bool indicating if all assertions in all test cases succeeded  
    """

    succeeded_asserts = 0
    passed_asserts = 0
    failed_asserts = 0 
    failed_tests = 0

    for report in test_suite_report:
        assert isinstance(report, Report), "entry in test_suite_summary is not type 'Report'"
        succeeded_asserts += report.succeededTests
        passed_asserts += report.passedTests
        failed_asserts += report.failedTests
        if report.failedTests > 0:
            failed_tests += 1 

    total_asserts = succeeded_asserts +  passed_asserts + failed_asserts
    total_tests = len(test_suite_report)
    succeeded_tests = total_tests - failed_tests

    print('')
    print('Summary:')
    print('Tests succeeded: %i/%i' % (succeeded_tests, total_tests))
    print('Assertions succeeded: %i/%i' % (succeeded_asserts, total_asserts))

    return failed_asserts == 0 


def timing_summary(timing:dict, verbose=False):
    """
    Timings summary for the test suite

    :param timing: Dict of timings for each test case 
    :param verbose: Give additional timings, useful for developers
    """

    times = [time for time in timing.values()]
    total_time = sum(times)
    print('Total test suite time (mins) : %.1f' % (total_time / 60.))

    if verbose:
        if len(times) > 0:
            avg_time = total_time / len(times)
        else:
            avg_time = 0.
        longest_time = 0. 
        longest_time_name = ''
        for name, time in timing.items():
            if time > longest_time:
                longest_time = time 
                longest_time_name = name 

        print('Average test time (s): %.1f' % avg_time)
        print('Longest test time (s): %.1f, taken by %s' % (longest_time, longest_time_name))
