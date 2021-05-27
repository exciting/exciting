import xml.etree.ElementTree as ET

from ..termcolor_wrapper import print_color
from .compare import compare_data, compare_info, compare_eigval, compare_INFO
from .path import Path
from .failure import Failure_code


def addTestXML(root, tstatus, tname, tdescription, tdirectory):
    test = ET.SubElement(root, 'test')

    status = ET.SubElement(test, 'status')
    status.text = tstatus

    name = ET.SubElement(test, 'name')
    name.text = tname

    description = ET.SubElement(test, 'description')
    for d in tdescription:
        line = ET.SubElement(description, 'line')
        line.text = d    
    
    directory = ET.SubElement(test, 'directory')
    directory.text = tdirectory


class Test(list):
    """
    Sub class of list. Stores Fails of a test.
    In:
        name            string          name of the test
        conditions      dict            Specifies how many tests may fail to pass the test.
    """
    def __init__(self, **kwargs):
        # TODO(A/B/H) Issue 67 Remove Defaults in favour of explicit tolerances files
        kwargsDefault = dict({'tolFloat' : 0.0,
                              'tolMSE' : 0.0,
                              'conditions' : dict({'float':0, 'array':0}),
                              'ignore' : [],
                              'tolValuePair' :{},
                              'tolIter' : 0,
                              'maxIter'  : 0,
                              'eigval' : {}})
        kwargs = {**kwargsDefault, **kwargs} 
        self.succeeded = False
        self.passed = False
        self.name = kwargs['name']
        self.tolFloat = kwargs['tolFloat']
        self.tolMSE = kwargs['tolMSE']
        self.conditions = kwargs['conditions']
        self.path = Path(self.name)
        self.ignore = kwargs['ignore']
        self.tolValuePair = kwargs['tolValuePair']
        self.tolIter = kwargs['tolIter']
        self.maxIter = kwargs['maxIter']
        self.eigval = kwargs['eigval']
        super(Test, self).__init__()

    def countFailures(self):
        """
        Counts the number of occurred failures.
        """
        failureCount = {'total' : len(self),
                         'format':0,
                         'float':0,
                         'array':0,
                         'string':0,
                         'reference':0}

        for failure in self:
            if failure.code == Failure_code.FORMAT:
                failureCount['format'] += 1
            elif failure.code == Failure_code.FLOAT:
                failureCount['float'] += 1
            elif failure.code == Failure_code.ARRAY:
                failureCount['array'] += 1
            elif failure.code == Failure_code.STRING:
                failureCount['string'] += 1

        return failureCount

    def evaluateFailureCount(self, keyList=['float', 'array']):
        failureCount = self.countFailures()
        passedList = []
        if failureCount['total']>0:
            self.succeeded = False
            for key in keyList:
                if failureCount[key]>0:
                    if failureCount[key]>self.conditions[key]:
                        passedList.append(False)
                    else:
                        passedList.append(True)
            self.passed = sum(passedList)
        else:
            self.succeeded = True
            self.passed = True
        return self

    def evaluate(self, ref_data, test_data):
        """
        Evaluates the test for two data sets with the given condition.
        In:
            ref_data       dict or list        test data
            test_data       dict or list        reference data
        """
        if 'info.xml' in self.name:
            compare_info(ref_data, test_data, self, self.tolFloat, self.tolMSE, self.path, self.tolIter, self.maxIter, self.ignore, self.tolValuePair)
        elif 'INFO.OUT' in self.name and 'WANNIER' not in self.name:
            compare_INFO(ref_data, test_data, self, self.tolFloat, self.tolMSE, self.path, self.tolIter, self.maxIter, self.ignore, self.tolValuePair)
        elif 'eigval.xml' in self.name:
            compare_eigval(ref_data, test_data, self, self.tolFloat, self.tolMSE, self.path, self.ignore, self.tolValuePair, self.eigval)
        else:                                                                               
            compare_data(ref_data, test_data, self, self.tolFloat, self.tolMSE, self.path, self.ignore, self.tolValuePair)

        self = self.evaluateFailureCount()
        return self

    def printFailList(self):
        """
        Prints out if the test has succeeded and if not, outputs the specific failures to stdout.
        """
        if len(self)>0:
            failureCount = self.countFailures()
            decription_tuple = (failureCount['total'],
                                failureCount['format'], 
                                failureCount['float'], self.conditions['float'], #is / allowed
                                failureCount['array'], self.conditions['array'], #is / allowed
                                failureCount['string'])
            if self.passed:
                print_color('    %s PASS'%self.name, 'yellow')
                print_color('     Failures TOTAL: %i, FORMAT: %i, FLOAT: %i/%i, ARRAY: %i/%i, STRING %i'%decription_tuple, 'yellow')
                for f in self:
                    f.printFailure(self.passed)
            else:
                print_color('    %s FAIL'%self.name, 'red')
                print_color('      Failures TOTAL: %i, FORMAT: %i, FLOAT: %i/%i, ARRAY: %i/%i, STRING %i'%decription_tuple, 'red')
                for f in self:
                    f.printFailure(self.passed)
        else:
            print_color('    %s SUCCESS'%self.name, 'green')

    def XMLoutput(self, root, tdirectory):
        if self.succeeded:
            tstatus = 'succeeded'
        elif self.passed:
            tstatus = 'passed'
        else:
            tstatus = 'failed'
        tname = self.name
        tdescription = []
        if len(self)>0:
            if self.passed:
                tdescription.append('%s passed.'%self.name)
            else:
                tdescription.append('%s failed.'%self.name)
            failureCount = self.countFailures()
            decription_tuple = (failureCount['total'],
                                failureCount['format'], 
                                failureCount['float'], self.conditions['float'], #is / allowed
                                failureCount['array'], self.conditions['array'], #is / allowed
                                failureCount['string'],
                                failureCount['reference'])
            tdescription.append('Failures TOTAL: %i, FORMAT: %i, FLOAT: %i/%i, ARRAY: %i/%i, STRING %i, REFERENCE: %i'%decription_tuple)
            for f in self:
                tdescription.append(str(f))
        else:
            tdescription.append('Test %s succeeded.'%self.name)
            tdescription.append('No failure occured.')

        addTestXML(root, tstatus, tname, tdescription, tdirectory)

def fromInit(testInit):
    condFloat  = int(testInit['condition'].split()[0])
    condMSE    = int(testInit['condition'].split()[1])
    return Test( name       = testInit['file'],
                 tolFloat   = float(testInit['tolFloat']),
                 tolMSE     = float(testInit['tolMSE']),
                 conditions = dict({'float':condFloat, 'array':condMSE}),
                 ignore     = testInit['ignore'].split(';'),
                 tolValuePair = testInit['tolValuePair'],
                 tolIter    = testInit['tolIter'],
                 maxIter    = testInit['maxIter'],
                 eigval     = testInit['eigval'])
