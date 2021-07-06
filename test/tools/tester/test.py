from typing import Union
import xml.etree.ElementTree as ET

from ..termcolor_wrapper import print_color
from .compare import compare_data, compare_info, compare_eigval, compare_INFO
from .path import Path
from .failure import Failure_code


class Test(list):
    """
    Sub class of list. Stores Fails of a test.
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
        self.floatFails = []
        self.intFails = []
        self.arrayFails = []
        self.stringFails = []
        self.otherFails = []
        super(Test, self).__init__()

    def count_failures(self) -> dict:
        """
        Counts the number of occurred failures.
        """
        failure_count = {'total' : len(self),
                         'format':0,
                         'float':0,
                         'array':0,
                         'string':0,
                         'integer':0,
                         'reference':0}

        failure_code_to_str = {Failure_code.FORMAT: 'format',
                                  Failure_code.FLOAT: 'float',
                                  Failure_code.ARRAY: 'array',
                                  Failure_code.STRING: 'string',
                                  Failure_code.INT: 'integer'
                                    }

        for failure in self:
            key = failure_code_to_str[failure.code]
            failure_count[key] += 1
        return failure_count

     

    def evaluate_failure_count(self, keyList=['float', 'array']):
        """
        Checks if failures occur in a test.        
        :param keyList: list of error keys, defaults to ['float', 'array']
        """
        failureCount = self.count_failures()
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


    def evaluate(self, test_data:Union[dict, list], ref_data:Union[dict, list], test_dir:str):
        """
        Evaluates the test for two data sets with the given condition.
        :param test_data:    output data from test
        :param ref_data:     reference output data
        :param test_dir:     path to test directory
        """


        if 'info.xml' == self.name:
            compare_info(test_data, 
                         ref_data, 
                         test_dir, 
                         self, 
                         self.tolIter, 
                         self.maxIter)
        elif 'INFO.OUT' == self.name:
            compare_INFO(test_data, 
                         ref_data, 
                         test_dir, 
                         self,   
                         self.tolIter, 
                         self.maxIter)
        elif 'eigval.xml' == self.name:
            compare_eigval(test_data, 
                           ref_data, 
                           test_dir, 
                           self, 
                           self.eigval)
        else:
            compare_data(test_data, 
                         ref_data, 
                         test_dir, 
                         self)

        self = self.evaluate_failure_count()
        return self

    def print_test_result(self):
        """
        Print the result of the test to stdout.
        """
        if len(self)>0:
            failureCount = self.count_failures()
            description_tuple = (failureCount['total'],
                                failureCount['format'], 
                                failureCount['float'], self.conditions['float'], # == is, allowed
                                failureCount['integer'],
                                failureCount['array'], self.conditions['array'], # == is, allowed
                                failureCount['string'])
            if self.passed:
                print_color('    %s PASS'%self.name, 'yellow')
                color = 'yellow'
            else:
                print_color('    %s FAIL'%self.name, 'red')
                color = 'red'

            print_color('      Failures TOTAL: %i, FORMAT: %i, FLOAT: %i/%i, INTEGER: %i, ARRAY: %i/%i, STRING %i'%description_tuple, color)

            sort_fail_list(self)

            failure_types = {'other': self.otherFails,
                             'float': self.floatFails, 
                             'integer': self.intFails, 
                             'array': self.arrayFails,
                             'string': self.stringFails}

            heading = {'float': float_heading(self), 
                       'integer': '         Result  Reference  Error  Tolerance  Key',
                       'array': '         MSE     Tolerance   Path',
                       'string': None,
                       'other': None}

            for key in failure_types:
                if key == 'other' or failureCount[key] > 0:
                    caps_key = key.upper()
                    print_color('        ' + caps_key + ' FAILURES', color)
                    if heading[key] != None:
                        print_color(heading[key], color) 
                    for f in failure_types[key]:
                        print_color(f.message, color)

        else:
            print_color('    %s SUCCESS'%self.name, 'green')




def sort_fail_list(self):
    """
    Sorts the failures by type.
    """
    for f in self:
        if f.code == Failure_code.FLOAT:
            self.floatFails.append(f)
        elif f.code == Failure_code.INT:
            self.intFails.append(f)
        elif f.code == Failure_code.ARRAY:
            self.arrayFails.append(f)
        elif f.code == Failure_code.STRING:
            self.stringFails.append(f)
        else:
            self.otherFails.append(f)


def float_heading(self):
    """
    Generates the heading for the list of float failures.
    """
    if 'eigval' in self.name:
        heading = '         (kpt, state)  Result        Reference     Error   Tolerance'
    else:
        heading = '         Line  Result        Reference     Error   Tolerance'
    return heading


def from_init(testInit:dict):
    """
    Reads data from the init dictionary.
    :param  testInit:   dictionary containing important information on the test case
  
    """

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
