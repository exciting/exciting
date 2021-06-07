'''
Unit tests for the test suite. Run the tests:
$python3 runselftest.py
'''
import unittest
import sys

from tools.tester.compare import *
from tools.tester.report import *
from tools.tester.path import *

tolValuePair = {}
tolValuePair['tol'] = '0'
tolValuePair['value'] = ''
      
class TestCompare(unittest.TestCase):
    def test_compare_data_NoFails(self):
        data1 = {'a':1.0, 'b':[1.0,2.0], 'c':{'a1':1.0}, 'd':"bla", 'e':[[1,2],3], 'f':['1','2']}
        data2 = {'a':1.0, 'b':[1.0,2.0], 'c':{'a1':1.0}, 'd':"bla", 'e':[[1,2],3], 'f':['1','2']}
        test = Test( name = 'NoFails', tolFloat = 1.0 , tolMSE = 1.0, tolValuePair= tolValuePair)
        test.evaluate(data1, data2)  
        self.assertEqual(test.countFailures(), dict({'total':0, 'format':0, 'float':0, 'array':0, 'string':0, 'reference':0}))

    def test_compare_data_FormatFails(self):
        data1 = {'a':1.0, 'b':[1.0,2.0], 'c':{'a1':1.0}, 'd':"bla", 'e':[[1,2],3], 'f':['1','2']}
        data2 = {'b':1.0, 'b':[1.0,2.0], 'c':{'a1':1.0}, 'd':"bla", 'e':[[1,2],3], 'f':['1','2']}
        test = Test( name = 'FormatFail1', tolFloat = 1.0 , tolMSE = 1.0, tolValuePair = tolValuePair)
        test.evaluate(data1, data2)
        self.assertEqual(test.countFailures(), dict({'total':1, 'format':1, 'float':0, 'array':0, 'string':0, 'reference':0}))

        data1 = {'a':1.0, 'b':[1.0,2.0,4.0], 'c':{'a1':1.0}, 'd':"bla", 'e':[[1,2],3], 'f':['1','2']}
        data2 = {'a':1.0, 'b':[1.0,2.0], 'c':{'a1':1.0}, 'd':"bla", 'e':[[1,2],3], 'f':['1','2']}
        test = Test( name = 'FormatFail2', tolFloat = 1.0 , tolMSE = 1.0, tolValuePair = tolValuePair)
        test.evaluate(data1, data2)
        self.assertEqual(test.countFailures(), dict({'total':1, 'format':1, 'float':0, 'array':0, 'string':0, 'reference':0}))

        data1 = {'a':1.0, 'b':[1.0,2.0], 'c':{'a1':1.0}, 'd':"bla", 'e':[[1,2], "a" ], 'f':['1','2']}
        data2 = {'a':1.0, 'b':[1.0,2.0], 'c':{'a1':1.0}, 'd':"bla", 'e':[[1,2], 3 ], 'f':['1','2']}
        test = Test( name = 'FormatFail3', tolFloat = 1.0 , tolMSE = 1.0, tolValuePair = tolValuePair)
        test.evaluate(data1, data2)
        self.assertEqual(test.countFailures(), dict({'total':1, 'format':1, 'float':0, 'array':0, 'string':0, 'reference':0}))
        
        data1 = {'a':'a', 'b':[1.0,2.0], 'c':{'a1':1.0}, 'd':"bla", 'e':[[1,2],3], 'f':['1','2']}
        data2 = {'a':1.0, 'b':[1.0,2,0], 'c':{'a1':1.0}, 'd':"bla", 'e':[[1,2],3], 'f':['1','2']}
        test = Test( name = 'FormatFail4', tolFloat = 1.0 , tolMSE = 1.0, tolValuePair = tolValuePair)
        test.evaluate(data1, data2)
        self.assertEqual(test.countFailures(), dict({'total':1, 'format':1, 'float':0, 'array':0, 'string':0, 'reference':0}))

    def test_compare_data_ListFails(self):
        data1 = {'a':5.0, 'b':[1.0,10.0], 'c':{'a1':1.0}, 'd':"blu", 'e':[[1,2],3], 'f':['1','2']}
        data2 = {'a':5.0, 'b':[1.0,2.0],  'c':{'a1':1.0}, 'd':"blu", 'e':[[1,2],3], 'f':['1','2']}
        test = Test( name = 'ListFail', tolFloat = 1.0 , tolMSE = 1.0, tolValuePair = tolValuePair)
        test.evaluate(data1, data2)
        self.assertEqual(test.countFailures(), dict({'total':1, 'format':0, 'float':0, 'array':1, 'string':0, 'reference':0}))

    def test_compare_data_FloatFails(self):
        data1 = {'a':5.0,   'b':[1.0,10.0], 'c':{'a1':1.0}, 'd':"blu", 'e':[[1,2],3], 'f':['1','2']}
        data2 = {'a':'1.0', 'b':[1.0,10.0], 'c':{'a1':1.0}, 'd':"blu", 'e':[[1,2],3], 'f':['1','2']}
        test = Test( name = 'FloatFail', tolFloat = 1.0 , tolMSE = 1.0, tolValuePair = tolValuePair)
        test.evaluate(data1, data2)
        self.assertEqual(test.countFailures(), dict({'total':1, 'format':0, 'float':1, 'array':0, 'string':0, 'reference':0}))

    def test_compare_data_StringFails(self):
        data1 = {'a':5.0, 'b':[1.0,10.0], 'c':{'a1':1.0}, 'd':"blu", 'e':[[1,2],3], 'f':['1','2']}
        data2 = {'a':5.0, 'b':[1.0,10.0], 'c':{'a1':1.0}, 'd':"bla", 'e':[[1,2],3], 'f':['1','2']}
        test = Test( name = 'StringFail', tolFloat = 1.0 , tolMSE = 1.0, tolValuePair = tolValuePair)
        test.evaluate(data1, data2)
        self.assertEqual(test.countFailures(), dict({'total':1, 'format':0, 'float':0, 'array':0, 'string':1, 'reference':0}))


class TestPath(unittest.TestCase):

    def test_path_add(self):
        path1 = Path('some')
        path2 = Path('path')
        path_ref = Path('some/path')
        self.assertEqual(path1+path2, path_ref)

    def test_path_split(self):
        path = 'some/path/to/file'
        path_list = [Path('some'),Path('path'),Path('to'),Path('file')]
        self.assertEqual(Path(path).split(), path_list)


        

if __name__ == '__main__':
    unittest.main()
