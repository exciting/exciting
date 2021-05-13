"""
Functions that check the type of an object
"""
import numpy as np

from .failure import *
from .path import *

# TODO(A/B/H) Issue 67 This function should use isinstance
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False


def checkFloatFail(test, test_data, ref_data, tolF, path, tolValuePair):
    if path.lastElement() in tolValuePair['value']:
        tol = float(tolValuePair['tol'])
    else:
        tol = tolF
    err = abs(float(test_data)-float(ref_data))
    if err > tol:
        test.append(Failure(Failure_code.FLOAT, error=err, tolerance=tol, path=path))

def checkArrayFail(test, test_data, ref_data, tolL, path, tolValuePair):
    if path.lastElement() in tolValuePair['value']:
        tol = tolValuePair['tol']
    else:
        tol = tolL
    mse = np.sqrt(np.mean((np.array(test_data)-np.array(ref_data))**2))
    if mse>tol:
        test.append(Failure(Failure_code.ARRAY, error = mse, tolerance = tol, path = path))    

def checkStringFail(test, test_data, ref_data, path):
    if test_data != ref_data:
        test.append(Failure(Failure_code.STRING, path = path))

def compare_data(test_data, ref_data, test, tolFloat, tolList, path, ignore, tolValuePair):
    #Check if test_data and ref_data are both lists.
    if isinstance(test_data, list) and isinstance(ref_data, list):
        #Check if the lengths are not the same .                       
        if len(test_data) != len(ref_data):
            test.append(Failure(Failure_code.FORMAT ,path=path))
            return
        else:
            #Check if test_data and ref_data can both be converted in float lists.
            if (sum([isinstance(x, float) for x in test_data]) == len(test_data)) and (sum([isinstance(x, float) for x in ref_data]) == len(ref_data)):
                checkArrayFail(test, test_data, ref_data, tolList, path, tolValuePair)
            
            else:           
                for i in range(0,len(test_data)):
                    path = path.update(str(i))
                    #Check of the elements of test_data and ref_data are dicts.
                    if isinstance(test_data[i], dict) and isinstance(test_data[i], dict):
                        compare_data(test_data[i], ref_data[i], test, tolFloat, tolList, path, ignore, tolValuePair)
                    
                    #Check of the elements of test_data and ref_data are lists.
                    elif isinstance(test_data[i], list) and isinstance(test_data[i], list):
                        compare_data(test_data[i], ref_data[i], test, tolFloat, tolList, path, ignore, tolValuePair)
                    
                    #Check of the elements of test_data and ref_data are floats (If maybe not all elements are floats).
                    elif isfloat(test_data[i]) and isfloat(ref_data[i]):
                        checkFloatFail(test, test_data[i], ref_data[i], tolFloat, path, tolValuePair)

                    elif isinstance(test_data[i], str) and isinstance(ref_data[i], str):
                        checkStringFail(test, test_data[i], ref_data[i], path)

                    else:
                        test.append(Failure(Failure_code.FORMAT ,path=path))
                        return
                    path = path.removeLastElement()
    
    #Check if test_data and ref_data are both dicts.
    elif isinstance(test_data, dict) and isinstance(ref_data, dict):
        #Check if the keys of test_data and ref_data are not identical.
        if test_data.keys() != ref_data.keys():
            test.append(Failure(Failure_code.FORMAT ,path=path))
        for key in test_data.keys():
            if key not in ref_data.keys():
                ignore.append(key)
        for key	in ref_data.keys():
            if key not in test_data.keys():
                ignore.append(key)
        for key in ref_data.keys():
            if key in ignore: break
            path = path.update(key)
            #Check if the corresponding elements are dicts.
            if isinstance(test_data[key], dict) and isinstance(ref_data[key], dict):
                compare_data(test_data[key], ref_data[key], test, tolFloat, tolList, path, ignore, tolValuePair)
                
            #Check if the corresponding elements are lists.
            elif isinstance(test_data[key], list) and isinstance(ref_data[key], list):
                compare_data(test_data[key], ref_data[key], test, tolFloat, tolList, path, ignore, tolValuePair)
                
            #Check if the corresponding elements are floats.
            elif isfloat(test_data[key]) and isfloat(ref_data[key]):
                checkFloatFail(test, test_data[key], ref_data[key], tolFloat, path, tolValuePair)
                
            #Check if the corresponding elements are strings.
            elif isinstance(test_data[key], str) and isinstance(ref_data[key], str):
                checkStringFail(test, test_data[key], ref_data[key], path)
                
            else:
                test.append(Failure(Failure_code.FORMAT ,path=path))
                return
            path = path.removeLastElement()
                
    elif isfloat(test_data) and isfloat(ref_data):
          checkFloatFail(test, test_data, ref_data, tolFloat, path, tolValuePair)
          
    else:
        test.append(Failure(Failure_code.FORMAT ,path=path))

def compare_INFO(test_data, ref_data, test, tolFloat, tolList, path, tolIter, maxIter, ignore, tolValuePair):
    
    # finds the last iteration
    list1 = []
    list2 = []
    for key in test_data['scl'].keys():
        list1.append(int(key))
    for key in ref_data['scl'].keys():
        list2.append(int(key))

    # checks wether the number of itereations is in the tolerance
    if maxIter != 0:
      if max(list1) > int(maxIter):
        test.append(Failure(Failure_code.FLOAT, error = max(list1)-int(maxIter), tolerance = 0, path = "INFO.OUT/Iterations"))
    else:
      checkFloatFail(test, max(list1), max(list2), int(tolIter), Path("INFO.OUT/Iterations"), tolValuePair)

    # TODO(A/B/H)issue 67. Add integer failure

    # compares the last iteration
    compare_data(test_data['scl'][str(max(list1))], ref_data['scl'][str(max(list2))], test, tolFloat, tolList, path, ignore, tolValuePair)

    # compares the initialization
    compare_data(test_data['initialization'], ref_data['initialization'], test, tolFloat, tolList, path, ignore, tolValuePair)
    
def compare_info(test_data, ref_data, test, tolFloat, tolList, path, tolIter, maxIter, ignore, tolValuePair):

    # finds the last iteration
    list1 = []
    list2 = []
    for key in test_data['scl'].keys():
        list1.append(int(key))
    for key in ref_data['scl'].keys():
        list2.append(int(key))
        
    # checks whether the number of iterations is in the tolerance
    if maxIter != 0:
      if max(list1) > int(maxIter):
        test.append(Failure(Failure_code.FLOAT, error = max(list1)-int(maxIter), tolerance = 0, path = "info.xml/Iterations"))
    else:
      checkFloatFail(test, max(list1), max(list2), int(tolIter), Path("info.xml/Iterations"), tolValuePair)

    # TODO(A/B/H)issue 67. Add integer failure
    
    # compares the energies calculated in the last iteration
    compare_data(test_data['scl'][str(max(list1))]['energies'], ref_data['scl'][str(max(list2))]['energies'], test, tolFloat, tolList, path, ignore, tolValuePair)
    
    # if moments are calculated, they are compared too
    if 'moments' in ref_data['scl'][str(max(list2))]:
        compare_data(test_data['scl'][str(max(list1))]['moments'], ref_data['scl'][str(max(list2))]['moments'], test, tolFloat, tolList, path, ignore, tolValuePair)


def compare_eigval(test_data, ref_data, test, tolFloat, tolList, path, ignore, tolValuePair, eigval):

    compare_data(test_data, ref_data, test, tolFloat, tolList, path, ignore, tolValuePair)
    if 'eigvalPair' in eigval:
        for key in eigval['eigvalPair'].split():
            key1 = key.split(',')
            ref = float(ref_data['kpt']['1']['state'][key1[1]]['eigenvalue']) - float(ref_data['kpt']['1']['state'][key1[0]]['eigenvalue'])
            dat = float(test_data['kpt']['1']['state'][key1[1]]['eigenvalue']) - float(test_data['kpt']['1']['state'][key1[0]]['eigenvalue'])
            checkFloatFail(test, dat, ref, float(eigval['tolEigval']), path, tolValuePair)
