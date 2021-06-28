from typing import Union
import json 
import numpy as np
import copy
import sys
from collections.abc import Hashable

def common_iterable(obj: Union[dict, list]):
    """
    Create an iterable dict or tuple, such that one can iterate over a
    dictionary or list with the same syntax:

        for index in common_iterable(d):
            print('index = key or integer index:', index)
            print('d[index] = dictionary value or list element:', d[index])
    """

    if isinstance(obj, dict):
        return obj
    else:
        return (index for index, value in enumerate(obj))
    
    
def __container_converter(data: Union[list,dict]):
    """Mutates the input dictionary, converting string representations of numerical data into numerical data.

    :param dict data: Dictionary with string values.
    :return:dict new_data: Dictionary with all string literal values converted to numerical values.
    """
    np_convert = {np.float64: float,
                  np.int32: int}

    for element in common_iterable(data): 
        if isinstance(data[element], dict) or isinstance(data[element], list):
           data[element] = __container_converter(data[element])
        elif isinstance(data[element], np.ndarray):
            data[element] = data[element].tolist()
        elif isinstance(data[element], (np.float64, np.int32)): 
            value = data[element]
            data[element] = np_convert[type(value)](value)
        elif isinstance(data[element], Hashable): 
            try:
                data[element] = json.loads(data[element])
            except:
                pass
        else: 
            message = 'Type not converted by dict converter: ' + str(type(data[element]))
            sys.exit(message)
    
    return data 

   
def container_converter(data: dict):
    """Mutates the input dictionary, converting string representations of numerical data into numerical data.
    
    :param dict data: Dictionary to be copied containing string values.
    :return: dict new_data: Copied dictionary with all string literal values converted to numerical values.
    """    
    
    new_data = copy.deepcopy(data)
    return __container_converter(new_data)