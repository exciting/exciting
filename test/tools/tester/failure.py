import os
import re
from enum import Enum
from typing import List
import numpy as np 
import warnings
import numpy as np




class Failure_code(Enum):
    FORMAT = 1
    FLOAT = 2
    ARRAY = 3
    STRING = 4
    RUN = 5
    REFERENCE = 6
    FILENOTEXIST = 7
    ERRORFILE = 8
    INT = 9


set_failure_code = {str: Failure_code.STRING,
                    int: Failure_code.INT, 
                    float: Failure_code.FLOAT, 
                    list: Failure_code.ARRAY,
                    np.ndarray: Failure_code.ARRAY
                    }


class Failure:
    """
    Class for reporting a failure that occurs when comparing test data to reference data.
    """
    def __init__(self, 
                 error = None,
                 tolerance = None,
                 err_msg = None,
                 path = None,
                 test_data = None,
                 ref_data = None,
                 test_dir = None,
                 test_name = None,
                 failure_code = None):

        if failure_code == None:
            failure_code = set_failure_code[type(ref_data)]

        self.path = path
        self.test_dir = test_dir
        self.test_name = test_name
        self.ref_data = ref_data
        self.test_data = test_data
        self.error = error
        self.tolerance = tolerance
        self.code = failure_code
        self.err_msg = err_msg
        self.message = self.set_failure_message()


    def set_failure_message(self):
        if self.code == Failure_code.FORMAT:
            message = "         FORMAT FAILURE: File from calculation has not the same format as the reference file."
        elif self.code == Failure_code.FLOAT:
            message = self.generate_float_message()
        elif self.code == Failure_code.INT:
            message = self.generate_int_message()
        elif self.code == Failure_code.ARRAY:
            message =  f"         {self.error :3.0e}   {self.tolerance}       {self.path}"
        elif self.code == Failure_code.STRING:
            message = f"         STRING FAILURE: Strings in line {self.find_line() :d} are not the same."
        elif self.code == Failure_code.RUN:
            message = f"         RUN FAILURE: Exciting run failed. All tests in this directory are aborted. Error message: {self.err_msg}"
        elif self.code == Failure_code.REFERENCE:
            message = f"         REFERENCE FAILURE: Reference for {self.err_msg} does not exist."
        elif self.code == Failure_code.FILENOTEXIST:
            message = f"         FILENOTEXIST FAILURE: File {self.err_msg} does not exist."
        elif self.code == Failure_code.ERRORFILE:
            message = f"         ERRORFILE: File {self.err_msg} is errornous."
        return message
       


    def find_line(self)-> List[int]:
        """
        Returns line of the reference file in which the failure occurs. 
        """
        file_path = os.path.join(self.test_dir, "ref", self.test_name)
        line_nr = []
        try:
            file = open(file_path, "r").readlines()
            for i, line in enumerate(file):
                if re.search(str(self.ref_data), line):
                    line_nr.append(i+1)
            if len(line_nr) > 1:
                warnings.warn('Value appears in more than one line.')
            elif len(line_nr)==0:
                line_nr.append(0)
        except:
            line_nr.append(0)
        return line_nr[0]

    
    def generate_float_message(self):
        """
        Generates the failure message for float failures.
        """
        if "eigval.xml" in str(self.path):
                                  # kpt state)                         Result                               Reference                           Error                             Tolerance
            message = f"         ( {str(self.path).split('/')[2] :2s},{str(self.path).split('/')[4]:>5s})   {float(self.test_data) : 011.8f}   {float(self.ref_data) : 011.8f}   {self.error :3.0e}   {self.tolerance :3.0e}"
        else:
                                # Line                       Result                             Reference                         Error                Tolerance                
            message = f"         {self.find_line() :d}   {float(self.test_data) : 011.8f}   {float(self.ref_data) : 011.8f}   {self.error :3.0e}   {self.tolerance :3.0e}"
            if ".xml" in str(self.path):
                message = message+'   '+str(self.path)

        return message


    def generate_int_message(self):
        """
        Generates the failure message for integer failures.
        """  
                           # Result                      Reference                     Error              Tolerance                  Key
        message = f"         {int(self.test_data) }      {int(self.ref_data) }         {self.error }      {self.tolerance }          {str(self.path) :s}"

        return message