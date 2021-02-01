from path import Path
from enum import Enum

from termcolor_wrapper import print_color

class Failure_code(Enum):
    FORMAT = 1
    FLOAT = 2
    ARRAY = 3
    STRING = 4
    RUN = 5
    REFERENCE = 6
    FILENOTEXIST = 7
    ERRORFILE = 8

class Failure():
    '''
    Class for reporting a Fail that occure when comparing test data to reference dat.
    In:
        code        int         code for differentiating different fails. 0 to 4 is occupied
        message     string      message that discribes why the fail occured
        path        Path      tells where the fail occured
    '''
    def __init__(self, failure_code, **kwargs):
        self.code = failure_code
        kwargsDefault = dict({'error':None, 'tolerance':None, 'err_msg':None, 'path':None})
        kwargs = {**kwargsDefault, **kwargs}
        self.path = kwargs["path"]

        if failure_code == Failure_code.FORMAT:
            self.message = "FORMAT FAILURE: File from calculation has not the same format as the reference file."
        elif failure_code == Failure_code.FLOAT:
            self.message = "FLOAT FAILURE: Error (%.3e) is bigger than tolerance (%.3e)."%(kwargs["error"], kwargs["tolerance"])
        elif failure_code == Failure_code.ARRAY:
            self.message =  "ARRAY FAILURE: Error (%.3e) is bigger than tolerance (%.3e)."%(kwargs["error"], kwargs["tolerance"])
        elif failure_code == Failure_code.STRING:
            self.message = "STRING FAILURE: Strings are not the same."
        elif failure_code == Failure_code.RUN:
            self.message = "RUN FAILURE: Exciting run failed. All tests in this directory are aborted. Error message: %s"%kwargs["err_msg"]
        elif failure_code == Failure_code.REFERENCE:
            self.message = "REFERENCE FAILURE: Reference for %s does not exist."%kwargs["err_msg"]
        elif failure_code == Failure_code.FILENOTEXIST:
            self.message = "FILENOTEXIST FAILURE: File %s does not exist."%kwargs["err_msg"]
        elif failure_code == Failure_code.ERRORFILE:
            self.message = "ERRORFILE: File %s is errornous."%kwargs["err_msg"]
    
    def __str__(self):
        if self.path:
            return self.message+'\n'+str(self.path)
        else:
            return self.message

    def printFailure(self, passed):
        text_color = 'yellow' if passed else 'red' 
        print_color('        %s'%self.message, text_color)
        if self.path:
            print_color('        %s'%self.path, text_color)
