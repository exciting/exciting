import Scientific

try:
    if Scientific.__version__ < 2.8:
        print 'your ScientifPython version is: ',Scientific.__version__ 
        print 'ScientificPython 2.8 or greater required for numpy support in NetCDF'
        raise Exception
except AttributeError:
    print 'It appears your ScientificPython has no version.'
    print 'That probably means it is not 2.8, which is required'
    raise Exception

