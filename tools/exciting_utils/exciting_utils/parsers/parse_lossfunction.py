#!/usr/bin/python

def parse_lossfunction( fname ):
    """
    Function that parses files containing loss function
    e.g. LOSS_FXCRPA_OC11_QMT001.OUT
    :param str fname: name of the file
    """
    xdata = []
    ydata = []
    file = open(fname,'r')
    for lines in file:
        if 'Frequency' in lines:
            break
    for lines in file:
        data = lines.split()
        xdata.append(float(data[0]))
        ydata.append(float(data[1]))
    file.close()
    return xdata, ydata