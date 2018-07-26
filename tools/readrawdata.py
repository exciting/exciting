#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This function reads generic numeric data from a file.
# The data needs to be structured in several blocks separated by one or multiple
# empty lines. Each block of data consists of columns of data separated by one of
# of the characters " " (space), "," (comma), ";" (semicolon) or "|" (vertical bar).
# Lines starting with "#" are treated as comments and are not read.
#
# If nparray=True (default) the function will return a three dimensional numpy array
# of the dimensions N3 x N2 x N1 which are stored in the output array dim=(N1,N2,N3).
# N3 gives the number of data blocks.
# N2 gives the (maximum) number of columns per block.
# N1 gives the (maximum) number of rows per block.
# If the blocks have different dimensions, missing entries are filled with zeros.
# data[0,1,2] references the third element in the second column of the first block.
#
# If nparray=False the function will return a list of lists of lists.
# data[0][2][1] references the third element in the second column of the first block.
# dim[n-1][0] contains the number of rows in the n-th block.
# dim[n-1][1] contains the (maximum) number of columns in the n-th block.
#_______________________________________________________________________________
import numpy as np

def readrawdata( filename, nparray=True):
    delim = [',', ';', '|'] # no space
    newset = 1
    rawdata = []
    rawdim = []
    with open( filename) as f:
        for l in f:
            for d in delim:
                l = l.replace( d, ' ')
            l = l.strip()
            if( l and (l[0] != '#')):
                if( newset):
                    newset = 0
                    rawdata.append( [])
                    rawdim.append( [0, 0])
                values = [float(v) for v in l.split()]
                rawdim[-1][0] += 1
                rawdim[-1][1] = max( rawdim[-1][1], len( values))
                rawdata[-1].append( values)
            elif( not l):
                if( not newset):
                    newset = 1
    
    if( nparray):
        dim = np.zeros( 3, dtype=int)
        for z in range( len( rawdim)):
            dim[2] += 1
            dim[1] = max( dim[1], rawdim[z][1])
            dim[0] = max( dim[0], rawdim[z][0])
        data = np.empty( (dim[2], dim[1], dim[0]))
        for z in range( dim[2]):
            for x in range( dim[0]):
                if( x < rawdim[z][0]):
                    for y in range( dim[1]):
                        if( y < len( rawdata[z][x])):
                            data[ z, y, x] = rawdata[z][x][y]
                        else:
                            data[ z, y, x] = np.nan
                else:
                    data[ z, :, x] = np.nan

        return data, dim
    else:
        return rawdata, rawdim
