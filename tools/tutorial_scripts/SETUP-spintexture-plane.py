#!/usr/bin/python2

"""
Create the spin texture input entry for a square.
"""

import os
import numpy as np
import argparse as ap

from exciting_tools.structure.lattice import reciprocal_lattice_vectors
from exciting_tools.parsers.input import parse_lattice_vectors

def option_parser():
    """
    Parse command line inputs 

    Parse: 
        directory
        center of the plot
        direction 1
        direction 2
        length
        grid
        bands

    :return input_options: Dictionary of parsed command line arguments 
    """
    p = ap.ArgumentParser(description=\
        'Usage: SETUP-spintexture-plane.py -d <directory> -c <center> -spv1 <spanvector1> -spv2 <spanvector2> -sl <length> -g <grid> -b <bandrange>')

    help_directory = 'Defines the relative path to the directory \
                      where input file for the spin texture calculation is stored.'

    help_center = 'Defines the center of the sqare. \
                   Needs to be three floats.'

    help_spanvector1 = 'Defines the first span direction. \
                        Needs to be orthogonal to spanvector2. \
                        Needs to be three floats.'

    help_spanvector2 = 'Defines the second span direction. \
                        Needs to be orthogonal to spanvector1. \
                        Needs to be three floats.'

    help_length = 'Defines the side length of the square.'

    help_grid = 'Defines number of grid points in each direction. \
                 Needs to be two ints.'
    
    help_bands = 'Defines the band range for that the spin texture is calculated. \
                  Needs to be two ints.'
    
    p.add_argument('-d',
                   metavar = 'directory',
                   help = help_directory,
                   type = str,
                   default = '.')
    
    p.add_argument('-c',
                   metavar = 'center',
                   help = help_center,
                   type = float,
                   nargs='+',
                   default = [0.0, 0.0, 0.0])
    
    p.add_argument('-spv1',
                   metavar = 'spanvector1',
                   help = help_spanvector1,
                   type = float,
                   nargs='+',
                   default = [1.0, 0.0, 0.0])
    
    p.add_argument('-spv2',
                   metavar = 'spanvector2',
                   help = help_spanvector2,
                   type = float,
                   nargs='+',
                   default = [0.0, 1.0, 0.0])
    
    p.add_argument('-sl',
                   metavar = 'length',
                   help = help_length,
                   type = float,
                   default = 0.1)

    p.add_argument('-g',
                   metavar = 'grid',
                   help = help_length,
                   type = str,
                   nargs='+',
                   default = [10, 10])
    
    p.add_argument('-b',
                   metavar = 'grid',
                   help = help_length,
                   type = str,
                   nargs='+',
                   default = [0, 1])
   
    args = p.parse_args()
    input_options = {}

    if os.path.isfile(os.path.join(args.d,'input.xml')):
        input_options['directory'] = args.d
    else:
        raise p.error('could not find input.xml in directory')

    try:
        center = np.array(args.c, dtype=float)
        len(center) == 3
        input_options['center'] = center
    except Exception() as e:
        raise p.error("center is errornous: " + str(e))

    try:
        spv1 = np.array(args.spv1, dtype=float)
        len(spv1) == 3
    except Exception() as e:
        raise p.error("spanvector1 is errornous: " + str(e))
    
    try:
        spv2 = np.array(args.spv2, dtype=float)
        len(spv2) == 3
    except Exception() as e:
        raise p.error("spanvector2 is errornous: " + str(e))

    if np.dot(spv1,spv2)==0:
        input_options['spanvector1'] = spv1
        input_options['spanvector2'] = spv2
    else:
        raise p.error('spanvector1 and spanvector2 must be orthogonal')
    
    input_options['length'] = args.sl

    try:
        grid = np.array(args.g, dtype=int)
        len(grid) == 2
        input_options['grid'] = grid
    except Exception() as e:
        raise p.error('grid is errornous: ' + str(e))

    try:
        bands = np.array(args.b, dtype=int)
        len(bands) == 2
        input_options['bands'] = bands
    except Exception() as e:
        raise p.error('bands is errornous: ' + str(e))
    
    return input_options

def main(directory, center, span_vector_1, span_vector_2, length, grid, bands):
    """
    Takes a square defined in cartesian coordinates and generates exciting input for a spin texture calculation on this square.
    :param directory:         relative path to exciting calculation
    :param center:            center of the spin texture plot
    :param span_vector_1:     spans the square
    :param span_vector_2:     spans the square
    :param length:            lengths of the square
    :param grid:              k-grid in the suqare
    :param bands:             bands considered for the spin texture calculation
    """
    lat_vec = parse_lattice_vectors(os.path.join(directory, "input.xml"))
    rec_lat = reciprocal_lattice_vectors(lat_vec)
    rec_lat_inv = np.linalg.inv(rec_lat)

    span_vector_1 = length * span_vector_1 / np.linalg.norm(span_vector_1)
    span_vector_2 = length * span_vector_2 / np.linalg.norm(span_vector_2)
    origin = center - 0.5 * (span_vector_1 + span_vector_2)
    point1 = origin + grid[0] / (grid[0]-1) * span_vector_1
    point2 = origin + grid[1] / (grid[1]-1) * span_vector_2

    point1 = rec_lat_inv.dot(point1)
    point2 = rec_lat_inv.dot(point2)
    origin = rec_lat_inv.dot(origin)
    center = rec_lat_inv.dot(center)

    print('Copy and paste the following lines into the properties element in the input.xml:')
    print('<spintext bands="%i %i">'%(bands[0], bands[1]))
    print('   <plot2d>')
    print('      <parallelogram grid = "%i %i">'%(grid[0], grid[1]))
    print('         <origin coord = "%f  %f  %f"/>'%(origin[0], origin[1], origin[2]))
    print('         <point  coord = "%f  %f  %f"/>'%(point1[0], point1[1], point1[2]))
    print('         <point  coord = "%f  %f  %f"/>'%(point2[0], point2[1], point2[2]))
    print('      </parallelogram>')
    print('   </plot2d>')
    print('</spintext>')

if __name__ == "__main__":
    input_options = option_parser()
    main(**input_options)
