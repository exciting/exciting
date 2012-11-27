"""
Function-like objects that creates cubic clusters.
"""

import numpy as np

from ase.cluster.cubic import FaceCenteredCubic
from ase.cluster.compounds import L1_2

def Octahedron(symbol, length, cutoff=0, latticeconstant=None, alloy=False):
    """
    Returns Face Centered Cubic clusters of the octahedral class depending
    on the choice of cutoff.

    Type                            Condition
    ----                            ---------
    Regular octahedron              cutoff = 0
    Truncated octahedron            cutoff > 0
    Regular truncated octahedron    length = 3 * cutoff + 1
    Cuboctahedron                   length = 2 * cutoff + 1

    Parameters
    ----------
    symbol: The chemical symbol or atomic number of the element(s).

    length: Number of atoms on the square edges of the complete octahedron.

    cutoff (optional): Number of layers cut at each vertex.

    latticeconstant (optional): The lattice constant. If not given,
    then it is extracted form ase.data.

    alloy (optional): If true the L1_2 structure is used. Default is False.
    """

    # Check length and cutoff
    if length < 2:
        raise ValueError("The lenght must be greater than one.")

    if cutoff < 0 or length < 2 * cutoff + 1:
        raise ValueError("The cutoff must fullfill: > 0 and <= (length - 1) / 2.")

    # Create cluster
    surfaces = [(1,1,1), (1,0,0)]
    if length % 2 == 0:
        center = np.array([0.5, 0.5, 0.5])
        layers = [length/2, length - 1 - cutoff]
    else:
        center = np.array([0.0, 0.0, 0.0])
        layers = [(length - 1)/2, length - 1 - cutoff]

    if not alloy:
        return FaceCenteredCubic(symbol, surfaces, layers, latticeconstant, center)
    else:
        return L1_2(symbol, surfaces, layers, latticeconstant, center)

