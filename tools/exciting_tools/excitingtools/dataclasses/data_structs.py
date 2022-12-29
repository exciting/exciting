""" Data Structures.

Data structure is defined as a container for data.
Many of these classes could be @dataclass, however excitingtools
retains support for python 3.6.
"""

class PointIndex:
    """ Container for (point, index) pair
    """
    def __init__(self, point, index: int):
        self.point = point
        self.index = index

class BandIndices:
    """Indices of valence band maximum and conduction band minimum"""
    def __init__(self, VBM: int, CBm: int):
        self.VBM = VBM
        self.CBm = CBm

class NumberOfStates:
    """Number of states. Useful when indexing does not start at 0/1
    """
    def __init__(self, first_state: int, last_state: int):
        self.first_state = first_state
        self.last_state = last_state
        self.n_states = last_state - first_state + 1
