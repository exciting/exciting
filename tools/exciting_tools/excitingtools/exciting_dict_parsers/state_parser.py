""" Parser for STATE.OUT binary file.
"""
import numpy as np
import struct
import sys


def get_byteorder_format_char(byteorder) -> str:
    """ Converts the byteorder string to formatting symbol for struct.unpack function

    :param byteorder: endianness of the byteorder ("little" or "big")
    :return: the character representing the endianness ('<' or '>')
    """
    chars = {
        "little": '<',
        "big": '>',
    }
    return chars[byteorder]


def read_array(file, shape, byteorder) -> np.ndarray:
    """ Read in a sequence of double values and reshape them into the wanted shape.

    :param file: binary file object, which contains the number sequence
    :param shape: the desired shape of the number sequence (with column-major order)
    :param byteorder: the endianness of the bytes representation
    :return: a numpy array of double values with the desired shape
    """
    # get the correct format character for the byteorder
    byteorder = get_byteorder_format_char(byteorder)
    # get the number of bytes for the array
    num_bytes = np.prod(shape) * 8
    # get all the double numbers
    values = struct.unpack(f"{byteorder}{np.prod(shape)}d", file.read(num_bytes))
    # return the values as correctly shaped numpy array
    return np.reshape(values, shape, order='F')


def read_complex_array(file, shape, byteorder) -> np.ndarray:
    """ Read in a sequence of complex values and reshape them into the wanted shape.

    :param file: binary file object, which contains the number sequence
    :param shape: the desired shape of the number sequence (with column-major order)
    :param byteorder: the endianness of the bytes representation
    :return: a numpy array of complex values with the desired shape
    """
    byteorder = get_byteorder_format_char(byteorder)
    array = np.reshape(
        [complex(*struct.unpack(f"{byteorder}2d", file.read(16))) for _ in range(np.prod(shape))], shape, order='F'
    )
    return array


def read_int(file, byteorder, num_bytes=4, signed=False) -> int:
    """ Read in a single integer.

    :param file: binary file object, which contains the number
    :param byteorder: the endianness of the bytes representation
    :param num_bytes: the number of bytes, which constitute the integer
    :param signed: weather the bytes should be interpreted as a signed or unsigned integer
    :return: the read in integer
    """
    return int.from_bytes(file.read(num_bytes), byteorder, signed=signed)


def read_integers(file, integers, dest, byteorder):
    """ Read in a sequence of named integers and store them in a dictionary.

    :param file: binary file object, which contains the integers
    :param integers: the keys of the integers for the dictionary
    :param dest: the dictionary, which stores the integers and gets modified
    :param byteorder: the endianness of the bytes representation
    """
    for integer in integers:
        num_bytes = read_int(file, byteorder)
        dest[integer] = read_int(file, byteorder, num_bytes)
        assert num_bytes == read_int(file, byteorder)


def parse_state_out(path, byteorder=sys.byteorder) -> dict:
    """ Parser for: STATE.OUT

    STATE.OUT is a binary file. For every 'Write' (in Fortran) there are 4 leading bytes, which are an integer
    equal to the number of bytes of the Fortran objects written. After this the leading 4 bytes are
    repeated.
    This file contains information about the version of exciting, which was used and all the information about
    the density and potentials.
    The functions for the muffin-tin region are stored as an expansion of real spherical harmonics.

    :param path: the path to binary file
    :param byteorder: the endianness of the file (defaults to the endianness of the system)
    :return: a dictionary containing all the state information
    """
    state = {}
    radial_meshes = []
    number_of_atoms = []
    with open(path, "rb") as file:
        # first read the version tuple (3 * 4 bytes) and the versionhash (40 bytes)
        assert read_int(file, byteorder) == 52
        state["version"] = tuple(read_int(file, byteorder) for _ in range(3))
        state["versionhash"] = file.read(40).decode("utf-8")
        assert read_int(file, byteorder) == 52
        # next read whether the calculation was spin polarized or not (4 bytes)
        assert read_int(file, byteorder) == 4
        state["spinpol"] = read_int(file, byteorder) == 1
        assert read_int(file, byteorder) == 4
        # next we read different integers
        integers = ("nspecies", "lmmaxvr", "nrmtmax")
        read_integers(file, integers, state, byteorder)
        # next for every species there are the number of atoms for this species, the number of points
        # of the radial mesh and finally the radial mesh
        for _ in range(state["nspecies"]):
            num_bytes = read_int(file, byteorder)
            number_of_atoms.append(read_int(file, byteorder, num_bytes))
            assert num_bytes == read_int(file, byteorder)
            num_bytes = read_int(file, byteorder)
            nrmt = read_int(file, byteorder, num_bytes)
            assert num_bytes == read_int(file, byteorder)
            num_bytes = read_int(file, byteorder)
            assert num_bytes // 8 == nrmt
            radial_meshes.append(read_array(file, nrmt, byteorder))
            assert num_bytes == read_int(file, byteorder)
        # save the number of atoms per species
        state["number of atoms"] = number_of_atoms
        # save the radial meshes
        state["muffintin radial meshes"] = radial_meshes

        # next read the g vector grid size (3 * 4 bytes)
        num_bytes = read_int(file, byteorder)
        state["g vector grid"] = tuple(read_int(file, byteorder) for _ in range(3))
        assert num_bytes == read_int(file, byteorder)

        # next again read a few integers
        integers = ("ngvec", "ndmag", "nspinor", "ldapu", "lmmaxlu")
        read_integers(file, integers, state, byteorder)

        # next read a few pairs of arrays
        arrays = (
            ("muffintin density", "interstitial density"),
            ("muffintin coulomb potential", "interstitial coulomb potential"),
            ("muffintin exchange-correlation potential", "interstitial exchange-correlation potential"),
        )
        muffintin_shape = (state["lmmaxvr"], state["nrmtmax"], sum(number_of_atoms))
        interstitial_shape = np.prod(state["g vector grid"])
        for muffintin_array, interstitial_array in arrays:
            num_bytes = read_int(file, byteorder)
            assert num_bytes // 8 == (np.prod(muffintin_shape) + interstitial_shape)
            state[muffintin_array] = read_array(file, muffintin_shape, byteorder)
            state[interstitial_array] = read_array(file, interstitial_shape, byteorder)
            assert num_bytes == read_int(file, byteorder)
        # the next entry is an array triple
        num_bytes = read_int(file, byteorder)
        assert num_bytes // 8 == (np.prod(muffintin_shape) + interstitial_shape + 2 * state["ngvec"])
        state["muffintin effective potential"] = read_array(file, muffintin_shape, byteorder)
        state["interstitial effective potential"] = read_array(file, interstitial_shape, byteorder)
        state["reciprocal interstitial effective potential"] = read_complex_array(file, state["ngvec"], byteorder)
        assert num_bytes == read_int(file, byteorder)

        # the next two arrays are only present if the calculation was spin polarized
        if state["spinpol"]:
            num_bytes = read_int(file, byteorder)
            state["muffintin magnetization"] = read_array(file, muffintin_shape, byteorder)
            state["interstitial magnetization"] = read_array(file, interstitial_shape, byteorder)
            assert num_bytes == read_int(file, byteorder)

        if state["ldapu"] != 0:
            shape = (state["lmmaxlu"], state["lmmaxlu"], state["nspinor"], state["nspinor"], sum(number_of_atoms))
            num_bytes = read_int(file, byteorder)
            assert num_bytes // 8 == 2 * np.prod(shape)
            state["vmatlu"] = read_complex_array(file, shape, byteorder)
            assert num_bytes == read_int(file, byteorder)

    return state
