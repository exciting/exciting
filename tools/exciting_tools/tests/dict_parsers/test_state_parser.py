from excitingtools.exciting_dict_parsers.state_parser import parse_state_out
from pathlib import Path
import numpy as np
import sys


def write_int(file, i: int, byteorder):
    """ Writes an integer as bytes to the file.

    :param file: file object to write the int
    :param i: integer, which bytes should be written
    :param byteorder: endianness of the byteorder ("little" or "big")
    """
    file.write(int.to_bytes(4, 4, byteorder))
    file.write(i.to_bytes(4, byteorder))
    file.write(int.to_bytes(4, 4, byteorder))


def write_array(file, a: np.ndarray, byteorder):
    """ Writes a numpy array to the file as bytes.

    :param file: file object to write the int
    :param a: numpy array, which bytes should be written
    :param byteorder: endianness of the byteorder ("little" or "big")
    """
    num_bytes = a.size * 8
    file.write(num_bytes.to_bytes(4, byteorder))
    file.write(a.tobytes(order='F'))
    file.write(num_bytes.to_bytes(4, byteorder))


def write_state(path, state: dict):
    """ Creates a file at the specified path and writes the dictionary in the same structure as STATE.OUT

    :param path: path to the file to which should be written to
    :param state: state dictionary, which should be converted to a binary STATE file
    """
    byteorder = sys.byteorder
    with open(path, "wb+") as file:
        # write version information
        file.write(int.to_bytes(52, 4, byteorder))
        for i in state["version"]:
            i: int
            file.write(i.to_bytes(4, byteorder))
        file.write("YomNFS0t4s8uhGxs4FCNGsAOp5DapSm7HAQ58aVd".encode("utf-8"))
        file.write(int.to_bytes(52, 4, byteorder))
        # write other general information
        for key in ("spinpol", "nspecies", "lmmaxvr", "nrmtmax"):
            write_int(file, state[key], byteorder)
        # next write the species specific information
        for i in range(state["nspecies"]):
            write_int(file, state["number of atoms"][i], byteorder)
            write_int(file, len(state["muffintin radial meshes"][i]), byteorder)
            write_array(file, state["muffintin radial meshes"][i], byteorder)
        # next write g vector grid
        file.write(int.to_bytes(12, 4, byteorder))
        for i in range(3):
            file.write(state["g vector grid"][i].to_bytes(4, byteorder))
        file.write(int.to_bytes(12, 4, byteorder))
        # next are again a few integers
        for key in ("ngvec", "ndmag", "nspinor", "ldapu", "lmmaxlu"):
            write_int(file, state[key], byteorder)
        # next write the array pairs
        arrays = (
            ("muffintin density", "interstitial density"),
            ("muffintin coulomb potential", "interstitial coulomb potential"),
            ("muffintin exchange-correlation potential", "interstitial exchange-correlation potential"),
        )
        for muffintin_array, interstitial_array in arrays:
            muffintin_array: np.ndarray = state[muffintin_array]
            interstitial_array: np.ndarray = state[interstitial_array]
            num_bytes = (muffintin_array.size + interstitial_array.size) * 8
            file.write(num_bytes.to_bytes(4, byteorder))
            file.write(muffintin_array.tobytes(order='F'))
            file.write(interstitial_array.tobytes(order='F'))
            file.write(num_bytes.to_bytes(4, byteorder))
        # Lastly write the potential arrays
        num_bytes = 8 * (state["muffintin effective potential"].size + state["interstitial effective potential"].size +
                         2 * state["reciprocal interstitial effective potential"].size)
        file.write(num_bytes.to_bytes(4, byteorder))
        file.write(state["muffintin effective potential"].tobytes(order='F'))
        file.write(state["interstitial effective potential"].tobytes(order='F'))
        file.write(state["reciprocal interstitial effective potential"].tobytes(order='F'))
        file.write(num_bytes.to_bytes(4, byteorder))


def test_parse_state_out(tmp_path: Path):
    directory = tmp_path

    # setting seed to avoid randomness in tests
    np.random.seed(0)

    state_ref = {
        "version": (1, 2, 3),
        "versionhash": "YomNFS0t4s8uhGxs4FCNGsAOp5DapSm7HAQ58aVd",
        "spinpol": False,
        "nspecies": 2,
        "lmmaxvr": 3,
        "nrmtmax": 3,
        "number of atoms": [1, 1],
        "muffintin radial meshes": [np.array([1.0, 2.0]), np.array([0.5, 0.75, 1.5])],
        "g vector grid": (2, 2, 2),
        "ngvec": 5,
        "ndmag": 0,
        "nspinor": 1,
        "ldapu": 0,
        "lmmaxlu": 16,
        "muffintin density": np.random.random((3, 3, 2)),
        "interstitial density": np.random.random(8),
        "muffintin coulomb potential": np.random.random((3, 3, 2)),
        "interstitial coulomb potential": np.random.random(8),
        "muffintin exchange-correlation potential": np.random.random((3, 3, 2)),
        "interstitial exchange-correlation potential": np.random.random(8),
        "muffintin effective potential": np.random.random((3, 3, 2)),
        "interstitial effective potential": np.random.random(8),
        "reciprocal interstitial effective potential": np.random.random(5) + 1j * np.random.random(5),
    }
    write_state(directory / "STATE.OUT", state_ref)
    state_out = parse_state_out(directory / "STATE.OUT")

    assert state_out["version"] == state_ref["version"]
    assert state_out["versionhash"] == state_ref["versionhash"]
    assert state_out["spinpol"] == state_ref["spinpol"]
    assert state_out["nspecies"] == state_ref["nspecies"]
    assert state_out["lmmaxvr"] == state_ref["lmmaxvr"]
    assert state_out["nrmtmax"] == state_ref["nrmtmax"]
    assert state_out["number of atoms"] == state_ref["number of atoms"]
    assert np.allclose(state_out["muffintin radial meshes"][0], state_ref["muffintin radial meshes"][0])
    assert np.allclose(state_out["muffintin radial meshes"][1], state_ref["muffintin radial meshes"][1])
    assert state_out["g vector grid"] == state_ref["g vector grid"]
    assert state_out["ngvec"] == state_ref["ngvec"]
    assert state_out["ndmag"] == state_ref["ndmag"]
    assert state_out["nspinor"] == state_ref["nspinor"]
    assert state_out["ldapu"] == state_ref["ldapu"]
    assert state_out["lmmaxlu"] == state_ref["lmmaxlu"]
    assert np.allclose(state_out["muffintin density"], state_ref["muffintin density"])
    assert np.allclose(state_out["interstitial density"], state_ref["interstitial density"])
    assert np.allclose(state_out["muffintin coulomb potential"], state_ref["muffintin coulomb potential"])
    assert np.allclose(state_out["interstitial coulomb potential"], state_ref["interstitial coulomb potential"])
    assert np.allclose(state_out["muffintin exchange-correlation potential"],
                       state_ref["muffintin exchange-correlation potential"])
    assert np.allclose(state_out["interstitial exchange-correlation potential"],
                       state_ref["interstitial exchange-correlation potential"])
    assert np.allclose(state_out["muffintin effective potential"], state_ref["muffintin effective potential"])
    assert np.allclose(state_out["interstitial effective potential"], state_ref["interstitial effective potential"])
    assert np.allclose(state_out["reciprocal interstitial effective potential"],
                       state_ref["reciprocal interstitial effective potential"])
