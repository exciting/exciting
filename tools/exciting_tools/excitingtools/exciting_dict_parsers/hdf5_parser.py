"""Parsers for HDF5 files."""

import numpy as np


def parse_hdf5_file_as_dict(fname: str) -> dict:
    """Parse the content of an hdf5 file as dictionary.
    Use this function only for small files.

    :param str fname: path to the file.
    :return dict h5_file_content: Content of the hdf5 file.
    """

    try:
        import h5py
    except ImportError:
        raise ImportError("h5py module not installed, but is required")

    def recursive_unpack(hdfobject, datadict):
        """Unpack an HDF5 data object to a dictionary recursively.

        :param h5py.Group hdfobject: HDF5 object to unpack
        :param dict datadict: Dictionary to unpack to.
        """
        for key, value in hdfobject.items():
            if isinstance(value, h5py.Group):
                datadict[key] = {}
                datadict[key] = recursive_unpack(value, datadict[key])

            elif isinstance(value, h5py.Dataset):
                datadict[key] = value[()]

        return datadict

    def convert_one_element_arrays(obj):
        """Replace np.array([a]) with a."""

        for key, value in obj.items():
            if isinstance(value, np.ndarray) and len(value) == 1:
                obj.update({key: value[0]})

            elif isinstance(value, dict):
                convert_one_element_arrays(value)

    file = h5py.File(fname)
    h5_file_content = recursive_unpack(file, {})
    if file:
        file.close()

    convert_one_element_arrays(h5_file_content)

    return h5_file_content
