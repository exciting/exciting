"""
! **************************************************************************************************
!  Copyright (C) 2020-2023 GreenX library
!  This file is distributed under the terms of the APACHE2 License.
!
! **************************************************************************************************

Inputs and assertions for test_gx_minimax_grid.f90

Would be better if one had python API to GX, and avoid the need for this intermediate binary.
Indeed, MUCH of the complication regarding running tests comes from this.

Issues and solutions:
---------------------
* Where should the python test write the input data for the fortran binary?
  - This is solved using tmp_path, which gives an absolute path.
"""
import numpy as np
from pathlib import Path
import pytest
import os

from pygreenx.run import BinaryRunner, BuildType


@pytest.fixture()
@pytest.mark.usefixtures("get_binary", "greenx_build_root")
def fortran_binary(get_binary, greenx_build_root):
    name = 'test_gx_minimax_grid.exe'
    _binary = get_binary(name)
    assert _binary is not None, f'{name} cannot be found in {greenx_build_root}'
    print(f'Binary source: {_binary}')
    return _binary


def mock_file(path, inputs_str):
    """ Write an inputs string to file.
    :param path:
    :param inputs_str:
    :return: Path object
    """
    file: Path = Path(path) / Path("inputs.dat")
    file.write_text(inputs_str)
    assert file.exists(), f"{file.name} does not exist"
    return file


class ETrans:
    def __init__(self, min, max):
        self.min = min
        self.max = max


# Note. This could be a fixture
def get_grids_and_weights(fortran_binary, n_mesh_points, e_trans: ETrans):

    # Fortran input file
    inputs_template = """
{}
n_mesh_points     {}
e_transition_min  {}
e_transition_max  {}
"""

    working_dir = os.getcwd()

    inputs_str = inputs_template.format(working_dir, n_mesh_points, e_trans.min, e_trans.max)
    file = mock_file(working_dir, inputs_str)

    # Run test
    runner = BinaryRunner(fortran_binary, BuildType.serial, args=[file.as_posix()])
    results = runner.run()
    assert results.success, f"Execution of {fortran_binary} failed"

    # Parse results
    tau_file = "tau.dat"
    freq_file = "freq.dat"

    return np.genfromtxt(tau_file), np.genfromtxt(freq_file)


# @pytest.mark.usefixtures("testdir")
def test_get_grids_and_weights(fortran_binary):
    """
    TODO(Maryam) Issue 24. Extend the test inputs, and include fringe cases.
    """
    # Test inputs
    n_mesh_points = 10
    e_trans = ETrans(2., 30.)

    tau_grid_weights, freq_grid_weights = get_grids_and_weights(fortran_binary,
                                                                n_mesh_points,
                                                                e_trans)

    ref_tau_grid_weights = np.array([[0.005918, 0.01525399],
                                     [0.03177845, 0.03688566],
                                     [0.08092057, 0.0622535],
                                     [0.15862351, 0.09469842],
                                     [0.27439677, 0.13947465],
                                     [0.44433573, 0.20467764],
                                     [0.69464678, 0.30283334],
                                     [1.06797564, 0.45564411],
                                     [1.63956744, 0.7122558],
                                     [2.58159534, 1.25940329]])

    ref_freq_grid_weights = np.array([[0.37397426, 0.75729004],
                                      [1.17910655, 0.87356340],
                                      [2.16745299, 1.12983421],
                                      [3.50267933, 1.57927297],
                                      [5.42267468, 2.32153469],
                                      [8.30517692, 3.54877702],
                                      [12.81707194, 5.68714236],
                                      [20.35677335, 9.94234857],
                                      [34.86933737, 21.24030777],
                                      [75.79619716, 80.07961331]])

    assert np.allclose(tau_grid_weights, ref_tau_grid_weights)
    assert np.allclose(freq_grid_weights, ref_freq_grid_weights)
