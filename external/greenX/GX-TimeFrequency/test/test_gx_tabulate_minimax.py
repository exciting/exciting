"""Testing minimax grid calculation
   Copyright (C) 2020-2023 GreenX library
   This file is distributed under the terms of the APACHE2 License.
"""
import enum
import numpy as np
import pytest

from pygreenx.run import BinaryRunner, BuildType, ProcessResults


@pytest.fixture()
@pytest.mark.usefixtures("get_binary", "greenx_build_root")
def fortran_binary(get_binary, greenx_build_root):
    name = 'gx_tabulate_grids.exe'
    _binary = get_binary(name)
    assert _binary is not None, f'{name} cannot be found in {greenx_build_root}'
    print(f'Binary source: {_binary}')
    return _binary


class ETrans:
    def __init__(self, min, max):
        self.min = min
        self.max = max


def parse_std_out(results: ProcessResults):
    """ Parse std.out of gx_tabulate_grids.exe

    Note, if the output format changes in ANY way, this breaks.
    One should switch to structured data.

    :param results subprocess results container.
    :return: Numpy array of errors.
    """
    # Data lines in output
    header = 5
    n_grids = 15
    n_columns = 7

    # Parse std out - only interested in the numerical data
    output_list = results.stdout.decode("utf-8").split('\n')[header:header + n_grids]

    tabulated_errors = np.empty(shape=(n_grids, n_columns))
    for i, line in enumerate(output_list[:]):
        # print(line)
        tabulated_errors[i, :] = np.asarray([float(x) for x in line.split()])

    return tabulated_errors


def get_tabulated_errors(fortran_binary, e_trans: ETrans):
    runner = BinaryRunner(fortran_binary, BuildType.serial,
                          args=['table', '-emin', f'{e_trans.min}', '-emax', f'{e_trans.max}'])
    results = runner.run()
    assert results.success, f"Execution of {fortran_binary} failed"
    return parse_std_out(results)


class Column(enum.Enum):
    """Reference Data Column Labels.
    """
    NumPoints = 0
    IErr = 1
    CosFTDualityError = 2
    MaxErrCosFTTimeToFreq = 3
    MaxErrCosFTFreqToTime = 4
    MaxErrSinFTimeToFreq = 5
    ERatio = 6


def test_tabulate_gx_minimax_grid(fortran_binary):

    e_trans = ETrans(0.4, 100.0)

    ref_errors_small_grid = np.array([[6,  0, 2.51200821E-04, 1.88667240E-03, 6.86380679E-03, 5.56242714E-02, 2.50000000E+02],
                                      [8,  0, 9.26518238E-04, 7.02310800E-04, 8.97132915E-04, 1.49383520E-02, 2.50000000E+02],
                                      [10, 0, 2.93904398E-02, 4.15097815E-04, 9.85375147E-04, 3.87442807E-03, 2.50000000E+02],
                                      [12, 0, 1.22980779E-03, 5.40966369E-05, 3.23113268E-05, 9.86318209E-04, 2.50000000E+02],
                                      [14, 0, 1.23655331E-03, 1.47661838E-05, 6.49328520E-06, 2.46521564E-04, 2.50000000E+02],
                                      [16, 0, 3.23839478E-03, 4.75444516E-06, 7.88188993E-07, 7.70304736E-05, 2.50000000E+02],
                                      [18, 0, 6.23300174E-03, 1.34768466E-06, 3.45186001E-07, 1.94942501E-05, 2.50000000E+02],
                                      [20, 0, 2.62469251E-03, 2.86296043E-07, 3.67016068E-08, 4.88850816E-06, 2.50000000E+02],
                                      [22, 0, 1.90928278E-02, 6.87320737E-08, 2.59110779E-07, 2.68402248E-07, 2.50000000E+02]])

    ref_errors_big_grid = np.array([[24, 0, 2.49244180E-02, 1.41793908E-08, 1.44000650E-08, 1.31534820E-07, 2.50000000E+02],
                                    [26, 0, 2.64119033E-01, 4.12358436E-09, 4.74668548E-09, 3.78405299E-08, 2.50000000E+02],
                                    [28, 0, 4.13766945E+04, 3.93976513E-09, 5.94764298E-09, 4.30436032E-08, 2.50000000E+02],
                                    [30, 0, 1.24498840E+06, 4.50563448E-09, 3.15970462E-08, 5.14232632E-08, 2.50000000E+02],
                                    [32, 0, 1.29256892E+06, 4.17286994E-09, 6.30107694E-08, 4.94132290E-08, 2.50000000E+02],
                                    [34, 0, 1.32139161E+06, 8.10191277E-09, 3.86664983E-07, 6.66088501E-08, 2.50000000E+02]])

    tabulated_errors = get_tabulated_errors(fortran_binary, e_trans)
    tabulated_errors_small_grids = tabulated_errors[0:9, :]
    tabulated_errors_large_grids = tabulated_errors[9:, :]

    # Grids for which values do not get large or very small
    assert np.allclose(tabulated_errors_small_grids[:, Column.CosFTDualityError.value],
                       ref_errors_small_grid[:, Column.CosFTDualityError.value])

    assert np.allclose(tabulated_errors_small_grids[:, Column.MaxErrCosFTTimeToFreq.value],
                       ref_errors_small_grid[:, Column.MaxErrCosFTTimeToFreq.value])

    assert np.allclose(tabulated_errors_small_grids[:, Column.MaxErrCosFTFreqToTime.value],
                       ref_errors_small_grid[:, Column.MaxErrCosFTFreqToTime.value])

    assert np.allclose(tabulated_errors_small_grids[:, Column.MaxErrSinFTimeToFreq.value],
                       ref_errors_small_grid[:, Column.MaxErrSinFTimeToFreq.value])

    # Grids for which values do get large or very small
    assert np.allclose(tabulated_errors_large_grids[:, Column.MaxErrCosFTTimeToFreq.value],
                       ref_errors_big_grid[:, Column.MaxErrCosFTTimeToFreq.value],
                       atol=1.e-7)

    assert np.allclose(tabulated_errors_large_grids[:, Column.MaxErrCosFTFreqToTime.value],
                       ref_errors_big_grid[:, Column.MaxErrCosFTFreqToTime.value],
                       atol=5.e-7)

    assert np.allclose(tabulated_errors_large_grids[:, Column.MaxErrSinFTimeToFreq.value],
                       ref_errors_big_grid[:, Column.MaxErrSinFTimeToFreq.value],
                       atol=1.e-6)

    # Alex gets a massive difference (~ 10%) grid 28 CosFTDualityError between
    # his Mac with openblas, and the Ubuntu CI with blas/lapack
    # Not sure if this is quantity is worth testing
