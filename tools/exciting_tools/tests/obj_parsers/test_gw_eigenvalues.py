import numpy as np
import pytest
from typing import Tuple

from excitingtools.exciting_obj_parsers.gw_eigenvalues import gw_eigenvalue_parser, _file_name, OxygenEvalQPColumns, \
    NitrogenEvalQPColumns
from excitingtools.dataclasses.eigenvalues import EigenValues

from excitingtools.dataclasses.data_structs import BandIndices
from excitingtools.utils.test_utils import MockFile


@pytest.fixture
def evalqp_nitrogen(tmp_path) -> Tuple[MockFile, BandIndices]:
    evalqp_string = """k-point #     1:    0.000000    0.000000    0.000000    0.015625
state    E_KS      E_HF       E_GW       Sx         Sc         Vxc         DE_HF        DE_GW       Znk
 10     0.19173    0.13755   -0.03139   -0.54312    0.02775   -0.48894   -0.05418   -0.22312    0.77301
 11     0.19173    0.13755   -0.03144   -0.54312    0.02769   -0.48894   -0.05418   -0.22317    0.77295
 12     0.19173    0.13755   -0.03143   -0.54312    0.02770   -0.48894   -0.05418   -0.22316    0.77289
 13     0.28440    0.45295    0.08701   -0.25690   -0.16165   -0.42545    0.16855   -0.19739    0.76784
 14     0.28440    0.45295    0.08704   -0.25690   -0.16162   -0.42545    0.16855   -0.19736    0.76796

k-point #     2:    0.000000    0.000000    0.250000    0.125000
 state    E_KS      E_HF       E_GW       Sx         Sc         Vxc         DE_HF        DE_GW       Znk
  10     0.04648   -0.06027   -0.17533   -0.56955    0.08039   -0.46280   -0.10674   -0.22181    0.72570
  11     0.16360    0.10244   -0.05906   -0.53950    0.03515   -0.47834   -0.06116   -0.22266    0.76772
  12     0.16360    0.10244   -0.05905   -0.53950    0.03515   -0.47834   -0.06116   -0.22265    0.76758
  13     0.26382    0.43091    0.06505   -0.28523   -0.16204   -0.45232    0.16709   -0.19877    0.77585
  14     0.31900    0.49772    0.12124   -0.24567   -0.17225   -0.42439    0.17872   -0.19777    0.76106

k-point #     3:    0.000000    0.000000    0.500000    0.062500
 state    E_KS      E_HF       E_GW       Sx         Sc         Vxc         DE_HF        DE_GW       Znk
  10    -0.06604   -0.23046   -0.28494   -0.60031    0.14134   -0.43589   -0.16442   -0.21890    0.70271
  11     0.14748    0.07621   -0.07605   -0.54935    0.04402   -0.47807   -0.07127   -0.22353    0.76467
  12     0.14748    0.07621   -0.07607   -0.54935    0.04398   -0.47807   -0.07127   -0.22356    0.76464
  13     0.24560    0.40782    0.04622   -0.29167   -0.15798   -0.45389    0.16222   -0.19938    0.77910
  14     0.31230    0.49063    0.11503   -0.21667   -0.17126   -0.39500    0.17833   -0.19727    0.76701
    """
    file = tmp_path / _file_name
    file.write_text(evalqp_string)
    return MockFile(file, string=evalqp_string), BandIndices(12, 13)


@pytest.fixture()
def evalqp_oxygen(tmp_path) -> Tuple[MockFile, BandIndices]:
    evalqp_string = """k-point #     1:    0.000000    0.000000    0.000000    0.125000
 state   E_KS       E_HF       E_GW       Sx         Re(Sc)     Im(Sc)     Vxc        DE_HF      DE_GW      Znk
  19    -0.08099   -0.28150   -0.12905   -0.97072    0.13696   -0.00000   -0.77021   -0.20051   -0.04806    0.75633
  20    -0.08099   -0.28150   -0.12891   -0.97072    0.13710   -0.00000   -0.77021   -0.20051   -0.04792    0.75570
  21    -0.08099   -0.28150   -0.12896   -0.97072    0.13704   -0.00000   -0.77021   -0.20051   -0.04797    0.75582
  22     0.06097    0.32400    0.08958   -0.36222   -0.22727   -0.00000   -0.62524    0.26303    0.02861    0.80010
  23     0.06097    0.32400    0.08957   -0.36222   -0.22728   -0.00000   -0.62524    0.26303    0.02860    0.80006
  24     0.17630    0.38074    0.18396   -0.27546   -0.19505    0.00000   -0.47989    0.20443    0.00766    0.81561
 
k-point #     2:    0.000000    0.000000    0.500000    0.500000
 state   E_KS       E_HF       E_GW       Sx         Re(Sc)     Im(Sc)     Vxc        DE_HF      DE_GW      Znk
  19    -0.15356   -0.36518   -0.19801   -0.93999    0.15132   -0.00008   -0.72837   -0.21162   -0.04445    0.73712
  20    -0.12231   -0.33123   -0.17047   -0.96038    0.14470   -0.00001   -0.75145   -0.20893   -0.04816    0.74983
  21    -0.12231   -0.33123   -0.17035   -0.96038    0.14482   -0.00001   -0.75145   -0.20893   -0.04804    0.74940
  22     0.08561    0.35313    0.11429   -0.38438   -0.23136    0.00001   -0.65190    0.26752    0.02868    0.79318
  23     0.08561    0.35313    0.11430   -0.38438   -0.23135    0.00001   -0.65190    0.26752    0.02869    0.79323
  24     0.16667    0.43794    0.19514   -0.48319   -0.23505    0.00001   -0.75446    0.27128    0.02848    0.78602
 
k-point #     3:    0.000000    0.500000    0.500000    0.375000
 state   E_KS       E_HF       E_GW       Sx         Re(Sc)     Im(Sc)     Vxc        DE_HF      DE_GW      Znk
  19    -0.11749   -0.30231   -0.15828   -0.93575    0.13070   -0.00001   -0.75094   -0.18481   -0.04078    0.75361
  20    -0.11749   -0.30231   -0.15818   -0.93575    0.13079   -0.00001   -0.75094   -0.18481   -0.04068    0.75304
  21    -0.06097   -0.25700   -0.10809   -0.99057    0.13377    0.00000   -0.79454   -0.19603   -0.04712    0.75686
  22     0.07678    0.30493    0.09569   -0.40229   -0.20437    0.00000   -0.63044    0.22815    0.01891    0.79516
  23     0.11412    0.40512    0.14613   -0.38364   -0.24999    0.00001   -0.67464    0.29100    0.03201    0.78053
  24     0.15691    0.42130    0.18404   -0.48012   -0.23004    0.00000   -0.74450    0.26438    0.02712    0.78971
"""
    file = tmp_path / _file_name
    file.write_text(evalqp_string)
    return MockFile(file, string=evalqp_string), BandIndices(21, 22)


def test_parse_evalqp_oxygen(evalqp_oxygen):
    """ Test that the parser correctly returns to the EigenValues object.
    """
    file, band_indices = evalqp_oxygen
    eigen_values: EigenValues = gw_eigenvalue_parser(file.full_path)

    ref_gw_eigenvalues = np.array([[-0.12905, -0.12891, -0.12896,  0.08958,  0.08957,  0.18396],
                                   [-0.19801, -0.17047, -0.17035,  0.11429,  0.11430,  0.19514],
                                   [-0.15828, -0.15818, -0.10809,  0.09569,  0.14613,  0.18404]])

    ref_k_points = np.array([[0.000000, 0.000000, 0.000000],
                             [0.000000, 0.000000, 0.500000],
                             [0.000000, 0.500000, 0.500000]])

    assert eigen_values.state_range.first_state == 19
    assert eigen_values.state_range.last_state == 24
    assert np.allclose(eigen_values.k_points, ref_k_points)
    assert eigen_values.k_indices == [1, 2, 3]
    assert np.allclose(eigen_values.all_eigenvalues, ref_gw_eigenvalues), "GW column eigenvalues, for all k-points"
    assert eigen_values.weights == [0.125, 0.5, 0.375000]


def test_parse_evalqp_oxygen_Znk(evalqp_oxygen):
    """ Test that the parser correctly returns to the EigenValues object.
    """
    file, band_indices = evalqp_oxygen
    eigen_values: EigenValues = gw_eigenvalue_parser(file.full_path, columns=OxygenEvalQPColumns.Znk)

    ref_gw_eigenvalues = np.array([[0.75633, 0.75570, 0.75582, 0.80010, 0.80006, 0.81561],
                                   [0.73712, 0.74983, 0.74940, 0.79318, 0.79323, 0.78602],
                                   [0.75361, 0.75304, 0.75686, 0.79516, 0.78053, 0.78971]])

    assert np.allclose(eigen_values.all_eigenvalues, ref_gw_eigenvalues), "Znk values, for all k-points"


def test_parse_evalqp_nitrogen(evalqp_nitrogen):
    """ Test that the parser correctly returns to the EigenValues object.
    """
    file, band_indices = evalqp_nitrogen
    eigen_values: EigenValues = gw_eigenvalue_parser(file.full_path, NitrogenEvalQPColumns.E_GW)

    ref_gw_eigenvalues = np.array([[-0.03139, -0.03144, -0.03143, 0.08701, 0.08704],
                                   [-0.17533, -0.05906, -0.05905, 0.06505, 0.12124],
                                   [-0.28494, -0.07605, -0.07607, 0.04622, 0.11503]])

    ref_k_points = np.array([[0.000000, 0.000000, 0.000000],
                             [0.000000, 0.000000, 0.250000],
                             [0.000000, 0.000000, 0.500000]])

    assert eigen_values.state_range.first_state == 10
    assert eigen_values.state_range.last_state == 14
    assert np.allclose(eigen_values.k_points, ref_k_points)
    assert eigen_values.k_indices == [1, 2, 3]
    assert np.allclose(eigen_values.all_eigenvalues, ref_gw_eigenvalues), "GW column eigenvalues, for all k-points"
    assert eigen_values.weights == [0.015625, 0.125, 0.0625]


def test_parse_evalqp_incompatible(evalqp_nitrogen):
    """ Test for the exception.
    Trying to parse _filename generated with nitrogen, but requesting oxygen column indexing.
    """
    file, band_indices = evalqp_nitrogen

    with pytest.raises(ValueError) as error:
        eigen_values: EigenValues = gw_eigenvalue_parser(file.full_path,
                                                         OxygenEvalQPColumns.E_GW)
    assert error.value.args[
               0] == "The requested data column is indexed according to exciting version OxygenEvalQPColumns," \
                     "which is not consistent with the columns of the parsed data. " \
                     "Check that your data was produced with the same code version."
