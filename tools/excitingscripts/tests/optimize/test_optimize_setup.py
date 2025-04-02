import numpy as np
import pytest
from excitingscripts.optimize.setup import (
    setup_optimize_lattice,
)
from excitingtools.exciting_obj_parsers.input_xml import parse_input_xml
from excitingtools.utils.test_utils import MockFile
from numpy.testing import assert_allclose


@pytest.fixture
def input_xml_mock(tmp_path):
    """ Mock 'input.xml' data."""
    input_xml_str = """<?xml version="1.0" ?>
    <input>

       <title>Be: Lattice optimization</title>

       <structure speciespath="/home/exciting/species">

          <crystal scale="4.300">
             <basevect>  1.00000000  0.00000000  0.00000000 </basevect>
             <basevect> -0.50000000  0.86602540  0.00000000 </basevect>
             <basevect>  0.00000000  0.00000000  1.50000000 </basevect>
          </crystal>

          <species speciesfile="Be.xml" rmt="1.95">
             <atom coord="0.66666667 0.33333333 0.75000000"/>
             <atom coord="0.33333333 0.66666667 0.25000000"/>
          </species>

       </structure>

       <groundstate 
          ngridk="6 6 4"
          xctype="GGA_PBE_SOL">
       </groundstate>

       <relax/>

    </input>
    """
    input_xml_file = tmp_path / "Be_opt.xml"
    input_xml_file.write_text(input_xml_str)
    return MockFile(input_xml_file, input_xml_str)

@pytest.mark.skip(reason="depends on sgroup software, so fails in CI test.")
def test_setup_optimize_lattice(input_xml_mock, tmp_path):
    run_directory = tmp_path
    max_strain = 0.03
    num_dist_str = 5
    opt_index = 1

    setup_optimize_lattice(run_directory, max_strain, num_dist_str, input_xml_mock.full_path, opt_index)

    reference_base_vectors = [
        np.array([[0.97, 0, 0], [-0.485, 0.840044638, 0], [0, 0, 1.455]]),
        np.array([[0.985, 0, 0], [-0.4925, 0.8530350189999999, 0], [0, 0, 1.4775]]),
        np.array([[1.00001, 0, 0], [-0.500005, 0.866034060254, 0], [0, 0, 1.500015]]),
        np.array([[1.015, 0, 0], [-0.5075, 0.8790157809999999, 0], [0, 0, 1.5225]]),
        np.array([[1.03, 0, 0], [-0.515, 0.8920061619999999, 0], [0, 0, 1.545]])
    ]

    for i in range(5):
        parsed_input = parse_input_xml(f"{tmp_path}/VOL/vol_{i + 1}/input.xml")
        assert_allclose(np.array(parsed_input.structure.lattice), reference_base_vectors[i])
