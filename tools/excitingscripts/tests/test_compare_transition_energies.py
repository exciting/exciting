import sys; print(sys.executable)

import numpy as np
import pytest
from excitingscripts.compare_transition_energies import determine_transition_energies
from excitingtools.utils.test_utils import MockFile

@pytest.fixture
def eigval_xml_mock(tmp_path) -> MockFile:
    """ Mock 'eigval.xml' data, containing only two k-sampling points.
    """
    eigval_xml_str = """<?xml version="1.0" encoding="UTF-8"?>
    <eigval>
      <kpt ik="1" vkl="0.000000000000e0 0.000000000000e0 0.000000000000e0">
        <state ist="1" eigenvalue="-4.719487348563e-1" occupancy="2.000000000000e0"/>
        <state ist="2" eigenvalue="3.898020901187e-1" occupancy="2.000000000000e0"/>
        <state ist="3" eigenvalue="3.898092645211e-1" occupancy="2.000000000000e0"/>
        <state ist="4" eigenvalue="3.898272209373e-1" occupancy="2.000000000000e0"/>
        <state ist="5" eigenvalue="6.515326498206e-1" occupancy="0.000000000000e0"/>
        <state ist="6" eigenvalue="6.515432261724e-1" occupancy="0.000000000000e0"/>
        <state ist="7" eigenvalue="6.515470146363e-1" occupancy="0.000000000000e0"/>
        <state ist="8" eigenvalue="9.530499684716e-1" occupancy="0.000000000000e0"/>
        <state ist="9" eigenvalue="1.165753467181e0" occupancy="0.000000000000e0"/>
        <state ist="10" eigenvalue="1.435123611263e0" occupancy="0.000000000000e0"/>
      </kpt>
      <kpt ik="2" vkl="0.000000000000e0 5.000000000000e-1 5.000000000000e-1">
        <state ist="1" eigenvalue="-1.243547220790e-1" occupancy="2.000000000000e0"/>
        <state ist="2" eigenvalue="-1.243526098397e-1" occupancy="2.000000000000e0"/>
        <state ist="3" eigenvalue="1.438428832526e-1" occupancy="2.000000000000e0"/>
        <state ist="4" eigenvalue="1.438448580139e-1" occupancy="2.000000000000e0"/>
        <state ist="5" eigenvalue="6.147426980054e-1" occupancy="0.000000000000e0"/>
        <state ist="6" eigenvalue="6.147433158914e-1" occupancy="0.000000000000e0"/>
        <state ist="7" eigenvalue="1.108334812219e0" occupancy="0.000000000000e0"/>
        <state ist="8" eigenvalue="1.108338304192e0" occupancy="0.000000000000e0"/>
        <state ist="9" eigenvalue="1.343653569414e0" occupancy="0.000000000000e0"/>
        <state ist="10" eigenvalue="1.343653905367e0" occupancy="0.000000000000e0"/>
      </kpt>
    </eigval> 
    """
    eigval_xml_file = tmp_path / "eigval.xml"
    eigval_xml_file.write_text(eigval_xml_str)
    return MockFile(eigval_xml_file, eigval_xml_str)

def test_determine_transition_energies(eigval_xml_mock, tmp_path):
    assert np.allclose( np.array(determine_transition_energies(root_directory=tmp_path)), 
                        np.array([7.121367508015419, 6.120261919200727]))
