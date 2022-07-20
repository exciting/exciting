"""
Test regex parser wrappers
"""
import pytest

from excitingtools.parser_utils.regex_parser import parse_value_regex


def test_parse_value_regex():

    test_string = """
    --------------------------------------------------------------------------------
    -                               frequency grid                                 -
    --------------------------------------------------------------------------------
     
     Type: < fgrid > gauleg2                                 
     Frequency axis: < fconv > imfreq                                  
     Number of frequencies: < nomeg >           32
     Cutoff frequency: < freqmax >    1.00000000000000  
    
    --------------------------------------------------------------------------------
    -                          Kohn-Sham band structure                            -
    --------------------------------------------------------------------------------
     
     Fermi energy:     0.0000
     Energy range:   -14.6863 1030.7919
     Band index of VBM:  21
     Band index of CBm:  22
     
     Indirect BandGap (eV):                    3.3206
     at k(VBM) =    0.000   0.500   0.500 ik =     3
        k(CBm) =    0.000   0.000   0.000 ik =     1
     Direct Bandgap at k(VBM) (eV):            3.7482
     Direct Bandgap at k(CBm) (eV):            3.8653
     
     --------------------------------------------------------------------------------
     -                            G0W0 band structure                               -
     --------------------------------------------------------------------------------
      
      Fermi energy:    -0.0054
      Energy range:   -16.2632 1031.4090
      Band index of VBM:  21
      Band index of CBm:  22
      
      Indirect BandGap (eV):                    5.3920
      at k(VBM) =    0.000   0.500   0.500 ik =     3
         k(CBm) =    0.000   0.000   0.000 ik =     1
      Direct Bandgap at k(VBM) (eV):            5.5472
      Direct Bandgap at k(CBm) (eV):            5.9646
    """

    assert {'Fermi energy': 0.0} == parse_value_regex(test_string, 'Fermi energy:'), \
        "Expect to match first string instance and return {key:float}"

    assert {'Indirect BandGap (eV)': 3.3206} == parse_value_regex(test_string, 'Indirect BandGap \\(eV\\):'), \
       "Expect to match first string instance and return {key:float}"

    assert {'Energy range': [-14.6863, 1030.7919]} == parse_value_regex(test_string, 'Energy range:'), \
        "Expect to match first string instance and return {key:List[float]}"

    assert {'Type: < fgrid >': 'gauleg2'} == parse_value_regex(test_string, 'Type: < fgrid >'), \
        "Expect to match first string instance and return {key:str}"


@pytest.mark.xfail
def test_failures_with_parse_value_regex():
    """
    This fails for timings because matching 'Initialization' returns ':        15.46',
    (for example), and the colon cannot be interpreted by eval()
    """
    test_string = """
     Initialization                             :        15.46
         - init_scf                             :         8.38
    """

    data = parse_value_regex(test_string, 'Initialization')
