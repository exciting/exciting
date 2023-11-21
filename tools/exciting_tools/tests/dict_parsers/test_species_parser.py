from excitingtools.exciting_dict_parsers.species_parser import parse_species_xml


def test_parse_species_xml():
    """ Test parsing of species.xml files.
    """
    assert isinstance(species_str, str), (
        "Expect species parser to handle strings of XML data, "
        "due to use of decorator"
    )

    species_dict = parse_species_xml(species_str)

    assert set(species_dict) == {'species', 'muffin_tin', 'atomic_states', 'basis'}, 'Top level species file keys'
    assert species_dict['species'] == {
        'chemicalSymbol': 'Zn', 'name': 'zinc', 'z': -30.0, 'mass': 119198.678
    }
    assert species_dict['muffin_tin'] == {
        'rmin': 1e-06, 'radius': 2.0, 'rinf': 21.8982, 'radialmeshPoints': 600.0
    }

    assert isinstance(species_dict['atomic_states'], list), 'Atomic states stored as a list'
    atomic_states = [{'n': 1, 'l': 0, 'kappa': 1, 'occ': 2.00000,
                      'core': True}, {'n': 2, 'l': 0, 'kappa': 1, 'occ': 2.00000, 'core': True},
                     {'n': 2, 'l': 1, 'kappa': 1, 'occ': 2.00000,
                      'core': True}, {'n': 2, 'l': 1, 'kappa': 2, 'occ': 4.00000, 'core': True},
                     {'n': 3, 'l': 0, 'kappa': 1, 'occ': 2.00000,
                      'core': False}, {'n': 3, 'l': 1, 'kappa': 1, 'occ': 2.00000, 'core': False},
                     {'n': 3, 'l': 1, 'kappa': 2, 'occ': 4.00000,
                      'core': False}, {'n': 3, 'l': 2, 'kappa': 2, 'occ': 4.00000, 'core': False},
                     {'n': 3, 'l': 2, 'kappa': 3, 'occ': 6.00000,
                      'core': False}, {'n': 4, 'l': 0, 'kappa': 1, 'occ': 2.00000, 'core': False}]
    assert species_dict['atomic_states'] == atomic_states

    basis = species_dict['basis']
    assert set(basis) == {'default', 'custom', 'lo'}, 'Keys for basis'

    assert basis['default'] == [{'type': 'lapw', 'trialEnergy': 0.1500, 'searchE': True}]

    # Custom apw, lapw or apw+lo
    assert basis['custom'] == [
        {'l': 0, 'type': 'lapw', 'trialEnergy': 1.35670550183736, 'searchE': False},
        {'l': 1, 'type': 'lapw', 'trialEnergy': -2.69952312512447,
         'searchE': False}, {'l': 2, 'type': 'lapw', 'trialEnergy': 0.00, 'searchE': False},
        {'l': 3, 'type': 'lapw', 'trialEnergy': 1.000,
         'searchE': False}, {'l': 4, 'type': 'lapw', 'trialEnergy': 1.000, 'searchE': False},
        {'l': 5, 'type': 'lapw', 'trialEnergy': 1.000, 'searchE': False}
    ]

    # All explicitly specified LOs
    los = [{'l': 0, 'wf': [{'matchingOrder': 0, 'searchE': False, 'trialEnergy': -4.37848525995355},
                           {'matchingOrder': 1, 'searchE': False, 'trialEnergy': -4.37848525995355}]},
           {'l': 0, 'wf': [{'matchingOrder': 0, 'searchE': False, 'trialEnergy': 1.35670550183736},
                           {'matchingOrder': 1, 'searchE': False, 'trialEnergy': 1.35670550183736}]},
           {'l': 0, 'wf': [{'matchingOrder': 0, 'searchE': False, 'trialEnergy': 1.35670550183736},
                           {'matchingOrder': 0, 'searchE': False, 'trialEnergy': -4.37848525995355}]},
           {'l': 0, 'wf': [{'matchingOrder': 1, 'searchE': False, 'trialEnergy': 1.35670550183736},
                           {'matchingOrder': 2, 'searchE': False, 'trialEnergy': 1.35670550183736}]},
           {'l': 1, 'wf': [{'matchingOrder': 0, 'searchE': False, 'trialEnergy': -2.69952312512447},
                           {'matchingOrder': 1, 'searchE': False, 'trialEnergy': -2.69952312512447}]},
           {'l': 1, 'wf': [{'matchingOrder': 1, 'searchE': False, 'trialEnergy': -2.69952312512447},
                           {'matchingOrder': 2, 'searchE': False, 'trialEnergy': -2.69952312512447}]},
           {'l': 2, 'wf': [{'matchingOrder': 0, 'searchE': False, 'trialEnergy': 0.0},
                           {'matchingOrder': 1, 'searchE': False, 'trialEnergy': 0.0}]},
           {'l': 2, 'wf': [{'matchingOrder': 1, 'searchE': False, 'trialEnergy': 0.0},
                           {'matchingOrder': 2, 'searchE': False, 'trialEnergy': 0.0}]},
           {'l': 3, 'wf': [{'matchingOrder': 0, 'searchE': False, 'trialEnergy': 1.0},
                           {'matchingOrder': 1, 'searchE': False, 'trialEnergy': 1.0}]},
           {'l': 3, 'wf': [{'matchingOrder': 1, 'searchE': False, 'trialEnergy': 1.0},
                           {'matchingOrder': 2, 'searchE': False, 'trialEnergy': 1.0}]},
           {'l': 4, 'wf': [{'matchingOrder': 0, 'searchE': False, 'trialEnergy': 1.0},
                           {'matchingOrder': 1, 'searchE': False, 'trialEnergy': 1.0}]},
           {'l': 4, 'wf': [{'matchingOrder': 1, 'searchE': False, 'trialEnergy': 1.0},
                           {'matchingOrder': 2, 'searchE': False, 'trialEnergy': 1.0}]},
           {'l': 5, 'wf': [{'matchingOrder': 0, 'searchE': False, 'trialEnergy': 1.0},
                           {'matchingOrder': 1, 'searchE': False, 'n': 5}]},
           {'l': 5, 'wf': [{'matchingOrder': 1, 'n': 5, 'searchE': False},
                           {'matchingOrder': 2, 'n': 6, 'searchE': False}]}]
    
    assert len(basis['lo']) == 14, "Number of explicitly-defined local orbitals"
    assert set(basis['lo'][0]) == {'l', 'wf'}, "Attributes defining a local orbital"
    assert basis['lo'] == los


def test_parse_species_xml_different_ordering():
    species_dict = parse_species_xml(species_str_diff_order)

    assert set(species_dict) == {'species', 'muffin_tin', 'atomic_states', 'basis'}, 'Top level species file keys'

    # References
    ref_species = {'chemicalSymbol': 'C', 'name': 'carbon', 'z': -6.0, 'mass': 21894.16673}
    ref_muffin_tin = {'rmin': 1e-05, 'radius': 1.45, 'rinf': 21.0932, 'radialmeshPoints': 250.0}
    ref_atomic_states = [{'n': 1, 'l': 0, 'kappa': 1, 'occ': 2.0, 'core': False},
                         {'n': 2, 'l': 0, 'kappa': 1, 'occ': 2.0, 'core': False},
                         {'n': 2, 'l': 1, 'kappa': 1, 'occ': 1.0, 'core': False},
                         {'n': 2, 'l': 1, 'kappa': 2, 'occ': 1.0, 'core': False}]
    ref_basis = {'default': [{'type': 'lapw', 'trialEnergy': 0.15, 'searchE': False}],
                 'custom': [{'l': 0, 'type': 'apw+lo', 'trialEnergy': 0.15, 'searchE': True},
                            {'l': 1, 'type': 'apw+lo', 'trialEnergy': 0.15, 'searchE': True}],
                 'lo': [{'l': 0, 'wfproj': False,
                         'wf': [{'matchingOrder': 0, 'trialEnergy': 0.15, 'searchE': True},
                                {'matchingOrder': 1, 'trialEnergy': 0.15, 'searchE': True},
                                {'matchingOrder': 0, 'trialEnergy': -9.8, 'searchE': False}]}]}
    
    assert species_dict['species'] == ref_species, "species data disagrees"
    assert species_dict['muffin_tin'] == ref_muffin_tin, "muffin tin data disagrees"
    assert species_dict['atomic_states'] == ref_atomic_states, "atomic state data disagrees"
    assert species_dict['basis']['default'] == ref_basis['default'], "basis default data disagrees"
    assert species_dict['basis']['custom'] == ref_basis['custom'], "basis custom data disagrees"
    assert species_dict['basis']['lo'] == ref_basis['lo'], "basis lo data disagrees"


species_str = """<?xml version="1.0" encoding="utf-8"?>
<spdb xsi:noNamespaceSchemaLocation="../../xml/species.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <sp chemicalSymbol="Zn" name="zinc" z="-30.0000" mass="119198.6780">
    <muffinTin rmin="0.100000E-05" radius="2.0000" rinf="21.8982" radialmeshPoints="600"/>
    <atomicState n="1" l="0" kappa="1" occ="2.00000" core="true"/>
    <atomicState n="2" l="0" kappa="1" occ="2.00000" core="true"/>
    <atomicState n="2" l="1" kappa="1" occ="2.00000" core="true"/>
    <atomicState n="2" l="1" kappa="2" occ="4.00000" core="true"/>
    <atomicState n="3" l="0" kappa="1" occ="2.00000" core="false"/>
    <atomicState n="3" l="1" kappa="1" occ="2.00000" core="false"/>
    <atomicState n="3" l="1" kappa="2" occ="4.00000" core="false"/>
    <atomicState n="3" l="2" kappa="2" occ="4.00000" core="false"/>
    <atomicState n="3" l="2" kappa="3" occ="6.00000" core="false"/>
    <atomicState n="4" l="0" kappa="1" occ="2.00000" core="false"/>
    <basis>
      <default type="lapw" trialEnergy="0.1500" searchE="true"/>

      <custom l="0" type="lapw" trialEnergy="1.35670550183736" searchE="false"/>
      <lo l="0">
        <wf matchingOrder="0" trialEnergy="-4.37848525995355" searchE="false"/>
	<wf matchingOrder="1" trialEnergy="-4.37848525995355" searchE="false"/>
      </lo>
      <lo l="0">
        <wf matchingOrder="0" trialEnergy="1.35670550183736" searchE="false"/>
        <wf matchingOrder="1" trialEnergy="1.35670550183736" searchE="false"/>
      </lo>
      <lo l="0">
        <wf matchingOrder="0" trialEnergy="1.35670550183736" searchE="false"/>
	<wf matchingOrder="0" trialEnergy="-4.37848525995355" searchE="false"/>
      </lo>
      <lo l="0">
        <wf matchingOrder="1" trialEnergy="1.35670550183736" searchE="false"/>
        <wf matchingOrder="2" trialEnergy="1.35670550183736" searchE="false"/>
      </lo>

      <custom l="1" type="lapw" trialEnergy="-2.69952312512447" searchE="false"/>
      <lo l="1">
        <wf matchingOrder="0" trialEnergy="-2.69952312512447" searchE="false"/>
        <wf matchingOrder="1" trialEnergy="-2.69952312512447" searchE="false"/>
      </lo>
      <lo l="1">
	<wf matchingOrder="1" trialEnergy="-2.69952312512447" searchE="false"/>
	<wf matchingOrder="2" trialEnergy="-2.69952312512447" searchE="false"/>
      </lo>

      <custom l="2" type="lapw" trialEnergy="0.00" searchE="false"/>
      <lo l="2">
	<wf matchingOrder="0" trialEnergy="0.00" searchE="false"/>
        <wf matchingOrder="1" trialEnergy="0.00" searchE="false"/>
      </lo>
      <lo l="2">
	<wf matchingOrder="1" trialEnergy="0.00" searchE="false"/>
        <wf matchingOrder="2" trialEnergy="0.00" searchE="false"/>
      </lo>

      <custom l="3" type="lapw" trialEnergy="1.000" searchE="false"/>
      <lo l="3">
	<wf matchingOrder="0" trialEnergy="1.000" searchE="false"/>
	<wf matchingOrder="1" trialEnergy="1.000" searchE="false"/>
      </lo>
      <lo l="3">
        <wf matchingOrder="1" trialEnergy="1.000" searchE="false"/>
        <wf matchingOrder="2" trialEnergy="1.000" searchE="false"/>
      </lo>

      <custom l="4" type="lapw" trialEnergy="1.000" searchE="false"/>
      <lo l="4">
	<wf matchingOrder="0" trialEnergy="1.000" searchE="false"/>
	<wf matchingOrder="1" trialEnergy="1.000" searchE="false"/>
      </lo>
      <lo l="4">
        <wf matchingOrder="1" trialEnergy="1.000" searchE="false"/>
        <wf matchingOrder="2" trialEnergy="1.000" searchE="false"/>
      </lo>

      <custom l="5" type="lapw" trialEnergy="1.000" searchE="false"/>
      <lo l="5">
	<wf matchingOrder="0" trialEnergy="1.000" searchE="false"/>
  <wf matchingOrder="1" n="5" searchE="false"/>      
      </lo>
      <lo l="5">
        <wf matchingOrder="1" n="5" searchE="false"/>
        <wf matchingOrder="2" n="6" searchE="false"/>
      </lo>

    </basis>
  </sp>
</spdb>
"""

species_str_diff_order = """<?xml version="1.0" encoding="UTF-8"?>
<spdb xsi:noNamespaceSchemaLocation="../../xml/species.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <sp chemicalSymbol="C" name="carbon" z="-6.00000" mass="21894.16673">
    <atomicState n="1" l="0" kappa="1" occ="2.00000" core="false"/>
    <atomicState n="2" l="0" kappa="1" occ="2.00000" core="false"/>
    <atomicState n="2" l="1" kappa="1" occ="1.00000" core="false"/>
    <atomicState n="2" l="1" kappa="2" occ="1.00000" core="false"/>
    <basis>
      <default type="lapw" trialEnergy="0.1500" searchE="false"/>
      <custom l="0" type="apw+lo" trialEnergy="0.1500" searchE="true"/>
      <lo l="0" wfproj="false">
	<wf matchingOrder="0" trialEnergy="0.1500" searchE="true"/>
	<wf matchingOrder="1" trialEnergy="0.1500" searchE="true"/>
	<wf matchingOrder="0" trialEnergy="-9.8" searchE="false"/>
      </lo>
      <custom l="1" type="apw+lo" trialEnergy="0.1500" searchE="true"/>
    </basis>
    <muffinTin rmin="0.100000E-04" radius="1.4500" rinf="21.0932" radialmeshPoints="250"/>
  </sp>
</spdb>
"""
