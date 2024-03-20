""" Test for dynamic class functions. """

from excitingtools.input.dynamic_class import class_constructor_string, class_name_uppercase, get_all_valid_subtrees, \
    give_class_dictionary


def test_get_all_valid_subtrees():
    ref_subtrees = {'crystal', 'shell', 'dfthalfparam', 'structure', 'symmetries',
                    'basevect', 'species', 'atom', 'LDAplusU'}
    assert set(get_all_valid_subtrees(["structure"])) == ref_subtrees


def test_class_dir():
    ref_dir = {'attributes': {'__doc__': 'Class for exciting structure input.',
                              '__module__': 'excitingtools.input.input_classes',
                              'name': 'structure'},
               'bases': '(ExcitingXMLInput, )'}
    assert give_class_dictionary("structure") == ref_dir


def test_class_constructor_string():
    test_input = {'Spin': {"bases": '(ExcitingXMLInput, )', "attributes": {
        "__doc__": "Class for exciting spin input.", 'name': "spin"}}}
    ref_string = ("from excitingtools.input.base_class import ExcitingXMLInput \n"
                  "ExcitingSpinInput = type('ExcitingSpinInput', (ExcitingXMLInput, ), "
                  "{'__doc__': 'Class for exciting spin input.', 'name': 'spin'}) \n")
    assert class_constructor_string(test_input) == ref_string


def test_class_name_uppercase():
    assert class_name_uppercase("groundstate") == "GroundState"
    assert class_name_uppercase("xs") == "XS"
    assert class_name_uppercase("structure") == "Structure"
    assert class_name_uppercase("LDAplusU") == "LDAplusU"
    assert class_name_uppercase("HartreeFock") == "HartreeFock"
    assert class_name_uppercase("EFG") == "EFG"
    assert class_name_uppercase("etCoeffComponents") == "EtCoeffComponents"
