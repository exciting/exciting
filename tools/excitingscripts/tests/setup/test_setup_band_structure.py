from excitingscripts.setup.band_structure import setup_band_structure

from excitingtools.exciting_dict_parsers.input_parser import parse_element_xml
from excitingtools.exciting_obj_parsers.input_xml import parse_input_xml


def test_setup_band_structure(input_xml_mock, tmp_path):
    properties_ref = {
        "bandstructure": {"plot1d": {"path": {
            "steps": 100,
            "point":
                [
                    {"coord": [0.0, 0.0, 0.0], "label": "G"},
                    {"coord": [0.5, 0.0, 0.5], "label": "X"},
                    {"coord": [0.5, 0.25, 0.75], "label": "W"},
                    {"coord": [0.375, 0.375, 0.75], "label": "K"},
                    {"coord": [0.0, 0.0, 0.0], "label": "G"},
                    {"coord": [0.5, 0.5, 0.5], "label": "L"},
                    {"coord": [0.625, 0.25, 0.625], "label": "U"},
                    {"coord": [0.5, 0.25, 0.75], "label": "W"},
                    {"coord": [0.5, 0.5, 0.5], "label": "L"},
                    {"coord": [0.375, 0.375, 0.75], "label": "K", "breakafter": True},
                    {"coord": [0.625, 0.25, 0.625], "label": "U"},
                    {"coord": [0.5, 0.0, 0.5], "label": "X"}]}}}}

    setup_band_structure(input_xml_mock.full_path, root_directory=tmp_path, overwrite=True)
    parsed_input = parse_input_xml(tmp_path / "input.xml")
    parsed_input_properties = parsed_input.properties.to_xml()  # pylint: disable=no-member

    assert properties_ref == parse_element_xml(parsed_input_properties)
