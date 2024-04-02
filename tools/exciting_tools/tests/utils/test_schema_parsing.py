"""Test the schema parsing and whether the file valid_attributes is up-to-date."""

from excitingtools.utils.schema_parsing import get_excitingtools_root, list_string_line_limit


def test_set_string_line_limit():
    properties_valid_subtrees = [
        "DFTD2",
        "EFG",
        "LSJ",
        "TSvdW",
        "bandstructure",
        "boltzequ",
        "chargedensityplot",
        "momentummatrix",
        "mossbauer",
        "mvecfield",
        "polarization",
        "raman",
        "shg",
        "spintext",
        "stm",
        "wannier",
        "wanniergap",
        "wannierplot",
        "wfplot",
        "xcmvecfield",
    ]
    reference_string = (
        'properties_valid_subtrees = ["DFTD2", "EFG", "LSJ", "TSvdW", '
        '"bandstructure", "boltzequ", "chargedensityplot", \n'
        '                             "momentummatrix", "mossbauer", "mvecfield", '
        '"polarization", "raman", "shg", "spintext", \n'
        '                             "stm", "wannier", "wanniergap", "wannierplot", '
        '"wfplot", "xcmvecfield"]'
    )
    assert list_string_line_limit("properties_valid_subtrees", properties_valid_subtrees) == reference_string


def test_get_exciting_root():
    assert get_excitingtools_root().name == "exciting_tools"
