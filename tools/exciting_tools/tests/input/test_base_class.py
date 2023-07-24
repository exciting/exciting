"""Test exiting input base class and its methods.

NOTE:
All attribute tests should assert on the XML tree content's as the attribute
order is not preserved by the ElementTree.tostring method. Elements appear to
be fine.
"""

import pathlib

from excitingtools.input.base_class import query_exciting_version


def test_query_exciting_version(tmp_path):
    """
    Test querying the exciting version.
    """
    version_inc_contents = """#define GITHASH "1a2087b0775a87059d53"
#define GITHASH2 "5d01a5475a10f00d0ad7"
#define COMPILERVERSION "GNU Fortran (MacPorts gcc9 9.3.0_4) 9.3.0"
#define VERSIONFROMDATE /21,12,01/
    """

    # Mock the version.inc file, and prepended path
    exciting_root = pathlib.Path(tmp_path)
    src = exciting_root / "src"
    src.mkdir()
    assert exciting_root.is_dir(), "exciting_root tmp_path directory does not exist"

    version_inc = exciting_root / "src" / "version.inc"
    version_inc.write_text(version_inc_contents)
    assert version_inc.exists(), "version.inc does not exist"

    versioning: dict = query_exciting_version(exciting_root)
    assert set(versioning.keys()) == {'compiler', 'git_hash'}, (
        "Expect `query_exciting_version` to return compiler used "
        "for last build, and exciting git hash")

    assert versioning['compiler'] == 'GNU Fortran (MacPorts gcc9 9.3.0_4) 9.3.0'
    assert versioning['git_hash'] == '1a2087b0775a87059d535d01a5475a10f00d0ad7'
