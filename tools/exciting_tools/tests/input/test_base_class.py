"""Test exiting input base class and its methods.

NOTE:
All attribute tests should assert on the XML tree content's as the attribute
order is not preserved by the ElementTree.tostring method. Elements appear to
be fine.
"""

from pathlib import Path

from excitingtools.input.base_class import query_exciting_version


def test_query_exciting_version(tmp_path: Path) -> None:
    """
    Test querying the exciting version.
    """
    version_inc_contents = """#define GITHASH "1a2087b0775a87059d53"
#define GITHASH2 "5d01a5475a10f00d0ad7"
#define COMPILERVERSION "GNU Fortran (MacPorts gcc9 9.3.0_4) 9.3.0"
#define VERSIONFROMDATE /21,12,01/
    """
    # Mock the version.inc file, and prepended path
    src = tmp_path / "src"
    src.mkdir()

    version_inc = tmp_path / "src" / "version.inc"
    version_inc.write_text(version_inc_contents)

    mod_misc = src / "mod_misc.F90"
    mod_misc.write_text("  !> Code version\n  character(40) :: versionname = 'NEON'\n")

    versioning: dict = query_exciting_version(tmp_path)
    assert set(versioning.keys()) == {"compiler", "git_hash", "major"}, (
        "Expect `query_exciting_version` to return compiler used "
        "for last build, exciting git hash and major version name."
    )

    assert versioning["compiler"] == "GNU Fortran (MacPorts gcc9 9.3.0_4) 9.3.0"
    assert versioning["git_hash"] == "1a2087b0775a87059d535d01a5475a10f00d0ad7"
    assert versioning["major"] == "NEON"
