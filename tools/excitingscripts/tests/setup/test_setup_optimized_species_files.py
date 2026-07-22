import pytest
from pathlib import Path
import sys
from excitingtools.species.species_file import SpeciesFile

TEST_DIR = Path(__file__).resolve().parent
SCRIPT_DIR = TEST_DIR.parent.parent / "excitingscripts" / "setup"

sys.path.insert(0, str(SCRIPT_DIR))

from optimized_species_files import optimize_species_from_folder

@pytest.fixture
def mock_files(tmp_path: Path) -> None:
    """Create an easy mock-file for the sake of testing"""
    species_name = "Ag"
    
    # 1. Maximized species file 
    xml_content = """<?xml version="1.0" encoding="UTF-8"?>
<spdb>
  <sp chemicalSymbol="Ag" name="silver" z="-47.0" mass="196631.7">
    <muffinTin rmin="0.00001" radius="2.0" points="400"/>
    <atomicState n="1" l="0" kappa="1" occ="2.0" core="true"/>
    <basis>
      <lo l="0">
        <wf matchingOrder="0" n="4"/>
        <wf matchingOrder="1" n="4"/>
      </lo>
      <lo l="0">
        <wf matchingOrder="1" n="5"/>
        <wf matchingOrder="2" n="5"/>
      </lo>
      <lo l="1">
        <wf matchingOrder="2" n="4"/>
        <wf matchingOrder="3" n="4"/>
      </lo>
    </basis>
  </sp>
</spdb>
"""
    (tmp_path / f"{species_name}.xml").write_text(xml_content, encoding="utf-8")
    # 2. LO-hierarchy
    out_content = """DBSV test precision for maximized species file: 6.291e-05

species    l   ns    mOs   precision    max iterations
Ag         0   5,5   1,2   1.000e-08    14
Ag         1   4,4   2,3   2.203e-03    16

Essential local orbitals:
species    l   ns    mOs
Ag         0   4,4   0,1
"""
    (tmp_path / f"lo_hierarchy_{species_name}.out").write_text(out_content, encoding="utf-8")


def test_optimize_species_from_folder(mock_files: None, tmp_path: Path) -> None:
    
    species_name = "Ag"
    precision_threshold= 1.0e-4
    
    optimize_species_from_folder(
        species_name=species_name,
        precision_threshold=precision_threshold,
        path_lo=tmp_path,
        path_xml=tmp_path,
        path_out=tmp_path
    )

    output_file = tmp_path / f"{species_name}_customized.xml"
    
    assert output_file.exists(), "no xml generated"

    species_obj =SpeciesFile.from_file(output_file)
    lo_elements = species_obj.basis["lo"]    
 
    assert len(lo_elements) == 2, f"Two remaining LO expected, but {len(lo_elements)} found."

    l_values = [lo["l"] for lo in lo_elements]
    
    assert l_values == [0, 1], f"Expected exactly LOs with l=[0, 1], but got {l_values}"   
 
    for lo in lo_elements:
        if lo["l"] == 0:
            wfs = lo["wf"]
            assert wfs[0]["n"] == 4, "wrong l=0 orbital deleted"
