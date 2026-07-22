from pathlib import Path
from xml.etree import ElementTree


RUN_DIRECTORY = Path(__file__).resolve().parent.parent / "run_excited_state_species_file_generation"
SPECIES_FILE = RUN_DIRECTORY / "Si_excited.xml"


def main() -> None:
    assert SPECIES_FILE.is_file(), f"Generated species file not found: {SPECIES_FILE}"

    species = ElementTree.parse(SPECIES_FILE).getroot().find(".//sp")
    assert species is not None, "Generated species file does not contain an sp element"
    assert species.get("chemicalSymbol") == "Si", "Generated species file is not for silicon"

    custom_lapws = {
        (int(custom.get("l")), int(custom.get("n")))
        for custom in species.findall("./basis/custom")
        if custom.get("type") == "lapw"
    }
    expected_custom_lapws = {(l, l + 1) for l in range(2, 5)}
    assert expected_custom_lapws <= custom_lapws, "Generated species file is missing custom LAPWs"

    lo_signatures = [
        (
            int(lo.get("l")),
            tuple(sorted((int(wf.get("n")), int(wf.get("matchingOrder"))) for wf in lo.findall("wf"))),
        )
        for lo in species.findall("./basis/lo")
    ]
    assert len(lo_signatures) == len(set(lo_signatures)), "Generated species file contains duplicate local orbitals"

    highest_n = {0: 9, 1: 9, 2: 10, 3: 10, 4: 11}
    for l, n in highest_n.items():
        assert (l, ((n, 0), (n, 1))) in lo_signatures, f"Missing highest-energy local orbital for l={l}"


if __name__ == "__main__":
    main()
