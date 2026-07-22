"""Generate species files for excited-state calculations.

This module contains tools to extend *exciting* ground-state species files with the additional
basis functions (local orbitals (LOs) and custom LAPWs) required for excited-state calculations.

Purpose
-------------------------
A ground-state species file only contains the basis functions (LOs, custom LAPWs) needed to
converge total-energy/ground-state properties. Excited-state calculations (e.g. BSE, GW)
require a much richer basis, extending both to higher principal quantum
numbers 'n' within already-present angular momentum channels 'l', as well as to entirely new,
higher 'l'-channels that are not part of the ground-state basis at all.

This tool cross-references an already existing ground-state `<species>.xml` file with a
pre-calculated `LO_RECOMMENDATION.OUT` file (containing Wigner-Seitz recommended linearization
energies for every (l, n) combination) and, given a user-supplied energy threshold [Ha], adds all
recommended basis functions whose trial energy lies at or below that threshold:

- For an angular momentum channel 'l' that is not yet part of the species file's valence basis:
    * one custom LAPW at the lowest possible principal quantum number (n = l + 1)
    * one local orbital with matchingOrders (0, 0) bridging each pair of adjacent required
      'n' states (skipped for the very first, lowest state, since there is no lower state to
      bridge from)
    * local orbitals with matchingOrders [(0, 1), (1, 2), (2, 3)] up to a given maximum
      matching order for every required (l, n) combination
- For an angular momentum channel 'l' that is already part of the species file's valence basis:
    * the same (0, 0) bridging LOs and higher matching order LOs, continuing on from the highest
      'n' already present in the file

Local orbitals (and custom LAPWs) that are already present in the species file are never
duplicated.

Command Line Interface (CLI) Options
------------------------------------
The script is executed via the terminal and accepts several arguments to customize the generation
process.

**Basic Usage:**

python3 -m excitingscripts.setup.excited_state_species_files [SPECIES...] --energy ENERGY [OPTIONS]

Available Arguments & Options:
- **Species**
	- **Type:** String
	- **Description:** One or more species names to optimize (e.g. O, Ag, Rh)

- **-e, --energy**
	- **Type:** Float
	- **Required:** Yes
	- **Description:** Energy threshold [Ha]. Recommended basis functions with a trial energy at or
	  below this threshold are added to the species file.

- **-mo, --max-matching-order**
	- **Type:** Integer
	- **Default Value:** 1
	- **Accepted Values:** 0, 1, 2, or 3
	- **Description:** Highest matching order used for the additional local orbitals.

- **--path-xml**
	- **Type:** String
	- **Default Value:** `./ground_state_species/`
	- **Description:** Root path where the script looks for the existing ground-state species file
	  `<SPECIES>.xml`

- **--path-lo**
	- **Type:** String
	- **Default Value:** `./lo_recommendations/`
	- **Description:** Root path where the script looks for the `LO_RECOMMENDATION.OUT` file

- **-o, --output**
	- **Type:** String
	- **Default Value:** `./excited_state_species/`
	- **Description:** Target directory path where the extended output files `<SPECIES>_excited.xml`
	  will be written

- **--search-e**
	- **Type:** Flag
	- **Default Value:** False
	- **Description:** If set, newly added basis functions have `searchE="true"` instead of `"false"`

- **--keep-lin-dep**
	- **Type:** Flag
	- **Default Value:** False
	- **Description:** If set, disables the automatic skipping of linearly dependent high-energy
	  local orbitals (see `SpeciesFile.get_first_helo_n`)
"""

from __future__ import annotations

import os
import re
from argparse import ArgumentParser
from pathlib import Path
from typing import Dict, List, Set, Tuple

from excitingtools.species.species_file import SpeciesFile

# Directory where this script is located
SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR.parent.parent / "data"

DEFAULT_XML_PATH = DATA_DIR / "ground_state_species"
DEFAULT_LO_PATH = DATA_DIR / "lo_recommendations"
DEFAULT_OUTPUT_PATH = DATA_DIR / "excited_state_species"

LO_RECOMMENDATION_FILENAME = "LO_RECOMMENDATION.OUT"
MAX_MATCHING_ORDER = 3

# Matches header lines of the form "# species: Si, l :  0"
_SPECIES_HEADER_RE = re.compile(r"species:\s*(?P<species>\S+?)\s*,\s*l\s*:\s*(?P<l>-?\d+)")


def parse_lo_recommendation(file_path: str | Path, species_name: str) -> Dict[int, List[Tuple[int, int, float]]]:
    """Read the LO_RECOMMENDATION.OUT file and extract, per angular momentum channel 'l', the
    recommended (node, n, trial_energy) triplets for a single species.

    The file may contain blocks for multiple species (see the `n_species` header); only blocks
    whose header matches `species_name` are extracted.

    :param file_path: Path to the LO_RECOMMENDATION.OUT file.
    :param species_name: Chemical symbol of the species to extract (e.g. 'Si').
    :return: Dictionary mapping l-channel to a list of (node, n, trial_energy) tuples.
    """
    recommendations: Dict[int, List[Tuple[int, int, float]]] = {}
    current_l = None
    is_relevant_block = False

    content = Path(file_path).read_text(encoding="utf-8")
    for raw_line in content.splitlines():
        line = raw_line.strip()

        if not line:
            continue

        if line.startswith("#"):
            match = _SPECIES_HEADER_RE.search(line)
            if match:
                current_l = int(match.group("l"))
                is_relevant_block = match.group("species") == species_name
                if is_relevant_block:
                    recommendations.setdefault(current_l, [])
            continue

        if not is_relevant_block or current_l is None:
            continue

        parts = line.split()
        if len(parts) != 3:
            continue

        try:
            node, n, energy = int(parts[0]), int(parts[1]), float(parts[2])
        except ValueError:
            continue

        recommendations[current_l].append((node, n, energy))

    return recommendations


def get_existing_l_channels(species_obj: SpeciesFile) -> Set[int]:
    """Collect the angular momentum channels 'l' already represented in the valence basis of a
    species file, considering non-core atomic states, custom basis functions, and local orbitals.

    :param species_obj: The parsed SpeciesFile object.
    :return: Set of l-channels already present in the species file's valence basis.
    """
    l_channels: Set[int] = set()

    for state in species_obj.atomic_states:
        if not state.get("core", False):
            l_channels.add(state["l"])

    for custom in species_obj.basis["custom"]:
        l_channels.add(custom["l"])

    for lo in species_obj.basis["lo"]:
        l_channels.add(lo["l"])

    return l_channels


def get_first_required_n(species_obj: SpeciesFile, l_channel: int, skip_lin_dep: bool = True) -> int:
    """Return the first principal quantum number not represented in the species basis.

    ``SpeciesFile.get_first_helo_n`` accounts for atomic states and local orbitals. Ordinary custom
    APW and LAPW basis functions do not themselves provide a high-energy local orbital and must not
    advance this number.

    :param species_obj: The parsed SpeciesFile object.
    :param l_channel: Angular momentum quantum number.
    :param skip_lin_dep: Forwarded to :meth:`SpeciesFile.get_first_helo_n`.
    :return: First principal quantum number for which a basis function should be added.
    """
    return species_obj.get_first_helo_n(l_channel, skip_lin_dep)


def determine_required_states(
    species_obj: SpeciesFile,
    recommendations: Dict[int, List[Tuple[int, int, float]]],
    energy_threshold: float,
    skip_lin_dep: bool = True,
) -> Dict[int, int]:
    """Determine, for every l-channel present in the LO recommendation data, how many additional
    consecutive 'n' states (starting right after whatever is already present in the species file)
    have a recommended trial energy at or below `energy_threshold`.

    :param species_obj: The parsed SpeciesFile object of the existing ground-state species file.
    :param recommendations: Mapping l -> list of (node, n, trial_energy) tuples, as returned by
        :func:`parse_lo_recommendation`.
    :param energy_threshold: Energy threshold [Ha]; states with trial_energy <= threshold are added.
    :param skip_lin_dep: Forwarded to :meth:`SpeciesFile.get_first_helo_n`.
    :return: Mapping l -> number of additional consecutive 'n' states required. l-channels that do
        not require any additional state are omitted.
    """
    required: Dict[int, int] = {}

    for l, entries in recommendations.items():
        first_n = get_first_required_n(species_obj, l, skip_lin_dep)
        sorted_entries = sorted(entries, key=lambda entry: entry[1])

        count = 0
        exceeded_threshold = False
        for _node, n, energy in sorted_entries:
            if n < first_n:
                continue
            if energy > energy_threshold:
                exceeded_threshold = True
                break
            count += 1

        if not exceeded_threshold:
            print(
                f"\nWarning: l={l} never reached an entry with trial_energy > "
                f"{energy_threshold} Ha.\n"
                "LO recommendation data may not extend far enough to determine the true cutoff n.\n"
                "Increase nodesmaxlo and rerun the LO recommendation.\n"
            )

        if count > 0:
            required[l] = count

    return required


def add_states_for_l_channel(
    species_obj: SpeciesFile,
    l: int,
    n_high_n: int,
    is_new_l_channel: bool,
    *,
    max_matching_order: int = 1,
    skip_lin_dep: bool = True,
    search_e: bool = False,
) -> List[str]:
    """Add the basis functions required to cover `n_high_n` additional states of angular momentum
    channel `l` to a species file object, following the *exciting* excited-state basis recipe:

    - If `l` is not yet part of the species file's valence basis, add a custom LAPW at the lowest
      required principal quantum number.
    - Add a bridging local orbital (matchingOrder 0, 0) between each pair of adjacent required
      states, skipped for the very first state of a brand-new l-channel (there is no lower state
      to bridge from in that case).
    - Add up to `max_matching_order` additional local orbitals per required state, with
      matchingOrders (0, 1), (1, 2), ..., up to (`max_matching_order` - 1, `max_matching_order`).

    Local orbitals and custom LAPWs already present in the species file are never duplicated.

    :param species_obj: The SpeciesFile object to modify in place.
    :param l: Angular momentum quantum number.
    :param n_high_n: Number of additional consecutive 'n' states to add for this l-channel.
    :param is_new_l_channel: Whether `l` is not yet part of the species file's valence basis.
    :param max_matching_order: Highest matching order to add per state. 0 adds only the bridging
        (0, 0) LO; 1 additionally adds a (0, 1) LO; 2 additionally adds a (1, 2) LO; 3 additionally
        adds a (2, 3) LO.
    :param skip_lin_dep: Forwarded to :meth:`SpeciesFile.get_first_helo_n`.
    :param search_e: Value of `searchE` for all newly added basis functions.
    :return: List of human-readable log messages describing what was added or skipped.
    """
    if not 0 <= max_matching_order <= MAX_MATCHING_ORDER:
        raise ValueError(f"max_matching_order must be between 0 and {MAX_MATCHING_ORDER}")

    log: List[str] = []
    first_n = get_first_required_n(species_obj, l, skip_lin_dep)
    if is_new_l_channel:
        already_custom = any(
            custom["l"] == l and custom["type"] == "lapw"
            for custom in species_obj.basis["custom"]
        )
        if not already_custom:
            species_obj.basis["custom"].append({"l": l, "type": "lapw", "n": first_n, "searchE": search_e})
            log.append(f"  [ADDED CUSTOM] lapw l={l} n={first_n}")
        else:
            log.append(f"  [SKIPPED] custom lapw l={l} n={first_n} already present")

    existing_los = {
        (lo["l"], tuple(sorted(wf["n"] for wf in lo["wf"])), tuple(sorted(wf["matchingOrder"] for wf in lo["wf"])))
        for lo in species_obj.basis["lo"]
    }
    existing_los.update(
        (custom["l"], (custom["n"], custom["n"]), (0, 1))
        for custom in species_obj.basis["custom"]
        if custom["type"] == "apw+lo" and "n" in custom
    )

    for nr_lo in range(n_high_n):
        n = first_n + nr_lo

        # There is no physical state below the lowest principal quantum number n = l + 1, so no
        # bridging LO can be formed there. This also applies when a custom basis function already
        # represents the l-channel.
        if n != l + 1:
            bridge_ns = tuple(sorted((n - 1, n)))
            bridge_key = (l, bridge_ns, (0, 0))
            if bridge_key not in existing_los:
                species_obj.add_lo(l, (n - 1, n), (0, 0), search_e=search_e)
                existing_los.add(bridge_key)
                log.append(f"  [ADDED LO] l={l} n=({n - 1},{n}) mO=(0,0)")
            else:
                log.append(f"  [SKIPPED] LO l={l} n=({n - 1},{n}) mO=(0,0) already present")

        # High-energy LOs (matchingOrders (0,1), (1,2), ..., up to max_matching_order)
        for mo in range(max_matching_order):
            helo_key = (l, (n, n), (mo, mo + 1))
            if helo_key not in existing_los:
                species_obj.add_lo(l, (n, n), (mo, mo + 1), search_e=search_e)
                existing_los.add(helo_key)
                log.append(f"  [ADDED LO] l={l} n=({n},{n}) mO=({mo},{mo + 1})")
            else:
                log.append(f"  [SKIPPED] LO l={l} n=({n},{n}) mO=({mo},{mo + 1}) already present")

    return log


def validate_species_file(species_obj: SpeciesFile) -> None:
    """Ensure that the species file uses principal-quantum-number-based basis definitions.

    Species files using explicit ``trialEnergy`` values for local orbitals or custom LAPWs are not
    supported.

    :param species_obj: The SpeciesFile object to validate.
    """

    # Check custom LAPWs
    for custom in species_obj.basis["custom"]:
        if "trialEnergy" in custom:
            raise ValueError(
                "Ground-state species files that use explicit trialEnergy values for custom APWs "
                "cannot be processed by this script. Transform the trial energies into principal "
                "quantum numbers or add high-energy basis functions by hand."
            )

    # Check local orbitals
    for lo in species_obj.basis["lo"]:
        for wf in lo["wf"]:
            if "trialEnergy" in wf:
                raise ValueError(
                    "Ground-state species files that use explicit trialEnergy values for local "
                    "orbitals cannot be processed by this script. Transform the trial energies "
                    "into principal quantum numbers or add high-energy basis functions by hand."
                )


def optimize_species_for_excited_states(
    species_name: str,
    energy_threshold: float,
    path_xml: str | Path,
    path_lo: str | Path,
    path_out: str | Path,
    max_matching_order: int,
    *,
    skip_lin_dep: bool = True,
    search_e: bool = False,
) -> None:
    """Extend a ground-state species file with the excited-state basis functions (LOs and custom
    LAPWs) recommended by an LO_RECOMMENDATION.OUT file, up to a given energy threshold, and write
    the result to disk.

    :param species_name: Name of the chemical species (e.g. 'Si').
    :param energy_threshold: Energy threshold [Ha]. Recommended states with trial energy <=
        threshold are added.
    :param path_xml: Root directory containing the existing ground-state `<species>.xml` file.
    :param path_lo: Root directory containing the LO_RECOMMENDATION.OUT file.
    :param path_out: Target directory where the extended XML file will be saved.
    :param max_matching_order: Highest matching order to add per state (0-3). 0 adds only the
        bridging (0, 0) LO; 1 additionally adds a (0, 1) LO; 2 additionally adds a (1, 2) LO; 3
        additionally adds a (2, 3) LO.
    :param skip_lin_dep: Forwarded to :meth:`SpeciesFile.get_first_helo_n`.
    :param search_e: Value of `searchE` for all newly added basis functions.
    """
    if not 0 <= max_matching_order <= MAX_MATCHING_ORDER:
        raise ValueError(f"max_matching_order must be between 0 and {MAX_MATCHING_ORDER}")

    xml_path = Path(path_xml) / f"{species_name}.xml"
    lo_recommendation_path = Path(path_lo) / LO_RECOMMENDATION_FILENAME
    output_dir = Path(path_out)

    if not xml_path.exists():
        raise FileNotFoundError(f"Ground-state species file for {species_name} could not be found at: {xml_path}")

    if not lo_recommendation_path.exists():
        raise FileNotFoundError(f"LO recommendation file could not be found at: {lo_recommendation_path}")

    output_dir.mkdir(parents=True, exist_ok=True)

    species_obj = SpeciesFile.from_file(xml_path)
    validate_species_file(species_obj)
    recommendations = parse_lo_recommendation(lo_recommendation_path, species_name)

    if not recommendations:
        raise ValueError(
            f"No LO recommendation data found for species '{species_name}' in {lo_recommendation_path}"
        )

    existing_l_channels = get_existing_l_channels(species_obj)
    required_states = determine_required_states(species_obj, recommendations, energy_threshold, skip_lin_dep)

    if not required_states:
        print(f"  [INFO] No additional basis functions required below {energy_threshold} Ha.")

    added_count = 0
    for l in sorted(required_states):
        n_high_n = required_states[l]
        is_new_l_channel = l not in existing_l_channels
        print(f"  Processing l={l} (new l-channel: {is_new_l_channel}), {n_high_n} state(s) below threshold")

        log = add_states_for_l_channel(
            species_obj,
            l,
            n_high_n,
            is_new_l_channel,
            skip_lin_dep=skip_lin_dep,
            search_e=search_e,
            max_matching_order=max_matching_order,
        )
        for entry in log:
            print(entry)
            if "[ADDED" in entry:
                added_count += 1

    output_file_path = output_dir / f"{species_name}_excited.xml"
    species_obj.write(output_file_path)

    print("-" * 40)
    print(f"  [SUCCESS] {added_count} basis functions added. file: {output_file_path}")


def main() -> None:
    parser = ArgumentParser(
        description="Extend exciting ground-state species files with the local orbitals and custom "
        "LAPWs required for excited-state calculations."
    )

    parser.add_argument(
        "species",
        type=str,
        nargs="+",
        help="One or more species names to process (e.g., Si O Ag)",
    )

    parser.add_argument(
        "--energy",
        "-e",
        type=float,
        required=True,
        dest="energy",
        help="Energy threshold [Ha]. Recommended basis functions with a trial energy at or below "
        "this threshold are added.",
    )

    parser.add_argument(
        "--max-matching-order",
        "-mo",
        type=int,
        choices=range(MAX_MATCHING_ORDER + 1),
        default=1,
        dest="max_mo",
        help="Highest matching order used for additional local orbitals (default: 1).",
    )

    parser.add_argument(
        "--path-xml",
        type=str,
        default=str(DEFAULT_XML_PATH),
        dest="path_xml",
        help="Root path for the existing ground-state species .xml files",
    )

    parser.add_argument(
        "--path-lo",
        type=str,
        default=str(DEFAULT_LO_PATH),
        dest="path_lo",
        help="Root path for the LO_RECOMMENDATION.OUT file",
    )

    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default=str(DEFAULT_OUTPUT_PATH),
        dest="path_out",
        help="Output directory for the extended excited-state species files",
    )

    parser.add_argument(
        "--search-e",
        action="store_true",
        dest="search_e",
        help='If set, newly added basis functions have searchE="true" instead of "false"',
    )

    parser.add_argument(
        "--keep-lin-dep",
        action="store_true",
        dest="keep_lin_dep",
        help="If set, disables automatic skipping of linearly dependent high-energy local orbitals",
    )

    args = parser.parse_args()

    for element in args.species:
        print(f"\nProcessing species: {element}")
        optimize_species_for_excited_states(
            species_name=element,
            energy_threshold=args.energy,
            max_matching_order=args.max_mo,
            path_xml=args.path_xml,
            path_lo=args.path_lo,
            path_out=args.path_out,
            skip_lin_dep=not args.keep_lin_dep,
            search_e=args.search_e,
        )
        print("." * 60)


if __name__ == "__main__":
    main()
