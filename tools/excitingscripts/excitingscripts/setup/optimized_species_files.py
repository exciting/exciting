"""Exciting Species File Optimizer.

This module contains tools and baseline data for optimizing *exciting* species files. 

The Purpose of Elemental Crystal Species Files
----------------------------------------------
In *exciting*, calculations require basis sets (including local orbitals (LOs)) tailored 
to each element. The files provided here were initially generated for elemental crystals 
to establish a "maximized," highly converged reference basis. 

Depending on your specific target precision, many local orbitals can be pruned to save time. 
This includes high-energy local orbitals that have negligible impact on ground-state properties, 
as well as standard local orbitals that can be safely removed if a lower, more lightweight 
calculation precision is sufficient for your application. 

The goal of this optimization tool is to cross-reference these maximized species files 
with pre-calculated benchmark calculations to dynamically prune non-essential orbitals 
based on their energy precision penalty, while strictly protecting those flagged in the 
benchmark results "Essential" whitelist to preserve the core physics. This yields a 
mathematically efficient, lightweight basis set that retains the required accuracy 
and reduces computational time.

Periodic Table Precision Map
----------------------------
The baseline convergence quality for all available elements is compiled in the benchmark 
dataset. The figure illustrating the initial test precision for each maximized species 
file can be found in the repository under:
`../_data/lo_hierarchies/PSE_DBSV_precision_map.pdf`

Command Line Interface (CLI) Options
------------------------------------
The script is executed via the terminal and accepts several arguments to customize the optimization process.

**Basic Usage:**

python3 `optimized_species_files.py` [SPECIES...] [OPTIONS]

Available Arguments & Options:
- **Species**
	- **Type:** String
	- **Description:** One or more species names to optimize (e.g. O, Ag, Rh)

- **-p, --precision**
	- **Type:** Float
	- **Default Value:** 1.0e-4
	- **Description:** Precision threshold for removing LOs

- **--path-lo**
	- **Type:** String
	- **Default Value:** `<package_install_dir>/_data/lo_hierarchies/`
	- **Description:** Root path where the script looks for the `lo_hierarchy_<SPECIES>.out`

- **--path-xml**
	- **Type:** String
	- **Default Value:** `<package_install_dir>/_data/max_species_files/`
	- **Description:** Root path where the script looks for the maximized species file `<SPECIES>.xml`

- **-o, --output**
	- **Type:** String
	- **Default Value:** `./customized_species/`
	- **Description:** Target directory path where the optimized output files `<SPECIES>_customized` will be written
"""

import os

from argparse import ArgumentParser
from typing import List, Tuple, Dict, Any
from pathlib import Path
from importlib.resources import files

from excitingtools.species.species_file import SpeciesFile

# Directory where this script is located
SCRIPT_DIR = Path(files("excitingscripts")).resolve()
DATA_DIR = SCRIPT_DIR / "_data"

DEFAULT_LO_PATH = DATA_DIR / "lo_hierarchies"
DEFAULT_XML_PATH = DATA_DIR / "max_species_files"
DEFAULT_OUTPUT_PATH = Path(os.getcwd()).resolve() / "customized_species"


def parse_lo_hierarchy(file_path: str | Path) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Read the lo_hierarchy_<species>.out file and extract non-essential and essential LOs.

    :param file_path: Path to the lo_hierarchy_<species>.out file.
    :return: A tuple containing two lists of dictionaries: (non_essential_los, essential_los).
    """
    non_essential = []
    essential_los = []
    mode = "non_essential"
    
    content = Path(file_path).read_text(encoding='utf-8')
    lines = content.splitlines()
    for raw_line in lines:
        line = raw_line.strip()
    
        if not line or line.startswith("species") or line.startswith("DBSV"):
            continue
    
        if "Essential local orbitals:" in line:
            mode = "essential"
            continue

        parts = line.split()
    
        if len(parts) < 4:
            continue

        l_val = int(parts[1])
        ns_val = tuple(sorted([int(x) for x in parts[2].split(',')]))
        mos_val = tuple(sorted([int(x) for x in parts[3].split(',')]))

        if mode == "non_essential" and len(parts) >= 5:
            prec_val = float(parts[4])
            non_essential.append({
                'precision': prec_val,
                'l_value': l_val,
                'n_values': ns_val,
                'matching_orders': mos_val
            })
        elif mode == "essential": 
            essential_los.append({
                'l_value': l_val,
                'n_values': ns_val,
                'matching_orders': mos_val
            })	

    return non_essential, essential_los


def remove_los(species_obj: SpeciesFile, l: int, ns: Tuple[int, ...], mOs: Tuple[int, ...]) -> bool:
    """Detect and remove specific local orbitals in the <species>.xml file object.

    :param species_obj: The parsed SpeciesFile object from excitingtools.
    :param l: Orbital quantum number l.
    :param ns: Tuple of principal quantum numbers.
    :param mOs: Tuple of matching orders.
    :return: True if the LO was found and removed, False otherwise.
    """
    target_ns = tuple(sorted(ns))
    target_mos = tuple(sorted(mOs))
    
    for lo in species_obj.basis["lo"]:
        if lo["l"] == l:
            current_ns = tuple(sorted([int(wf["n"]) for wf in lo["wf"] if "n" in wf]))
            current_mos = tuple(sorted([int(wf["matchingOrder"]) for wf in lo["wf"]]))
            
            if current_ns == target_ns and current_mos == target_mos:
                species_obj.basis["lo"].remove(lo)
                return True
    return False


def optimize_species_from_folder(species_name: str, precision_threshold: float, path_lo: str | Path, path_xml: str | Path, path_out: str | Path) -> None:
    """Optimize species files by filtering out non-essential LOs below a given precision threshold and writing the results to disk.

    :param species_name: Name of the chemical species (e.g. 'Yb').
    :param precision_threshold: Threshold for precision [Ha]. LOs with precision < precision_threshold are removed.
    :param path_lo: Root directory containing the lo_hierarchies output files.
    :param path_xml: Root directory containing the maximum species XML files.
    :param path_out: Target directory where the optimized XML file will be saved.
    """
    xml_path = Path(path_xml) / f"{species_name}.xml"
    hierarchy_path = Path(path_lo) / f"lo_hierarchy_{species_name}.out"
    output_dir = Path(path_out) 
    
    if not xml_path.exists():
        raise FileNotFoundError(f"Baseline XML file for {species_name} could not be found at: {xml_path}")

    if not hierarchy_path.exists():
        raise FileNotFoundError(f"Hierarchy log file (.out) for {species_name} could not be found at: {hierarchy_path}")

    output_dir.mkdir(parents=True, exist_ok=True)   

    non_essential, essential_los = parse_lo_hierarchy(hierarchy_path)
    species_obj = SpeciesFile.from_file(xml_path)

    protected: set[tuple[int, tuple[int, ...], tuple[int, ...]]] = set()
    for e in essential_los:
        protected.add((e['l_value'], e['n_values'], e['matching_orders']))

    removed_count = 0
    for cand in non_essential:
        prec = cand['precision']
        l = cand['l_value']
        ns = cand['n_values']
        mOs = cand['matching_orders']

        if prec < precision_threshold:
            if (l, ns, mOs) not in protected:
                if remove_los(species_obj, l, ns, mOs):
                    print(f"  [DELETED] l={l} ns={ns} mOs={mOs} (Precision: {prec:.2e})")
                    removed_count += 1
            else:
                print(f"  [PROTECTED] l={l} ns={ns} is essential!")

    output_file_path = output_dir / f"{species_name}_customized.xml"
    species_obj.write(output_file_path)

    print("-" * 40)
    print(f"  [SUCCESS] {removed_count} LOs removed. file: {output_file_path}")


def main() -> None:
    parser = ArgumentParser(description="Optimize exciting species files by removing non-essential local orbitals.")

    parser.add_argument("species",
                        type=str,
                        nargs="+",
                        help="One or more species names to optimize (e.g., Yb He Ag)")

    parser.add_argument("--precision", "-p",
                        type=float,
                        default=1.0e-4,
                        dest="precision",
                        help="Total energy precision threshold [Ha] for removing LOs (default: 1.0e-4 Ha)")

    parser.add_argument("--path-lo",
                        type=str,
                        default=str(DEFAULT_LO_PATH),
                        dest="path_lo",
                        help="Root path for lo_hierarchies .out files")

    parser.add_argument("--path-xml",
                        type=str,
                        default=str(DEFAULT_XML_PATH),
                        dest="path_xml",
                        help="Root path for maximized species .xml files")

    parser.add_argument("--output", "-o",
                        type=str,
                        default=str(DEFAULT_OUTPUT_PATH),
                        dest="path_out",
                        help="Output directory for customized species files")

    args = parser.parse_args()

    for element in args.species:
        print(f"\nProcessing species: {element}")
        optimize_species_from_folder(
            species_name=element,
            precision_threshold=args.precision,
            path_lo=args.path_lo,
            path_xml=args.path_xml,
            path_out=args.path_out )
        print("." * 60)


if __name__ == "__main__":
    main()
