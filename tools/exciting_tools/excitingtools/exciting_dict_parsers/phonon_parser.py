"""Parsers for exciting phonon files."""

import re
from pathlib import Path
from typing import Any, Dict, Union

import numpy as np

from excitingtools.parser_utils.parser_utils import numpy_gen_from_txt

path_type = Union[Path, str]


def parse_phonon_out(filename: path_type) -> dict:
    """Parse the phonon output file to extract phonon data for each q-point and mode.

    :param filename: Path to PHONON.OUT file.
    :return: Dictionary containing phonon data indexed by q-point.
    """

    phonon_data, current_q_point, current_mode = {}, {}, {}

    with open(filename) as f:
        lines = [line.strip() for line in f if line.strip()]

    for line in lines:
        if "q-point" in line:
            q_point_info = line.split(":")[0].split()
            q_point_index = int(q_point_info[0])
            q_vector = list(map(float, q_point_info[1:]))

            current_q_point = {"q_vector": q_vector, "modes": []}

            phonon_data[str(q_point_index)] = current_q_point

        elif "mode" in line:
            mode_info = line.split(":")[0].split()
            mode_index = int(mode_info[0])
            frequency = float(mode_info[1])

            current_mode = {"mode_index": str(mode_index), "frequency": frequency, "eigenvector_info": []}
            current_q_point["modes"].append(current_mode)

        else:
            eigenvector_info = line.split(":")[0].split() if ":" in line else line.split()
            species, atom, polarisation = map(int, eigenvector_info[:3])
            eigenvector_component_real, eigenvector_component_imag = map(float, eigenvector_info[-2:])

            current_mode["eigenvector_info"].append(
                {
                    "species": species,
                    "atom": atom,
                    "polarisation": polarisation,
                    "eigenvector_component_real": eigenvector_component_real,
                    "eigenvector_component_imag": eigenvector_component_imag,
                }
            )

    return phonon_data


def parse_dyn_out(filename: path_type) -> dict:
    """Parse the dynamical matrix output file.

    :param filename: Path to DYN_Q????_????_????.OUT file.
    :return: Dictionary containing dynamical matrix entry for each species, atom and polarization.
    """

    dyn_data = {}

    with open(filename) as f:
        lines = [line.strip() for line in f if line.strip()]

    i = 0
    for line in lines:
        row = line.split(":")
        if len(row) == 2:
            i += 1
            dyn_real, dyn_imag = map(float, row[0].split())
            species, atom, polarisation = map(int, [s.split("=")[1] for s in row[1].split(",")])
            dyn_data[str(i)] = {
                "species": species,
                "atom": atom,
                "polarisation": polarisation,
                "dynmat_real": dyn_real,
                "dynmat_imag": dyn_imag,
            }

    return dyn_data


def parse_epsinf_out(filename: path_type) -> dict:
    """Parse the high-frequency dielectric constant output file.

    :param filename: Path to EPSINF.OUT file.
    :return: Dictionary containing high frequency dielectric constant.
    """
    return {"epsinf": numpy_gen_from_txt(filename)}


def parse_zstar_out(filename: path_type) -> dict:
    """Parse the Born-effective charges output file.

    :param filename: Path to ZSTAR.OUT file.
    :return: Dictionary containing the Born-effective charge tensor for each atom.
    """
    header_re = re.compile(
        r"# species\s+(\d+)\s+atom\s+(\d+)\s+\((\w+)\s+(\d+)\)\s*:\s*([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)"
    )

    zstar_data: Dict[str, Any] = {"atoms": {}, "acoustic_sum_rule_correction": None}
    lines = Path(filename).read_text().splitlines()
    iat = 0
    i = 0

    while i < len(lines):
        line = lines[i]

        # Atom header
        m = header_re.match(line)
        if m:
            iat += 1
            species, atom, symbol, species_atom = map(lambda x: int(x) if x.isdigit() else x, m.groups()[:4])
            position = np.array(list(map(float, m.groups()[4:])), dtype=float)
            tensor = np.loadtxt(lines[i + 1 : i + 4])
            zstar_data["atoms"][str(iat)] = {
                "species": species,
                "atom": atom,
                "symbol": symbol,
                "species_atom": species_atom,
                "position": position,
                "tensor": tensor,
            }
            i += 4
            continue

        # Acoustic sum rule correction
        if line.startswith("# Acoustic sum rule correction"):
            zstar_data["acoustic_sum_rule_correction"] = np.loadtxt(lines[i + 1 : i + 4])
            break

        i += 1

    return zstar_data
