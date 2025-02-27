"""Parsers for exciting phonon files."""

from pathlib import Path
from typing import Union


def parse_phonon_out(filename: Union[str, Path]) -> dict:
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
