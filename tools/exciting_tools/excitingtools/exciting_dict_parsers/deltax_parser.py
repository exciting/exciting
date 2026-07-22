"""Parser for DELTAX.OUT file"""

from pathlib import Path
from typing import Union

path_type = Union[Path, str]


def parse_deltax_out(name: path_type) -> dict:
    """Parser for: DELTAX.OUT

    :param name: the path of the file
    :return: a dictionary containing all the derivative discontinuities
    """
    lines = Path(name).read_text().splitlines()

    # first read the number of k-points and number of states
    n_k_points, n_states = map(int, lines[0].strip().split())

    # then we can check if the number of lines is correct
    assert (n_states + 1) * n_k_points + 1 == len(lines)

    data = {"n_k_points": n_k_points, "n_states": n_states, "derivative discontinuities": []}
    for i in range(n_k_points):
        k_point = list(map(float, lines[(n_states + 1) * i + 1].strip().split()))
        discontinuities = [
            float(line.strip()) for line in lines[(n_states + 1) * i + 2 : (n_states + 1) * i + 2 + n_states]
        ]
        data["derivative discontinuities"].append({"k_point": k_point, "discontinuities": discontinuities})

    return data
