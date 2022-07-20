"""Parsers for GW's EPS00_GW.OUT output.
"""
import numpy as np

from excitingtools.parser_utils.parser_decorators import accept_file_name
from excitingtools.utils.utils import get_new_line_indices

# Output file name
_file_name = 'EPS00_GW.OUT'


def parse_eps00_frequencies(file_string: str) -> dict:
    """Parse frequencies from EPS00_GW.OUT.

    :param str file_string: Input string
    :return dict frequencies: Dictionary of frequencues {index: frequency}
    """
    initial_header_size = 3
    block_size = 6

    # Start file_lines on first block header
    file_lines = file_string.splitlines()
    file = file_lines[initial_header_size:]
    n_freq = int(len(file) / block_size)

    frequencies = {}
    j = 0
    for i in range(0, n_freq):
        j = i * block_size
        frequencies[i + 1] = float(file[j].split()[-1])

    assert len(frequencies) == n_freq

    return frequencies


def parse_eps00_blocks(file_string: str, n_freq: int) -> dict:
    """Parser for epsilon blocks in EPS00_GW.OUT.

    File of the form:

    (dielectric tensor, random phase approximation)

    frequency index and value:      1    0.01985507
    real part, imaginary part below
     152.69882148    0.00000000    0.00000000         0.00000000    0.00000000    0.00000000
       0.00000000  152.69882148    0.00000000         0.00000000    0.00000000    0.00000000
       0.00000000    0.00000000  152.69882148         0.00000000    0.00000000    0.00000000

    frequency index and value:      2    0.02771249
     real part, imaginary part below
        8.22189228    0.00000000    0.00000000        -0.00000000    0.00000000    0.00000000
        0.00000000    8.22189228    0.00000000         0.00000000   -0.00000000    0.00000000
        0.00000000    0.00000000    8.22189228         0.00000000    0.00000000   -0.00000000

    :param str file_string: Input string
    :param int n_freq: Number of frequency points
    :return dict data: Dict containing real and imaginary parts of the dielectric function,
     eps00, for each frequency point. Frequency indexing (keys) start at 1.
    """
    line = get_new_line_indices(file_string)
    assert file_string[line[0]:line[1]].isspace(), (
        "First line of EPS00_GW.OUT must be a whiteline")

    initial_header = 3
    offset = initial_header - 1
    # Header lines + eps00
    repeat_block_size = 6

    data = {}
    file_string = file_string.splitlines()[initial_header:]

    def extract_eps_i(line: str) -> tuple:
        eps_i_re = np.array(line.split()[0:3], dtype=float)
        eps_i_img = np.array(line.split()[3:6], dtype=float)
        return eps_i_re, eps_i_img

    for i_freq in range(0, n_freq):
        i = offset + i_freq * repeat_block_size

        eps_x_re, eps_x_img = extract_eps_i(file_string[i])
        eps_y_re, eps_y_img = extract_eps_i(file_string[i + 1])
        eps_z_re, eps_z_img = extract_eps_i(file_string[i + 2])

        data[i_freq + 1] = {
            're': np.array([eps_x_re, eps_y_re, eps_z_re]),
            'img': np.array([eps_x_img, eps_y_img, eps_z_img])
        }

    return data


@accept_file_name
def parse_eps00_gw(file_string: str) -> dict:
    """ Parser frequency grid and epsilon00 from EPS00_GW.OUT.

    :param str file_string: Input string
    :return dict data: Dict containing the frequency, and real and imaginary parts
     of the dielectric function, eps00, for each frequency point.
    """
    frequencies = parse_eps00_frequencies(file_string)
    eps00 = parse_eps00_blocks(file_string, len(frequencies))
    assert len(frequencies) == len(eps00), \
        "Expect eps00 to have n frequencies consistent with the frequency grid"

    # Repackage frequency points and eps00 together
    data = {}
    for i_freq in range(1, len(frequencies) + 1):
        data[i_freq] = {'frequency': frequencies[i_freq], 'eps00': eps00[i_freq]}

    return data
