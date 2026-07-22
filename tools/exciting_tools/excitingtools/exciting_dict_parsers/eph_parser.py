"""Parsers for exciting electron-phonon files."""

import re
from os import walk
from pathlib import Path
from typing import Union

import numpy as np

path_type = Union[Path, str]


def parse_evalqp_out(filename: path_type) -> dict:
    """Parse the quasi-particle energies and self-energy corrections from an electron-phonon calculation.

    :param filename: Path to EVALQP_T####.dat file.
    :return: Dictionary containing quasi-particle energies and self-energy corrections per k-point.
    """

    header = re.compile(r"k-point #\s*(\d+):\s*([-\d.Ee+]+)\s*([-\d.Ee+]+)\s*([-\d.Ee+]+)\s*([-\d.Ee+]+)")
    lines = Path(filename).read_text().splitlines()

    eval_qp = {"points": {}}

    i = 0
    while i < len(lines):
        m = header.match(lines[i])
        words = lines[i].split()
        if len(words) == 0:
            i += 1
            continue
        if m:
            ip = int(m.groups()[0])
            eval_qp["points"][ip] = {"vector": list(map(float, m.groups()[1:4]))}
            keys = lines[i + 1].split()
            for key in keys:
                eval_qp["points"][ip][key] = []
            i += 1
        elif words[0].isdigit():
            eval_qp["points"][ip][keys[0]].append(int(words[0]))
            for key, word in zip(keys[1:], words[1:]):
                eval_qp["points"][ip][key].append(float(word))
        i += 1

    return eval_qp


def parse_bandgap_renormalization(directories: list[path_type], prefixes: list[str] = ["EVALQP"]) -> dict:
    """Parse the quasi-particle energies of band egdes for all files in a directory.

    :param directory: Path to directory that contains EVALQP_T####.dat files.
    :return: Dictionary containing quasi-particle energies of band edges per k-point and temperature.
    """

    # extract all matching files from directory
    files = []
    for directory in directories:
        files_ = [f for f in next(walk(directory), (None, None, []))[2] if f.startswith((*prefixes,))]
        files += [str(Path(directory) / f) for f in files_]

    # find temperatures and bring files in correct order
    temps = [re.search(r"T\d{4}", f) for f in files]
    temps = np.array([int(t[0][1:]) for t in temps if t])
    idx = np.argsort(temps)
    files = [files[i] for i in idx]
    temps = [int(temps[i]) for i in idx]

    # read files
    result = {}
    for file, temp in zip(files, temps):
        data = parse_evalqp_out(file)
        for point in data["points"].values():
            key = ", ".join([f"{x:20.12g}" for x in point["vector"]])
            (occupied,) = np.where(np.array(point["E_KS[Ha]"]) < 0)
            ivb = -1
            if len(occupied):
                ivb = point["state"][occupied[-1]]
                if key not in result:
                    result[key] = {"vector": point["vector"], "temperatures": {}}
                if temp not in result[key]["temperatures"]:
                    result[key]["temperatures"][temp] = {}
                result[key]["valence band index"] = ivb
                result[key]["valence band KS energy"] = point["E_KS[Ha]"][point["state"].index(ivb)]
                result[key]["temperatures"][temp]["valence band QP energy"] = (
                    point["Re(E_EPH)[Ha]"][point["state"].index(ivb)]
                    + 1j * point["Im(E_EPH)[Ha]"][point["state"].index(ivb)]
                )
            if ivb + 1 in point["state"]:
                icb = ivb + 1
                if key not in result:
                    result[key] = {"vector": point["vector"], "temperatures": {}}
                if temp not in result[key]["temperatures"]:
                    result[key]["temperatures"][temp] = {}
                result[key]["conduction band index"] = icb
                result[key]["conduction band KS energy"] = point["E_KS[Ha]"][point["state"].index(icb)]
                result[key]["temperatures"][temp]["conduction band QP energy"] = (
                    point["Re(E_EPH)[Ha]"][point["state"].index(icb)]
                    + 1j * point["Im(E_EPH)[Ha]"][point["state"].index(icb)]
                )

    # restructure results
    keys = list(result.keys())
    result["points"] = list(result.values())
    for key in keys:
        result.pop(key)

    return result
