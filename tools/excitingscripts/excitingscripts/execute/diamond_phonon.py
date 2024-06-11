"""Execute phonon diamond calculations."""

import json
import os
import time
from argparse import ArgumentParser
from pathlib import Path
from typing import Union, Dict
from xml.etree import ElementTree

from excitingscripts.execute.single import run_exciting


def execute_diamond_phonon(
    work_dir: Union[Path, str] = "workdir", excitingroot: str = os.getenv("EXCITINGROOT")
) -> Dict[str, Union[str, Dict[float, Dict[str, float]]]]:
    """Executes a series of exciting diamond calculations to get phonons.

    :param work_dir: Working directory containing the input files
    :param excitingroot: root directory of exciting
    :return: phonon results from the calculations
    """
    work_dir = Path(work_dir)
    assert work_dir.exists(), f"{work_dir} does not exist"

    with open(work_dir / "INFO-diamond-phonon.json") as fid:
        phonon_info: dict = json.load(fid)
    x_calc = "X-phonon-calculation mode" in phonon_info

    run_directories = work_dir.glob("displ*")
    results = {}

    for run_dir in run_directories:
        dir_name = run_dir.name
        print(f"Now running exciting for run_dir: {dir_name}")

        start_time = time.time()
        run_exciting(run_dir.as_posix(), excitingroot)
        print(f" Finished exciting run. Elapsed time: {time.time() - start_time:.2} seconds.")

        # can't use excitingtools parser here because it doesn't capture forces
        info_xml = run_dir / "info.xml"
        info = ElementTree.fromstring(info_xml.read_text())
        scl_info = info.find("groundstate").find("scl")

        energy = float(scl_info.findall("iter")[-1].find("energies").get("totalEnergy"))

        force = scl_info.find("structure").find("species").findall("atom")[1].find("forces").find("totalforce")
        force = round(float(force.get("x")), 10)

        displacement = float(dir_name[6:])
        results[displacement] = {"energy": energy, "force": force}
        if x_calc and abs(displacement) > 2e-6:
            results[-displacement] = results[displacement]

    return {"exciting_version": excitingroot + "/bin/exciting_smp", "results": dict(sorted(results.items()))}


def main() -> None:
    parser = ArgumentParser(
        description="Executing a series of exciting diamond calculations to get phonons."
    )

    parser.add_argument(
        "--work-directory",
        "-w",
        type=str,
        default="workdir",
        dest="work_directory",
        help="Working directory",
    )

    args = parser.parse_args()

    results = execute_diamond_phonon(args.work_directory)

    with open(f"{args.work_directory}/phonon_results.json", "w") as f:
        json.dump(results, f, indent=4)


if __name__ == "__main__":
    main()
