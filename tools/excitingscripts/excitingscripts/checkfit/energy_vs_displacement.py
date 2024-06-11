"""Python script for extracting the derivatives at zero displacement of the energy-vs-displacement curves."""

import json
from pathlib import Path

from excitingscripts.checkfit.checkfit import quantity_specific_checkfit, arg_parser


def main() -> None:
    quantity = "energy"
    factor = 2 / 3

    with open("INFO-diamond-phonon.json") as fid:
        phonon_info: dict = json.load(fid)
    if "X-phonon-calculation mode" in phonon_info:
        factor = 1 / 4
        print("This is the X-phonon-calculation:")
    if Path("quantum-espresso").exists():
        factor *= 0.5

    args = arg_parser(quantity).parse_args()

    check_fit_func = quantity_specific_checkfit(quantity, factor, 2)
    check_fit_func(args.maximum_displacement_fit, args.order_of_derivative, args.atomic_mass)


if __name__ == "__main__":
    main()
