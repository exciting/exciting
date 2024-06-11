"""Python script for extracting the derivatives at zero displacement of the force-vs-displacement curves."""

from excitingscripts.checkfit.checkfit import quantity_specific_checkfit, arg_parser


def main() -> None:
    quantity = "force"
    factor = 2
    args = arg_parser(quantity).parse_args()

    check_fit_func = quantity_specific_checkfit(quantity, factor, 1)
    check_fit_func(args.maximum_displacement_fit, args.order_of_derivative, args.atomic_mass)


if __name__ == "__main__":
    main()
