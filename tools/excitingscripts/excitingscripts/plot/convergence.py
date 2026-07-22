"""Visualize convergence results.

Located at `excitingscripts/plot/convergence.py`.

Call as:

```bash
python3 -m excitingscripts.plot.convergence plot_mode
```
Where <code>plot_mode</code> is either <code>k</code> for plotting energy curves with varying values of the <code><span style="color:green">groundstate</span></code> attribute <code><span style="color:MediumBlue">ngridk</span></code>, <code>r</code> for varying values of the <code><span style="color:green">groundstate</span></code> attribute  <code><span style="color:MediumBlue">rgkmax</span></code> or <code>rk</code> for a 3D plot with varying values of both attributes.
"""

from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np


def plot_convergence_k() -> None:
    """Plot energy curves for varying values of the groundstate attribute ngridk.
    """
    convergence_data = np.loadtxt("./convergence-test")

    rgkmax = convergence_data[0, 1]
    k_values = convergence_data[:, 0]
    energy_values = convergence_data[:, 2]

    plt.xlabel("ngridk", fontsize=18)
    plt.ylabel("Energy [Ha]", fontsize=18)

    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tick_params(length=5, width=2)
    plt.rcParams["axes.linewidth"] = 4

    plt.plot(k_values, energy_values, color="red", linewidth=2)
    plt.scatter(k_values, energy_values, color="green", label=f"rgkmax = {rgkmax}", linewidth=3, zorder=2)

    plt.legend(prop={"size": 14})
    plt.grid(linestyle='--', linewidth=0.5)

def plot_convergence_r() -> None:
    """Plot energy curves for varying values of the groundstate attribute rgkmax.
    """
    convergence_data = np.loadtxt("./convergence-test")

    ngridk = convergence_data[0, 0]
    rgkmax_values = convergence_data[:, 1]
    energy_values = convergence_data[:, 2]

    plt.xlabel("rgkmax", fontsize=18)
    plt.ylabel("Energy [Ha]", fontsize=18)

    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tick_params(length=5, width=2)
    plt.rcParams["axes.linewidth"] = 4

    plt.plot(rgkmax_values, energy_values, color="red", linewidth=2)
    plt.scatter(rgkmax_values, energy_values, color="green", label=f"ngridk = {ngridk}", linewidth=3, zorder=2)

    plt.legend(prop={"size": 14})
    plt.grid(linestyle='--', linewidth=0.5)

def plot_convergence_rk() -> None:
    """Plot energy curves for varying values of the groundstate attributes ngridk and rgkmax.
    """
    convergence_data = np.loadtxt("./convergence-test")

    k_values = convergence_data[:, 0]
    rgkmax_values = convergence_data[:, 1]
    energy_values = convergence_data[:, 2]

    k_grids = list(sorted(np.unique(k_values).tolist()))
    rgkmax_range = list(sorted(np.unique(rgkmax_values).tolist()))

    grid_map = {
        k_value: {
            rgkmax_value: (k_idx, rgkmax_idx)
            for rgkmax_idx, rgkmax_value in enumerate(rgkmax_range)
        }
        for k_idx, k_value in enumerate(k_grids)
    }

    energy_values = np.empty((len(k_grids), len(rgkmax_range)))
    for k_val, rgkmax_val, energy in convergence_data:
        array_index = grid_map[k_val][rgkmax_val]
        energy_values[array_index[0], array_index[1]] = energy

    ax = plt.axes(projection='3d')

    X, Y = np.meshgrid(rgkmax_range, k_grids)
    ax.plot_wireframe(X, Y, energy_values, color='darkviolet')

    ax.set_xlabel('rgkmax', fontsize=12)
    ax.set_ylabel('ngridk', fontsize=12)
    ax.set_zlabel("Energy [Ha]", fontsize=12)
    ax.zaxis.labelpad = 18
    ax.zaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    ax.tick_params(axis='z', which='major', pad=10)

    ax.grid(False)

def main() -> None:
    parser = ArgumentParser(description="""Plot energy curves for varying values for the groundstate attributes ngridk
                                        and rgkmax.""")

    parser.add_argument("plot_mode",
                        type=str,
                        nargs=1,
                        help="convergence plot mode, expected values: k, r or rk")

    parser.add_argument("-sh", "--show",
                        action="store_true",
                        help="show plot")

    args = parser.parse_args()

    if args.plot_mode[0] == "k":
        plot_convergence_k()
    if args.plot_mode[0] == "r":
        plot_convergence_r()
    if args.plot_mode[0] == "rk":
        plot_convergence_rk()

    plt.tight_layout()
    plt.savefig('PLOT.png', orientation='portrait', format='png', dpi=300)

    if args.show:
        plt.show()

if __name__ == "__main__":
    main()
