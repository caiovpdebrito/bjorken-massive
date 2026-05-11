#!/usr/bin/env python3
from aesthetic import *

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_metadata(input_file):
    mode = None
    n_values = None
    l_values = None

    with open(input_file) as file:
        for line in file:
            if not line.startswith("#"):
                continue

            tokens = line[1:].split()

            if len(tokens) == 0:
                continue

            if tokens[0] == "mode":
                mode = tokens[1]

            if tokens[0] == "n_values":
                n_values = np.array([int(value) for value in tokens[1:]])

            if tokens[0] == "l_values":
                l_values = np.array([int(value) for value in tokens[1:]])

    return mode, n_values, l_values


def interpolate_curve(x, y, num_points=400):
    x_unique, indices = np.unique(x, return_index=True)
    y_unique = y[indices]

    if len(x_unique) < 2:
        return x_unique, y_unique

    try:
        from scipy.interpolate import PchipInterpolator

        interpolator = PchipInterpolator(x_unique, y_unique)
        x_dense = np.linspace(x_unique.min(), x_unique.max(), num_points)
        y_dense = interpolator(x_dense)

    except ImportError:
        x_dense = np.linspace(x_unique.min(), x_unique.max(), num_points)
        y_dense = np.interp(x_dense, x_unique, y_unique)

    return x_dense, y_dense


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/mass_scan_n.dat")
    parser.add_argument("--tau-label", default=None)
    args = parser.parse_args()

    input_file = Path(args.input)
    output_dir = Path("figures")
    output_file = output_dir / f"{input_file.stem}.pdf"

    mode, n_values, l_values = read_metadata(input_file)

    if mode is None or n_values is None or l_values is None:
        raise RuntimeError("ERROR: input file is missing metadata header.")

    data = np.loadtxt(input_file, comments="#")
    data = np.atleast_2d(data)

    z = data[:, 0]
    moments = data[:, 1:]

    order = np.argsort(z)
    z = z[order]
    moments = moments[order, :]

    if mode == "n":
        color_values = n_values
        colorbar_label = r"$n$"
        fixed_l = l_values[0]
        ylabel = rf"$|\varrho_{{n,{fixed_l}}}|/\varrho^{{\mathrm{{eq}}}}_{{n,0}}$"

    elif mode == "l":
        color_values = l_values
        colorbar_label = r"$l$"
        fixed_n = n_values[0]
        ylabel = rf"$|\varrho_{{{fixed_n},l}}|/\varrho^{{\mathrm{{eq}}}}_{{{fixed_n},0}}$"

    else:
        raise RuntimeError("ERROR: unknown scan mode.")

    if len(color_values) != moments.shape[1]:
        raise RuntimeError("ERROR: metadata does not match number of data columns.")

    cmap = plt.colormaps["coolwarm"]

    if color_values.min() == color_values.max():
        norm = plt.Normalize(vmin=color_values.min() - 0.5,
                             vmax=color_values.max() + 0.5)
    else:
        norm = plt.Normalize(vmin=color_values.min(), vmax=color_values.max())

    fig, ax = plt.subplots()

    for column, color_value in enumerate(color_values):
        color = cmap(norm(color_value))

        z_dense, moment_dense = interpolate_curve(z, moments[:, column])

        ax.plot(z_dense, moment_dense, color=color, linewidth=1.8)
        ax.scatter(z, moments[:, column], color=color, s=32, zorder=3)

    # Aesthetic options for plots
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label(colorbar_label)

    ax.set_xlabel(r"$m$ [GeV]")
    ax.set_ylabel(ylabel)

    ax.set_xlim(z[0], z[-1])

    ax.set_title(rf"$\tau = {args.tau_label}$ fm")

    fig.tight_layout()

    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file, bbox_inches="tight")

    print(f"Plot written to {output_file}")

    plt.show()


if __name__ == "__main__":
    main()