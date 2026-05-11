#!/usr/bin/env python3
from aesthetic import *
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("l", type=int)
    parser.add_argument("--nmin", type=int, default=0)
    args = parser.parse_args()

    input_file = Path("data") / f"moments_l{args.l}.dat"
    output_dir = Path("figures")
    output_file = output_dir / f"moments_l{args.l}.pdf"

    data = np.loadtxt(input_file)

    x = data[:, 0]
    y = data[:, 1:]

    n_values = args.nmin + np.arange(y.shape[1])

    cmap = plt.colormaps["coolwarm"]
    norm = plt.Normalize(vmin=n_values.min(), vmax=n_values.max())

    fig, ax = plt.subplots()

    for column, n in zip(range(0, y.shape[1], 2), n_values[::2]):
        ax.loglog(x, y[:, column], color=cmap(norm(n)), linewidth=1.8)

    # Add colorbar
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label(r"$n$")

    ax.set_xlabel(r"$\tau/\tau_R$")
    ax.set_ylabel(rf'$|\varrho_{{n,{args.l}}}|/\varrho^{{\mathrm{{eq}}}}_{{n,0}}$')
    ax.set_title(rf"$l = {args.l}$")

    fig.tight_layout()

    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file, bbox_inches="tight")

    print(f"Plot written to {output_file}")

    plt.show()

###
##
#

if __name__ == "__main__":
    main()
