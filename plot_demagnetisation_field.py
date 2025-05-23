
#!/usr/bin/env python3
"""
plot_demag_cross_section.py

Read a CSV file with columns x,y,z,fx,fy,fz and plot a 2D cross‑section
halfway through the Z plane, color‑coded by the fz component, with all arrows
normalized to unit length.
"""
import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm

def plot_cross_section(
    csv_file: str,
    output_file: str = None,
    atol_fraction: float = 0.1,
    arrow_scale: float = 15.0
):
    # --- 1) Load and verify ---
    df = pd.read_csv(csv_file)
    required = {'x','y','z','fx','fy','fz'}
    if not required.issubset(df.columns):
        missing = required - set(df.columns)
        raise ValueError(f"Missing columns in CSV: {missing}")

    x  = df['x'].values
    y  = df['y'].values
    z  = df['z'].values
    fx = df['fx'].values
    fy = df['fy'].values
    fz = df['fz'].values

    # --- 2) pick mid‑plane in Z ---
    unique_z = np.unique(z)
    mid_idx  = len(unique_z) // 2
    z0       = unique_z[mid_idx]
    # thickness = fraction of grid spacing
    if len(unique_z) > 1:
        dz = abs(unique_z[1] - unique_z[0]) * atol_fraction
    else:
        dz = atol_fraction

    mask = np.isclose(z, z0, atol=dz)
    x_slice  = x[mask]
    y_slice  = y[mask]
    fx_slice = fx[mask]
    fy_slice = fy[mask]
    fz_slice = fz[mask]

    # --- 3) normalize in‑plane vectors to unit length ---
    lengths = np.hypot(fx_slice, fy_slice)
    # avoid division by zero
    lengths[lengths == 0] = 1.0
    u = fx_slice / (2*lengths)
    v = fy_slice / (2*lengths)

    # --- 4) plot quiver with color = fz ---
    norm = Normalize(vmin=fz_slice.min(), vmax=fz_slice.max())
    cmap = cm.viridis

    fig, ax = plt.subplots(figsize=(8,8))
    q = ax.quiver(
        x_slice, y_slice,
        u, v,
        fz_slice,
        cmap=cmap, norm=norm,
        scale=arrow_scale,
        pivot='mid',
        angles='xy',
        width=0.005
    )
    cbar = fig.colorbar(q, ax=ax)
    cbar.set_label(r'$S_z$')

    ax.set_xlabel(r'x $(\AA)$')
    ax.set_ylabel(r'y $(\AA)$')
    ax.set_aspect('equal', 'box')
    plt.tight_layout()

    # --- 5) save or show ---
    if output_file:
        fig.savefig(output_file, dpi=300)
        print(f"Saved figure to {output_file}")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Plot normalized‐length demag field cross‑section from CSV"
    )
    parser.add_argument('csv_file',
                        help="Path to input CSV (must have x,y,z,fx,fy,fz)")
    parser.add_argument('-o','--output',
                        help="Output image file (PNG, PDF, etc.). If omitted, shows interactively.")
    parser.add_argument('--atol-fraction', type=float, default=0.1,
                        help="Fraction of z‑spacing for slice thickness (default: 0.1)")
    parser.add_argument('--arrow-scale', type=float, default=15.0,
                        help="Quiver ‘scale’ parameter to size arrows (default: 15)")

    args = parser.parse_args()
    try:
        plot_cross_section(
            args.csv_file,
            args.output,
            args.atol_fraction,
            args.arrow_scale
        )
    except Exception as e:
        print("Error:", e, file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
