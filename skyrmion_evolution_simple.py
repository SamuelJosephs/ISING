
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm
import matplotlib.animation as animation
import pandas as pd

def save_heatmap(input_dir, output_dir, output_prefix, lattice_param, frame):
    # Read data from the specified CSV file in the input directory
    filename = f"{input_dir}/frame_{frame:05d}.csv"
    data = pd.read_csv(filename)
    
    # Extract coordinates and Sz components
    x = data['x'].values
    y = data['y'].values
    z = data['z'].values
    Sz = data['Sz'].values
    
    # Select the middle z-slice for the heatmap
    unique_z = np.unique(z)
    middle_z = unique_z[len(unique_z) // 2]
    z_slice_mask = np.isclose(z, middle_z, atol=0.1 * lattice_param)
    
    x_slice = x[z_slice_mask]
    y_slice = y[z_slice_mask]
    Sz_slice = Sz[z_slice_mask]
    
    # Set up colormap and normalize based on min and max Sz values in this slice
    cmap_z = cm.viridis
    norm_heat = Normalize(vmin=Sz_slice.min(), vmax=Sz_slice.max())
    
    # Create a new figure for the heatmap
    fig_heatmap, ax_heatmap = plt.subplots(figsize=(8, 8))
    
    # Try to reshape the data into a regular grid if possible
    x_unique = np.sort(np.unique(x_slice))
    y_unique = np.sort(np.unique(y_slice))
    
    if len(x_unique) * len(y_unique) == len(x_slice):
        # Data is on a regular grid. Build the heatmap array.
        heatmap = np.zeros((len(y_unique), len(x_unique)))
        for i, y_val in enumerate(y_unique):
            for j, x_val in enumerate(x_unique):
                idx = np.where((np.isclose(x_slice, x_val, atol=1e-8)) & 
                               (np.isclose(y_slice, y_val, atol=1e-8)))[0]
                if idx.size > 0:
                    heatmap[i, j] = Sz_slice[idx[0]]
        cax = ax_heatmap.imshow(
            heatmap,
            extent=[x_unique[0], x_unique[-1], y_unique[0], y_unique[-1]],
            origin='lower', cmap=cmap_z, norm=norm_heat
        )
    else:
        # If data is not on a regular grid, fall back to a scatter plot
        cax = ax_heatmap.scatter(x_slice, y_slice, c=Sz_slice, cmap=cmap_z, norm=norm_heat)
    
    ax_heatmap.set_xlabel(r"x $(\AA)$")
    ax_heatmap.set_ylabel(r"y $(\AA)$")
    fig_heatmap.colorbar(cax, ax=ax_heatmap, label="Sz")
    
    # Save the heatmap figure as PNG and PDF
    heatmap_png = f"{output_dir}/{output_prefix}_sz_heatmap.png"
    heatmap_pdf = f"{output_dir}/{output_prefix}_sz_heatmap.pdf"
    fig_heatmap.savefig(heatmap_png, dpi=300)
    fig_heatmap.savefig(heatmap_pdf, dpi=300)
    print(f"Heatmap saved as {heatmap_png} and {heatmap_pdf}")
    plt.close(fig_heatmap)

def main():
    # Check for proper command line arguments
    if len(sys.argv) < 3:
        print("Usage: python script.py <input_dir> <output_dir> [output_prefix]")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    output_prefix = sys.argv[3] if len(sys.argv) > 3 else "skyrmion_evolution"

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Read simulation parameters from info.txt in the input directory
    with open(f"{input_dir}/info.txt", 'r') as f:
        lines = f.readlines()
        num_frames = int(lines[0].strip())
        mesh_params = lines[1].strip().split()
        num_cells = int(mesh_params[0])
        lattice_param = float(mesh_params[1])
    
    # Create figure for the cross-sectional view (quiver plot)
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Set up colormap for Sz (z component) and a ScalarMappable for the colorbar
    cmap_z = cm.viridis
    norm = Normalize(vmin=-1, vmax=1)  # Dummy initialization; will update per frame
    sm = cm.ScalarMappable(cmap=cmap_z, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label("Sz")
    
    # Simulation parameters for time display
    total_time = 10.0  # Should match the value in the simulation program
    
    def update_plot(frame):
        # Read data from CSV file in the input directory
        filename = f"{input_dir}/frame_{frame:05d}.csv"
        data = pd.read_csv(filename)
        
        # Extract coordinates and spin components
        x = data['x'].values
        y = data['y'].values
        z = data['z'].values
        Sx = data['Sx'].values
        Sy = data['Sy'].values
        Sz = data['Sz'].values
        
        # Select the middle z-slice for the cross-sectional view
        unique_z = np.unique(z)
        middle_z = unique_z[len(unique_z) // 2]
        z_slice_mask = np.isclose(z, middle_z, atol=0.1 * lattice_param)
        
        x_slice = x[z_slice_mask]
        y_slice = y[z_slice_mask]
        Sx_slice = Sx[z_slice_mask]
        Sy_slice = Sy[z_slice_mask]
        Sz_slice = Sz[z_slice_mask]
        
        # Update normalization based on current frame's Sz min and max values
        norm = Normalize(vmin=Sz_slice.min(), vmax=Sz_slice.max())
        sm.set_norm(norm)
        cbar.update_normal(sm)
        
        # Determine arrow colors based on the updated normalization
        arrow_colors = cmap_z(norm(Sz_slice))
        
        # Clear the axis and plot the in-plane spin components as arrows
        ax.clear()
        ax.quiver(
            x_slice, y_slice, Sx_slice, Sy_slice,
            color=arrow_colors, scale=15, pivot='mid', angles='xy'
        )
        
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_aspect('equal')
        ax.set_xlim(0, num_cells * lattice_param)
        ax.set_ylim(0, num_cells * lattice_param)
        
        # Optionally, show time in title
        time_value = frame * (total_time / (num_frames - 1))
        ax.set_title(f"t = {time_value:.2f}")
        
        return ax,
    
    # Create the animation for the cross-sectional view
    ani = animation.FuncAnimation(
        fig, update_plot, frames=range(num_frames),
        interval=200, blit=False
    )
    
    # Save the animation to an MP4 file using the specified output prefix
    animation_mp4 = f"{output_dir}/{output_prefix}_cross_section.mp4"
    ani.save(animation_mp4, writer='ffmpeg', fps=10, dpi=200)
    print(f"Animation saved as {animation_mp4}")
    
    # --- NEW: save the *final* quiver frame as a PDF ---
    last_frame = num_frames - 1
    update_plot(last_frame)  # draws onto `fig`/`ax`
    quiver_pdf = f"{output_dir}/{output_prefix}_cross_section_final.pdf"
    fig.savefig(quiver_pdf, dpi=300)
    print(f"Final cross‑section frame saved as {quiver_pdf}")
    plt.close(fig)
    
    # Generate and save the heatmap for the Sz components using the final frame
    save_heatmap(input_dir, output_dir, output_prefix, lattice_param, frame=last_frame)

if __name__ == "__main__":
    main()
