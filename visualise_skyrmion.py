import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import os
import pandas as pd

def main():
    # Directory containing the output files
    output_dir = "skyrmion_evolution"
    
    # Read simulation parameters
    with open(f"{output_dir}/info.txt", 'r') as f:
        lines = f.readlines()
        num_frames = int(lines[0].strip())
        mesh_params = lines[1].strip().split()
        num_cells = int(mesh_params[0])
        lattice_param = float(mesh_params[1])
    
    # Create figure for animation
    fig = plt.figure(figsize=(16, 9))
    fig.suptitle("Skyrmion Evolution in 3D", fontsize=16)
    
    # 2D views
    ax1 = fig.add_subplot(2, 2, 1)  # Top-left: xy-plane (Sx, Sy)
    ax2 = fig.add_subplot(2, 2, 2)  # Top-right: Sz heatmap
    
    # 3D view
    ax3 = fig.add_subplot(2, 1, 2, projection='3d')  # Bottom: 3D view
    
    # Set up colormaps
    cmap = cm.coolwarm
    cmap_z = cm.viridis
    norm = Normalize(vmin=-1, vmax=1)
    
    # Function to update the plot for each frame
    def update_plot(frame):
        print(f"Processing frame {frame}/{num_frames-1}")
        
        # Read data from CSV file
        filename = f"{output_dir}/frame_{frame:05d}.csv"
        data = pd.read_csv(filename)
        
        # Extract coordinates and spin components
        x = data['x'].values
        y = data['y'].values
        z = data['z'].values
        Sx = data['Sx'].values
        Sy = data['Sy'].values
        Sz = data['Sz'].values
        
        # Clear axes
        ax1.clear()
        ax2.clear()
        ax3.clear()
        
        # Select a z-slice for 2D views
        unique_z = np.unique(z)
        middle_z = unique_z[len(unique_z) // 2]
        z_slice_mask = np.isclose(z, middle_z, atol=0.1*lattice_param)
        
        x_slice = x[z_slice_mask]
        y_slice = y[z_slice_mask]
        Sx_slice = Sx[z_slice_mask]
        Sy_slice = Sy[z_slice_mask]
        Sz_slice = Sz[z_slice_mask]
        
        # 1. Top-left: In-plane components (Sx, Sy) as arrows
        ax1.quiver(x_slice, y_slice, Sx_slice, Sy_slice, scale=15, 
                  pivot='mid', color='black', angles='xy')
        ax1.set_title(f"In-plane spin components at z={middle_z:.1f}")
        ax1.set_xlabel("x")
        ax1.set_ylabel("y")
        ax1.set_aspect('equal')
        ax1.set_xlim(0, num_cells * lattice_param)
        ax1.set_ylim(0, num_cells * lattice_param)
        
        # 2. Top-right: Out-of-plane component (Sz) as a colormap
        sc = ax2.scatter(x_slice, y_slice, c=Sz_slice, cmap=cmap_z, 
                        norm=norm, s=50)
        ax2.set_title(f"Sz component at z={middle_z:.1f}")
        ax2.set_xlabel("x")
        ax2.set_ylabel("y")
        ax2.set_aspect('equal')
        ax2.set_xlim(0, num_cells * lattice_param)
        ax2.set_ylim(0, num_cells * lattice_param)
        
        # Add colorbar if it's the first frame
        if frame == 0 and not hasattr(fig, 'cbar'):
            fig.cbar = plt.colorbar(sc, ax=ax2)
            fig.cbar.set_label("Sz")
        
        # 3. Full 3D visualization of spins
        # Down-sample for cleaner visualization in 3D
        sample_rate = 4  # Adjust this value based on your mesh density
        sample_mask = np.zeros(len(x), dtype=bool)
        
        # Create a grid of indices to sample
        unique_x = np.unique(x)
        unique_y = np.unique(y)
        unique_z = np.unique(z)
        
        # Every sample_rate points in each direction
        x_indices = np.arange(0, len(unique_x), sample_rate)
        y_indices = np.arange(0, len(unique_y), sample_rate)
        z_indices = np.arange(0, len(unique_z), sample_rate)
        
        # Combine to get 3D indices
        for i in x_indices:
            x_val = unique_x[i]
            for j in y_indices:
                y_val = unique_y[j]
                for k in z_indices:
                    z_val = unique_z[k]
                    # Find atoms close to this grid point
                    idx = np.where(
                        (np.abs(x - x_val) < 0.1*lattice_param) & 
                        (np.abs(y - y_val) < 0.1*lattice_param) & 
                        (np.abs(z - z_val) < 0.1*lattice_param)
                    )[0]
                    if len(idx) > 0:
                        sample_mask[idx[0]] = True
        
        # Get sampled data
        x_sampled = x[sample_mask]
        y_sampled = y[sample_mask]
        z_sampled = z[sample_mask]
        Sx_sampled = Sx[sample_mask]
        Sy_sampled = Sy[sample_mask]
        Sz_sampled = Sz[sample_mask]
        
        # Plot 3D arrows for spins
        ax3.quiver(x_sampled, y_sampled, z_sampled, 
                  Sx_sampled, Sy_sampled, Sz_sampled,
                  color=plt.cm.viridis(norm(Sz_sampled)), 
                  arrow_length_ratio=0.2, 
                  normalize=True)
        
        # Add colored sphere markers for atom positions, colored by Sz
        ax3.scatter(x_sampled, y_sampled, z_sampled, 
                   c=Sz_sampled, cmap=cmap_z, norm=norm, s=20)
        
        ax3.set_title("3D Visualization of Skyrmion")
        ax3.set_xlabel("x")
        ax3.set_ylabel("y")
        ax3.set_zlabel("z")
        ax3.set_xlim(0, num_cells * lattice_param)
        ax3.set_ylim(0, num_cells * lattice_param)
        ax3.set_zlim(0, num_cells * lattice_param)
        
        # Rotate the view for better visualization
        ax3.view_init(elev=30, azim=frame % 360)  # Rotate the view
        
        # Update frame title
        time_value = frame * (total_time / (num_frames - 1))
        fig.suptitle(f"Skyrmion Evolution in 3D (t = {time_value:.2f})", fontsize=16)
        
        return ax1, ax2, ax3

    # Simulation parameters for time display
    total_time = 10.0  # Should match the value in the Fortran program
    
    # Create animation
    ani = animation.FuncAnimation(fig, update_plot, frames=range(num_frames), 
                                 interval=200, blit=False)
    
    # Save animation
    ani.save('skyrmion_evolution.mp4', writer='ffmpeg', fps=10, dpi=200)
    
    print("Animation saved as skyrmion_evolution.mp4")

if __name__ == "__main__":
    main()
