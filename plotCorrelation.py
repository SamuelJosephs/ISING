import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

def load_correlation_data(num_files=60):
    """
    Load correlation data files.
    Returns a list of numpy arrays containing the data.
    """
    data_arrays = []
    for i in range(1, num_files + 1):
        filename = f'./correlation_results/correlation_{i}'
        try:
            data = np.loadtxt(filename)
            data_arrays.append(data)
        except Exception as e:
            print(f"Error loading file {filename}: {e}")
    
    return data_arrays

def create_correlation_gif(data_arrays, output_file='correlation_animation.gif', 
                         interval=500, fps=2):
    """
    Create a GIF animation of correlation plots.
    """
    # Set up the figure and axis
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Find global min and max values for consistent y-axis
    all_data = np.concatenate(data_arrays)
    y_min, y_max = all_data.min(), all_data.max()
    y_padding = (y_max - y_min) * 0.1
    y_min -= y_padding
    y_max += y_padding
    
    def animate(frame):
        ax.clear()
        
        # Get current data
        data = data_arrays[frame]
        x = np.arange(len(data))
        
        # Create the plot
        ax.plot(x, data, 'b-', lw=2)
        ax.set_ylim(y_min, y_max)
        
        # Add labels and title
        ax.set_xlabel('Index')
        ax.set_ylabel('Correlation Value')
        ax.set_title(f'Correlation Plot {frame + 1}/{len(data_arrays)}')
        
        # Add grid
        ax.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
    
    # Create the animation
    anim = FuncAnimation(fig, animate, frames=len(data_arrays), 
                        interval=interval, blit=False)
    
    # Save as GIF
    writer = PillowWriter(fps=fps)
    anim.save(output_file, writer=writer)
    plt.close()

def main():
    # Load all correlation data
    print("Loading correlation data...")
    data_arrays = load_correlation_data()
    
    if not data_arrays:
        print("No correlation data files found!")
        return
    
    print(f"Found {len(data_arrays)} correlation files")
    
    # Create the GIF
    print("Creating correlation animation...")
    create_correlation_gif(data_arrays)
    print("Animation saved as 'correlation_animation.gif'")

if __name__ == "__main__":
    main()
