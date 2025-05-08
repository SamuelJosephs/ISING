import numpy as np
import matplotlib.pyplot as plt 

def plot_heatmap_from_csv(filepath: str, save_name: str, cmap: str = 'viridis') -> None:
    """
    Reads a CSV file at the given filepath, plots a heatmap of the data with a colorbar,
    and saves the figure to a file with the given save_name.

    Parameters:
    - filepath: Path to the CSV input file (string).
    - save_name: Filename (with extension) for saving the heatmap image (string).
    - cmap: Colormap for the heatmap (optional, default 'viridis').

    Example:
        plot_heatmap_from_csv('data/my_data.csv', 'output/heatmap.png')
    """
    # Load the data from CSV using numpy
    data = np.loadtxt(filepath, delimiter=',')
    data = np.transpose(data)
    # Create the plot
    fig, ax = plt.subplots()
    cax = ax.imshow(data, aspect='auto', cmap=cmap)
    fig.colorbar(cax, ax=ax)

    # Add labels and title if desired
    ax.set_xlabel('Columns')
    ax.set_ylabel('Rows')
    ax.set_title('Heatmap of CSV Data')

    # Save the figure
    plt.tight_layout()
    plt.savefig(save_name)
    plt.close(fig)


plot_heatmap_from_csv("./density_mask.csv", "density_mask.png")
plot_heatmap_from_csv("./density_matrix.csv", "density_matrix.png")
plot_heatmap_from_csv("./visited_array.csv", "visited_array.png")

