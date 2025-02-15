import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import glob
import re
import pandas as pd
from scipy.optimize import curve_fit

def exp_decay(x, A, xi):
    """Exponential decay function: A * exp(-x/xi)"""
    return A * np.exp(-x/xi)

def extract_iteration(filename):
    """Extract iteration number from filename."""
    match = re.search(r'correlation_(\d+)$', filename)
    if match:
        return int(match.group(1))
    return 0

def read_correlation_data():
    """Read all correlation files and return sorted data."""
    files = glob.glob('./correlation_results/correlation_*')
    data_dict = {}
    
    for file in files:
        iteration = extract_iteration(file)
        try:
            data = np.loadtxt(file)
            df = pd.DataFrame(data, columns=['position', 'correlation'])
            df = df.sort_values('position')
            data_dict[iteration] = df
        except Exception as e:
            print(f"Error reading file {file}: {e}")
    
    sorted_iterations = sorted(data_dict.keys())
    return [data_dict[i] for i in sorted_iterations], sorted_iterations

def fit_exponential(df):
    """Fit exponential decay to data, excluding R=0 point for fitting only."""
    df_fit = df[df['position'] > 0].copy()
    
    try:
        popt, _ = curve_fit(exp_decay, 
                           df_fit['position'], 
                           df_fit['correlation'],
                           p0=[1.0, 1.0],
                           bounds=([0, 0], [np.inf, np.inf]))
        return popt
    except RuntimeError:
        print("Fitting failed")
        return None

def create_animation():
    """Create and save the correlation animation."""
    data_list, iterations = read_correlation_data()
    
    if not data_list:
        print("No data files found!")
        return
    
    # Create figure with two subplots using gridspec_kw for height ratios
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), 
                                  gridspec_kw={'height_ratios': [2, 1]})
    fig.tight_layout(pad=3.0)
    
    # First subplot - Correlation vs Position
    data_line, = ax1.plot([], [], 'b.', label='Data')
    fit_line, = ax1.plot([], [], 'r-', label='Exp Fit')
    
    ax1.set_xlabel('Position')
    ax1.set_ylabel('Correlation')
    ax1.set_title('Ising Model Correlation vs Position')
    
    # Set axis limits for first subplot
    all_positions = np.concatenate([df['position'] for df in data_list])
    all_correlations = np.concatenate([df['correlation'] for df in data_list])
    
    ax1.set_xlim(min(all_positions), max(all_positions))
    ax1.set_ylim(min(all_correlations), max(all_correlations))
    ax1.grid(True)
    ax1.legend()
    
    # Text for fit parameters
    fit_text = ax1.text(0.02, 0.98, '', transform=ax1.transAxes, 
                       verticalalignment='top')
    
    # Second subplot - Xi vs Iteration with log scale
    xi_line, = ax2.plot([], [], 'g.-', label='ξ')
    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('ξ')
    ax2.set_title('Correlation Length (ξ) vs Iteration')
    ax2.set_yscale('log')
    ax2.grid(True, which="both")  # Grid lines for both major and minor ticks
    ax2.legend()
    
    # Initialize lists to store xi values
    xi_values = []
    iteration_values = []
    
    def init():
        """Initialize animation."""
        data_line.set_data([], [])
        fit_line.set_data([], [])
        xi_line.set_data([], [])
        fit_text.set_text('')
        return data_line, fit_line, xi_line, fit_text
    
    def animate(frame):
        """Animation function called for each frame."""
        df = data_list[frame]
        current_iteration = iterations[frame]
        
        # Plot correlation data
        data_line.set_data(df['position'], df['correlation'])
        
        # Fit exponential and plot
        popt = fit_exponential(df)
        if popt is not None:
            A, xi = popt
            x_fit = np.linspace(min(df['position'][df['position'] > 0]), 
                              max(df['position']), 100)
            y_fit = exp_decay(x_fit, A, xi)
            fit_line.set_data(x_fit, y_fit)
            
            # Update fit parameters text
            fit_text.set_text(f'A = {A:.3f}\nξ = {xi:.3f}')
            
            # Update xi vs iteration plot
            xi_values.append(xi)
            iteration_values.append(current_iteration)
            xi_line.set_data(iteration_values, xi_values)
            
            # Update xi plot limits for log scale
            ax2.set_xlim(min(iterations), max(iterations))
            if len(xi_values) > 0:
                y_min = max(min(xi_values) * 0.5, 1e-10)  # Avoid zero in log plot
                y_max = max(xi_values) * 2
                ax2.set_ylim(y_min, y_max)
        
        ax1.set_title(f'Ising Model Correlation vs Position (Iteration {current_iteration})')
        return data_line, fit_line, xi_line, fit_text
    
    # Create animation
    anim = FuncAnimation(fig, animate, init_func=init,
                        frames=len(data_list), 
                        interval=200,
                        blit=True)
    
    # Save animation
    anim.save('correlation_animation.gif', writer='pillow')
    plt.close()

if __name__ == "__main__":
    create_animation()
    print("Animation saved as 'correlation_animation.gif'")
