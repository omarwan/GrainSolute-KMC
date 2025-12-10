import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import ipywidgets as widgets
from IPython.display import display
from matplotlib.colors import LinearSegmentedColormap, Normalize

# Define custom color map for values between 10 and 20
colors = [(0, 0, 0), (1, 1, 1)]  # Gradient from blue to green
custom_cmap = LinearSegmentedColormap.from_list('custom_cmap', colors, N=50)
# Normalize to restrict the color map to values 10 to 20
norm = Normalize(vmin=1, vmax=50)

q = 3
steps = np.arange(0, 1000000, 100000)
for step in steps:
    directory = "./output/"
    grains_file_path = f'{directory}/grains_grid_step_{step}.csv'
    solutes_file_path = f'{directory}/solutes_grid_step_{step}.csv'
    
    grains_df = pd.read_csv(grains_file_path, header=None)
    solutes_df = pd.read_csv(solutes_file_path, header=None)
    
    grains_grid = grains_df.to_numpy()
    solutes_grid = solutes_df.to_numpy()
    
    
    grid_size_x, grid_size_y = grains_grid.shape[0:2]
    
    solute_positions = np.argwhere(solutes_grid > 0)

    

    plt.figure(figsize=(5,5),dpi=200)
    plt.imshow(grains_grid, cmap='tab20_r')
    plt.scatter(solute_positions[:, 1], solute_positions[:, 0], color='k', s=10)
    plt.xticks([])
    plt.yticks([])
    print(len(solute_positions))

    plt.tight_layout()
    plt.show()
    
    # fig, ax = plt.subplots(2, 1, figsize=(8, 6))  # Adjusted figsize for better proportions
    
    # # Panel 1: Grains grid with solute positions
    # ax[0].imshow(grains_grid, cmap='gray', norm=norm)
    # ax[0].scatter(solute_positions[:, 1], solute_positions[:, 0], color='red', s=0.4, label="Solutes")
    # ax[0].set_xticks([])
    # ax[0].set_yticks([])
    # ax[0].invert_yaxis()  
    
    # # Panel 2: Average number of solutes per line (X-axis)
    # rows, cols = grains_grid.shape
    # line_solute_counts = np.bincount(solute_positions[:, 1], minlength=cols)  # Count solutes per column
    # average_solutes = line_solute_counts / rows  # Compute average per column along X-axis
    
    # y_positions = np.arange(cols)
    # ax[1].plot(y_positions, average_solutes, color='blue', linewidth=2, label='Average Solutes per Line')
    # ax[1].set_xticks([])
    # ax[1].set_yticks([])
    # ax[1].set_xlim([0,511])
    # ax[1].set_ylim([0,0.5])
    # # Adjust layout and show the plot
    # plt.tight_layout()
    # plt.show()
