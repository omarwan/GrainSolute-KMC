import numpy as np
import pandas as pd

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import ipywidgets as widgets
from IPython.display import display
steps = np.arange(0, 100000000, 10000)

for step in steps:
    directory = "./output2/"
    grains_file_path = f'{directory}/grains_grid_step_{step}.csv'
    solutes_file_path = f'{directory}/solutes_grid_step_{step}.csv'

    
    grains_df = pd.read_csv(grains_file_path, header=None)
    solutes_df = pd.read_csv(solutes_file_path, header=None)
    
    grains_grid = grains_df.to_numpy()
    solutes_grid = solutes_df.to_numpy()
    
    grid_size_x, grid_size_y = solutes_grid.shape[0:2]
    # Get the dimensions of the grid
    m, n = solutes_grid.shape[0:2]
    
    # Open a file to write the dump data
    with open(f'./dumps/solute_grid_{step}.dump', 'w') as file:
        # Write the header
        file.write(f"ITEM: TIMESTEP\n")
        file.write(f"0\n")
        file.write(f"ITEM: NUMBER OF ATOMS\n")
        file.write(f"{ m * n}\n")  # Updated to account for both grids
        file.write(f"ITEM: BOX BOUNDS pp pp pp\n")
        file.write(f"0.0 {n}.0\n0.0 {m}.0\n0.0 2.0\n")  # Updated to account for the added height
        file.write(f"ITEM: ATOMS id x y z\n")
    
        # Write the atom data
        id_counter = 1
        for i in range(m):
            for j in range(n):
                solute_value = solutes_grid[i, j]
                # Assuming z-coordinate is 0 for the first grid
                file.write(f"{solute_value} {j} {i} 0.0\n")
                id_counter += 1
    
    # Open a file to write the dump data
    with open(f'./dumps/grains_grid_{step}.dump', 'w') as file:
        # Write the header
        file.write(f"ITEM: TIMESTEP\n")
        file.write(f"0\n")
        file.write(f"ITEM: NUMBER OF ATOMS\n")
        file.write(f"{ m * n}\n")  # Updated to account for both grids
        file.write(f"ITEM: BOX BOUNDS pp pp pp\n")
        file.write(f"0.0 {n}.0\n0.0 {m}.0\n0.0 2.0\n")  # Updated to account for the added height
        file.write(f"ITEM: ATOMS id x y z\n")
    
        id_counter = 1
        for i in range(m):
            for j in range(n):
                grain_value = grains_grid[i, j]
                # Assuming z-coordinate is 0 for the first grid
                file.write(f"{grain_value} {j} {i} 0.0\n")
                id_counter += 1
    