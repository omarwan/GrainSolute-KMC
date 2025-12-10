import numpy as np
import os
import matplotlib.pyplot as plt

# Parameters
grid_size_x = 256
grid_size_y = 1024
circle_label = 1
matrix_label = 25
solute_label = 1
density_step = 0.05
densities = np.arange(0.00,1.01, density_step)
folder_name = "generated_data_256x1024"
boundary_width = 32

# Create folder
os.makedirs(folder_name, exist_ok=True)

# Generate structure_grains.csv
structure_grains = np.full((grid_size_x, grid_size_y), matrix_label, dtype=int)

for i in range(grid_size_x):
    for j in range(grid_size_y):
        structure_grains[i, j] = 1
        if j < boundary_width or j > grid_size_y - boundary_width:
            structure_grains[i, j] = 2

np.savetxt(f"{folder_name}/structure_grains.csv", structure_grains, fmt="%d", delimiter=",")




# Generate structure_solutes_{density}.csv for each density
for density in densities:
    structure_solutes = np.zeros((grid_size_x, grid_size_y), dtype=int)
    num_solutes = int(density * grid_size_y * grid_size_x)
    
    # Flatten indices and select top-left to bottom-right order
    solute_indices = np.arange(grid_size_x * grid_size_y)
    np.random.shuffle(solute_indices)  # Shuffle to ensure random selection
    selected_indices = np.sort(solute_indices[:num_solutes])  # Sort to ensure top-left to bottom-right order
    
    # Assign unique labels to the selected solutes
    for i, idx in enumerate(selected_indices):
        structure_solutes.flat[idx] = i + 1  # Unique label starting from 1
    plt.figure()
    plt.imshow(structure_grains)
    solute_positions = np.argwhere(structure_solutes > 0)
    
    plt.scatter(solute_positions[:, 1], solute_positions[:, 0],  color='r',  s=0.1)

    # Save the solute structure
    np.savetxt(f"{folder_name}/structure_solutes_c{density:.2f}.csv", structure_solutes, fmt="%d", delimiter=",")



# Confirm creation
os.listdir(folder_name)
plt.imshow(structure_grains)
solute_positions = np.argwhere(structure_solutes > 0)

plt.scatter(solute_positions[:, 1], solute_positions[:, 0],  color='r',  s=0.1)
