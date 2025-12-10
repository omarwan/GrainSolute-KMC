#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h> // For getpid()
#include "rng.h"

#define INIT_FULLY_BULK 1
#define INIT_RANDOM 2
#define INIT_FROM_FILE 3

void initialize_solutes_grid(int **solutes_grid, int grid_size_x, int grid_size_y, int init_type, double density, const char *file_path, RngState *rng) {
    FILE *file = NULL;
    char line[100000];
    int i, j;
    int solute_counter = 1; // Counter for assigning unique numbers to solutes

    // Initialize all grid values to a default value (e.g., 0)
    for (int i = 0; i < grid_size_x; i++) {
        for (int j = 0; j < grid_size_y; j++) {
            solutes_grid[i][j] = 0; // Default value
        }
    }

    switch (init_type) {
        case INIT_FROM_FILE:
            file = fopen(file_path, "r");
            if (!file) {
                perror("Error opening solute structure file");
                return; // Use return here, consider modifying function signature to return an error code.
            }

            i = 0;
            while (fgets(line, sizeof(line), file) && i < grid_size_x) {
                char *token = strtok(line, ",");
                j = 0;
                while (token != NULL && j < grid_size_y) {
                    solutes_grid[i][j] = atoi(token);
                    token = strtok(NULL, ",");
                    j++;
                }
                i++;
            }
            fclose(file); // Ensure file is closed here
            break;

        default: // Treat any other case as INIT_RANDOM
            {
                int total_cells = grid_size_x * grid_size_y;
                int num_solutes = (int)(density * total_cells + 0.5); // Round to nearest integer
                int *values = (int *)malloc(total_cells * sizeof(int));

                // Step 1: Fill the array with exact number of 1's and 0's
                for (i = 0; i < total_cells; i++) {
                    values[i] = i < num_solutes ? 1 : 0;
                }

                // Step 2: Shuffle the array using Fisher-Yates shuffle
                for (i = total_cells - 1; i > 0; i--) {
                    int j = (int)(rng_next_u32(rng) % (uint32_t)(i + 1));
                    int temp = values[i];
                    values[i] = values[j];
                    values[j] = temp;
                }

                // Step 3: Copy the shuffled values back into the 2D solutes_grid
                for (i = 0; i < grid_size_x; i++) {
                    for (j = 0; j < grid_size_y; j++) {
                        if (values[i * grid_size_y + j] == 1) {
                            solutes_grid[i][j] = solute_counter++;
                        }
                    }
                }
                free(values); // Free the 1D array
            }
            break;
    }
}
