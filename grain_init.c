#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "grain_init.h"
#include "rng.h"

#define INIT_FULLY_BULK 1
#define INIT_RANDOM 2
#define INIT_FROM_FILE 3

void initialize_grains_grid(int **grains_grid, int grid_size_x, int grid_size_y, int init_type, int q, const char *file_path, RngState *rng) {
    // Initialize all grid values to a default value (e.g., 0)
    for (int i = 0; i < grid_size_x; i++) {
        for (int j = 0; j < grid_size_y; j++) {
            grains_grid[i][j] = 0; // Default value
        }
    }

    FILE *file = NULL; 

    switch (init_type) {
        case INIT_FULLY_BULK:
            // Fully bulk initialization
            for (int i = 0; i < grid_size_x; i++) {
                for (int j = 0; j < grid_size_y; j++) {
                    grains_grid[i][j] = 1;
                }
            }
            break;

        case INIT_RANDOM:
            // Random initialization
            if (q <= 0) {
                fprintf(stderr, "Error: q must be greater than 0 for random initialization.\n");
                exit(EXIT_FAILURE);
            }
            for (int i = 0; i < grid_size_x; i++) {
                for (int j = 0; j < grid_size_y; j++) {
                    grains_grid[i][j] = (int)(rng_next_u32(rng) % q) + 1;  // Random values in [1, q]
                }
            }
            break;

        case INIT_FROM_FILE:
            // File-based initialization
            file = fopen(file_path, "r");
            if (!file) {
                perror("Error opening grain structure file");
                exit(EXIT_FAILURE);
            }

            char line[100000];
            int i = 0;
            while (fgets(line, sizeof(line), file) && i < grid_size_x) {
                char *token = strtok(line, ",");
                int j = 0;
                while (token != NULL && j < grid_size_y) {
                    grains_grid[i][j] = atoi(token);
                    token = strtok(NULL, ",");
                    j++;
                }
                if (j < grid_size_y) {
                    fprintf(stderr, "Error: File row %d has fewer columns than expected.\n", i);
                    fclose(file);
                    exit(EXIT_FAILURE);
                }
                i++;
            }
            if (i < grid_size_x) {
                fprintf(stderr, "Error: File has fewer rows than expected.\n");
                fclose(file);
                exit(EXIT_FAILURE);
            }

            fclose(file);
            break;

        default:
            // Invalid initialization type
            fprintf(stderr, "Unknown initialization type\n");
            exit(EXIT_FAILURE);
    }

    // Validation: Ensure no uninitialized cells remain
    for (int i = 0; i < grid_size_x; i++) {
        for (int j = 0; j < grid_size_y; j++) {
            if (grains_grid[i][j] == 0) {
                fprintf(stderr, "Error: Uninitialized cell found at (%d, %d).\n", i, j);
                exit(EXIT_FAILURE);
            }
        }
    }
}
