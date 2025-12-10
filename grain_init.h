#ifndef GRAIN_INIT_H
#define GRAIN_INIT_H

#include <stdio.h>
#include "rng.h"

// Definition of initialization types
#define INIT_FULLY_BULK 1
#define INIT_RANDOM 2
#define INIT_FROM_FILE 3

// Function declaration
void initialize_grains_grid(int **grains_grid, int grid_size_x, int grid_size_y, int init_type, int q, const char *file_path, RngState *rng);

#endif
